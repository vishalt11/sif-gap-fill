"""
Predict monthly crop-pure winter-wheat SIF for selected Bavarian NUTS3 regions.

The trained model expects 19 predictor channels on a 20 m grid and a
200 x 200 pixel (4 km x 4 km) input window. This script:

1. Uses the existing monthly Sentinel-2 WASP composites.
2. Arithmetic-averages GLASS FAPAR composites whose nominal dates fall in
   each calendar month.
3. Arithmetic-averages every available QA-accepted VIIRS daily PAR raster
   in each calendar month.
4. Builds the annual crop channels and a strict winter-wheat mask. A 20 m
   output pixel is strict winter wheat only when its four underlying 10 m
   crop pixels are all winter wheat.
5. Runs batched CNN inference only for output blocks containing strict
   winter-wheat pixels inside the selected NUTS3 regions.
6. Writes monthly NUTS3 means without saving full-tile 20 m prediction maps.

PAR note:
The local PAR archive may contain only selected OCO-2 dates rather than every
calendar day. The script records par_days_available so incomplete monthly
means remain explicit.

Before running this file once, run:

    source("export_nuts3_regions_for_wheat_inference.R")

That converts the existing sf RDS into a GeoPackage readable by GeoPandas.
"""

from __future__ import annotations

from contextlib import ExitStack
from dataclasses import dataclass
import calendar
import math
from pathlib import Path
import re
from typing import Iterable

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from rasterio.enums import Resampling
from rasterio.features import rasterize
from rasterio.transform import Affine, array_bounds
from rasterio.warp import reproject, transform_bounds
from rasterio.windows import Window, from_bounds
import torch
from torch import nn

import prepare_sentinel2_multisif_cnn_chips as sentinel


# ---------------------------------------------------------------------------
# Configuration

NUTS3_PATH = Path("data/nuts3_regions_80pct_in_mgrs.gpkg")
NUTS3_LAYER = "nuts3_regions_80pct_in_mgrs"

SENTINEL_DIR = Path("data/geodes_wasp_zips")
FAPAR_DIR = Path("data/glass_geotiff/fapar")
PAR_DIR = Path("data/viirs_vnp18a2_daily_mean_par_germany_native")
CROP_DIR = Path("data/crop_type_tif")

CHECKPOINT_PATH = Path(
    "data/cnn_sentinel2_chips/results/"
    "sentinel2_spatial_aggregate_4km_unet.pt"
)

OUTPUT_DIR = Path("data/winter_wheat_yield_model/monthly_crop_pure_sif")
LONG_OUTPUT_PATH = OUTPUT_DIR / "nuts3_monthly_strict_wheat_sif_long.csv"
WIDE_OUTPUT_PATH = OUTPUT_DIR / "nuts3_monthly_strict_wheat_sif_wide.csv"
SCENE_LOG_PATH = OUTPUT_DIR / "scene_processing_log.csv"

TARGET_TILES = ("32UNA", "32UPU", "32UQV")
YEARS = (2019, 2020, 2021, 2022, 2024)
MONTHS = (3, 4, 5, 6, 7)

MODEL_RES_M = 20.0
MODEL_SIZE = 200
OUTPUT_BLOCK_SIZE = 100
MODEL_EDGE_MARGIN = (MODEL_SIZE - OUTPUT_BLOCK_SIZE) // 2

# Reduce this if GPU memory is limited. CPU inference also uses this batch size.
INFERENCE_BATCH_SIZE = 8

# Use None for all scenes. Set an integer for a smoke test.
MAX_SCENES: int | None = None

# Resume completed tile/year/month scenes from the long output table.
RESUME = True

# Full prediction maps are intentionally not written.
SAVE_FULL_MAPS = False

WINTER_WHEAT_CODE = sentinel.CROP_CODES["winter_wheat"]
STRICT_WHEAT_THRESHOLD = 1.0 - 1e-6
PAR_ACCEPTED_QA_CODES = sentinel.PAR_ACCEPTED_QA_CODES
CHANNEL_NAMES = sentinel.CHANNEL_NAMES

# The shared formulas use these globals when allocating 4 km arrays.
sentinel.CHIP_SIZE_M = MODEL_SIZE * MODEL_RES_M
sentinel.CHIP_RES_M = MODEL_RES_M
sentinel.CHIP_SIZE = MODEL_SIZE
sentinel.FAPAR_DIR = FAPAR_DIR


# ---------------------------------------------------------------------------
# Model architecture


def group_count(channels: int) -> int:
    for groups in (8, 4, 2, 1):
        if channels % groups == 0:
            return groups
    return 1


class ConvBlock(nn.Module):
    def __init__(self, in_channels: int, out_channels: int):
        super().__init__()
        groups = group_count(out_channels)
        self.block = nn.Sequential(
            nn.Conv2d(in_channels, out_channels, 3, padding=1, bias=False),
            nn.GroupNorm(groups, out_channels),
            nn.SiLU(inplace=True),
            nn.Conv2d(out_channels, out_channels, 3, padding=1, bias=False),
            nn.GroupNorm(groups, out_channels),
            nn.SiLU(inplace=True),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.block(x)


class ThreeLevelUNet(nn.Module):
    def __init__(self, in_channels: int, base_channels: int = 16):
        super().__init__()
        self.enc1 = ConvBlock(in_channels, base_channels)
        self.pool1 = nn.MaxPool2d(2)
        self.enc2 = ConvBlock(base_channels, base_channels * 2)
        self.pool2 = nn.MaxPool2d(2)
        self.enc3 = ConvBlock(base_channels * 2, base_channels * 4)
        self.pool3 = nn.MaxPool2d(2)

        self.bottleneck = ConvBlock(base_channels * 4, base_channels * 8)

        self.up3 = nn.ConvTranspose2d(
            base_channels * 8, base_channels * 4, 2, stride=2
        )
        self.dec3 = ConvBlock(base_channels * 8, base_channels * 4)
        self.up2 = nn.ConvTranspose2d(
            base_channels * 4, base_channels * 2, 2, stride=2
        )
        self.dec2 = ConvBlock(base_channels * 4, base_channels * 2)
        self.up1 = nn.ConvTranspose2d(
            base_channels * 2, base_channels, 2, stride=2
        )
        self.dec1 = ConvBlock(base_channels * 2, base_channels)
        self.out = nn.Conv2d(base_channels, 1, kernel_size=1)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        e1 = self.enc1(x)
        e2 = self.enc2(self.pool1(e1))
        e3 = self.enc3(self.pool2(e2))
        bottleneck = self.bottleneck(self.pool3(e3))

        d3 = self.dec3(torch.cat([self.up3(bottleneck), e3], dim=1))
        d2 = self.dec2(torch.cat([self.up2(d3), e2], dim=1))
        d1 = self.dec1(torch.cat([self.up1(d2), e1], dim=1))
        return self.out(d1)


@dataclass(frozen=True)
class ModelBundle:
    model: nn.Module
    device: torch.device
    channel_mean: np.ndarray
    channel_std: np.ndarray
    target_mean: float
    target_std: float
    normalize_target: bool
    calibration_intercept: float
    calibration_slope: float


def load_checkpoint() -> ModelBundle:
    if not CHECKPOINT_PATH.exists():
        raise FileNotFoundError(CHECKPOINT_PATH)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    try:
        checkpoint = torch.load(
            CHECKPOINT_PATH,
            map_location=device,
            weights_only=False,
        )
    except TypeError:
        # Compatibility with older PyTorch releases that predate weights_only.
        checkpoint = torch.load(CHECKPOINT_PATH, map_location=device)

    checkpoint_channels = [str(value) for value in checkpoint["channel_names"]]
    if checkpoint_channels != CHANNEL_NAMES:
        raise ValueError(
            "Checkpoint channel order differs from the preparation helper.\n"
            f"Checkpoint: {checkpoint_channels}\n"
            f"Expected:   {CHANNEL_NAMES}"
        )

    model = ThreeLevelUNet(
        in_channels=len(checkpoint_channels),
        base_channels=int(checkpoint["base_channels"]),
    ).to(device)
    model.load_state_dict(checkpoint["model_state_dict"])
    model.eval()

    channel_mean = np.asarray(checkpoint["channel_mean"], dtype=np.float32)
    channel_std = np.asarray(checkpoint["channel_std"], dtype=np.float32)
    if channel_mean.shape != (len(CHANNEL_NAMES),):
        raise ValueError(f"Unexpected channel mean shape: {channel_mean.shape}")
    if channel_std.shape != (len(CHANNEL_NAMES),):
        raise ValueError(f"Unexpected channel std shape: {channel_std.shape}")
    if not np.isfinite(channel_std).all() or (channel_std <= 0).any():
        raise ValueError("Checkpoint contains invalid channel standard deviations")

    calibration = checkpoint.get("calibration", {})
    return ModelBundle(
        model=model,
        device=device,
        channel_mean=channel_mean,
        channel_std=channel_std,
        target_mean=float(checkpoint.get("target_mean", 0.0)),
        target_std=float(checkpoint.get("target_std", 1.0)),
        normalize_target=bool(checkpoint.get("normalize_target", False)),
        calibration_intercept=float(calibration.get("intercept", 0.0)),
        calibration_slope=float(calibration.get("slope", 1.0)),
    )


# ---------------------------------------------------------------------------
# Grid and file helpers


@dataclass(frozen=True)
class GridSpec:
    row_off: int
    col_off: int
    height: int
    width: int
    transform: Affine
    crs: object

    @property
    def bounds(self) -> tuple[float, float, float, float]:
        return array_bounds(self.height, self.width, self.transform)


@dataclass(frozen=True)
class OutputBlock:
    output_row_start: int
    output_row_end: int
    output_col_start: int
    output_col_end: int
    input_row_start: int
    input_col_start: int

    @property
    def input_window(self) -> Window:
        return Window(
            self.input_col_start,
            self.input_row_start,
            MODEL_SIZE,
            MODEL_SIZE,
        )

    @property
    def output_slice_in_prediction(self) -> tuple[slice, slice]:
        row0 = self.output_row_start - self.input_row_start
        row1 = self.output_row_end - self.input_row_start
        col0 = self.output_col_start - self.input_col_start
        col1 = self.output_col_end - self.input_col_start
        return slice(row0, row1), slice(col0, col1)


def normalize_mgrs_tile(value: object) -> str:
    tile = str(value).strip().upper()
    if tile.startswith("T"):
        tile = tile[1:]
    if not re.fullmatch(r"\d{2}[A-Z]{3}", tile):
        raise ValueError(f"Invalid MGRS tile: {value}")
    return tile


def require_columns(
    table: pd.DataFrame,
    columns: Iterable[str],
    label: str,
) -> None:
    missing = [column for column in columns if column not in table.columns]
    if missing:
        raise ValueError(f"{label} is missing columns: {missing}")


def sentinel_product_path(tile: str, year: int, month: int) -> Path | None:
    year_dir = SENTINEL_DIR / tile / str(year)
    if not year_dir.exists():
        return None

    pattern = f"SENTINEL2*_{year}{month:02d}*-*_L3A_T{tile}_*"
    matches = sorted(path for path in year_dir.glob(pattern) if path.is_dir())
    if not matches:
        return None
    if len(matches) != 1:
        raise ValueError(
            f"Expected one Sentinel product for {tile}, {year}-{month:02d}; "
            f"found {len(matches)}: {matches}"
        )
    return matches[0]


def product_date(product_path: Path) -> pd.Timestamp:
    match = re.search(r"_(\d{8})-", product_path.name)
    if match is None:
        raise ValueError(f"Cannot parse product date from {product_path.name}")
    return pd.to_datetime(match.group(1), format="%Y%m%d")


def crop_path(year: int) -> Path:
    path = CROP_DIR / f"croptypes_{year}.tif"
    if not path.exists():
        raise FileNotFoundError(path)
    return path


def fapar_doys_for_month(year: int, month: int) -> list[int]:
    return [
        doy
        for doy in sentinel.FAPAR_COMPOSITE_DOYS
        if (
            pd.Timestamp(year=year, month=1, day=1)
            + pd.Timedelta(days=doy - 1)
        ).month
        == month
    ]


def par_files_for_month(year: int, month: int) -> list[tuple[pd.Timestamp, Path, Path]]:
    year_dir = PAR_DIR / str(year)
    pattern = (
        f"VNP18A2.002_{year}-{month:02d}-??_"
        "Daily_Mean_PAR_VIIRS_Sinusoidal_native.tif"
    )
    rows = []
    for par_path in sorted(year_dir.glob(pattern)):
        match = re.search(r"_(\d{4}-\d{2}-\d{2})_Daily_Mean_PAR_", par_path.name)
        if match is None:
            continue
        date = pd.Timestamp(match.group(1))
        qa_path = year_dir / (
            f"VNP18A2.002_{date:%Y-%m-%d}_"
            "PAR_Quality_VIIRS_Sinusoidal_native.tif"
        )
        if not qa_path.exists():
            raise FileNotFoundError(qa_path)
        rows.append((date, par_path, qa_path))
    return rows


def expanded_window_for_bounds(
    bounds: tuple[float, float, float, float],
    transform: Affine,
    pad_pixels: int,
) -> Window:
    window = from_bounds(*bounds, transform=transform)
    col_start = math.floor(window.col_off) - pad_pixels
    row_start = math.floor(window.row_off) - pad_pixels
    col_stop = math.ceil(window.col_off + window.width) + pad_pixels
    row_stop = math.ceil(window.row_off + window.height) + pad_pixels
    return Window(
        col_start,
        row_start,
        max(1, col_stop - col_start),
        max(1, row_stop - row_start),
    )


def windows_intersect(
    left: float,
    bottom: float,
    right: float,
    top: float,
    other_left: float,
    other_bottom: float,
    other_right: float,
    other_top: float,
) -> bool:
    return not (
        right <= other_left
        or left >= other_right
        or top <= other_bottom
        or bottom >= other_top
    )


def roi_grid_from_regions(
    reference: rasterio.io.DatasetReader,
    regions: gpd.GeoDataFrame,
) -> GridSpec:
    projected = regions.to_crs(reference.crs)
    xmin, ymin, xmax, ymax = projected.total_bounds
    region_window = from_bounds(xmin, ymin, xmax, ymax, reference.transform)

    # Include enough context for every 200 x 200 model input assembled around
    # the 100 x 100 non-overlapping output blocks.
    pad = MODEL_SIZE
    col0 = max(0, math.floor(region_window.col_off) - pad)
    row0 = max(0, math.floor(region_window.row_off) - pad)
    col1 = min(
        reference.width,
        math.ceil(region_window.col_off + region_window.width) + pad,
    )
    row1 = min(
        reference.height,
        math.ceil(region_window.row_off + region_window.height) + pad,
    )
    window = Window(col0, row0, col1 - col0, row1 - row0)
    return GridSpec(
        row_off=int(row0),
        col_off=int(col0),
        height=int(row1 - row0),
        width=int(col1 - col0),
        transform=reference.window_transform(window),
        crs=reference.crs,
    )


# ---------------------------------------------------------------------------
# Persistent Sentinel and crop readers


class SceneRasterReader:
    def __init__(self, product_path: Path, crop_raster_path: Path):
        self.product_path = product_path
        self.crop_raster_path = crop_raster_path
        self.stack = ExitStack()
        self.band_sources: dict[str, rasterio.io.DatasetReader] = {}
        self.flag_sources: dict[str, rasterio.io.DatasetReader] = {}
        self.crop_source: rasterio.io.DatasetReader | None = None

    def __enter__(self) -> "SceneRasterReader":
        for band_id, flag_resolution in sentinel.SENTINEL_BANDS.items():
            band_path = sentinel.sentinel_band_path(self.product_path, band_id)
            flag_path = sentinel.sentinel_flag_path(
                self.product_path,
                flag_resolution,
            )
            band = self.stack.enter_context(rasterio.open(band_path))
            flag = self.stack.enter_context(rasterio.open(flag_path))
            if (
                band.crs != flag.crs
                or band.transform != flag.transform
                or band.width != flag.width
                or band.height != flag.height
            ):
                raise ValueError(
                    f"Band and flag grids differ: {band_path.name}, "
                    f"{flag_path.name}"
                )
            self.band_sources[band_id] = band
            self.flag_sources[band_id] = flag

        self.crop_source = self.stack.enter_context(
            rasterio.open(self.crop_raster_path)
        )
        return self

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        self.stack.close()

    def read_reflectance(
        self,
        band_id: str,
        dst_transform: Affine,
        dst_crs,
    ) -> np.ndarray:
        band_src = self.band_sources[band_id]
        flag_src = self.flag_sources[band_id]
        dst_bounds = array_bounds(MODEL_SIZE, MODEL_SIZE, dst_transform)
        source_bounds = transform_bounds(
            dst_crs,
            band_src.crs,
            *dst_bounds,
            densify_pts=21,
        )
        window = expanded_window_for_bounds(
            source_bounds,
            band_src.transform,
            pad_pixels=1,
        )
        band = band_src.read(
            1,
            window=window,
            boundless=True,
            fill_value=sentinel.REFLECTANCE_NODATA,
        )
        flag = flag_src.read(1, window=window, boundless=True, fill_value=0)
        source_transform = band_src.window_transform(window)

        valid = (
            (flag == sentinel.LAND_FLAG_VALUE)
            & (band != sentinel.REFLECTANCE_NODATA)
            & np.isfinite(band)
        )
        source = np.full(band.shape, sentinel.WARP_NODATA, dtype=np.float32)
        source[valid] = (
            band[valid].astype(np.float32)
            / sentinel.REFLECTANCE_QUANTIFICATION_VALUE
        )

        destination = np.full(
            (MODEL_SIZE, MODEL_SIZE),
            np.nan,
            dtype=np.float32,
        )
        reproject(
            source=source,
            destination=destination,
            src_transform=source_transform,
            src_crs=band_src.crs,
            src_nodata=sentinel.WARP_NODATA,
            dst_transform=dst_transform,
            dst_crs=dst_crs,
            dst_nodata=np.nan,
            resampling=Resampling.average,
        )
        return destination

    def sentinel_indices(
        self,
        dst_transform: Affine,
        dst_crs,
    ) -> list[np.ndarray]:
        bands = {
            band_id: self.read_reflectance(band_id, dst_transform, dst_crs)
            for band_id in sentinel.SENTINEL_BANDS
        }
        ndmi = sentinel.normalized_difference(bands["B8"], bands["B11"])
        ndvi = sentinel.normalized_difference(bands["B8"], bands["B4"])
        evi = sentinel.enhanced_vegetation_index(
            bands["B8"],
            bands["B4"],
            bands["B2"],
        )
        nirv = bands["B8"] * ndvi
        ndre = sentinel.normalized_difference(bands["B8A"], bands["B5"])
        return [ndmi, ndvi, evi, nirv, ndre]

    def read_crop_values(
        self,
        dst_transform: Affine,
        dst_crs,
        dst_height: int,
        dst_width: int,
    ) -> tuple[np.ndarray, Affine, object]:
        if self.crop_source is None:
            raise RuntimeError("Crop reader is not open")
        src = self.crop_source
        dst_bounds = array_bounds(dst_height, dst_width, dst_transform)
        source_bounds = transform_bounds(
            dst_crs,
            src.crs,
            *dst_bounds,
            densify_pts=21,
        )
        window = expanded_window_for_bounds(
            source_bounds,
            src.transform,
            pad_pixels=2,
        )
        values = src.read(
            1,
            window=window,
            boundless=True,
            fill_value=sentinel.NON_CROP_CODE,
        )
        return values, src.window_transform(window), src.crs

    def crop_fraction(
        self,
        crop_values: np.ndarray,
        crop_transform: Affine,
        crop_crs,
        codes: Iterable[int],
        dst_transform: Affine,
        dst_crs,
        dst_height: int,
        dst_width: int,
    ) -> np.ndarray:
        source = np.isin(crop_values, list(codes)).astype(np.float32)
        destination = np.zeros((dst_height, dst_width), dtype=np.float32)
        reproject(
            source=source,
            destination=destination,
            src_transform=crop_transform,
            src_crs=crop_crs,
            dst_transform=dst_transform,
            dst_crs=dst_crs,
            dst_nodata=0.0,
            resampling=Resampling.average,
        )
        return np.clip(destination, 0.0, 1.0)

    def crop_channels(
        self,
        month: int,
        dst_transform: Affine,
        dst_crs,
    ) -> list[np.ndarray]:
        crop_values, crop_transform, crop_crs = self.read_crop_values(
            dst_transform,
            dst_crs,
            MODEL_SIZE,
            MODEL_SIZE,
        )
        channels = [
            self.crop_fraction(
                crop_values,
                crop_transform,
                crop_crs,
                codes,
                dst_transform,
                dst_crs,
                MODEL_SIZE,
                MODEL_SIZE,
            )
            for codes in sentinel.CROP_GROUPS.values()
        ]
        channels.append(
            self.crop_fraction(
                crop_values,
                crop_transform,
                crop_crs,
                sentinel.active_codes_for_month(month),
                dst_transform,
                dst_crs,
                MODEL_SIZE,
                MODEL_SIZE,
            )
        )

        non_crop_source = (
            ~np.isin(crop_values, sentinel.KNOWN_CROP_CODES)
        ).astype(np.float32)
        non_crop = np.zeros((MODEL_SIZE, MODEL_SIZE), dtype=np.float32)
        reproject(
            source=non_crop_source,
            destination=non_crop,
            src_transform=crop_transform,
            src_crs=crop_crs,
            dst_transform=dst_transform,
            dst_crs=dst_crs,
            dst_nodata=1.0,
            resampling=Resampling.average,
        )
        channels.append(np.clip(non_crop, 0.0, 1.0))
        return [channel.astype(np.float32) for channel in channels]

    def strict_wheat_mask(self, grid: GridSpec) -> np.ndarray:
        crop_values, crop_transform, crop_crs = self.read_crop_values(
            grid.transform,
            grid.crs,
            grid.height,
            grid.width,
        )
        wheat_fraction = self.crop_fraction(
            crop_values,
            crop_transform,
            crop_crs,
            [WINTER_WHEAT_CODE],
            grid.transform,
            grid.crs,
            grid.height,
            grid.width,
        )
        return wheat_fraction >= STRICT_WHEAT_THRESHOLD


# ---------------------------------------------------------------------------
# Region and monthly coarse predictor grids


def load_regions() -> gpd.GeoDataFrame:
    if not NUTS3_PATH.exists():
        raise FileNotFoundError(
            f"{NUTS3_PATH} does not exist. Run "
            "export_nuts3_regions_for_wheat_inference.R first."
        )
    regions = gpd.read_file(NUTS3_PATH, layer=NUTS3_LAYER)
    require_columns(
        regions,
        [
            "nuts_id",
            "nuts3",
            "mgrs_tile",
            "nuts3_area_m2",
            "overlap_area_m2",
            "nuts3_area_in_mgrs_pct",
            "geometry",
        ],
        "NUTS3 regions",
    )
    regions = regions.loc[~regions.geometry.is_empty & regions.geometry.notna()].copy()
    regions["mgrs_tile"] = regions["mgrs_tile"].map(normalize_mgrs_tile)
    regions = regions.loc[regions["mgrs_tile"].isin(TARGET_TILES)].copy()
    if regions.empty:
        raise ValueError("No target NUTS3 regions remain after tile filtering")
    if regions["nuts_id"].duplicated().any():
        duplicates = regions.loc[
            regions["nuts_id"].duplicated(keep=False),
            ["nuts_id", "mgrs_tile"],
        ]
        raise ValueError(f"Duplicate NUTS3 IDs:\n{duplicates}")
    return regions.reset_index(drop=True)


def rasterize_regions(
    regions: gpd.GeoDataFrame,
    grid: GridSpec,
) -> tuple[np.ndarray, dict[int, dict]]:
    projected = regions.to_crs(grid.crs).reset_index(drop=True)
    region_lookup: dict[int, dict] = {}
    shapes = []
    for index, row in projected.iterrows():
        code = index + 1
        region_lookup[code] = {
            "nuts_id": row["nuts_id"],
            "nuts3": row["nuts3"],
            "mgrs_tile": row["mgrs_tile"],
            "nuts3_area_m2": float(row["nuts3_area_m2"]),
            "overlap_area_m2": float(row["overlap_area_m2"]),
            "nuts3_area_in_mgrs_pct": float(row["nuts3_area_in_mgrs_pct"]),
        }
        shapes.append((row.geometry, code))

    region_ids = rasterize(
        shapes,
        out_shape=(grid.height, grid.width),
        transform=grid.transform,
        fill=0,
        all_touched=False,
        dtype=np.int16,
    )
    return region_ids, region_lookup


def warp_continuous_to_grid(
    path: Path,
    grid: GridSpec,
    valid_min: float,
    valid_max: float,
    qa_path: Path | None = None,
    accepted_qa: Iterable[int] | None = None,
) -> np.ndarray | None:
    with ExitStack() as stack:
        src = stack.enter_context(rasterio.open(path))
        qa_src = (
            stack.enter_context(rasterio.open(qa_path))
            if qa_path is not None
            else None
        )
        if src.crs is None:
            raise ValueError(f"Raster has no CRS: {path}")
        if qa_src is not None and (
            src.crs != qa_src.crs
            or src.transform != qa_src.transform
            or src.width != qa_src.width
            or src.height != qa_src.height
        ):
            raise ValueError(f"PAR and QA grids differ: {path}, {qa_path}")

        source_bounds = transform_bounds(
            grid.crs,
            src.crs,
            *grid.bounds,
            densify_pts=21,
        )
        if not windows_intersect(
            *source_bounds,
            src.bounds.left,
            src.bounds.bottom,
            src.bounds.right,
            src.bounds.top,
        ):
            return None

        window = expanded_window_for_bounds(
            source_bounds,
            src.transform,
            pad_pixels=2,
        )
        fill_value = (
            src.nodata if src.nodata is not None else sentinel.WARP_NODATA
        )
        values = src.read(
            1,
            window=window,
            boundless=True,
            fill_value=fill_value,
        ).astype(np.float32)
        source_transform = src.window_transform(window)

        valid = (
            np.isfinite(values)
            & (values >= valid_min)
            & (values <= valid_max)
        )
        if src.nodata is not None:
            valid &= values != src.nodata

        if qa_src is not None:
            qa = qa_src.read(
                1,
                window=window,
                boundless=True,
                fill_value=0,
            )
            valid &= np.isin(qa, list(accepted_qa or []))

        source = np.full(values.shape, sentinel.WARP_NODATA, dtype=np.float32)
        source[valid] = values[valid]
        destination = np.full(
            (grid.height, grid.width),
            np.nan,
            dtype=np.float32,
        )
        reproject(
            source=source,
            destination=destination,
            src_transform=source_transform,
            src_crs=src.crs,
            src_nodata=sentinel.WARP_NODATA,
            dst_transform=grid.transform,
            dst_crs=grid.crs,
            dst_nodata=np.nan,
            resampling=Resampling.bilinear,
        )
        return destination


def running_nanmean(
    arrays: Iterable[np.ndarray | None],
    shape: tuple[int, int],
) -> tuple[np.ndarray, int]:
    total = np.zeros(shape, dtype=np.float32)
    count = np.zeros(shape, dtype=np.uint16)
    n_arrays = 0
    for array in arrays:
        if array is None:
            continue
        finite = np.isfinite(array)
        total[finite] += array[finite].astype(np.float32, copy=False)
        count[finite] += 1
        n_arrays += 1
    result = np.full(shape, np.nan, dtype=np.float32)
    valid = count > 0
    result[valid] = (total[valid] / count[valid]).astype(np.float32)
    return result, n_arrays


def monthly_fapar_grid(
    year: int,
    month: int,
    grid: GridSpec,
) -> tuple[np.ndarray, list[int]]:
    doys = fapar_doys_for_month(year, month)
    if not doys:
        raise ValueError(f"No configured FAPAR composites for {year}-{month:02d}")

    used_doys = []

    def composite_arrays() -> Iterable[np.ndarray]:
        for doy in doys:
            tile_arrays = []
            for tile in sentinel.MODIS_TILES:
                path = sentinel.fapar_path(year, doy, tile)
                tile_arrays.append(
                    warp_continuous_to_grid(
                        path,
                        grid,
                        valid_min=0.0,
                        valid_max=1.0,
                    )
                )
            composite, n_tiles = running_nanmean(
                tile_arrays,
                (grid.height, grid.width),
            )
            if n_tiles > 0:
                used_doys.append(doy)
                yield composite

    monthly, _ = running_nanmean(
        composite_arrays(),
        (grid.height, grid.width),
    )
    return monthly, used_doys


def monthly_par_grid(
    year: int,
    month: int,
    grid: GridSpec,
) -> tuple[np.ndarray, list[pd.Timestamp]]:
    files = par_files_for_month(year, month)
    arrays = (
        warp_continuous_to_grid(
            par_path,
            grid,
            valid_min=sentinel.PAR_VALID_MIN,
            valid_max=sentinel.PAR_VALID_MAX,
            qa_path=qa_path,
            accepted_qa=PAR_ACCEPTED_QA_CODES,
        )
        for _, par_path, qa_path in files
    )
    monthly, n_arrays = running_nanmean(
        arrays,
        (grid.height, grid.width),
    )
    dates = [date for date, _, _ in files]
    if n_arrays != len(dates):
        raise ValueError(
            f"{len(dates) - n_arrays} PAR rasters did not intersect the scene grid"
        )
    return monthly, dates


# ---------------------------------------------------------------------------
# Window construction and inference


def partition_edges(start: int, stop: int, block_size: int) -> list[int]:
    if stop <= start:
        return []
    edges = [start]
    while edges[-1] + block_size < stop:
        edges.append(edges[-1] + block_size)
    if edges[-1] != stop:
        edges.append(stop)
    return edges


def centred_input_start(
    output_start: int,
    output_stop: int,
    raster_size: int,
) -> int:
    output_size = output_stop - output_start
    start = output_start - (MODEL_SIZE - output_size) // 2
    return min(max(start, 0), raster_size - MODEL_SIZE)


def candidate_output_blocks(
    reference: rasterio.io.DatasetReader,
    grid: GridSpec,
    target_mask: np.ndarray,
) -> list[OutputBlock]:
    target_rows, target_cols = np.nonzero(target_mask)
    if target_rows.size == 0:
        return []

    global_row_min = grid.row_off + int(target_rows.min())
    global_row_max = grid.row_off + int(target_rows.max()) + 1
    global_col_min = grid.col_off + int(target_cols.min())
    global_col_max = grid.col_off + int(target_cols.max()) + 1

    valid_row_start = MODEL_EDGE_MARGIN
    valid_row_stop = reference.height - MODEL_EDGE_MARGIN
    valid_col_start = MODEL_EDGE_MARGIN
    valid_col_stop = reference.width - MODEL_EDGE_MARGIN

    row_edges = partition_edges(
        valid_row_start,
        valid_row_stop,
        OUTPUT_BLOCK_SIZE,
    )
    col_edges = partition_edges(
        valid_col_start,
        valid_col_stop,
        OUTPUT_BLOCK_SIZE,
    )

    blocks = []
    for row0, row1 in zip(row_edges[:-1], row_edges[1:]):
        if row1 <= global_row_min or row0 >= global_row_max:
            continue
        for col0, col1 in zip(col_edges[:-1], col_edges[1:]):
            if col1 <= global_col_min or col0 >= global_col_max:
                continue

            intersect_row0 = max(row0, grid.row_off)
            intersect_row1 = min(row1, grid.row_off + grid.height)
            intersect_col0 = max(col0, grid.col_off)
            intersect_col1 = min(col1, grid.col_off + grid.width)
            if (
                intersect_row0 >= intersect_row1
                or intersect_col0 >= intersect_col1
            ):
                continue

            local_rows = slice(
                intersect_row0 - grid.row_off,
                intersect_row1 - grid.row_off,
            )
            local_cols = slice(
                intersect_col0 - grid.col_off,
                intersect_col1 - grid.col_off,
            )
            if not target_mask[local_rows, local_cols].any():
                continue

            blocks.append(
                OutputBlock(
                    output_row_start=row0,
                    output_row_end=row1,
                    output_col_start=col0,
                    output_col_end=col1,
                    input_row_start=centred_input_start(
                        row0,
                        row1,
                        reference.height,
                    ),
                    input_col_start=centred_input_start(
                        col0,
                        col1,
                        reference.width,
                    ),
                )
            )
    return blocks


def coarse_grid_slice(
    values: np.ndarray,
    grid: GridSpec,
    block: OutputBlock,
) -> np.ndarray:
    row0 = block.input_row_start - grid.row_off
    col0 = block.input_col_start - grid.col_off
    row1 = row0 + MODEL_SIZE
    col1 = col0 + MODEL_SIZE
    if row0 < 0 or col0 < 0 or row1 > grid.height or col1 > grid.width:
        raise ValueError(
            "Inference input window falls outside the prepared regional grid. "
            "Increase the ROI padding."
        )
    return values[row0:row1, col0:col1]


def build_input_chip(
    reader: SceneRasterReader,
    reference: rasterio.io.DatasetReader,
    grid: GridSpec,
    block: OutputBlock,
    month: int,
    monthly_fapar: np.ndarray,
    monthly_par: np.ndarray,
    bundle: ModelBundle,
) -> np.ndarray:
    transform = reference.window_transform(block.input_window)
    spectral_indices = reader.sentinel_indices(transform, reference.crs)
    crop_channels = reader.crop_channels(month, transform, reference.crs)

    fapar = coarse_grid_slice(monthly_fapar, grid, block)
    par = coarse_grid_slice(monthly_par, grid, block)
    apar = fapar * par
    month_angle = 2.0 * np.pi * month / 12.0
    constant_shape = (MODEL_SIZE, MODEL_SIZE)

    channels = [
        *spectral_indices,
        fapar,
        par,
        apar,
        *crop_channels,
        np.full(constant_shape, np.sin(month_angle), dtype=np.float32),
        np.full(constant_shape, np.cos(month_angle), dtype=np.float32),
    ]
    if len(channels) != len(CHANNEL_NAMES):
        raise RuntimeError(
            f"Built {len(channels)} channels; expected {len(CHANNEL_NAMES)}"
        )

    x = np.stack(channels, axis=0).astype(np.float32)
    x = (
        x - bundle.channel_mean[:, None, None]
    ) / bundle.channel_std[:, None, None]
    return np.nan_to_num(x, nan=0.0, posinf=0.0, neginf=0.0)


def denormalize_prediction(
    prediction: np.ndarray,
    bundle: ModelBundle,
) -> np.ndarray:
    if bundle.normalize_target:
        return prediction * bundle.target_std + bundle.target_mean
    return prediction


def block_target_arrays(
    block: OutputBlock,
    grid: GridSpec,
    region_ids: np.ndarray,
    strict_wheat: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, tuple[slice, slice]]:
    global_row0 = max(block.output_row_start, grid.row_off)
    global_row1 = min(block.output_row_end, grid.row_off + grid.height)
    global_col0 = max(block.output_col_start, grid.col_off)
    global_col1 = min(block.output_col_end, grid.col_off + grid.width)

    region_rows = slice(global_row0 - grid.row_off, global_row1 - grid.row_off)
    region_cols = slice(global_col0 - grid.col_off, global_col1 - grid.col_off)
    prediction_rows = slice(
        global_row0 - block.input_row_start,
        global_row1 - block.input_row_start,
    )
    prediction_cols = slice(
        global_col0 - block.input_col_start,
        global_col1 - block.input_col_start,
    )

    ids = region_ids[region_rows, region_cols]
    wheat = strict_wheat[region_rows, region_cols]
    return ids, wheat, (prediction_rows, prediction_cols)


def accumulate_predictions(
    predictions: np.ndarray,
    blocks: list[OutputBlock],
    grid: GridSpec,
    region_ids: np.ndarray,
    strict_wheat: np.ndarray,
    raw_sum: np.ndarray,
    calibrated_sum: np.ndarray,
    prediction_count: np.ndarray,
    bundle: ModelBundle,
) -> None:
    for prediction, block in zip(predictions, blocks):
        raw_map = denormalize_prediction(prediction, bundle)
        calibrated_map = (
            bundle.calibration_intercept
            + bundle.calibration_slope * raw_map
        )
        ids, wheat, prediction_slice = block_target_arrays(
            block,
            grid,
            region_ids,
            strict_wheat,
        )
        selected_ids = ids[wheat & (ids > 0)]
        if selected_ids.size == 0:
            continue

        raw_values = raw_map[prediction_slice][wheat & (ids > 0)]
        calibrated_values = calibrated_map[prediction_slice][wheat & (ids > 0)]
        n_codes = len(raw_sum)
        raw_sum += np.bincount(
            selected_ids,
            weights=raw_values,
            minlength=n_codes,
        )[:n_codes]
        calibrated_sum += np.bincount(
            selected_ids,
            weights=calibrated_values,
            minlength=n_codes,
        )[:n_codes]
        prediction_count += np.bincount(
            selected_ids,
            minlength=n_codes,
        )[:n_codes]


@torch.no_grad()
def infer_scene(
    bundle: ModelBundle,
    reader: SceneRasterReader,
    reference: rasterio.io.DatasetReader,
    grid: GridSpec,
    region_ids: np.ndarray,
    strict_wheat: np.ndarray,
    month: int,
    monthly_fapar: np.ndarray,
    monthly_par: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    target_mask = strict_wheat & (region_ids > 0)
    blocks = candidate_output_blocks(reference, grid, target_mask)
    n_codes = int(region_ids.max()) + 1
    raw_sum = np.zeros(n_codes, dtype=np.float64)
    calibrated_sum = np.zeros(n_codes, dtype=np.float64)
    prediction_count = np.zeros(n_codes, dtype=np.int64)

    pending_x: list[np.ndarray] = []
    pending_blocks: list[OutputBlock] = []

    def flush() -> None:
        if not pending_x:
            return
        x = torch.from_numpy(np.stack(pending_x)).to(
            bundle.device,
            non_blocking=True,
        )
        with torch.autocast(
            device_type=bundle.device.type,
            dtype=torch.float16,
            enabled=bundle.device.type == "cuda",
        ):
            prediction = bundle.model(x)[:, 0]
        predictions = prediction.float().cpu().numpy()
        accumulate_predictions(
            predictions,
            pending_blocks,
            grid,
            region_ids,
            strict_wheat,
            raw_sum,
            calibrated_sum,
            prediction_count,
            bundle,
        )
        pending_x.clear()
        pending_blocks.clear()

    for index, block in enumerate(blocks, start=1):
        pending_x.append(
            build_input_chip(
                reader,
                reference,
                grid,
                block,
                month,
                monthly_fapar,
                monthly_par,
                bundle,
            )
        )
        pending_blocks.append(block)
        if len(pending_x) >= INFERENCE_BATCH_SIZE:
            flush()
        if index % 25 == 0 or index == len(blocks):
            print(f"    inference blocks: {index:,} / {len(blocks):,}")
    flush()
    return raw_sum, calibrated_sum, prediction_count, len(blocks)


# ---------------------------------------------------------------------------
# Output assembly


def existing_completed_scenes() -> set[tuple[str, int, int]]:
    if not RESUME or not LONG_OUTPUT_PATH.exists():
        return set()
    existing = pd.read_csv(LONG_OUTPUT_PATH)
    require_columns(
        existing,
        ["mgrs_tile", "year", "month", "scene_status"],
        "Existing output",
    )
    complete = existing.loc[
        existing["scene_status"].isin(
            {
                "success",
                "missing_sentinel_product",
                "missing_monthly_par",
                "no_strict_wheat",
            }
        )
    ]
    return {
        (str(row.mgrs_tile), int(row.year), int(row.month))
        for row in complete.itertuples()
    }


def replace_scene_rows(rows: list[dict]) -> None:
    table = pd.DataFrame(rows)
    scene = table.iloc[0]
    if LONG_OUTPUT_PATH.exists():
        existing = pd.read_csv(LONG_OUTPUT_PATH)
        same_scene = (
            (existing["mgrs_tile"].astype(str) == str(scene["mgrs_tile"]))
            & (existing["year"].astype(int) == int(scene["year"]))
            & (existing["month"].astype(int) == int(scene["month"]))
        )
        table = pd.concat(
            [existing.loc[~same_scene], table],
            ignore_index=True,
        )
    table = table.sort_values(
        ["mgrs_tile", "year", "month", "nuts_id"],
        kind="stable",
    )
    table.to_csv(LONG_OUTPUT_PATH, index=False)


def missing_scene_rows(
    regions: gpd.GeoDataFrame,
    year: int,
    month: int,
    status: str,
    error: str = "",
) -> list[dict]:
    return [
        {
            "nuts_id": row.nuts_id,
            "nuts3": row.nuts3,
            "mgrs_tile": row.mgrs_tile,
            "year": year,
            "month": month,
            "month_name": calendar.month_name[month],
            "sentinel_product": "",
            "sentinel_product_date": "",
            "fapar_composite_doys": "",
            "fapar_composites_available": 0,
            "par_dates_available": "",
            "par_days_available": 0,
            "sif_wheat_strict_raw": np.nan,
            "sif_wheat_strict_calibrated": np.nan,
            "strict_wheat_pixels_in_tile": np.nan,
            "strict_wheat_pixels_predicted": 0,
            "strict_wheat_area_predicted_km2": 0.0,
            "prediction_coverage_pct": np.nan,
            "nuts3_area_in_mgrs_pct": row.nuts3_area_in_mgrs_pct,
            "n_inference_blocks": 0,
            "scene_status": status,
            "error": error,
        }
        for row in regions.itertuples()
    ]


def successful_scene_rows(
    region_lookup: dict[int, dict],
    region_ids: np.ndarray,
    strict_wheat: np.ndarray,
    raw_sum: np.ndarray,
    calibrated_sum: np.ndarray,
    prediction_count: np.ndarray,
    year: int,
    month: int,
    product_path: Path,
    fapar_doys: list[int],
    par_dates: list[pd.Timestamp],
    n_blocks: int,
) -> list[dict]:
    strict_ids = region_ids[strict_wheat & (region_ids > 0)]
    total_counts = np.bincount(
        strict_ids,
        minlength=len(raw_sum),
    )[: len(raw_sum)]

    rows = []
    for code, metadata in region_lookup.items():
        predicted_n = int(prediction_count[code])
        total_n = int(total_counts[code])
        if predicted_n > 0:
            raw_mean = float(raw_sum[code] / predicted_n)
            calibrated_mean = float(calibrated_sum[code] / predicted_n)
        else:
            raw_mean = np.nan
            calibrated_mean = np.nan

        coverage = (
            100.0 * predicted_n / total_n if total_n > 0 else np.nan
        )
        rows.append(
            {
                **metadata,
                "year": year,
                "month": month,
                "month_name": calendar.month_name[month],
                "sentinel_product": str(product_path),
                "sentinel_product_date": product_date(product_path).date(),
                "fapar_composite_doys": ";".join(map(str, fapar_doys)),
                "fapar_composites_available": len(fapar_doys),
                "par_dates_available": ";".join(
                    date.strftime("%Y-%m-%d") for date in par_dates
                ),
                "par_days_available": len(par_dates),
                "sif_wheat_strict_raw": raw_mean,
                "sif_wheat_strict_calibrated": calibrated_mean,
                "strict_wheat_pixels_in_tile": total_n,
                "strict_wheat_pixels_predicted": predicted_n,
                "strict_wheat_area_predicted_km2": (
                    predicted_n * MODEL_RES_M * MODEL_RES_M / 1_000_000.0
                ),
                "prediction_coverage_pct": coverage,
                "n_inference_blocks": n_blocks,
                "scene_status": (
                    "success" if predicted_n > 0 else "no_strict_wheat"
                ),
                "error": "",
            }
        )
    return rows


def write_wide_output() -> None:
    long = pd.read_csv(LONG_OUTPUT_PATH)
    successful = long.loc[long["scene_status"] == "success"].copy()
    successful["month_label"] = successful["month"].map(
        {
            3: "March",
            4: "April",
            5: "May",
            6: "June",
            7: "July",
        }
    )

    identifiers = [
        "nuts_id",
        "nuts3",
        "mgrs_tile",
        "year",
        "nuts3_area_in_mgrs_pct",
    ]
    base = long[identifiers].drop_duplicates().set_index(identifiers)
    raw = successful.pivot_table(
        index=identifiers,
        columns="month_label",
        values="sif_wheat_strict_raw",
        aggfunc="first",
    ).rename(columns=lambda value: f"SIF_{value}")
    calibrated = successful.pivot_table(
        index=identifiers,
        columns="month_label",
        values="sif_wheat_strict_calibrated",
        aggfunc="first",
    ).rename(columns=lambda value: f"SIF_{value}_calibrated")
    pixels = successful.pivot_table(
        index=identifiers,
        columns="month_label",
        values="strict_wheat_pixels_predicted",
        aggfunc="first",
    ).rename(columns=lambda value: f"wheat_pixels_{value}")

    wide = base.join(raw, how="left").join(calibrated, how="left")
    wide = wide.join(pixels, how="left")
    desired_months = ["March", "April", "May", "June", "July"]
    for month in desired_months:
        for column in (
            f"SIF_{month}",
            f"SIF_{month}_calibrated",
            f"wheat_pixels_{month}",
        ):
            if column not in wide.columns:
                wide[column] = np.nan

    ordered = [
        *(f"SIF_{month}" for month in desired_months),
        *(f"SIF_{month}_calibrated" for month in desired_months),
        *(f"wheat_pixels_{month}" for month in desired_months),
    ]
    wide = wide.reset_index()[identifiers + ordered]
    wide.to_csv(WIDE_OUTPUT_PATH, index=False)


def write_scene_log(scene_rows: list[dict]) -> None:
    scene = scene_rows[0]
    log_row = pd.DataFrame(
        [
            {
                "mgrs_tile": scene["mgrs_tile"],
                "year": scene["year"],
                "month": scene["month"],
                "scene_statuses": ";".join(
                    sorted({row["scene_status"] for row in scene_rows})
                ),
                "n_regions": len(scene_rows),
                "n_successful_regions": sum(
                    row["scene_status"] == "success" for row in scene_rows
                ),
                "par_days_available": scene["par_days_available"],
                "fapar_composites_available": scene[
                    "fapar_composites_available"
                ],
                "error": scene["error"],
            }
        ]
    )
    if SCENE_LOG_PATH.exists():
        existing = pd.read_csv(SCENE_LOG_PATH)
        same_scene = (
            (existing["mgrs_tile"].astype(str) == str(scene["mgrs_tile"]))
            & (existing["year"].astype(int) == int(scene["year"]))
            & (existing["month"].astype(int) == int(scene["month"]))
        )
        log_row = pd.concat(
            [existing.loc[~same_scene], log_row],
            ignore_index=True,
        )
    log_row = log_row.sort_values(
        ["mgrs_tile", "year", "month"],
        kind="stable",
    )
    log_row.to_csv(SCENE_LOG_PATH, index=False)


# ---------------------------------------------------------------------------
# Main workflow


def process_scene(
    tile: str,
    year: int,
    month: int,
    tile_regions: gpd.GeoDataFrame,
    bundle: ModelBundle,
) -> list[dict]:
    product_path = sentinel_product_path(tile, year, month)
    if product_path is None:
        print("  missing Sentinel product")
        return missing_scene_rows(
            tile_regions,
            year,
            month,
            "missing_sentinel_product",
        )

    b5_path = sentinel.sentinel_band_path(product_path, "B5")
    with rasterio.open(b5_path) as reference:
        if reference.crs is None:
            raise ValueError(f"Sentinel raster has no CRS: {b5_path}")
        if not (
            math.isclose(abs(reference.transform.a), MODEL_RES_M, abs_tol=1e-6)
            and math.isclose(abs(reference.transform.e), MODEL_RES_M, abs_tol=1e-6)
        ):
            raise ValueError(
                f"Expected a {MODEL_RES_M:g} m B5 grid: {b5_path}"
            )

        grid = roi_grid_from_regions(reference, tile_regions)
        region_ids, region_lookup = rasterize_regions(tile_regions, grid)

        print(
            f"  regional grid: {grid.height:,} x {grid.width:,} pixels; "
            f"{len(region_lookup)} NUTS3 regions"
        )
        monthly_fapar, used_fapar_doys = monthly_fapar_grid(
            year,
            month,
            grid,
        )
        monthly_par, par_dates = monthly_par_grid(year, month, grid)
        if not par_dates:
            return missing_scene_rows(
                tile_regions,
                year,
                month,
                "missing_monthly_par",
                "No available daily PAR rasters for this month",
            )

        with SceneRasterReader(product_path, crop_path(year)) as reader:
            strict_wheat = reader.strict_wheat_mask(grid)
            target_pixels = int(
                (strict_wheat & (region_ids > 0)).sum(dtype=np.int64)
            )
            print(f"  strict winter-wheat target pixels: {target_pixels:,}")
            if target_pixels == 0:
                return missing_scene_rows(
                    tile_regions,
                    year,
                    month,
                    "no_strict_wheat",
                )

            raw_sum, calibrated_sum, prediction_count, n_blocks = infer_scene(
                bundle,
                reader,
                reference,
                grid,
                region_ids,
                strict_wheat,
                month,
                monthly_fapar,
                monthly_par,
            )

    return successful_scene_rows(
        region_lookup,
        region_ids,
        strict_wheat,
        raw_sum,
        calibrated_sum,
        prediction_count,
        year,
        month,
        product_path,
        used_fapar_doys,
        par_dates,
        n_blocks,
    )


def main() -> None:
    if SAVE_FULL_MAPS:
        raise NotImplementedError(
            "Full-map output is intentionally disabled in this regional workflow"
        )

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    regions = load_regions()
    bundle = load_checkpoint()
    completed = existing_completed_scenes()

    print(f"Device: {bundle.device}")
    print(f"Loaded {len(regions)} NUTS3 regions")
    print(f"Model channels ({len(CHANNEL_NAMES)}): {CHANNEL_NAMES}")
    print(
        "PAR monthly means use all QA-valid daily rasters that are actually "
        "available in each local year/month folder."
    )

    scenes = [
        (tile, year, month)
        for tile in TARGET_TILES
        for year in YEARS
        for month in MONTHS
    ]
    if MAX_SCENES is not None:
        scenes = scenes[:MAX_SCENES]

    processed = 0
    for scene_index, (tile, year, month) in enumerate(scenes, start=1):
        scene_key = (tile, year, month)
        if scene_key in completed:
            print(
                f"[{scene_index}/{len(scenes)}] {tile} {year}-{month:02d}: "
                "already complete"
            )
            continue

        tile_regions = regions.loc[regions["mgrs_tile"] == tile].copy()
        if tile_regions.empty:
            print(
                f"[{scene_index}/{len(scenes)}] {tile} {year}-{month:02d}: "
                "no assigned regions"
            )
            continue

        print(f"[{scene_index}/{len(scenes)}] {tile} {year}-{month:02d}")
        try:
            rows = process_scene(
                tile,
                year,
                month,
                tile_regions,
                bundle,
            )
        except Exception as error:
            rows = missing_scene_rows(
                tile_regions,
                year,
                month,
                "failed",
                f"{type(error).__name__}: {error}",
            )
            replace_scene_rows(rows)
            write_scene_log(rows)
            raise

        replace_scene_rows(rows)
        write_scene_log(rows)
        processed += 1
        print(
            f"  wrote {len(rows)} region rows; "
            f"{sum(row['scene_status'] == 'success' for row in rows)} successful"
        )

    if LONG_OUTPUT_PATH.exists():
        write_wide_output()
    print(f"Processed {processed} new scenes")
    print(f"Long output: {LONG_OUTPUT_PATH}")
    print(f"Wide output: {WIDE_OUTPUT_PATH}")
    print(f"Scene log: {SCENE_LOG_PATH}")


if __name__ == "__main__":
    main()
