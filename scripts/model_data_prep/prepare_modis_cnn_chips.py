"""
Prepare 16x16 MODIS-scale raster chips for CNN/U-Net SIF modelling.

This script builds one multi-channel chip per OCO-2 SIF footprint:
  - FAPAR, EVI, NDVI from converted GLASS MODIS 250 m GeoTIFF tile files
  - PAR from one converted global GLASS MODIS 0.05 degree GeoTIFF file
  - APAR = FAPAR * PAR
  - crop-type fractions aggregated from 10 m crop rasters to the chip grid
  - month sin/cos and latitude/longitude constant channels
  - fractional OCO-2 footprint mask on the same chip grid

The temporal match uses the GLASS 8-day composite interval that contains the
SIF observation date. For example, DOY 033 covers DOY 033-040, DOY 041 covers
DOY 041-048, and so on.
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from functools import lru_cache
import glob
from pathlib import Path
from typing import Iterable

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from rasterio.enums import Resampling
from rasterio.features import rasterize
from rasterio.transform import Affine, from_origin
from rasterio.warp import reproject, transform_geom
from rasterio.windows import Window, from_bounds
from shapely import wkb
from shapely.geometry import MultiPoint, Polygon, box, mapping, shape


# ---------------------------------------------------------------------------
# Config

SIF_CSV_PATH = Path("data/338k_cnn_sif.csv")

GEOTIFF_ROOT = Path("data/glass_geotiff")
FAPAR_DIR = GEOTIFF_ROOT / "fapar"
EVI_DIR = GEOTIFF_ROOT / "evi"
NDVI_DIR = GEOTIFF_ROOT / "ndvi"
PAR_DIR = GEOTIFF_ROOT / "par"
CROP_DIR = Path("data/crop_type_tif")

OUTPUT_DIR = Path("data/cnn_modis_chips/modis_8day_250m_16x16")

YEARS = set(range(2019, 2025))
GLASS_DOYS = np.arange(33, 210, 8, dtype=np.int16)
MODIS_TILES = ("h18v03", "h18v04")

CHIP_SIZE = 16
CHIP_RES_M = 250.0
MASK_OVERSAMPLE = 8
SHARD_SIZE = 512

# Keep a stratified development set first. Set SAMPLE_MODE = "none" later to
# generate chips for all rows.
SAMPLE_MODE = "stratified"  # "stratified" or "none"
SAMPLES_PER_STATE_YEAR_MONTH = 500
RANDOM_SEED = 42

# Optional hard cap after stratified sampling. Keep None for the 50k-ish
# stratified development run.
MAX_ROWS: int | None = None
# Keep this at 1 for the first test run. Increase to 2-4 after the output looks
# sane. On Windows, each process keeps its own lightweight raster metadata cache.
N_WORKERS = 4

MODIS_SINUSOIDAL_CRS = (
    "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 "
    "+units=m +no_defs"
)

PRODUCTS = {
    "fapar": {
        "directory": FAPAR_DIR,
        "valid_min": 0.0,
        "valid_max": 1.0,
    },
    "evi": {
        "directory": EVI_DIR,
        "valid_min": -1.0,
        "valid_max": 1.0,
    },
    "ndvi": {
        "directory": NDVI_DIR,
        "valid_min": -1.0,
        "valid_max": 1.0,
    },
}


CROP_GROUPS = {
    "winter_wheat_fraction": [11],
    "winter_barley_fraction": [12],
    "rapeseed_fraction": [71],
    "maize_fraction": [30],
    "grassland_fraction": [82, 83],
    "other_crop_fraction": [13, 14, 21, 22, 23, 40, 50, 60, 81, 90, 100, 110, 111],
}

CHANNEL_NAMES = [
    "fapar",
    "evi",
    "ndvi",
    "par",
    "apar",
    *CROP_GROUPS.keys(),
    "non_crop_fraction",
    "month_sin",
    "month_cos",
    "latitude",
    "longitude",
]


# ---------------------------------------------------------------------------
# Small raster containers


@dataclass(frozen=True)
class RasterSource:
    path: Path
    transform: Affine
    crs: object
    bounds: tuple[float, float, float, float]
    nodata: float | None


def raster_bounds(transform: Affine, width: int, height: int) -> tuple[float, float, float, float]:
    west = transform.c
    north = transform.f
    east = west + width * transform.a
    south = north + height * transform.e
    return (min(west, east), min(south, north), max(west, east), max(south, north))


def chip_bounds(transform: Affine) -> tuple[float, float, float, float]:
    return raster_bounds(transform, CHIP_SIZE, CHIP_SIZE)


def intersects(a: tuple[float, float, float, float], b: tuple[float, float, float, float]) -> bool:
    return a[0] < b[2] and a[2] > b[0] and a[1] < b[3] and a[3] > b[1]


# ---------------------------------------------------------------------------
# GeoTIFF reading


def read_geotiff_source(path: Path) -> RasterSource:
    with rasterio.open(path) as src:
        bounds = (src.bounds.left, src.bounds.bottom, src.bounds.right, src.bounds.top)
        return RasterSource(
            path=path,
            transform=src.transform,
            crs=src.crs,
            bounds=bounds,
            nodata=src.nodata,
        )


def clean_product_array(arr: np.ndarray, product: str) -> np.ndarray:
    cfg = PRODUCTS[product]
    out = arr.astype(np.float32, copy=True)
    out[(out < cfg["valid_min"]) | (out > cfg["valid_max"])] = np.nan
    return out


def clean_par_array(arr: np.ndarray) -> np.ndarray:
    out = arr.astype(np.float32, copy=True)
    out[out < 0] = np.nan
    return out


def find_one_file(pattern: str) -> Path:
    matches = sorted(Path(match) for match in glob.glob(pattern))
    if len(matches) != 1:
        raise FileNotFoundError(f"Expected 1 file for pattern {pattern}, found {len(matches)}")
    return matches[0]


@lru_cache(maxsize=512)
def load_modis_tile(product: str, year: int, doy: int, tile: str) -> RasterSource:
    cfg = PRODUCTS[product]
    pattern = str(cfg["directory"] / tile / str(year) / f"*A{year}{doy:03d}.{tile}.*.tif")
    path = find_one_file(pattern)
    return read_geotiff_source(path)


@lru_cache(maxsize=256)
def load_par(year: int, doy: int) -> RasterSource:
    path = find_one_file(str(PAR_DIR / str(year) / f"*A{year}{doy:03d}.*.tif"))
    return read_geotiff_source(path)


# ---------------------------------------------------------------------------
# SIF loading and temporal matching


def build_sif_polygon(row: pd.Series) -> Polygon:
    coords = [
        (row["Lon_corner1"], row["Lat_corner1"]),
        (row["Lon_corner2"], row["Lat_corner2"]),
        (row["Lon_corner3"], row["Lat_corner3"]),
        (row["Lon_corner4"], row["Lat_corner4"]),
    ]
    polygon = Polygon(coords)
    if not polygon.is_valid:
        polygon = MultiPoint(coords).convex_hull
    return polygon


def containing_glass_doy(sif_doy: int) -> int:
    """Return the GLASS composite start DOY whose 8-day window contains sif_doy."""
    candidates = GLASS_DOYS[GLASS_DOYS <= sif_doy]
    if len(candidates) == 0:
        return int(GLASS_DOYS[0])
    if sif_doy > int(GLASS_DOYS[-1]) + 7:
        return int(GLASS_DOYS[-1])
    return int(candidates[-1])


def stratified_sample_sif(df: pd.DataFrame) -> pd.DataFrame:
    if SAMPLE_MODE == "none":
        return df.copy()

    if SAMPLE_MODE != "stratified":
        raise ValueError(f"Unknown SAMPLE_MODE: {SAMPLE_MODE}")

    sampled_parts = []
    for _, group in df.groupby(["state", "sif_year", "month"], sort=False, observed=True):
        n = min(len(group), SAMPLES_PER_STATE_YEAR_MONTH)
        sampled_parts.append(group.sample(n=n, random_state=RANDOM_SEED))

    if not sampled_parts:
        return df.head(0).copy()

    sampled = pd.concat(sampled_parts, axis=0)
    return sampled.reset_index(drop=True)


def load_sif_footprints() -> gpd.GeoDataFrame:
    df = pd.read_csv(SIF_CSV_PATH)
    df = df.copy()

    df["Delta_Date"] = pd.to_datetime(df["Delta_Date"], errors="coerce")
    df["sif_year"] = df["Delta_Date"].dt.year
    df["sif_doy"] = df["Delta_Date"].dt.dayofyear
    df["month"] = df["Delta_Date"].dt.month

    df = df[
        df["sif_year"].isin(YEARS)
        & df["month"].between(2, 7)
        & df["Daily_SIF_740nm"].notna()
        & df["Delta_Date"].notna()
    ].copy()

    if "sif_row_id" not in df.columns:
        df.insert(0, "sif_row_id", np.arange(1, len(df) + 1, dtype=np.int64))

    df = stratified_sample_sif(df)

    if MAX_ROWS is not None:
        df = df.head(MAX_ROWS).copy()

    df["sif_year"] = df["sif_year"].astype(int)
    df["sif_doy"] = df["sif_doy"].astype(int)
    df["month"] = df["month"].astype(int)
    df["Delta_Date"] = df["Delta_Date"].dt.date

    df["composite_doy"] = df["sif_doy"].map(containing_glass_doy).astype(np.int16)
    df["composite_start_date"] = pd.to_datetime(
        df["sif_year"].astype(str) + df["composite_doy"].astype(str).str.zfill(3),
        format="%Y%j",
    ).dt.date
    df["composite_end_date"] = (
        pd.to_datetime(df["composite_start_date"]) + pd.Timedelta(days=7)
    ).dt.date
    df["composite_day_offset"] = df["sif_doy"] - df["composite_doy"]

    df = df.sort_values(["sif_year", "composite_doy", "state", "month", "sif_row_id"])

    geometry = df.apply(build_sif_polygon, axis=1)
    gdf = gpd.GeoDataFrame(
        df.drop(columns=["geometry"], errors="ignore"),
        geometry=geometry,
        crs="EPSG:4326",
    )
    gdf = gdf.reset_index(drop=True)

    print(f"Loaded {len(df):,} SIF rows after filtering/sampling")
    return gdf


# ---------------------------------------------------------------------------
# Chip generation


def make_chip_transform(geom_modis) -> Affine:
    centroid = geom_modis.centroid
    chip_width_m = CHIP_SIZE * CHIP_RES_M
    west = centroid.x - chip_width_m / 2
    north = centroid.y + chip_width_m / 2
    return from_origin(west, north, CHIP_RES_M, CHIP_RES_M)


def reproject_to_chip(
    source: RasterSource,
    dst_transform: Affine,
    dst_crs: str = MODIS_SINUSOIDAL_CRS,
    resampling: Resampling = Resampling.bilinear,
) -> np.ndarray:
    dst = np.full((CHIP_SIZE, CHIP_SIZE), np.nan, dtype=np.float32)
    with rasterio.open(source.path) as src:
        reproject(
            source=rasterio.band(src, 1),
            destination=dst,
            src_transform=src.transform,
            src_crs=src.crs,
            src_nodata=src.nodata,
            dst_transform=dst_transform,
            dst_crs=dst_crs,
            dst_nodata=np.nan,
            resampling=resampling,
        )
    return dst


def mosaic_modis_product_chip(product: str, year: int, doy: int, dst_transform: Affine) -> np.ndarray:
    bounds = chip_bounds(dst_transform)
    chip_arrays = []

    for tile in MODIS_TILES:
        src = load_modis_tile(product, year, doy, tile)
        if not intersects(bounds, src.bounds):
            continue
        chip_arrays.append(clean_product_array(reproject_to_chip(src, dst_transform), product))

    if not chip_arrays:
        return np.full((CHIP_SIZE, CHIP_SIZE), np.nan, dtype=np.float32)

    stacked = np.stack(chip_arrays, axis=0)
    with np.errstate(invalid="ignore"):
        return np.nanmean(stacked, axis=0).astype(np.float32)


def par_chip(year: int, doy: int, dst_transform: Affine) -> np.ndarray:
    return clean_par_array(reproject_to_chip(load_par(year, doy), dst_transform))


def fractional_footprint_mask(geom_modis, dst_transform: Affine) -> np.ndarray:
    hi_size = CHIP_SIZE * MASK_OVERSAMPLE
    hi_res = CHIP_RES_M / MASK_OVERSAMPLE
    hi_transform = from_origin(dst_transform.c, dst_transform.f, hi_res, hi_res)

    hi_mask = rasterize(
        [(mapping(geom_modis), 1.0)],
        out_shape=(hi_size, hi_size),
        transform=hi_transform,
        fill=0.0,
        dtype="float32",
        all_touched=False,
    )

    return hi_mask.reshape(
        CHIP_SIZE,
        MASK_OVERSAMPLE,
        CHIP_SIZE,
        MASK_OVERSAMPLE,
    ).mean(axis=(1, 3)).astype(np.float32)


def _window_from_bounds_safe(bounds: tuple[float, float, float, float], transform: Affine) -> Window:
    window = from_bounds(*bounds, transform=transform)
    return window.round_offsets().round_lengths()


@lru_cache(maxsize=16)
def crop_path_for_year(year: int) -> Path:
    path = CROP_DIR / f"croptypes_{year}.tif"
    if not path.exists():
        raise FileNotFoundError(path)
    return path


def crop_group_fraction_chip(
    crop_values: np.ndarray,
    crop_transform: Affine,
    crop_crs,
    group_codes: Iterable[int],
    dst_transform: Affine,
) -> np.ndarray:
    group_mask = np.isin(crop_values, list(group_codes)).astype(np.float32)
    dst = np.zeros((CHIP_SIZE, CHIP_SIZE), dtype=np.float32)

    reproject(
        source=group_mask,
        destination=dst,
        src_transform=crop_transform,
        src_crs=crop_crs,
        dst_transform=dst_transform,
        dst_crs=MODIS_SINUSOIDAL_CRS,
        dst_nodata=0.0,
        resampling=Resampling.average,
    )

    return np.clip(dst, 0.0, 1.0).astype(np.float32)


def crop_fractions_chip(year: int, dst_transform: Affine) -> list[np.ndarray]:
    chip_poly_modis = box(*chip_bounds(dst_transform))

    with rasterio.open(crop_path_for_year(year)) as src:
        chip_poly_crop = transform_geom(
            MODIS_SINUSOIDAL_CRS,
            src.crs,
            mapping(chip_poly_modis),
            precision=6,
        )
        crop_bounds = tuple(shape(chip_poly_crop).bounds)

        window = _window_from_bounds_safe(crop_bounds, src.transform)
        crop_values = src.read(1, window=window, boundless=True, fill_value=0)
        crop_transform = src.window_transform(window)

        crop_channels = [
            crop_group_fraction_chip(
                crop_values,
                crop_transform,
                src.crs,
                codes,
                dst_transform,
            )
            for codes in CROP_GROUPS.values()
        ]

        crop_code_union = sorted({code for codes in CROP_GROUPS.values() for code in codes})
        non_crop_values = (~np.isin(crop_values, crop_code_union)).astype(np.float32)
        non_crop = np.zeros((CHIP_SIZE, CHIP_SIZE), dtype=np.float32)

        reproject(
            source=non_crop_values,
            destination=non_crop,
            src_transform=crop_transform,
            src_crs=src.crs,
            dst_transform=dst_transform,
            dst_crs=MODIS_SINUSOIDAL_CRS,
            dst_nodata=1.0,
            resampling=Resampling.average,
        )

    crop_channels.append(np.clip(non_crop, 0.0, 1.0).astype(np.float32))
    return crop_channels


def constant_channel(value: float) -> np.ndarray:
    return np.full((CHIP_SIZE, CHIP_SIZE), value, dtype=np.float32)


def build_chip(row: pd.Series, geom_modis) -> tuple[np.ndarray, np.ndarray, dict]:
    year = int(row["sif_year"])
    doy = int(row["composite_doy"])
    month = int(row["month"])

    dst_transform = make_chip_transform(geom_modis)

    fapar = mosaic_modis_product_chip("fapar", year, doy, dst_transform)
    evi = mosaic_modis_product_chip("evi", year, doy, dst_transform)
    ndvi = mosaic_modis_product_chip("ndvi", year, doy, dst_transform)
    par = par_chip(year, doy, dst_transform)
    apar = fapar * par

    crop_channels = crop_fractions_chip(year, dst_transform)
    mask = fractional_footprint_mask(geom_modis, dst_transform)

    month_angle = 2.0 * np.pi * month / 12.0
    month_sin = constant_channel(np.sin(month_angle))
    month_cos = constant_channel(np.cos(month_angle))
    latitude = constant_channel(float(row["Latitude"]))
    longitude = constant_channel(float(row["Longitude"]))

    channels = [
        fapar,
        evi,
        ndvi,
        par,
        apar,
        *crop_channels,
        month_sin,
        month_cos,
        latitude,
        longitude,
    ]

    x = np.stack(channels, axis=0).astype(np.float32)

    valid_fraction = {
        f"{name}_valid_fraction": float(np.isfinite(channel).mean())
        for name, channel in zip(CHANNEL_NAMES, channels)
        if name not in {"month_sin", "month_cos", "latitude", "longitude"}
    }

    # Most neural-network training code does not accept NaN tensors. Keep the
    # valid fractions in metadata so you can filter samples later if needed.
    x = np.nan_to_num(x, nan=0.0, posinf=0.0, neginf=0.0)

    metadata = {
        "sif_row_id": row["sif_row_id"],
        "Delta_Date": row["Delta_Date"],
        "state": row["state"],
        "sif_year": year,
        "sif_doy": int(row["sif_doy"]),
        "composite_doy": doy,
        "composite_start_date": row["composite_start_date"],
        "composite_end_date": row["composite_end_date"],
        "composite_day_offset": int(row["composite_day_offset"]),
        "Daily_SIF_740nm": float(row["Daily_SIF_740nm"]),
        "Latitude": float(row["Latitude"]),
        "Longitude": float(row["Longitude"]),
        "mask_sum": float(mask.sum()),
        **valid_fraction,
    }

    return x, mask, metadata


# ---------------------------------------------------------------------------
# Optional multiprocessing helpers


def make_payload(idx: int, row: pd.Series, geom_modis) -> tuple[int, dict, bytes]:
    row_dict = row.to_dict()
    row_dict.pop("geometry", None)
    return idx, row_dict, geom_modis.wkb


def build_chip_from_payload(payload: tuple[int, dict, bytes]) -> tuple[int, np.ndarray, np.ndarray, dict]:
    idx, row_dict, geom_wkb = payload
    row = pd.Series(row_dict)
    geom_modis = wkb.loads(geom_wkb)
    x, mask, metadata = build_chip(row, geom_modis)
    return idx, x, mask, metadata


# ---------------------------------------------------------------------------
# Writing shards


def write_shard(
    shard_id: int,
    xs: list[np.ndarray],
    masks: list[np.ndarray],
    ys: list[float],
    row_ids: list[int],
) -> Path:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    path = OUTPUT_DIR / f"chips_{shard_id:05d}.npz"
    np.savez_compressed(
        path,
        X=np.stack(xs, axis=0).astype(np.float32),
        footprint_mask=np.stack(masks, axis=0).astype(np.float32),
        y=np.asarray(ys, dtype=np.float32),
        sif_row_id=np.asarray(row_ids),
        channel_names=np.asarray(CHANNEL_NAMES),
    )
    return path


def prepare_chips() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    sif_gdf = load_sif_footprints()
    sif_modis = sif_gdf.to_crs(MODIS_SINUSOIDAL_CRS)

    metadata_rows = []
    xs: list[np.ndarray] = []
    masks: list[np.ndarray] = []
    ys: list[float] = []
    row_ids: list[int] = []
    shard_id = 0

    def handle_result(result: tuple[int, np.ndarray, np.ndarray, dict]) -> None:
        nonlocal shard_id, xs, masks, ys, row_ids
        idx, x, mask, metadata = result
        if idx % 100 == 0:
            print(f"Prepared chip {idx + 1:,} / {len(sif_gdf):,}")
        if metadata["mask_sum"] <= 0:
            print(f"Skipping sif_row_id={metadata['sif_row_id']} because footprint mask is empty")
            return

        xs.append(x)
        masks.append(mask)
        ys.append(metadata["Daily_SIF_740nm"])
        row_ids.append(metadata["sif_row_id"])
        metadata_rows.append(metadata)

        if len(xs) == SHARD_SIZE:
            written = write_shard(shard_id, xs, masks, ys, row_ids)
            print(f"Wrote {written}")
            shard_id += 1
            xs, masks, ys, row_ids = [], [], [], []

    payloads = (
        make_payload(idx, row, sif_modis.geometry.iloc[idx])
        for idx, row in sif_gdf.iterrows()
    )

    if N_WORKERS == 1:
        for payload in payloads:
            handle_result(build_chip_from_payload(payload))
    else:
        with ProcessPoolExecutor(max_workers=N_WORKERS) as executor:
            for result in executor.map(build_chip_from_payload, payloads, chunksize=8):
                handle_result(result)

    if xs:
        written = write_shard(shard_id, xs, masks, ys, row_ids)
        print(f"Wrote {written}")

    metadata = pd.DataFrame(metadata_rows)
    metadata.to_csv(OUTPUT_DIR / "chip_metadata.csv", index=False)

    channel_table = pd.DataFrame(
        {"channel_index": np.arange(len(CHANNEL_NAMES)), "channel_name": CHANNEL_NAMES}
    )
    channel_table.to_csv(OUTPUT_DIR / "channel_names.csv", index=False)

    print(f"Done. Wrote {len(metadata):,} chips to {OUTPUT_DIR}")


if __name__ == "__main__":
    prepare_chips()
