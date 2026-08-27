"""Extract annual wheat composition and monthly crop-pure NIRv by NUTS3.

This is a lightweight companion to
``predict_monthly_winter_wheat_sif_nuts3.py``. It does not load the CNN or
repeat SIF inference. It uses the same NUTS3 regions, assigned MGRS tiles,
20 m grid, Sentinel validity rules and strict winter-wheat definition.

Outputs:
  - nuts3_wheat_composition_nirv_long.csv
  - nuts3_wheat_composition_nirv_wide.csv
  - nuts3_monthly_strict_wheat_sif_nirv_wide.csv

``ww_pct`` is stored as a fraction from 0 to 1. It is the winter-wheat area
divided by all known-crop area inside the NUTS3 region's assigned tile.
Monthly NIRv is averaged over strict 20 m winter-wheat pixels, where all four
underlying 10 m crop pixels are winter wheat.
"""

from __future__ import annotations

import calendar
import math
from contextlib import ExitStack
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from rasterio.enums import Resampling
from rasterio.transform import Affine, array_bounds
from rasterio.warp import reproject, transform_bounds
from rasterio.windows import Window, from_bounds

import predict_monthly_winter_wheat_sif_nuts3 as wheat_sif
import prepare_sentinel2_multisif_cnn_chips as sentinel


# ---------------------------------------------------------------------------
# Configuration

OUTPUT_DIR = wheat_sif.OUTPUT_DIR
SIF_WIDE_PATH = wheat_sif.WIDE_OUTPUT_PATH

LONG_OUTPUT_PATH = (
    OUTPUT_DIR / "nuts3_wheat_composition_nirv_long.csv"
)
WIDE_OUTPUT_PATH = (
    OUTPUT_DIR / "nuts3_wheat_composition_nirv_wide.csv"
)
MERGED_OUTPUT_PATH = (
    OUTPUT_DIR / "nuts3_monthly_strict_wheat_sif_nirv_wide.csv"
)

MGRS_REFERENCE_DIR = Path("data/temp_data/mgrs_tifs")

TARGET_TILES = wheat_sif.TARGET_TILES
YEARS = wheat_sif.YEARS
MONTHS = wheat_sif.MONTHS

MODEL_RES_M = wheat_sif.MODEL_RES_M
BLOCK_SIZE = 1024
TEN_METRE_PIXELS_PER_20M_PIXEL = 4.0

WINTER_WHEAT_CODE = wheat_sif.WINTER_WHEAT_CODE
STRICT_WHEAT_THRESHOLD = wheat_sif.STRICT_WHEAT_THRESHOLD


# ---------------------------------------------------------------------------
# Small grid helpers

@dataclass(frozen=True)
class RasterBlock:
    row_start: int
    row_stop: int
    col_start: int
    col_stop: int
    transform: Affine

    @property
    def height(self) -> int:
        return self.row_stop - self.row_start

    @property
    def width(self) -> int:
        return self.col_stop - self.col_start

    @property
    def rows(self) -> slice:
        return slice(self.row_start, self.row_stop)

    @property
    def cols(self) -> slice:
        return slice(self.col_start, self.col_stop)


def regional_grid(
    reference: rasterio.io.DatasetReader,
    regions,
) -> wheat_sif.GridSpec:
    """Create a 20 m grid around the regions without CNN context padding."""
    projected = regions.to_crs(reference.crs)
    xmin, ymin, xmax, ymax = projected.total_bounds
    region_window = from_bounds(
        xmin,
        ymin,
        xmax,
        ymax,
        transform=reference.transform,
    )

    col0 = max(0, math.floor(region_window.col_off))
    row0 = max(0, math.floor(region_window.row_off))
    col1 = min(
        reference.width,
        math.ceil(region_window.col_off + region_window.width),
    )
    row1 = min(
        reference.height,
        math.ceil(region_window.row_off + region_window.height),
    )
    if col1 <= col0 or row1 <= row0:
        raise ValueError("NUTS3 regions do not overlap the reference raster")

    window = Window(col0, row0, col1 - col0, row1 - row0)
    return wheat_sif.GridSpec(
        row_off=int(row0),
        col_off=int(col0),
        height=int(row1 - row0),
        width=int(col1 - col0),
        transform=reference.window_transform(window),
        crs=reference.crs,
    )


def iter_region_blocks(
    grid: wheat_sif.GridSpec,
    region_ids: np.ndarray,
) -> list[RasterBlock]:
    blocks = []
    for row_start in range(0, grid.height, BLOCK_SIZE):
        row_stop = min(row_start + BLOCK_SIZE, grid.height)
        for col_start in range(0, grid.width, BLOCK_SIZE):
            col_stop = min(col_start + BLOCK_SIZE, grid.width)
            ids = region_ids[row_start:row_stop, col_start:col_stop]
            if not np.any(ids > 0):
                continue
            blocks.append(
                RasterBlock(
                    row_start=row_start,
                    row_stop=row_stop,
                    col_start=col_start,
                    col_stop=col_stop,
                    transform=(
                        grid.transform
                        * Affine.translation(col_start, row_start)
                    ),
                )
            )
    return blocks


def source_window_for_block(
    source: rasterio.io.DatasetReader,
    block: RasterBlock,
    destination_crs,
    pad_pixels: int,
) -> Window:
    bounds = array_bounds(block.height, block.width, block.transform)
    source_bounds = transform_bounds(
        destination_crs,
        source.crs,
        *bounds,
        densify_pts=21,
    )
    return wheat_sif.expanded_window_for_bounds(
        source_bounds,
        source.transform,
        pad_pixels=pad_pixels,
    )


def reference_b5_path(tile: str, year: int) -> Path:
    for month in MONTHS:
        product_path = wheat_sif.sentinel_product_path(tile, year, month)
        if product_path is not None:
            return sentinel.sentinel_band_path(product_path, "B5")

    matches = sorted(
        path
        for path in MGRS_REFERENCE_DIR.rglob("*.tif")
        if f"T{tile}" in path.name and "_B5" in path.stem
    )
    if len(matches) != 1:
        raise FileNotFoundError(
            f"Expected one B5 reference for {tile}, {year}; "
            f"found {len(matches)}"
        )
    return matches[0]


# ---------------------------------------------------------------------------
# Crop composition

class CropReader:
    def __init__(self, path: Path):
        self.path = path
        self.source: rasterio.io.DatasetReader | None = None

    def __enter__(self) -> "CropReader":
        self.source = rasterio.open(self.path)
        if self.source.crs is None:
            raise ValueError(f"Crop raster has no CRS: {self.path}")
        return self

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        if self.source is not None:
            self.source.close()

    def fractions(
        self,
        block: RasterBlock,
        destination_crs,
    ) -> tuple[np.ndarray, np.ndarray]:
        if self.source is None:
            raise RuntimeError("Crop reader is not open")

        window = source_window_for_block(
            self.source,
            block,
            destination_crs,
            pad_pixels=2,
        )
        crop_values = self.source.read(
            1,
            window=window,
            boundless=True,
            fill_value=sentinel.NON_CROP_CODE,
        )
        if self.source.nodata is not None:
            crop_values[crop_values == self.source.nodata] = (
                sentinel.NON_CROP_CODE
            )
        source_transform = self.source.window_transform(window)

        crop_source = np.isin(
            crop_values,
            sentinel.KNOWN_CROP_CODES,
        ).astype(np.float32)
        wheat_source = (
            crop_values == WINTER_WHEAT_CODE
        ).astype(np.float32)

        crop_fraction = np.zeros(
            (block.height, block.width),
            dtype=np.float32,
        )
        wheat_fraction = np.zeros_like(crop_fraction)
        for source_values, destination in (
            (crop_source, crop_fraction),
            (wheat_source, wheat_fraction),
        ):
            reproject(
                source=source_values,
                destination=destination,
                src_transform=source_transform,
                src_crs=self.source.crs,
                dst_transform=block.transform,
                dst_crs=destination_crs,
                dst_nodata=0.0,
                resampling=Resampling.average,
            )

        return (
            np.clip(crop_fraction, 0.0, 1.0),
            np.clip(wheat_fraction, 0.0, 1.0),
        )


def annual_crop_composition(
    crop_raster_path: Path,
    grid: wheat_sif.GridSpec,
    region_ids: np.ndarray,
    region_lookup: dict[int, dict],
    blocks: list[RasterBlock],
    year: int,
) -> tuple[pd.DataFrame, np.ndarray]:
    n_codes = max(region_lookup) + 1
    crop_20m_area_sum = np.zeros(n_codes, dtype=np.float64)
    wheat_20m_area_sum = np.zeros(n_codes, dtype=np.float64)
    strict_wheat_count = np.zeros(n_codes, dtype=np.int64)
    strict_wheat = np.zeros(region_ids.shape, dtype=bool)

    with CropReader(crop_raster_path) as reader:
        for block in blocks:
            ids = region_ids[block.rows, block.cols]
            crop_fraction, wheat_fraction = reader.fractions(
                block,
                grid.crs,
            )
            valid_region = ids > 0

            crop_20m_area_sum += np.bincount(
                ids[valid_region],
                weights=crop_fraction[valid_region],
                minlength=n_codes,
            )[:n_codes]
            wheat_20m_area_sum += np.bincount(
                ids[valid_region],
                weights=wheat_fraction[valid_region],
                minlength=n_codes,
            )[:n_codes]

            strict_block = (
                wheat_fraction >= STRICT_WHEAT_THRESHOLD
            ) & valid_region
            strict_wheat[block.rows, block.cols] = strict_block
            strict_wheat_count += np.bincount(
                ids[strict_block],
                minlength=n_codes,
            )[:n_codes]

    rows = []
    for code, metadata in region_lookup.items():
        crop_equivalent_10m = (
            crop_20m_area_sum[code]
            * TEN_METRE_PIXELS_PER_20M_PIXEL
        )
        wheat_equivalent_10m = (
            wheat_20m_area_sum[code]
            * TEN_METRE_PIXELS_PER_20M_PIXEL
        )
        ww_pct = (
            wheat_equivalent_10m / crop_equivalent_10m
            if crop_equivalent_10m > 0
            else np.nan
        )
        rows.append(
            {
                **metadata,
                "year": year,
                "winter_wheat_pixel_equivalents_10m": (
                    wheat_equivalent_10m
                ),
                "total_crop_pixel_equivalents_10m": (
                    crop_equivalent_10m
                ),
                "winter_wheat_area_km2": (
                    wheat_equivalent_10m * 100.0 / 1_000_000.0
                ),
                "total_crop_area_km2": (
                    crop_equivalent_10m * 100.0 / 1_000_000.0
                ),
                "ww_pct": ww_pct,
                "strict_wheat_pixels_20m": int(
                    strict_wheat_count[code]
                ),
            }
        )

    return pd.DataFrame(rows), strict_wheat


# ---------------------------------------------------------------------------
# Monthly NIRv

class NirvReader:
    def __init__(self, product_path: Path):
        self.product_path = product_path
        self.stack = ExitStack()
        self.b4 = None
        self.b8 = None
        self.flag = None

    def __enter__(self) -> "NirvReader":
        b4_path = sentinel.sentinel_band_path(self.product_path, "B4")
        b8_path = sentinel.sentinel_band_path(self.product_path, "B8")
        flag_path = sentinel.sentinel_flag_path(
            self.product_path,
            "R1",
        )
        self.b4 = self.stack.enter_context(rasterio.open(b4_path))
        self.b8 = self.stack.enter_context(rasterio.open(b8_path))
        self.flag = self.stack.enter_context(rasterio.open(flag_path))

        sources = (self.b4, self.b8, self.flag)
        if any(source.crs is None for source in sources):
            raise ValueError(
                f"Sentinel raster has no CRS: {self.product_path}"
            )
        first = self.b4
        for source in sources[1:]:
            if (
                source.crs != first.crs
                or source.transform != first.transform
                or source.width != first.width
                or source.height != first.height
            ):
                raise ValueError(
                    "B4, B8 and FLG_R1 grids do not match in "
                    f"{self.product_path}"
                )
        return self

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        self.stack.close()

    def read_nirv(
        self,
        block: RasterBlock,
        destination_crs,
    ) -> np.ndarray:
        if self.b4 is None or self.b8 is None or self.flag is None:
            raise RuntimeError("NIRv reader is not open")

        window = source_window_for_block(
            self.b4,
            block,
            destination_crs,
            pad_pixels=1,
        )
        b4 = self.b4.read(
            1,
            window=window,
            boundless=True,
            fill_value=sentinel.REFLECTANCE_NODATA,
        )
        b8 = self.b8.read(
            1,
            window=window,
            boundless=True,
            fill_value=sentinel.REFLECTANCE_NODATA,
        )
        flag = self.flag.read(
            1,
            window=window,
            boundless=True,
            fill_value=0,
        )
        source_transform = self.b4.window_transform(window)

        destinations = []
        for band in (b4, b8):
            valid = (
                (flag == sentinel.LAND_FLAG_VALUE)
                & (band != sentinel.REFLECTANCE_NODATA)
                & np.isfinite(band)
            )
            source = np.full(
                band.shape,
                sentinel.WARP_NODATA,
                dtype=np.float32,
            )
            source[valid] = (
                band[valid].astype(np.float32)
                / sentinel.REFLECTANCE_QUANTIFICATION_VALUE
            )
            destination = np.full(
                (block.height, block.width),
                np.nan,
                dtype=np.float32,
            )
            reproject(
                source=source,
                destination=destination,
                src_transform=source_transform,
                src_crs=self.b4.crs,
                src_nodata=sentinel.WARP_NODATA,
                dst_transform=block.transform,
                dst_crs=destination_crs,
                dst_nodata=np.nan,
                resampling=Resampling.average,
            )
            destinations.append(destination)

        b4_20m, b8_20m = destinations
        ndvi = sentinel.normalized_difference(b8_20m, b4_20m)
        return (b8_20m * ndvi).astype(np.float32)


def monthly_nirv_rows(
    product_path: Path,
    grid: wheat_sif.GridSpec,
    region_ids: np.ndarray,
    region_lookup: dict[int, dict],
    strict_wheat: np.ndarray,
    blocks: list[RasterBlock],
    year: int,
    month: int,
) -> list[dict]:
    n_codes = max(region_lookup) + 1
    value_sum = np.zeros(n_codes, dtype=np.float64)
    value_count = np.zeros(n_codes, dtype=np.int64)
    strict_count = np.bincount(
        region_ids[strict_wheat & (region_ids > 0)],
        minlength=n_codes,
    )[:n_codes]

    with NirvReader(product_path) as reader:
        for block in blocks:
            ids = region_ids[block.rows, block.cols]
            wheat = strict_wheat[block.rows, block.cols]
            target = wheat & (ids > 0)
            if not np.any(target):
                continue

            nirv = reader.read_nirv(block, grid.crs)
            valid = target & np.isfinite(nirv)
            if not np.any(valid):
                continue
            value_sum += np.bincount(
                ids[valid],
                weights=nirv[valid],
                minlength=n_codes,
            )[:n_codes]
            value_count += np.bincount(
                ids[valid],
                minlength=n_codes,
            )[:n_codes]

    rows = []
    for code, metadata in region_lookup.items():
        count = int(value_count[code])
        total = int(strict_count[code])
        rows.append(
            {
                **metadata,
                "year": year,
                "month": month,
                "month_name": calendar.month_name[month],
                "sentinel_product": str(product_path),
                "sentinel_product_date": (
                    wheat_sif.product_date(product_path).date()
                ),
                "mean_nirv_strict_wheat": (
                    float(value_sum[code] / count)
                    if count > 0
                    else np.nan
                ),
                "nirv_valid_strict_wheat_pixels_20m": count,
                "strict_wheat_pixels_20m": total,
                "nirv_coverage_pct": (
                    100.0 * count / total
                    if total > 0
                    else np.nan
                ),
                "scene_status": (
                    "success"
                    if count > 0
                    else "no_valid_strict_wheat_nirv"
                ),
                "error": "",
            }
        )
    return rows


def missing_nirv_rows(
    region_lookup: dict[int, dict],
    year: int,
    month: int,
    status: str,
    error: str = "",
) -> list[dict]:
    return [
        {
            **metadata,
            "year": year,
            "month": month,
            "month_name": calendar.month_name[month],
            "sentinel_product": "",
            "sentinel_product_date": pd.NaT,
            "mean_nirv_strict_wheat": np.nan,
            "nirv_valid_strict_wheat_pixels_20m": 0,
            "strict_wheat_pixels_20m": np.nan,
            "nirv_coverage_pct": np.nan,
            "scene_status": status,
            "error": error,
        }
        for metadata in region_lookup.values()
    ]


# ---------------------------------------------------------------------------
# Output tables

def make_wide_table(
    composition: pd.DataFrame,
    nirv_long: pd.DataFrame,
) -> pd.DataFrame:
    keys = [
        "nuts_id",
        "nuts3",
        "mgrs_tile",
        "year",
        "nuts3_area_in_mgrs_pct",
    ]
    composition_columns = [
        "winter_wheat_pixel_equivalents_10m",
        "total_crop_pixel_equivalents_10m",
        "winter_wheat_area_km2",
        "total_crop_area_km2",
        "ww_pct",
        "strict_wheat_pixels_20m",
    ]
    base = composition[keys + composition_columns].copy()

    successful = nirv_long.loc[
        nirv_long["scene_status"].eq("success")
    ].copy()
    nirv = successful.pivot_table(
        index=keys,
        columns="month_name",
        values="mean_nirv_strict_wheat",
        aggfunc="first",
    ).rename(columns=lambda value: f"mean_nirv_{value}")
    valid_pixels = successful.pivot_table(
        index=keys,
        columns="month_name",
        values="nirv_valid_strict_wheat_pixels_20m",
        aggfunc="first",
    ).rename(columns=lambda value: f"nirv_valid_pixels_{value}")

    wide = (
        base.set_index(keys)
        .join(nirv, how="left")
        .join(valid_pixels, how="left")
        .reset_index()
    )
    month_names = [calendar.month_name[month] for month in MONTHS]
    for month_name in month_names:
        for column in (
            f"mean_nirv_{month_name}",
            f"nirv_valid_pixels_{month_name}",
        ):
            if column not in wide.columns:
                wide[column] = np.nan

    ordered = [
        *keys,
        *composition_columns,
        *(f"mean_nirv_{month}" for month in month_names),
        *(f"nirv_valid_pixels_{month}" for month in month_names),
    ]
    return wide[ordered]


def merge_with_sif_wide(features_wide: pd.DataFrame) -> pd.DataFrame:
    if not SIF_WIDE_PATH.exists():
        raise FileNotFoundError(
            f"Completed SIF wide table not found: {SIF_WIDE_PATH}"
        )
    sif_wide = pd.read_csv(SIF_WIDE_PATH, low_memory=False)
    keys = ["nuts_id", "mgrs_tile", "year"]
    wheat_sif.require_columns(sif_wide, keys, "SIF wide table")
    wheat_sif.require_columns(features_wide, keys, "NIRv feature table")

    sif_wide["mgrs_tile"] = sif_wide["mgrs_tile"].map(
        wheat_sif.normalize_mgrs_tile
    )
    features_wide = features_wide.copy()
    features_wide["mgrs_tile"] = features_wide["mgrs_tile"].map(
        wheat_sif.normalize_mgrs_tile
    )

    new_columns = [
        "winter_wheat_pixel_equivalents_10m",
        "total_crop_pixel_equivalents_10m",
        "winter_wheat_area_km2",
        "total_crop_area_km2",
        "ww_pct",
        "strict_wheat_pixels_20m",
        *(f"mean_nirv_{calendar.month_name[m]}" for m in MONTHS),
        *(f"nirv_valid_pixels_{calendar.month_name[m]}" for m in MONTHS),
    ]
    return sif_wide.merge(
        features_wide[keys + new_columns],
        on=keys,
        how="left",
        validate="one_to_one",
    )


# ---------------------------------------------------------------------------
# Main workflow

def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    regions = wheat_sif.load_regions()

    composition_parts = []
    nirv_rows = []
    total_year_groups = len(TARGET_TILES) * len(YEARS)
    group_index = 0

    for tile in TARGET_TILES:
        tile_regions = regions.loc[
            regions["mgrs_tile"].eq(tile)
        ].copy()
        if tile_regions.empty:
            continue

        for year in YEARS:
            group_index += 1
            print(
                f"[{group_index}/{total_year_groups}] "
                f"{tile} {year}: crop composition"
            )
            reference_path = reference_b5_path(tile, year)
            with rasterio.open(reference_path) as reference:
                if reference.crs is None:
                    raise ValueError(
                        f"Reference raster has no CRS: {reference_path}"
                    )
                if not (
                    math.isclose(
                        abs(reference.transform.a),
                        MODEL_RES_M,
                        abs_tol=1e-6,
                    )
                    and math.isclose(
                        abs(reference.transform.e),
                        MODEL_RES_M,
                        abs_tol=1e-6,
                    )
                ):
                    raise ValueError(
                        f"Expected a {MODEL_RES_M:g} m reference: "
                        f"{reference_path}"
                    )

                grid = regional_grid(reference, tile_regions)
                region_ids, region_lookup = wheat_sif.rasterize_regions(
                    tile_regions,
                    grid,
                )
                blocks = iter_region_blocks(grid, region_ids)

            composition, strict_wheat = annual_crop_composition(
                wheat_sif.crop_path(year),
                grid,
                region_ids,
                region_lookup,
                blocks,
                year,
            )
            composition_parts.append(composition)
            print(
                f"  regional grid: {grid.height:,} x {grid.width:,}; "
                f"{len(blocks)} occupied blocks"
            )

            for month in MONTHS:
                product_path = wheat_sif.sentinel_product_path(
                    tile,
                    year,
                    month,
                )
                print(
                    f"  {calendar.month_name[month]}: ",
                    end="",
                )
                if product_path is None:
                    print("missing Sentinel product")
                    nirv_rows.extend(
                        missing_nirv_rows(
                            region_lookup,
                            year,
                            month,
                            "missing_sentinel_product",
                        )
                    )
                    continue

                try:
                    rows = monthly_nirv_rows(
                        product_path,
                        grid,
                        region_ids,
                        region_lookup,
                        strict_wheat,
                        blocks,
                        year,
                        month,
                    )
                except Exception as error:
                    nirv_rows.extend(
                        missing_nirv_rows(
                            region_lookup,
                            year,
                            month,
                            "failed",
                            f"{type(error).__name__}: {error}",
                        )
                    )
                    raise
                nirv_rows.extend(rows)
                print(
                    f"{sum(row['scene_status'] == 'success' for row in rows)} "
                    f"/ {len(rows)} regions successful"
                )

    composition = pd.concat(composition_parts, ignore_index=True)
    nirv_long = pd.DataFrame(nirv_rows)
    long_output = nirv_long.merge(
        composition,
        on=[
            "nuts_id",
            "nuts3",
            "mgrs_tile",
            "year",
            "nuts3_area_m2",
            "overlap_area_m2",
            "nuts3_area_in_mgrs_pct",
        ],
        how="left",
        validate="many_to_one",
        suffixes=("", "_annual"),
    )
    wide_output = make_wide_table(composition, nirv_long)
    merged_output = merge_with_sif_wide(wide_output)

    long_output.to_csv(LONG_OUTPUT_PATH, index=False)
    wide_output.to_csv(WIDE_OUTPUT_PATH, index=False)
    merged_output.to_csv(MERGED_OUTPUT_PATH, index=False)

    print(f"Saved long table: {LONG_OUTPUT_PATH}")
    print(f"Saved companion wide table: {WIDE_OUTPUT_PATH}")
    print(f"Saved merged SIF + NIRv table: {MERGED_OUTPUT_PATH}")
    print(
        "Merged rows with ww_pct: "
        f"{merged_output['ww_pct'].notna().sum():,} / "
        f"{len(merged_output):,}"
    )


if __name__ == "__main__":
    main()
