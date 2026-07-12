"""
Prepare multi-footprint Sentinel-2 CNN chips on a common 20 m grid.

Inputs are produced by diagnose_sentinel2_multisif_density_cluster_chips.R:
  - one row per 6 km chip in sentinel2_multi_sif_6km_chip_manifest.csv
  - one row per assigned SIF footprint in
    sentinel2_multi_sif_6km_chip_assignments.csv

Each output sample contains:
  X:                 [19 channels, 300, 300], stored as float16
  footprint_masks:   [max_footprints, 300, 300], stored as float16
  y_targets:         [max_footprints], stored as float32
  target_accept:     [max_footprints]
  footprint_valid:   [max_footprints]

Predictor NaNs are intentionally preserved. During model training, channel
statistics should ignore NaNs; normalize first, then replace remaining NaNs
with zero so missing pixels map to the normalized channel mean.
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from functools import lru_cache
import glob
import json
import math
from pathlib import Path
import re
from typing import Iterable

import numpy as np
import pandas as pd
import rasterio
from rasterio.enums import Resampling
from rasterio.features import rasterize
from rasterio.transform import Affine, from_origin
from rasterio.warp import reproject, transform_bounds, transform_geom
from rasterio.windows import Window, from_bounds
from shapely.geometry import MultiPoint, Polygon, mapping, shape


# ---------------------------------------------------------------------------
# Config

CHIP_DIAGNOSTIC_DIR = Path("data/sentinel2_multisif_chip_diagnostics")
CHIP_MANIFEST_PATH = CHIP_DIAGNOSTIC_DIR / "sentinel2_multi_sif_6km_chip_manifest.csv"
CHIP_ASSIGNMENTS_PATH = (
    CHIP_DIAGNOSTIC_DIR / "sentinel2_multi_sif_6km_chip_assignments.csv"
)

FAPAR_DIR = Path("data/glass_geotiff/fapar")
PAR_DIR = Path("data/viirs_vnp18a2_daily_mean_par_germany_native")
CROP_DIR = Path("data/crop_type_tif")

OUTPUT_DIR = Path(
    "data/cnn_sentinel2_chips/multisif_6km_20m_indices_fapar_active_crop"
)

TARGET_COLUMN = "target_modis_sif"
FINAL_CHECK_COLUMN = "final_check_modis_sif"

CHIP_SIZE_M = 6000.0
CHIP_RES_M = 20.0
CHIP_SIZE = 300
MAX_FOOTPRINTS = 10
MIN_VALID_FOOTPRINTS = 4
MIN_FOOTPRINT_INSIDE_FRACTION = 0.90

# A 4x supersampled mask gives 5 m subpixels on the 20 m output grid.
MASK_OVERSAMPLE = 4

# Keep shards modest: one 19-channel chip plus ten masks is still large.
SHARD_SIZE = 16

# Start with a smoke test. Set to None only after checking its outputs.
#MAX_CHIPS: int | None = 10
MAX_CHIPS: int | None = None
# ProcessPoolExecutor requires at least 1. Use 1 for sequential processing.
# After the smoke test, 2 may help on fast storage; many workers can cause
# disk contention because the Sentinel-2 GeoTIFFs are large.
N_WORKERS = 1

REFLECTANCE_QUANTIFICATION_VALUE = 10000.0
REFLECTANCE_NODATA = -10000
LAND_FLAG_VALUE = 4
NON_CROP_CODE = 0
WARP_NODATA = -9999.0

PAR_VALID_MIN = 0.0
PAR_VALID_MAX = 700.0
PAR_FILL_VALUE = -1.0
PAR_ACCEPTED_QA_CODES = (1, 2)

MODIS_TILES = ("h18v03", "h18v04")
FAPAR_COMPOSITE_DOYS = tuple(range(33, 210, 8))

# Bands used by the five requested indices. R1 is the 10 m flag grid and R2
# is the 20 m flag grid.
SENTINEL_BANDS = {
    "B2": "R1",
    "B4": "R1",
    "B8": "R1",
    "B5": "R2",
    "B8A": "R2",
    "B11": "R2",
}


# Crop legend codes from data/crop_type_tif/LEGEND_CropTypes.txt
CROP_CODES = {
    "winter_wheat": 11,
    "winter_barley": 12,
    "winter_rye": 13,
    "other_winter_cereals": 14,
    "spring_wheat": 21,
    "spring_barley": 22,
    "spring_oat": 23,
    "maize": 30,
    "legumes": 40,
    "potato": 50,
    "sugar_beet": 60,
    "rapeseed": 71,
    "clover_alfalfa": 81,
    "arable_grass": 82,
    "permanent_grassland": 83,
    "vineyard": 90,
    "fruit_trees_and_other_woody_vegetation": 100,
    "hops": 110,
    "other_agricultural_use": 111,
}

KNOWN_CROP_CODES = sorted(CROP_CODES.values())

CROP_GROUPS = {
    "winter_wheat_fraction": [CROP_CODES["winter_wheat"]],
    "winter_barley_fraction": [CROP_CODES["winter_barley"]],
    "winter_rye_fraction": [CROP_CODES["winter_rye"]],
    "maize_fraction": [CROP_CODES["maize"]],
    "grass_forage_fraction": [
        CROP_CODES["clover_alfalfa"],
        CROP_CODES["arable_grass"],
        CROP_CODES["permanent_grassland"],
    ],
    "woody_fraction": [
        CROP_CODES["vineyard"],
        CROP_CODES["fruit_trees_and_other_woody_vegetation"],
    ],
}

ASSIGNED_CROP_CODES = sorted({code for codes in CROP_GROUPS.values() for code in codes})
CROP_GROUPS["other_crop_fraction"] = [
    code for code in KNOWN_CROP_CODES if code not in ASSIGNED_CROP_CODES
]

ACTIVE_GROWTH_MONTHS_BY_CODE = {
    CROP_CODES["winter_wheat"]: list(range(2, 8)) + [10, 11],
    CROP_CODES["winter_barley"]: list(range(2, 7)) + [10, 11],
    CROP_CODES["winter_rye"]: list(range(2, 8)) + [10, 11],
    CROP_CODES["other_winter_cereals"]: list(range(2, 8)) + [10, 11],
    CROP_CODES["spring_wheat"]: list(range(3, 9)),
    CROP_CODES["spring_barley"]: list(range(3, 9)),
    CROP_CODES["spring_oat"]: list(range(3, 9)),
    CROP_CODES["maize"]: list(range(5, 11)),
    CROP_CODES["legumes"]: list(range(4, 10)),
    CROP_CODES["potato"]: list(range(4, 10)),
    CROP_CODES["sugar_beet"]: list(range(4, 11)),
    CROP_CODES["rapeseed"]: list(range(2, 8)) + [9, 10, 11],
    CROP_CODES["clover_alfalfa"]: list(range(3, 11)),
    CROP_CODES["arable_grass"]: list(range(3, 12)),
    CROP_CODES["permanent_grassland"]: list(range(3, 12)),
    CROP_CODES["vineyard"]: list(range(4, 11)),
    CROP_CODES["fruit_trees_and_other_woody_vegetation"]: list(range(3, 11)),
    CROP_CODES["hops"]: list(range(4, 10)),
    CROP_CODES["other_agricultural_use"]: list(range(3, 11)),
}

CHANNEL_NAMES = [
    "ndmi",
    "ndvi",
    "evi",
    "nirv",
    "ndre",
    "fapar",
    "par",
    "apar",
    *CROP_GROUPS.keys(),
    "active_crop_fraction",
    "non_crop_fraction",
    "month_sin",
    "month_cos",
]


# ---------------------------------------------------------------------------
# Generic helpers


def require_columns(df: pd.DataFrame, columns: Iterable[str], label: str) -> None:
    missing = [column for column in columns if column not in df.columns]
    if missing:
        raise ValueError(f"{label} is missing required columns: {missing}")


def to_float(value) -> float:
    if pd.isna(value):
        return np.nan
    return float(value)


def accept_to_bool(value) -> bool:
    return str(value).strip().lower() == "accept"


def constant_channel(value: float) -> np.ndarray:
    return np.full((CHIP_SIZE, CHIP_SIZE), value, dtype=np.float32)


def chip_transform(chip_row: pd.Series) -> Affine:
    return from_origin(
        float(chip_row["chip_xmin"]),
        float(chip_row["chip_ymax"]),
        CHIP_RES_M,
        CHIP_RES_M,
    )


def chip_bounds(chip_row: pd.Series) -> tuple[float, float, float, float]:
    return (
        float(chip_row["chip_xmin"]),
        float(chip_row["chip_ymin"]),
        float(chip_row["chip_xmax"]),
        float(chip_row["chip_ymax"]),
    )


def expanded_window_for_bounds(
    bounds: tuple[float, float, float, float],
    transform: Affine,
    pad_pixels: int = 1,
) -> Window:
    raw = from_bounds(*bounds, transform=transform)
    col_start = math.floor(raw.col_off) - pad_pixels
    row_start = math.floor(raw.row_off) - pad_pixels
    col_stop = math.ceil(raw.col_off + raw.width) + pad_pixels
    row_stop = math.ceil(raw.row_off + raw.height) + pad_pixels
    return Window(
        col_start,
        row_start,
        col_stop - col_start,
        row_stop - row_start,
    )


def mosaic_nanmean(arrays: list[np.ndarray]) -> np.ndarray:
    if not arrays:
        return np.full((CHIP_SIZE, CHIP_SIZE), np.nan, dtype=np.float32)

    stacked = np.stack(arrays, axis=0)
    finite = np.isfinite(stacked)
    counts = finite.sum(axis=0)
    sums = np.where(finite, stacked, 0.0).sum(axis=0)
    output = np.full(stacked.shape[1:], np.nan, dtype=np.float32)
    np.divide(sums, counts, out=output, where=counts > 0)
    return output.astype(np.float32)


def normalized_difference(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    denominator = a + b
    output = np.full(a.shape, np.nan, dtype=np.float32)
    valid = np.isfinite(a) & np.isfinite(b) & (np.abs(denominator) > 1e-6)
    output[valid] = (a[valid] - b[valid]) / denominator[valid]
    return output


def fapar_composite_doy(sif_doy: int) -> int:
    composite_doy = 1 + 8 * ((int(sif_doy) - 1) // 8)
    if composite_doy not in FAPAR_COMPOSITE_DOYS:
        raise ValueError(
            f"SIF DOY {sif_doy} maps to FAPAR composite DOY {composite_doy}, "
            f"outside the downloaded range {FAPAR_COMPOSITE_DOYS[0]}-"
            f"{FAPAR_COMPOSITE_DOYS[-1]}"
        )
    return composite_doy


# ---------------------------------------------------------------------------
# Manifest and assignment loading


def load_chip_tables() -> tuple[pd.DataFrame, dict[str, pd.DataFrame]]:
    manifest = pd.read_csv(CHIP_MANIFEST_PATH)
    assignments = pd.read_csv(CHIP_ASSIGNMENTS_PATH)

    require_columns(
        manifest,
        [
            "chip_id",
            "mgrs_tile_t",
            "Delta_Date",
            "sif_year",
            "sif_doy",
            "chip_rows",
            "chip_cols",
            "chip_xmin",
            "chip_ymin",
            "chip_xmax",
            "chip_ymax",
            "chip_status",
            "n_sif",
        ],
        "chip manifest",
    )
    require_columns(
        assignments,
        [
            "chip_id",
            "sif_row_id",
            "mgrs_tile_t",
            "Delta_Date",
            "sif_year",
            "sif_doy",
            "product_path",
            "Lat_corner1",
            "Lat_corner2",
            "Lat_corner3",
            "Lat_corner4",
            "Lon_corner1",
            "Lon_corner2",
            "Lon_corner3",
            "Lon_corner4",
            TARGET_COLUMN,
            FINAL_CHECK_COLUMN,
        ],
        "chip assignments",
    )

    manifest = manifest.copy()
    assignments = assignments.copy()
    manifest["Delta_Date"] = pd.to_datetime(manifest["Delta_Date"], errors="raise").dt.date
    assignments["Delta_Date"] = pd.to_datetime(
        assignments["Delta_Date"], errors="raise"
    ).dt.date

    manifest_numeric = [
        "sif_year",
        "sif_doy",
        "chip_rows",
        "chip_cols",
        "chip_xmin",
        "chip_ymin",
        "chip_xmax",
        "chip_ymax",
        "n_sif",
    ]
    assignment_numeric = [
        "sif_row_id",
        "sif_year",
        "sif_doy",
        "track_score",
        "Latitude",
        "Longitude",
        "Lat_corner1",
        "Lat_corner2",
        "Lat_corner3",
        "Lat_corner4",
        "Lon_corner1",
        "Lon_corner2",
        "Lon_corner3",
        "Lon_corner4",
        TARGET_COLUMN,
        "Metadata.MeasurementMode",
    ]

    for column in manifest_numeric:
        manifest[column] = pd.to_numeric(manifest[column], errors="raise")
    for column in assignment_numeric:
        if column in assignments.columns:
            assignments[column] = pd.to_numeric(assignments[column], errors="coerce")

    manifest = manifest[
        (manifest["chip_status"] == "inside")
        & (manifest["chip_rows"].astype(int) == CHIP_SIZE)
        & (manifest["chip_cols"].astype(int) == CHIP_SIZE)
    ].copy()

    manifest = manifest.sort_values(
        ["sif_year", "sif_doy", "mgrs_tile_t", "chip_id"]
    ).reset_index(drop=True)

    if MAX_CHIPS is not None:
        manifest = manifest.head(MAX_CHIPS).copy()

    assignments = assignments[assignments["chip_id"].isin(manifest["chip_id"])].copy()
    sort_columns = ["chip_id"]
    if "track_score" in assignments.columns:
        sort_columns.append("track_score")
    sort_columns.append("sif_row_id")
    assignments = assignments.sort_values(sort_columns)

    groups: dict[str, pd.DataFrame] = {}
    for chip_id, group in assignments.groupby("chip_id", sort=False):
        group = group.reset_index(drop=True)
        manifest_row = manifest.loc[manifest["chip_id"] == chip_id].iloc[0]

        checks = {
            "Delta_Date": group["Delta_Date"].nunique(dropna=False),
            "mgrs_tile_t": group["mgrs_tile_t"].nunique(dropna=False),
            "product_path": group["product_path"].nunique(dropna=False),
        }
        if "Metadata.MeasurementMode" in group.columns:
            checks["Metadata.MeasurementMode"] = group[
                "Metadata.MeasurementMode"
            ].nunique(dropna=False)

        invalid = {name: count for name, count in checks.items() if count != 1}
        if invalid:
            raise ValueError(f"chip_id={chip_id} is not homogeneous: {invalid}")

        if group.iloc[0]["Delta_Date"] != manifest_row["Delta_Date"]:
            raise ValueError(f"Date mismatch between manifest and assignments for {chip_id}")
        if group.iloc[0]["mgrs_tile_t"] != manifest_row["mgrs_tile_t"]:
            raise ValueError(f"MGRS tile mismatch between manifest and assignments for {chip_id}")

        groups[chip_id] = group

    missing_groups = sorted(set(manifest["chip_id"]) - set(groups))
    if missing_groups:
        raise ValueError(f"Manifest chips have no footprint assignments: {missing_groups[:5]}")

    print(
        f"Loaded {len(manifest):,} inside chips and "
        f"{len(assignments):,} assigned SIF footprints"
    )
    return manifest, groups


# ---------------------------------------------------------------------------
# Sentinel-2 bands and spectral indices


def product_id_from_path(product_path: Path) -> str:
    return product_path.name


def sentinel_product_tile(product_path: Path) -> str:
    match = re.search(r"_T(\d{2}[A-Z]{3})_", product_id_from_path(product_path))
    if match is None:
        raise ValueError(f"Could not parse T-prefixed MGRS tile from {product_path}")
    return f"T{match.group(1)}"


def sentinel_band_path(product_path: Path, band_id: str) -> Path:
    product_id = product_id_from_path(product_path)
    path = product_path / f"{product_id}_FRC_{band_id}.tif"
    if not path.exists():
        raise FileNotFoundError(path)
    return path


def sentinel_flag_path(product_path: Path, flag_resolution: str) -> Path:
    product_id = product_id_from_path(product_path)
    path = product_path / "MASKS" / f"{product_id}_FLG_{flag_resolution}.tif"
    if not path.exists():
        raise FileNotFoundError(path)
    return path


def read_reflectance_band_to_chip(
    product_path: Path,
    band_id: str,
    bounds: tuple[float, float, float, float],
    dst_transform: Affine,
    dst_crs,
) -> np.ndarray:
    flag_resolution = SENTINEL_BANDS[band_id]
    band_path = sentinel_band_path(product_path, band_id)
    flag_path = sentinel_flag_path(product_path, flag_resolution)

    with rasterio.open(band_path) as band_src, rasterio.open(flag_path) as flag_src:
        if band_src.crs != flag_src.crs or band_src.transform != flag_src.transform:
            raise ValueError(
                f"Band and flag grids differ for {band_path.name} and {flag_path.name}"
            )
        if band_src.width != flag_src.width or band_src.height != flag_src.height:
            raise ValueError(
                f"Band and flag dimensions differ for {band_path.name} and {flag_path.name}"
            )

        window = expanded_window_for_bounds(bounds, band_src.transform, pad_pixels=1)
        band = band_src.read(
            1,
            window=window,
            boundless=True,
            fill_value=REFLECTANCE_NODATA,
        )
        flag = flag_src.read(1, window=window, boundless=True, fill_value=0)
        source_transform = band_src.window_transform(window)

        valid = (
            (flag == LAND_FLAG_VALUE)
            & (band != REFLECTANCE_NODATA)
            & np.isfinite(band)
        )
        source = np.full(band.shape, WARP_NODATA, dtype=np.float32)
        source[valid] = band[valid].astype(np.float32) / REFLECTANCE_QUANTIFICATION_VALUE

        destination = np.full((CHIP_SIZE, CHIP_SIZE), np.nan, dtype=np.float32)
        reproject(
            source=source,
            destination=destination,
            src_transform=source_transform,
            src_crs=band_src.crs,
            src_nodata=WARP_NODATA,
            dst_transform=dst_transform,
            dst_crs=dst_crs,
            dst_nodata=np.nan,
            resampling=Resampling.average,
        )

    return destination


def sentinel_indices_chip(
    product_path: Path,
    mgrs_tile_t: str,
    bounds: tuple[float, float, float, float],
    dst_transform: Affine,
) -> tuple[list[np.ndarray], object]:
    product_tile = sentinel_product_tile(product_path)
    if product_tile != mgrs_tile_t:
        raise ValueError(
            f"Product tile {product_tile} does not match assignment tile "
            f"{mgrs_tile_t}: {product_path}"
        )

    b5_path = sentinel_band_path(product_path, "B5")
    with rasterio.open(b5_path) as reference:
        dst_crs = reference.crs
        if dst_crs is None:
            raise ValueError(f"Sentinel reference raster has no CRS: {b5_path}")
        if not (
            math.isclose(abs(reference.transform.a), CHIP_RES_M, abs_tol=1e-6)
            and math.isclose(abs(reference.transform.e), CHIP_RES_M, abs_tol=1e-6)
        ):
            raise ValueError(
                f"Expected a {CHIP_RES_M:g} m B5 grid, found "
                f"{reference.transform.a}, {reference.transform.e}: {b5_path}"
            )

        expected_bounds = (
            dst_transform.c,
            dst_transform.f + CHIP_SIZE * dst_transform.e,
            dst_transform.c + CHIP_SIZE * dst_transform.a,
            dst_transform.f,
        )
        tolerance = 1e-6
        if any(abs(a - b) > tolerance for a, b in zip(bounds, expected_bounds)):
            raise ValueError("Chip bounds and transform are inconsistent")

        inside_reference = (
            bounds[0] >= reference.bounds.left - tolerance
            and bounds[1] >= reference.bounds.bottom - tolerance
            and bounds[2] <= reference.bounds.right + tolerance
            and bounds[3] <= reference.bounds.top + tolerance
        )
        if not inside_reference:
            raise ValueError(
                f"Inside chip falls outside its Sentinel product raster: {product_path}"
            )

        col_offset = (bounds[0] - reference.transform.c) / reference.transform.a
        row_offset = (reference.transform.f - bounds[3]) / abs(reference.transform.e)
        if not (
            math.isclose(col_offset, round(col_offset), abs_tol=1e-6)
            and math.isclose(row_offset, round(row_offset), abs_tol=1e-6)
        ):
            raise ValueError(
                f"Chip is not aligned to the Sentinel B5 pixel grid: {product_path}"
            )

    bands = {
        band_id: read_reflectance_band_to_chip(
            product_path,
            band_id,
            bounds,
            dst_transform,
            dst_crs,
        )
        for band_id in SENTINEL_BANDS
    }

    ndmi = normalized_difference(bands["B8"], bands["B11"])
    ndvi = normalized_difference(bands["B8"], bands["B4"])

    evi_denominator = bands["B8"] + 6.0 * bands["B4"] - 7.5 * bands["B2"] + 1.0
    evi = np.full((CHIP_SIZE, CHIP_SIZE), np.nan, dtype=np.float32)
    evi_valid = (
        np.isfinite(bands["B8"])
        & np.isfinite(bands["B4"])
        & np.isfinite(bands["B2"])
        & (np.abs(evi_denominator) > 1e-6)
    )
    evi[evi_valid] = 2.5 * (
        (bands["B8"][evi_valid] - bands["B4"][evi_valid])
        / evi_denominator[evi_valid]
    )

    nirv = bands["B8"] * ndvi
    ndre = normalized_difference(bands["B8A"], bands["B5"])

    return [ndmi, ndvi, evi, nirv, ndre], dst_crs


# ---------------------------------------------------------------------------
# FAPAR


@lru_cache(maxsize=512)
def fapar_path(year: int, doy: int, tile: str) -> Path:
    pattern = str(FAPAR_DIR / tile / str(year) / f"*A{year}{doy:03d}.{tile}.*.tif")
    matches = sorted(Path(path) for path in glob.glob(pattern))
    if len(matches) != 1:
        raise FileNotFoundError(
            f"Expected one FAPAR file for year={year}, doy={doy}, tile={tile}; "
            f"found {len(matches)}"
        )
    return matches[0]


def reproject_fapar_tile_to_chip(
    path: Path,
    bounds: tuple[float, float, float, float],
    dst_transform: Affine,
    dst_crs,
) -> np.ndarray | None:
    with rasterio.open(path) as src:
        source_bounds = transform_bounds(
            dst_crs,
            src.crs,
            *bounds,
            densify_pts=21,
        )

        intersects = not (
            source_bounds[2] <= src.bounds.left
            or source_bounds[0] >= src.bounds.right
            or source_bounds[3] <= src.bounds.bottom
            or source_bounds[1] >= src.bounds.top
        )
        if not intersects:
            return None

        window = expanded_window_for_bounds(source_bounds, src.transform, pad_pixels=2)
        values = src.read(
            1,
            window=window,
            boundless=True,
            fill_value=src.nodata if src.nodata is not None else WARP_NODATA,
        ).astype(np.float32)
        source_transform = src.window_transform(window)

        invalid = ~np.isfinite(values) | (values < 0.0) | (values > 1.0)
        if src.nodata is not None:
            invalid |= values == src.nodata
        values[invalid] = WARP_NODATA

        destination = np.full((CHIP_SIZE, CHIP_SIZE), np.nan, dtype=np.float32)
        reproject(
            source=values,
            destination=destination,
            src_transform=source_transform,
            src_crs=src.crs,
            src_nodata=WARP_NODATA,
            dst_transform=dst_transform,
            dst_crs=dst_crs,
            dst_nodata=np.nan,
            resampling=Resampling.bilinear,
        )

    return destination


def fapar_chip(
    year: int,
    composite_doy: int,
    bounds: tuple[float, float, float, float],
    dst_transform: Affine,
    dst_crs,
) -> np.ndarray:
    arrays = []
    for tile in MODIS_TILES:
        values = reproject_fapar_tile_to_chip(
            fapar_path(year, composite_doy, tile),
            bounds,
            dst_transform,
            dst_crs,
        )
        if values is not None:
            arrays.append(values)
    return mosaic_nanmean(arrays)


# ---------------------------------------------------------------------------
# Daily VIIRS PAR


@lru_cache(maxsize=2048)
def par_paths(date_iso: str) -> tuple[Path, Path]:
    year = date_iso[:4]
    year_dir = PAR_DIR / year
    par_path = (
        year_dir
        / f"VNP18A2.002_{date_iso}_Daily_Mean_PAR_VIIRS_Sinusoidal_native.tif"
    )
    qa_path = (
        year_dir
        / f"VNP18A2.002_{date_iso}_PAR_Quality_VIIRS_Sinusoidal_native.tif"
    )

    if not par_path.exists():
        raise FileNotFoundError(par_path)
    if not qa_path.exists():
        raise FileNotFoundError(qa_path)
    return par_path, qa_path


def par_chip(
    delta_date,
    bounds: tuple[float, float, float, float],
    dst_transform: Affine,
    dst_crs,
) -> np.ndarray:
    date_iso = pd.Timestamp(delta_date).strftime("%Y-%m-%d")
    par_path, qa_path = par_paths(date_iso)

    with rasterio.open(par_path) as par_src, rasterio.open(qa_path) as qa_src:
        if par_src.crs is None or qa_src.crs is None:
            raise ValueError(f"PAR or QA raster has no CRS for {date_iso}")
        if par_src.crs != qa_src.crs or par_src.transform != qa_src.transform:
            raise ValueError(
                f"PAR and QA grids differ for {date_iso}: "
                f"{par_path.name}, {qa_path.name}"
            )
        if par_src.width != qa_src.width or par_src.height != qa_src.height:
            raise ValueError(f"PAR and QA dimensions differ for {date_iso}")

        # EPSG:9001 in the source WKT identifies the metre unit. Rasterio uses
        # the complete embedded Clarke-1866 sinusoidal CRS for this transform.
        source_bounds = transform_bounds(
            dst_crs,
            par_src.crs,
            *bounds,
            densify_pts=21,
        )

        intersects = not (
            source_bounds[2] <= par_src.bounds.left
            or source_bounds[0] >= par_src.bounds.right
            or source_bounds[3] <= par_src.bounds.bottom
            or source_bounds[1] >= par_src.bounds.top
        )
        if not intersects:
            raise ValueError(
                f"Sentinel chip does not intersect the Germany PAR raster for {date_iso}"
            )

        window = expanded_window_for_bounds(source_bounds, par_src.transform, pad_pixels=2)
        par_values = par_src.read(
            1,
            window=window,
            boundless=True,
            fill_value=PAR_FILL_VALUE,
        ).astype(np.float32)
        qa_values = qa_src.read(
            1,
            window=window,
            boundless=True,
            fill_value=0,
        )
        source_transform = par_src.window_transform(window)

        valid = (
            np.isin(qa_values, PAR_ACCEPTED_QA_CODES)
            & np.isfinite(par_values)
            & (par_values != PAR_FILL_VALUE)
            & (par_values >= PAR_VALID_MIN)
            & (par_values <= PAR_VALID_MAX)
        )
        if par_src.nodata is not None:
            valid &= par_values != par_src.nodata

        # Apply categorical QA on the native grid before any interpolation.
        source = np.full(par_values.shape, WARP_NODATA, dtype=np.float32)
        source[valid] = par_values[valid]

        destination = np.full((CHIP_SIZE, CHIP_SIZE), np.nan, dtype=np.float32)
        reproject(
            source=source,
            destination=destination,
            src_transform=source_transform,
            src_crs=par_src.crs,
            src_nodata=WARP_NODATA,
            dst_transform=dst_transform,
            dst_crs=dst_crs,
            dst_nodata=np.nan,
            resampling=Resampling.bilinear,
        )

    return destination


# ---------------------------------------------------------------------------
# Crop fractions


@lru_cache(maxsize=16)
def crop_path_for_year(year: int) -> Path:
    path = CROP_DIR / f"croptypes_{year}.tif"
    if not path.exists():
        raise FileNotFoundError(path)
    return path


def active_codes_for_month(month: int) -> list[int]:
    return [
        code
        for code, active_months in ACTIVE_GROWTH_MONTHS_BY_CODE.items()
        if month in active_months
    ]


def fraction_from_crop_codes(
    crop_values: np.ndarray,
    crop_transform: Affine,
    crop_crs,
    group_codes: Iterable[int],
    dst_transform: Affine,
    dst_crs,
) -> np.ndarray:
    codes = list(group_codes)
    if not codes:
        return np.zeros((CHIP_SIZE, CHIP_SIZE), dtype=np.float32)

    source = np.isin(crop_values, codes).astype(np.float32)
    destination = np.zeros((CHIP_SIZE, CHIP_SIZE), dtype=np.float32)
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
    return np.clip(destination, 0.0, 1.0).astype(np.float32)


def crop_fraction_chips(
    year: int,
    month: int,
    bounds: tuple[float, float, float, float],
    dst_transform: Affine,
    dst_crs,
) -> list[np.ndarray]:
    with rasterio.open(crop_path_for_year(year)) as src:
        source_bounds = transform_bounds(
            dst_crs,
            src.crs,
            *bounds,
            densify_pts=21,
        )
        window = expanded_window_for_bounds(source_bounds, src.transform, pad_pixels=2)
        crop_values = src.read(1, window=window, boundless=True, fill_value=NON_CROP_CODE)
        crop_transform = src.window_transform(window)

        channels = [
            fraction_from_crop_codes(
                crop_values,
                crop_transform,
                src.crs,
                codes,
                dst_transform,
                dst_crs,
            )
            for codes in CROP_GROUPS.values()
        ]

        active_crop = fraction_from_crop_codes(
            crop_values,
            crop_transform,
            src.crs,
            active_codes_for_month(month),
            dst_transform,
            dst_crs,
        )

        non_crop_source = (~np.isin(crop_values, KNOWN_CROP_CODES)).astype(np.float32)
        non_crop = np.zeros((CHIP_SIZE, CHIP_SIZE), dtype=np.float32)
        reproject(
            source=non_crop_source,
            destination=non_crop,
            src_transform=crop_transform,
            src_crs=src.crs,
            dst_transform=dst_transform,
            dst_crs=dst_crs,
            dst_nodata=1.0,
            resampling=Resampling.average,
        )

    channels.append(np.clip(active_crop, 0.0, 1.0).astype(np.float32))
    channels.append(np.clip(non_crop, 0.0, 1.0).astype(np.float32))
    return channels


# ---------------------------------------------------------------------------
# Footprint masks


def build_sif_polygon_wgs84(row: pd.Series) -> Polygon:
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


def build_sif_polygon_projected(row: pd.Series, dst_crs):
    polygon = build_sif_polygon_wgs84(row)
    projected = transform_geom(
        "EPSG:4326",
        dst_crs,
        mapping(polygon),
        precision=6,
    )
    return shape(projected)


def fractional_footprint_mask(geometry, dst_transform: Affine) -> np.ndarray:
    high_size = CHIP_SIZE * MASK_OVERSAMPLE
    high_res = CHIP_RES_M / MASK_OVERSAMPLE
    high_transform = from_origin(dst_transform.c, dst_transform.f, high_res, high_res)

    high_mask = rasterize(
        [(mapping(geometry), 1.0)],
        out_shape=(high_size, high_size),
        transform=high_transform,
        fill=0.0,
        dtype="float32",
        all_touched=False,
    )

    return high_mask.reshape(
        CHIP_SIZE,
        MASK_OVERSAMPLE,
        CHIP_SIZE,
        MASK_OVERSAMPLE,
    ).mean(axis=(1, 3)).astype(np.float32)


def mask_inside_fraction(mask: np.ndarray, geometry) -> float:
    expected_mask_sum = float(geometry.area) / (CHIP_RES_M * CHIP_RES_M)
    if expected_mask_sum <= 0:
        return 0.0
    return float(min(mask.sum() / expected_mask_sum, 1.0))


def build_multisif_masks(
    footprints: pd.DataFrame,
    dst_transform: Affine,
    dst_crs,
) -> tuple[dict, list[dict]]:
    masks = np.zeros((MAX_FOOTPRINTS, CHIP_SIZE, CHIP_SIZE), dtype=np.float32)
    y_targets = np.full(MAX_FOOTPRINTS, np.nan, dtype=np.float32)
    target_accept = np.zeros(MAX_FOOTPRINTS, dtype=np.uint8)
    footprint_valid = np.zeros(MAX_FOOTPRINTS, dtype=np.uint8)
    sif_row_ids = np.full(MAX_FOOTPRINTS, -1, dtype=np.int64)
    mask_sums = np.zeros(MAX_FOOTPRINTS, dtype=np.float32)
    inside_fractions = np.zeros(MAX_FOOTPRINTS, dtype=np.float32)

    kept_rows: list[dict] = []
    skipped_rows: list[dict] = []

    for _, footprint in footprints.iterrows():
        geometry = build_sif_polygon_projected(footprint, dst_crs)
        mask = fractional_footprint_mask(geometry, dst_transform)
        mask_sum = float(mask.sum())
        inside_fraction = mask_inside_fraction(mask, geometry)

        skip_reason = ""
        if mask_sum <= 0:
            skip_reason = "empty_mask"
        elif inside_fraction < MIN_FOOTPRINT_INSIDE_FRACTION:
            skip_reason = "clipped_footprint"

        footprint_meta = {
            "sif_row_id": int(footprint["sif_row_id"]),
            TARGET_COLUMN: to_float(footprint[TARGET_COLUMN]),
            FINAL_CHECK_COLUMN: footprint[FINAL_CHECK_COLUMN],
            "mask_sum": mask_sum,
            "mask_inside_fraction": inside_fraction,
            "state": footprint.get("state", ""),
            "hzs": footprint.get("hzs", ""),
            "measurement_mode": to_float(
                footprint.get("Metadata.MeasurementMode", np.nan)
            ),
        }

        if skip_reason:
            skipped_rows.append(
                {**footprint_meta, "slot": -1, "skip_reason": skip_reason}
            )
            continue

        slot = len(kept_rows)
        masks[slot] = mask
        y_targets[slot] = to_float(footprint[TARGET_COLUMN])
        target_accept[slot] = accept_to_bool(footprint[FINAL_CHECK_COLUMN])
        footprint_valid[slot] = 1
        sif_row_ids[slot] = int(footprint["sif_row_id"])
        mask_sums[slot] = mask_sum
        inside_fractions[slot] = inside_fraction
        kept_rows.append({**footprint_meta, "slot": slot, "skip_reason": ""})

        if len(kept_rows) == MAX_FOOTPRINTS:
            break

    payload = {
        "footprint_masks": masks.astype(np.float16),
        "y_targets": y_targets,
        "target_accept": target_accept,
        "footprint_valid": footprint_valid,
        "sif_row_ids": sif_row_ids,
        "mask_sums": mask_sums,
        "mask_inside_fractions": inside_fractions,
        "n_valid_footprints": int(footprint_valid.sum()),
        "n_skipped_footprints": len(skipped_rows),
    }
    return payload, kept_rows + skipped_rows


# ---------------------------------------------------------------------------
# Complete chip assembly


def build_predictor_chip(
    chip_row: pd.Series,
    footprints: pd.DataFrame,
    dst_transform: Affine,
) -> tuple[np.ndarray, dict, object, Path, int]:
    first = footprints.iloc[0]
    product_path = Path(str(first["product_path"]))
    if not product_path.is_dir():
        raise FileNotFoundError(product_path)

    mgrs_tile_t = str(chip_row["mgrs_tile_t"])
    bounds = chip_bounds(chip_row)
    spectral_indices, dst_crs = sentinel_indices_chip(
        product_path,
        mgrs_tile_t,
        bounds,
        dst_transform,
    )

    year = int(chip_row["sif_year"])
    sif_doy = int(chip_row["sif_doy"])
    composite_doy = fapar_composite_doy(sif_doy)
    month = pd.Timestamp(chip_row["Delta_Date"]).month

    fapar = fapar_chip(
        year,
        composite_doy,
        bounds,
        dst_transform,
        dst_crs,
    )
    par = par_chip(
        chip_row["Delta_Date"],
        bounds,
        dst_transform,
        dst_crs,
    )
    apar = fapar * par
    crop_channels = crop_fraction_chips(
        year,
        month,
        bounds,
        dst_transform,
        dst_crs,
    )

    month_angle = 2.0 * np.pi * month / 12.0
    channels = [
        *spectral_indices,
        fapar,
        par,
        apar,
        *crop_channels,
        constant_channel(np.sin(month_angle)),
        constant_channel(np.cos(month_angle)),
    ]

    if len(channels) != len(CHANNEL_NAMES):
        raise RuntimeError(
            f"Built {len(channels)} channels but expected {len(CHANNEL_NAMES)}"
        )

    # Avoid float16 overflow without imposing additional ecological clipping.
    float16_max = np.finfo(np.float16).max
    for channel in channels:
        channel[~np.isfinite(channel)] = np.nan
        channel[np.abs(channel) > float16_max] = np.nan

    valid_fraction = {
        f"{name}_valid_fraction": float(np.isfinite(channel).mean())
        for name, channel in zip(CHANNEL_NAMES, channels)
    }

    # NaNs are intentionally retained for nan-aware training normalization.
    x = np.stack(channels, axis=0).astype(np.float16)
    return x, valid_fraction, dst_crs, product_path, composite_doy


def build_chip(
    chip_row: pd.Series,
    footprints: pd.DataFrame,
) -> tuple[dict, dict, list[dict]]:
    dst_transform = chip_transform(chip_row)
    x, valid_fraction, dst_crs, product_path, composite_doy = build_predictor_chip(
        chip_row,
        footprints,
        dst_transform,
    )
    mask_payload, footprint_metadata = build_multisif_masks(
        footprints,
        dst_transform,
        dst_crs,
    )

    measurement_mode = to_float(
        footprints.iloc[0].get("Metadata.MeasurementMode", np.nan)
    )
    chip_metadata = {
        "chip_id": chip_row["chip_id"],
        "Delta_Date": chip_row["Delta_Date"],
        "sif_year": int(chip_row["sif_year"]),
        "sif_doy": int(chip_row["sif_doy"]),
        "fapar_composite_doy": composite_doy,
        "par_date": chip_row["Delta_Date"],
        "mgrs_tile_t": chip_row["mgrs_tile_t"],
        "product_path": str(product_path),
        "measurement_mode": measurement_mode,
        "chip_status": chip_row["chip_status"],
        "chip_rows": int(chip_row["chip_rows"]),
        "chip_cols": int(chip_row["chip_cols"]),
        "chip_size_m": float(chip_row.get("chip_size_m", CHIP_SIZE_M)),
        "chip_xmin": float(chip_row["chip_xmin"]),
        "chip_ymin": float(chip_row["chip_ymin"]),
        "chip_xmax": float(chip_row["chip_xmax"]),
        "chip_ymax": float(chip_row["chip_ymax"]),
        "n_sif_assigned": int(len(footprints)),
        "n_valid_footprints": mask_payload["n_valid_footprints"],
        "n_skipped_footprints": mask_payload["n_skipped_footprints"],
        "target_accepted_footprints": int(
            (mask_payload["target_accept"] * mask_payload["footprint_valid"]).sum()
        ),
        "states": chip_row.get("states", ""),
        "hzs_values": chip_row.get("hzs_values", ""),
        "mean_target_modis_sif": to_float(
            chip_row.get("mean_target_modis_sif", np.nan)
        ),
        **valid_fraction,
    }

    return {"X": x, **mask_payload}, chip_metadata, footprint_metadata


# ---------------------------------------------------------------------------
# Optional multiprocessing


def make_payload(
    index: int,
    chip_row: pd.Series,
    footprints: pd.DataFrame,
) -> tuple[int, dict, list[dict]]:
    return index, chip_row.to_dict(), footprints.to_dict(orient="records")


def build_chip_from_payload(
    payload: tuple[int, dict, list[dict]],
) -> tuple[int, dict, dict, list[dict]]:
    index, chip_dict, footprint_dicts = payload
    chip_payload, chip_metadata, footprint_metadata = build_chip(
        pd.Series(chip_dict),
        pd.DataFrame(footprint_dicts),
    )
    return index, chip_payload, chip_metadata, footprint_metadata


# ---------------------------------------------------------------------------
# Shard writing


def write_shard(
    shard_id: int,
    xs: list[np.ndarray],
    masks: list[np.ndarray],
    y_targets: list[np.ndarray],
    target_accepts: list[np.ndarray],
    footprint_valids: list[np.ndarray],
    sif_row_ids: list[np.ndarray],
    mask_sums: list[np.ndarray],
    mask_inside_fractions: list[np.ndarray],
    chip_ids: list[str],
) -> Path:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    path = OUTPUT_DIR / f"chips_{shard_id:05d}.npz"
    np.savez_compressed(
        path,
        X=np.stack(xs, axis=0).astype(np.float16),
        footprint_masks=np.stack(masks, axis=0).astype(np.float16),
        y_targets=np.stack(y_targets, axis=0).astype(np.float32),
        target_accept=np.stack(target_accepts, axis=0).astype(np.uint8),
        footprint_valid=np.stack(footprint_valids, axis=0).astype(np.uint8),
        sif_row_ids=np.stack(sif_row_ids, axis=0).astype(np.int64),
        mask_sums=np.stack(mask_sums, axis=0).astype(np.float32),
        mask_inside_fractions=np.stack(mask_inside_fractions, axis=0).astype(np.float32),
        chip_id=np.asarray(chip_ids),
        target_names=np.asarray([TARGET_COLUMN]),
        final_check_names=np.asarray([FINAL_CHECK_COLUMN]),
        channel_names=np.asarray(CHANNEL_NAMES),
    )
    return path


def prepare_chips() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    manifest, assignment_groups = load_chip_tables()

    metadata_rows: list[dict] = []
    footprint_metadata_rows: list[dict] = []

    xs: list[np.ndarray] = []
    masks: list[np.ndarray] = []
    y_targets: list[np.ndarray] = []
    target_accepts: list[np.ndarray] = []
    footprint_valids: list[np.ndarray] = []
    sif_row_ids: list[np.ndarray] = []
    mask_sums: list[np.ndarray] = []
    mask_inside_fractions: list[np.ndarray] = []
    chip_ids: list[str] = []
    shard_id = 0

    def flush_shard() -> None:
        nonlocal shard_id
        nonlocal xs, masks, y_targets, target_accepts, footprint_valids
        nonlocal sif_row_ids, mask_sums, mask_inside_fractions, chip_ids

        if not xs:
            return

        written = write_shard(
            shard_id,
            xs,
            masks,
            y_targets,
            target_accepts,
            footprint_valids,
            sif_row_ids,
            mask_sums,
            mask_inside_fractions,
            chip_ids,
        )
        print(f"Wrote {written}")
        shard_id += 1
        xs = []
        masks = []
        y_targets = []
        target_accepts = []
        footprint_valids = []
        sif_row_ids = []
        mask_sums = []
        mask_inside_fractions = []
        chip_ids = []

    def handle_result(result: tuple[int, dict, dict, list[dict]]) -> None:
        index, chip_payload, chip_metadata, footprint_metadata = result
        if index % 10 == 0:
            print(f"Prepared chip {index + 1:,} / {len(manifest):,}")

        chip_id = chip_metadata["chip_id"]
        for row in footprint_metadata:
            footprint_metadata_rows.append({**row, "chip_id": chip_id})

        if chip_metadata["n_valid_footprints"] < MIN_VALID_FOOTPRINTS:
            print(
                f"Skipping chip_id={chip_id}: only "
                f"{chip_metadata['n_valid_footprints']} footprint masks remain"
            )
            return

        xs.append(chip_payload["X"])
        masks.append(chip_payload["footprint_masks"])
        y_targets.append(chip_payload["y_targets"])
        target_accepts.append(chip_payload["target_accept"])
        footprint_valids.append(chip_payload["footprint_valid"])
        sif_row_ids.append(chip_payload["sif_row_ids"])
        mask_sums.append(chip_payload["mask_sums"])
        mask_inside_fractions.append(chip_payload["mask_inside_fractions"])
        chip_ids.append(chip_id)
        metadata_rows.append(chip_metadata)

        if len(xs) == SHARD_SIZE:
            flush_shard()

    payloads = (
        make_payload(index, chip_row, assignment_groups[chip_row["chip_id"]])
        for index, chip_row in manifest.iterrows()
    )

    if N_WORKERS == 1:
        for payload in payloads:
            handle_result(build_chip_from_payload(payload))
    else:
        with ProcessPoolExecutor(max_workers=N_WORKERS) as executor:
            for result in executor.map(
                build_chip_from_payload,
                payloads,
                chunksize=1,
            ):
                handle_result(result)

    flush_shard()

    metadata = pd.DataFrame(metadata_rows)
    metadata.to_csv(OUTPUT_DIR / "chip_metadata.csv", index=False)

    footprint_metadata = pd.DataFrame(footprint_metadata_rows)
    footprint_metadata.to_csv(OUTPUT_DIR / "footprint_metadata.csv", index=False)

    pd.DataFrame(
        {
            "channel_index": np.arange(len(CHANNEL_NAMES)),
            "channel_name": CHANNEL_NAMES,
        }
    ).to_csv(OUTPUT_DIR / "channel_names.csv", index=False)

    pd.DataFrame(
        {
            "target_index": [0],
            "target_name": [TARGET_COLUMN],
            "final_check_column": [FINAL_CHECK_COLUMN],
        }
    ).to_csv(OUTPUT_DIR / "target_names.csv", index=False)

    config = {
        "chip_size_m": CHIP_SIZE_M,
        "chip_resolution_m": CHIP_RES_M,
        "chip_rows": CHIP_SIZE,
        "chip_cols": CHIP_SIZE,
        "channels": CHANNEL_NAMES,
        "target": TARGET_COLUMN,
        "final_check": FINAL_CHECK_COLUMN,
        "max_footprints": MAX_FOOTPRINTS,
        "min_valid_footprints": MIN_VALID_FOOTPRINTS,
        "min_footprint_inside_fraction": MIN_FOOTPRINT_INSIDE_FRACTION,
        "mask_oversample": MASK_OVERSAMPLE,
        "predictor_storage_dtype": "float16",
        "mask_storage_dtype": "float16",
        "target_storage_dtype": "float32",
        "predictor_nan_policy": "preserve; normalize with nan-aware stats then fill with 0",
        "par_units": "W/m2",
        "par_valid_range": [PAR_VALID_MIN, PAR_VALID_MAX],
        "par_accepted_qa_codes": list(PAR_ACCEPTED_QA_CODES),
        "par_resampling": "bilinear after native-grid QA masking",
    }
    with (OUTPUT_DIR / "dataset_config.json").open("w", encoding="ascii") as file:
        json.dump(config, file, indent=2)

    print(f"Done. Wrote {len(metadata):,} chips to {OUTPUT_DIR}")
    print(
        "Training reminder: cast X/masks to float32, compute nan-aware training "
        "statistics, normalize, then replace normalized NaNs with zero. Average "
        "footprint losses within each chip before averaging across the batch."
    )


if __name__ == "__main__":
    prepare_chips()
