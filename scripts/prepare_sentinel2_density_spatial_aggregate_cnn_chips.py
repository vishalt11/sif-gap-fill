"""Build model-ready CNN chips for density-aggregated 4 km SIF windows.

The input tables are the consolidated density-window manifest and footprint
assignments retained after the Sentinel-2 valid-fraction filter. Each output
sample contains:

  X:                    [19, 200, 200] predictor channels, stored as float16
  aggregate_weight_map: [200, 200] equal-footprint weights, stored as float16
  y_aggregate:          mean target_modis_sif across assigned footprints

The downloaded Sentinel-2 L2A GeoTIFF already represents the exact 4 km
window at 20 m resolution. Its six named reflectance bands are used to derive
five spectral indices. The four Sentinel quality bands are deliberately not
included in X or copied into the model metadata.

No train/validation/test split is made here. Predictor normalization and all
partitioning must be fitted or defined later from the training partition only.
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
import json
import math
from pathlib import Path
import re
from typing import Iterable

import numpy as np
import pandas as pd
import rasterio

from model_data_prep import prepare_sentinel2_multisif_cnn_chips as predictors


# ---------------------------------------------------------------------------
# Configuration

PROJECT_ROOT = Path(__file__).resolve().parent

AGGREGATION_DIR = (
    PROJECT_ROOT
    / "data"
    / "density_aggregation"
    / "sentinel2_spatial_aggregation_density_4000m_landcover_combined_"
    "mask60_min4_s2valid95_par_available"
)
MANIFEST_PATH = AGGREGATION_DIR / "density_cluster_4000m_aggregate_manifest.csv"
ASSIGNMENTS_PATH = AGGREGATION_DIR / "density_cluster_4000m_sif_assignments.csv"

OUTPUT_DIR = (
    PROJECT_ROOT
    / "data"
    / "cnn_sentinel2_chips"
    / "spatial_aggregate_density_4km_20m_s2valid95"
)

TARGET_COLUMN = "target_modis_sif"
AGGREGATE_TARGET_COLUMN = "aggregated_target_modis_sif"

CHIP_SIZE_M = 4000.0
CHIP_RES_M = 20.0
CHIP_SIZE = 200
MIN_FOOTPRINTS = 4
MIN_ASSIGNED_FOOTPRINT_INSIDE_FRACTION = 0.60
MIN_SENTINEL_VALID_FRACTION = 0.95
MASK_OVERSAMPLE = 4

# The FAPAR coverage used here is intentionally restricted to the two h18
# MODIS tiles that cover the German study area.
FAPAR_MODIS_TILES = ("h18v03", "h18v04")

SENTINEL_REFLECTANCE_BANDS = (
    "B02",
    "B04",
    "B05",
    "B08",
    "B8A",
    "B11",
)

# Sixteen samples keep the compressed shard files at a manageable size.
SHARD_SIZE = 16

# Set to a small positive integer for a smoke test. Leave as None for all
# consolidated windows.
MAX_CHIPS: int | None = None

# Four workers are a reasonable starting point for the raster-heavy workload.
# Reduce this if simultaneous reads cause storage contention.
N_WORKERS = 4

# Prevent shards from different runs being mixed in one output directory.
FAIL_IF_OUTPUT_EXISTS = True

CHANNEL_NAMES = predictors.CHANNEL_NAMES
if len(CHANNEL_NAMES) != 19:
    raise RuntimeError(f"Expected 19 predictor channels, found {len(CHANNEL_NAMES)}")

# Reuse the established FAPAR, PAR, crop and footprint-mask routines at the
# dimensions of the new density-window dataset.
predictors.FAPAR_DIR = PROJECT_ROOT / "data" / "glass_geotiff" / "fapar"
predictors.PAR_DIR = (
    PROJECT_ROOT / "data" / "viirs_vnp18a2_daily_mean_par_germany_native"
)
predictors.CROP_DIR = PROJECT_ROOT / "data" / "crop_type_tif"
predictors.CHIP_SIZE_M = CHIP_SIZE_M
predictors.CHIP_RES_M = CHIP_RES_M
predictors.CHIP_SIZE = CHIP_SIZE
predictors.MASK_OVERSAMPLE = MASK_OVERSAMPLE
predictors.MODIS_TILES = FAPAR_MODIS_TILES


# ---------------------------------------------------------------------------
# Generic helpers

def require_columns(
    table: pd.DataFrame,
    columns: Iterable[str],
    table_name: str,
) -> None:
    missing = [column for column in columns if column not in table.columns]
    if missing:
        raise ValueError(f"{table_name} is missing required columns: {missing}")


def parse_bool(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False)
    return (
        series.astype("string")
        .str.strip()
        .str.lower()
        .isin({"true", "t", "1", "yes"})
    )


def normalize_mgrs_tile(value: object) -> str:
    tile = str(value).strip().upper()
    if not tile.startswith("T"):
        tile = f"T{tile}"
    if re.fullmatch(r"T[0-9]{2}[C-X][A-Z]{2}", tile) is None:
        raise ValueError(f"Invalid MGRS tile: {value}")
    return tile


def epsg_from_mgrs_tile(value: object) -> int:
    tile = normalize_mgrs_tile(value)
    zone = int(tile[1:3])
    latitude_band = tile[3]
    return (32600 if latitude_band >= "N" else 32700) + zone


def resolve_project_path(value: object) -> Path:
    if pd.isna(value) or not str(value).strip():
        raise ValueError("Empty Sentinel-2 raster path")
    path = Path(str(value).strip())
    if not path.is_absolute():
        path = PROJECT_ROOT / path
    return path.resolve()


def optional_value(row: pd.Series, column: str, default=""):
    value = row.get(column, default)
    if pd.isna(value):
        return default
    return value


def window_bounds(row: pd.Series) -> tuple[float, float, float, float]:
    return (
        float(row["cell_xmin"]),
        float(row["cell_ymin"]),
        float(row["cell_xmax"]),
        float(row["cell_ymax"]),
    )


def require_nonempty_channel(name: str, values: np.ndarray, sample_id: str) -> None:
    if values.shape != (CHIP_SIZE, CHIP_SIZE):
        raise ValueError(
            f"{sample_id}: channel {name} has shape {values.shape}, expected "
            f"{(CHIP_SIZE, CHIP_SIZE)}"
        )
    if not np.isfinite(values).any():
        raise ValueError(f"{sample_id}: channel {name} contains no finite pixels")


# ---------------------------------------------------------------------------
# Input tables

def load_aggregation_tables() -> tuple[pd.DataFrame, dict[str, pd.DataFrame]]:
    manifest = pd.read_csv(MANIFEST_PATH, low_memory=False)
    assignments = pd.read_csv(ASSIGNMENTS_PATH, low_memory=False)

    require_columns(
        manifest,
        [
            "aggregation_id",
            "cell_id",
            "aggregation_size_m",
            "cell_pixels",
            "cell_xmin",
            "cell_ymin",
            "cell_xmax",
            "cell_ymax",
            "mgrs_tile_t",
            "window_crs",
            "Delta_Date",
            "sif_year",
            "sif_month",
            "sif_doy",
            "measurement_mode",
            "n_footprints",
            AGGREGATE_TARGET_COLUMN,
            "sentinel2_raster_path",
            "sentinel2_valid_fraction_before_fill",
            "sentinel2_neighbour_fill_complete_bool",
            "passes_sentinel2_quality_filter",
        ],
        "aggregate manifest",
    )
    require_columns(
        assignments,
        [
            "aggregation_id",
            "sif_row_id",
            "source_csv_row",
            "Delta_Date",
            "sif_year",
            "sif_month",
            "sif_doy",
            "measurement_mode",
            "mgrs_tile_t",
            TARGET_COLUMN,
            "Lat_corner1",
            "Lat_corner2",
            "Lat_corner3",
            "Lat_corner4",
            "Lon_corner1",
            "Lon_corner2",
            "Lon_corner3",
            "Lon_corner4",
            "window_mask_inside_fraction",
        ],
        "SIF assignments",
    )

    manifest = manifest.copy()
    assignments = assignments.copy()
    manifest["aggregation_id"] = manifest["aggregation_id"].astype(str)
    assignments["aggregation_id"] = assignments["aggregation_id"].astype(str)

    if manifest["aggregation_id"].duplicated().any():
        duplicate = manifest.loc[
            manifest["aggregation_id"].duplicated(), "aggregation_id"
        ].iloc[0]
        raise ValueError(f"Duplicate manifest aggregation_id: {duplicate}")

    manifest["Delta_Date"] = pd.to_datetime(
        manifest["Delta_Date"], errors="raise"
    ).dt.date
    assignments["Delta_Date"] = pd.to_datetime(
        assignments["Delta_Date"], errors="raise"
    ).dt.date

    manifest_numeric = [
        "aggregation_size_m",
        "cell_pixels",
        "cell_xmin",
        "cell_ymin",
        "cell_xmax",
        "cell_ymax",
        "sif_year",
        "sif_month",
        "sif_doy",
        "measurement_mode",
        "n_footprints",
        AGGREGATE_TARGET_COLUMN,
        "sentinel2_valid_fraction_before_fill",
    ]
    assignment_numeric = [
        "sif_row_id",
        "source_csv_row",
        "sif_year",
        "sif_month",
        "sif_doy",
        "measurement_mode",
        TARGET_COLUMN,
        "Lat_corner1",
        "Lat_corner2",
        "Lat_corner3",
        "Lat_corner4",
        "Lon_corner1",
        "Lon_corner2",
        "Lon_corner3",
        "Lon_corner4",
        "window_mask_inside_fraction",
    ]
    for column in manifest_numeric:
        manifest[column] = pd.to_numeric(manifest[column], errors="raise")
    for column in assignment_numeric:
        assignments[column] = pd.to_numeric(assignments[column], errors="raise")

    manifest["mgrs_tile_t"] = manifest["mgrs_tile_t"].map(normalize_mgrs_tile)
    assignments["mgrs_tile_t"] = assignments["mgrs_tile_t"].map(
        normalize_mgrs_tile
    )

    accepted = (
        parse_bool(manifest["passes_sentinel2_quality_filter"])
        & parse_bool(manifest["sentinel2_neighbour_fill_complete_bool"])
        & (
            manifest["sentinel2_valid_fraction_before_fill"]
            >= MIN_SENTINEL_VALID_FRACTION
        )
    )
    if not accepted.all():
        raise ValueError(
            "The consolidated manifest contains rows that fail its Sentinel-2 "
            "quality criteria"
        )

    size_ok = np.isclose(
        manifest["aggregation_size_m"], CHIP_SIZE_M, rtol=0.0, atol=0.01
    )
    pixels_ok = manifest["cell_pixels"].astype(int) == CHIP_SIZE
    width_ok = np.isclose(
        manifest["cell_xmax"] - manifest["cell_xmin"],
        CHIP_SIZE_M,
        rtol=0.0,
        atol=0.01,
    )
    height_ok = np.isclose(
        manifest["cell_ymax"] - manifest["cell_ymin"],
        CHIP_SIZE_M,
        rtol=0.0,
        atol=0.01,
    )
    footprint_count_ok = manifest["n_footprints"].astype(int) >= MIN_FOOTPRINTS
    if not (size_ok & pixels_ok & width_ok & height_ok & footprint_count_ok).all():
        raise ValueError("Manifest contains an invalid 4 km density window")

    expected_crs = manifest["mgrs_tile_t"].map(
        lambda tile: f"EPSG:{epsg_from_mgrs_tile(tile)}"
    )
    observed_crs = manifest["window_crs"].astype(str).str.strip().str.upper()
    if not (expected_crs == observed_crs).all():
        raise ValueError("Manifest window_crs does not match its MGRS tile")

    manifest = manifest.sort_values(
        ["sif_year", "sif_doy", "mgrs_tile_t", "aggregation_id"]
    ).reset_index(drop=True)
    if MAX_CHIPS is not None:
        manifest = manifest.head(MAX_CHIPS).copy()

    selected_ids = set(manifest["aggregation_id"])
    assignments = assignments[
        assignments["aggregation_id"].isin(selected_ids)
    ].copy()
    assignments = assignments.sort_values(
        ["aggregation_id", "sif_row_id", "source_csv_row"]
    )

    if not np.isfinite(assignments[TARGET_COLUMN]).all():
        raise ValueError("Assignments contain non-finite SIF targets")
    if (
        assignments["window_mask_inside_fraction"]
        < MIN_ASSIGNED_FOOTPRINT_INSIDE_FRACTION - 1e-6
    ).any():
        raise ValueError(
            "Assignments contain a footprint below the 60% window-coverage filter"
        )

    manifest_lookup = manifest.set_index("aggregation_id", drop=False)
    groups: dict[str, pd.DataFrame] = {}
    for aggregation_id, group in assignments.groupby(
        "aggregation_id", sort=False
    ):
        group = group.reset_index(drop=True)
        manifest_row = manifest_lookup.loc[aggregation_id]

        homogeneous_columns = [
            "Delta_Date",
            "mgrs_tile_t",
            "measurement_mode",
        ]
        invalid = {
            column: int(group[column].nunique(dropna=False))
            for column in homogeneous_columns
            if group[column].nunique(dropna=False) != 1
        }
        if invalid:
            raise ValueError(
                f"{aggregation_id} has non-homogeneous assignments: {invalid}"
            )

        if group.iloc[0]["Delta_Date"] != manifest_row["Delta_Date"]:
            raise ValueError(f"Date mismatch for {aggregation_id}")
        if group.iloc[0]["mgrs_tile_t"] != manifest_row["mgrs_tile_t"]:
            raise ValueError(f"MGRS tile mismatch for {aggregation_id}")
        if int(group.iloc[0]["measurement_mode"]) != int(
            manifest_row["measurement_mode"]
        ):
            raise ValueError(f"Measurement-mode mismatch for {aggregation_id}")
        if len(group) != int(manifest_row["n_footprints"]):
            raise ValueError(
                f"Footprint count mismatch for {aggregation_id}: manifest="
                f"{int(manifest_row['n_footprints'])}, assignments={len(group)}"
            )

        assignment_target = float(group[TARGET_COLUMN].mean())
        manifest_target = float(manifest_row[AGGREGATE_TARGET_COLUMN])
        if not np.isclose(
            assignment_target, manifest_target, rtol=1e-6, atol=1e-7
        ):
            raise ValueError(
                f"Aggregate target mismatch for {aggregation_id}: "
                f"manifest={manifest_target}, assignments={assignment_target}"
            )
        groups[aggregation_id] = group

    missing = sorted(selected_ids - set(groups))
    if missing:
        raise ValueError(f"Manifest windows have no assignments: {missing[:5]}")

    print(
        f"Loaded {len(manifest):,} density windows and "
        f"{len(assignments):,} assigned SIF footprints"
    )
    return manifest, groups


def preflight_predictor_files(manifest: pd.DataFrame) -> None:
    """Check each required file once without reading its raster pixels."""
    missing_sentinel = []
    for value in manifest["sentinel2_raster_path"].drop_duplicates():
        path = resolve_project_path(value)
        if not path.is_file():
            missing_sentinel.append(str(path))
    if missing_sentinel:
        raise FileNotFoundError(
            "Missing Sentinel-2 rasters; first paths: "
            + ", ".join(missing_sentinel[:5])
        )

    unique_dates = sorted(set(manifest["Delta_Date"]))
    for delta_date in unique_dates:
        date_iso = pd.Timestamp(delta_date).strftime("%Y-%m-%d")
        predictors.par_paths(date_iso)

    fapar_keys = sorted(
        {
            (
                int(row.sif_year),
                predictors.fapar_composite_doy(int(row.sif_doy)),
            )
            for row in manifest[["sif_year", "sif_doy"]].itertuples(index=False)
        }
    )
    for year, composite_doy in fapar_keys:
        for tile in FAPAR_MODIS_TILES:
            predictors.fapar_path(year, composite_doy, tile)

    for year in sorted(set(manifest["sif_year"].astype(int))):
        predictors.crop_path_for_year(year)

    print(
        "Predictor file preflight passed for "
        f"{len(manifest):,} Sentinel rasters, {len(unique_dates):,} PAR dates, "
        f"{len(fapar_keys):,} FAPAR year/composites and "
        f"{manifest['sif_year'].nunique():,} crop-map years"
    )


# ---------------------------------------------------------------------------
# Sentinel-2 indices

def read_sentinel_indices(
    manifest_row: pd.Series,
) -> tuple[list[np.ndarray], object, object, Path]:
    aggregation_id = str(manifest_row["aggregation_id"])
    raster_path = resolve_project_path(manifest_row["sentinel2_raster_path"])
    expected_tile = normalize_mgrs_tile(manifest_row["mgrs_tile_t"])
    expected_epsg = epsg_from_mgrs_tile(expected_tile)
    expected_date = pd.Timestamp(manifest_row["Delta_Date"]).strftime("%Y-%m-%d")
    expected_bounds = window_bounds(manifest_row)

    with rasterio.open(raster_path) as src:
        if src.width != CHIP_SIZE or src.height != CHIP_SIZE:
            raise ValueError(
                f"{aggregation_id}: Sentinel raster is {src.width} x {src.height}, "
                f"expected {CHIP_SIZE} x {CHIP_SIZE}"
            )
        if src.crs is None or src.crs.to_epsg() != expected_epsg:
            raise ValueError(
                f"{aggregation_id}: Sentinel CRS {src.crs} does not match "
                f"EPSG:{expected_epsg}"
            )
        if not (
            math.isclose(abs(src.transform.a), CHIP_RES_M, abs_tol=1e-6)
            and math.isclose(abs(src.transform.e), CHIP_RES_M, abs_tol=1e-6)
            and math.isclose(src.transform.b, 0.0, abs_tol=1e-9)
            and math.isclose(src.transform.d, 0.0, abs_tol=1e-9)
        ):
            raise ValueError(
                f"{aggregation_id}: Sentinel raster is not an unrotated 20 m grid"
            )

        observed_bounds = (
            float(src.bounds.left),
            float(src.bounds.bottom),
            float(src.bounds.right),
            float(src.bounds.top),
        )
        if any(
            abs(observed - expected) > 0.01
            for observed, expected in zip(observed_bounds, expected_bounds)
        ):
            raise ValueError(
                f"{aggregation_id}: Sentinel bounds {observed_bounds} do not "
                f"match manifest bounds {expected_bounds}"
            )

        descriptions = tuple(src.descriptions)
        missing = [
            name for name in SENTINEL_REFLECTANCE_BANDS if name not in descriptions
        ]
        if missing:
            raise ValueError(
                f"{aggregation_id}: Sentinel raster lacks bands {missing}; "
                f"descriptions={descriptions}"
            )
        if len(set(descriptions)) != len(descriptions):
            raise ValueError(
                f"{aggregation_id}: Sentinel raster has duplicate band descriptions"
            )

        tags = src.tags()
        if tags.get("aggregation_id") != aggregation_id:
            raise ValueError(f"{aggregation_id}: Sentinel aggregation tag mismatch")
        if normalize_mgrs_tile(tags.get("mgrs_tile", "")) != expected_tile:
            raise ValueError(f"{aggregation_id}: Sentinel MGRS tag mismatch")
        if tags.get("target_date") != expected_date:
            raise ValueError(f"{aggregation_id}: Sentinel target-date tag mismatch")

        indexes = [descriptions.index(name) + 1 for name in SENTINEL_REFLECTANCE_BANDS]
        reflectance = src.read(indexes, masked=True).filled(np.nan).astype(np.float32)
        dst_transform = src.transform
        dst_crs = src.crs

    bands = {
        name: reflectance[index]
        for index, name in enumerate(SENTINEL_REFLECTANCE_BANDS)
    }
    for name, values in bands.items():
        require_nonempty_channel(name, values, aggregation_id)

    ndmi = predictors.normalized_difference(bands["B08"], bands["B11"])
    ndvi = predictors.normalized_difference(bands["B08"], bands["B04"])
    evi = predictors.enhanced_vegetation_index(
        bands["B08"], bands["B04"], bands["B02"]
    )
    nirv = bands["B08"] * ndvi
    ndre = predictors.normalized_difference(bands["B8A"], bands["B05"])

    indices = [ndmi, ndvi, evi, nirv, ndre]
    for name, values in zip(CHANNEL_NAMES[:5], indices):
        require_nonempty_channel(name, values, aggregation_id)
    return indices, dst_transform, dst_crs, raster_path


# ---------------------------------------------------------------------------
# Predictor and supervision assembly

def build_predictor_chip(
    manifest_row: pd.Series,
) -> tuple[np.ndarray, object, object, Path, int]:
    aggregation_id = str(manifest_row["aggregation_id"])
    spectral_indices, dst_transform, dst_crs, raster_path = read_sentinel_indices(
        manifest_row
    )
    bounds = window_bounds(manifest_row)
    year = int(manifest_row["sif_year"])
    month = int(manifest_row["sif_month"])
    sif_doy = int(manifest_row["sif_doy"])
    composite_doy = predictors.fapar_composite_doy(sif_doy)

    fapar = predictors.fapar_chip(
        year,
        composite_doy,
        bounds,
        dst_transform,
        dst_crs,
    )
    par = predictors.par_chip(
        manifest_row["Delta_Date"],
        bounds,
        dst_transform,
        dst_crs,
    )
    apar = fapar * par
    crop_channels = predictors.crop_fraction_chips(
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
        predictors.constant_channel(np.sin(month_angle)),
        predictors.constant_channel(np.cos(month_angle)),
    ]
    if len(channels) != len(CHANNEL_NAMES):
        raise RuntimeError(
            f"{aggregation_id}: built {len(channels)} channels, expected "
            f"{len(CHANNEL_NAMES)}"
        )

    float16_max = np.finfo(np.float16).max
    for name, channel in zip(CHANNEL_NAMES, channels):
        channel[~np.isfinite(channel)] = np.nan
        channel[np.abs(channel) > float16_max] = np.nan
        require_nonempty_channel(name, channel, aggregation_id)

    # NaNs are preserved for nan-aware normalization after the dataset split.
    x = np.stack(channels, axis=0).astype(np.float16)
    return x, dst_transform, dst_crs, raster_path, composite_doy


def build_aggregate_weight_map(
    footprints: pd.DataFrame,
    dst_transform,
    dst_crs,
) -> tuple[np.ndarray, list[dict], dict]:
    normalized_masks: list[np.ndarray] = []
    footprint_rows: list[dict] = []

    for slot, (_, footprint) in enumerate(footprints.iterrows()):
        geometry = predictors.build_sif_polygon_projected(footprint, dst_crs)
        mask = predictors.fractional_footprint_mask(geometry, dst_transform)
        mask_sum = float(mask.sum(dtype=np.float64))
        rasterized_inside_fraction = predictors.mask_inside_fraction(mask, geometry)

        if not np.isfinite(mask_sum) or mask_sum <= 0:
            raise ValueError(
                "Empty footprint mask for aggregation_id="
                f"{footprint['aggregation_id']}, sif_row_id="
                f"{int(footprint['sif_row_id'])}"
            )

        normalized_masks.append(mask / mask_sum)
        footprint_rows.append(
            {
                "slot": slot,
                "sif_row_id": int(footprint["sif_row_id"]),
                "source_csv_row": int(footprint["source_csv_row"]),
                "source_dataset": optional_value(
                    footprint, "source_dataset", ""
                ),
                "observed_sif": float(footprint[TARGET_COLUMN]),
                "assigned_window_inside_fraction": float(
                    footprint["window_mask_inside_fraction"]
                ),
                "rasterized_mask_inside_fraction": float(
                    rasterized_inside_fraction
                ),
                "mask_sum_pixels": mask_sum,
                "mask_area_km2_on_chip": (
                    mask_sum * CHIP_RES_M**2 / 1_000_000.0
                ),
                "land_cover": optional_value(footprint, "land_cover", ""),
                "land_cover_class": optional_value(
                    footprint, "land_cover_class", ""
                ),
                "BKR10_ID": optional_value(footprint, "BKR10_ID", ""),
                "BKR_NAME": optional_value(footprint, "BKR_NAME", ""),
            }
        )

    # Each footprint mask is normalized separately before averaging, so every
    # accepted footprint contributes equally regardless of its polygon area.
    aggregate_weight = np.mean(
        np.stack(normalized_masks, axis=0), axis=0
    ).astype(np.float32)
    weight_sum = float(aggregate_weight.sum(dtype=np.float64))
    if not np.isclose(weight_sum, 1.0, rtol=1e-6, atol=1e-6):
        raise ValueError(f"Aggregate weight map sums to {weight_sum}, expected 1")
    aggregate_weight /= aggregate_weight.sum(dtype=np.float64)

    positive = aggregate_weight > 0
    mask_inside_values = [
        row["rasterized_mask_inside_fraction"] for row in footprint_rows
    ]
    diagnostics = {
        "aggregate_weight_sum_float32": float(
            aggregate_weight.sum(dtype=np.float64)
        ),
        "aggregate_support_pixels": int(positive.sum()),
        "aggregate_support_fraction": float(positive.mean()),
        "effective_weighted_pixels": float(
            1.0 / np.square(aggregate_weight.astype(np.float64)).sum()
        ),
        "mean_rasterized_mask_inside_fraction": float(
            np.mean(mask_inside_values)
        ),
        "min_rasterized_mask_inside_fraction": float(
            np.min(mask_inside_values)
        ),
        "max_rasterized_mask_inside_fraction": float(
            np.max(mask_inside_values)
        ),
    }
    return aggregate_weight.astype(np.float16), footprint_rows, diagnostics


def build_sample(
    manifest_row: pd.Series,
    footprints: pd.DataFrame,
) -> tuple[dict, dict, list[dict]]:
    x, dst_transform, dst_crs, raster_path, fapar_doy = build_predictor_chip(
        manifest_row
    )
    weight_map, footprint_metadata, weight_diagnostics = (
        build_aggregate_weight_map(footprints, dst_transform, dst_crs)
    )

    assignment_target = float(footprints[TARGET_COLUMN].mean())
    manifest_target = float(manifest_row[AGGREGATE_TARGET_COLUMN])
    if not np.isclose(assignment_target, manifest_target, rtol=1e-6, atol=1e-7):
        raise ValueError(
            f"Target changed while building {manifest_row['aggregation_id']}"
        )

    metadata = {
        "aggregation_id": str(manifest_row["aggregation_id"]),
        "cell_id": str(manifest_row["cell_id"]),
        "source_dataset": optional_value(manifest_row, "source_dataset", ""),
        "Delta_Date": manifest_row["Delta_Date"],
        "sif_year": int(manifest_row["sif_year"]),
        "sif_month": int(manifest_row["sif_month"]),
        "sif_doy": int(manifest_row["sif_doy"]),
        "measurement_mode": int(manifest_row["measurement_mode"]),
        "mgrs_tile_t": str(manifest_row["mgrs_tile_t"]),
        "window_crs": str(manifest_row["window_crs"]),
        "sentinel2_raster_path": str(raster_path),
        "fapar_composite_doy": int(fapar_doy),
        "par_date": manifest_row["Delta_Date"],
        "n_footprints": int(len(footprints)),
        "y_aggregate": assignment_target,
        "target_median": predictors.to_float(
            manifest_row.get("median_target_modis_sif", np.nan)
        ),
        "target_min": predictors.to_float(
            manifest_row.get("min_target_modis_sif", np.nan)
        ),
        "target_max": predictors.to_float(
            manifest_row.get("max_target_modis_sif", np.nan)
        ),
        "target_sd": predictors.to_float(
            manifest_row.get("sd_target_modis_sif", np.nan)
        ),
        "target_se": predictors.to_float(
            manifest_row.get("se_target_modis_sif", np.nan)
        ),
        "majority_land_cover": optional_value(
            manifest_row, "majority_land_cover", ""
        ),
        "majority_land_cover_class": optional_value(
            manifest_row, "majority_land_cover_class", ""
        ),
        "land_cover_majority_fraction": predictors.to_float(
            manifest_row.get("land_cover_majority_fraction", np.nan)
        ),
        "majority_BKR10_ID": optional_value(
            manifest_row, "majority_BKR10_ID", ""
        ),
        "majority_BKR_NAME": optional_value(
            manifest_row, "majority_BKR_NAME", ""
        ),
        "states": optional_value(manifest_row, "states", ""),
        "hzs_values": optional_value(manifest_row, "hzs_values", ""),
        "cell_xmin": float(manifest_row["cell_xmin"]),
        "cell_ymin": float(manifest_row["cell_ymin"]),
        "cell_xmax": float(manifest_row["cell_xmax"]),
        "cell_ymax": float(manifest_row["cell_ymax"]),
        **weight_diagnostics,
    }

    for row in footprint_metadata:
        row["aggregation_id"] = str(manifest_row["aggregation_id"])

    sample = {
        "X": x,
        "aggregate_weight_map": weight_map,
        "y_aggregate": np.float32(assignment_target),
        "n_footprints": np.int16(len(footprints)),
    }
    return sample, metadata, footprint_metadata


# ---------------------------------------------------------------------------
# Optional multiprocessing

def make_payload(
    index: int,
    manifest_row: pd.Series,
    footprints: pd.DataFrame,
) -> tuple[int, dict, list[dict]]:
    return index, manifest_row.to_dict(), footprints.to_dict(orient="records")


def build_sample_from_payload(
    payload: tuple[int, dict, list[dict]],
) -> tuple[int, dict, dict, list[dict]]:
    index, manifest_dict, footprint_dicts = payload
    sample, metadata, footprint_metadata = build_sample(
        pd.Series(manifest_dict),
        pd.DataFrame(footprint_dicts),
    )
    return index, sample, metadata, footprint_metadata


# ---------------------------------------------------------------------------
# Shards and output tables

def check_output_directory() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    generated_names = [
        "chip_metadata.csv",
        "footprint_metadata.csv",
        "channel_names.csv",
        "target_names.csv",
        "dataset_config.json",
    ]
    existing = list(OUTPUT_DIR.glob("chips_*.npz"))
    existing.extend(
        OUTPUT_DIR / name
        for name in generated_names
        if (OUTPUT_DIR / name).exists()
    )
    if FAIL_IF_OUTPUT_EXISTS and existing:
        names = ", ".join(path.name for path in existing[:5])
        raise FileExistsError(
            f"Output directory already contains generated files ({names}). "
            "Use an empty output directory before running this script."
        )


def write_shard(
    shard_id: int,
    xs: list[np.ndarray],
    weights: list[np.ndarray],
    targets: list[np.float32],
    footprint_counts: list[np.int16],
    aggregation_ids: list[str],
) -> Path:
    path = OUTPUT_DIR / f"chips_{shard_id:05d}.npz"
    np.savez_compressed(
        path,
        X=np.stack(xs, axis=0).astype(np.float16),
        aggregate_weight_map=np.stack(weights, axis=0).astype(np.float16),
        y_aggregate=np.asarray(targets, dtype=np.float32),
        n_footprints=np.asarray(footprint_counts, dtype=np.int16),
        aggregation_id=np.asarray(aggregation_ids),
        channel_names=np.asarray(CHANNEL_NAMES),
        target_name=np.asarray([AGGREGATE_TARGET_COLUMN]),
    )
    return path


def prepare_chips() -> None:
    check_output_directory()
    manifest, assignment_groups = load_aggregation_tables()
    preflight_predictor_files(manifest)

    metadata_rows: list[dict] = []
    footprint_metadata_rows: list[dict] = []
    xs: list[np.ndarray] = []
    weights: list[np.ndarray] = []
    targets: list[np.float32] = []
    footprint_counts: list[np.int16] = []
    aggregation_ids: list[str] = []
    shard_id = 0

    def flush_shard() -> None:
        nonlocal shard_id, xs, weights, targets, footprint_counts, aggregation_ids
        if not xs:
            return
        written = write_shard(
            shard_id,
            xs,
            weights,
            targets,
            footprint_counts,
            aggregation_ids,
        )
        print(f"Wrote {written}")
        shard_id += 1
        xs = []
        weights = []
        targets = []
        footprint_counts = []
        aggregation_ids = []

    def handle_result(result: tuple[int, dict, dict, list[dict]]) -> None:
        index, sample, metadata, footprint_metadata = result
        if index % 10 == 0:
            print(f"Prepared window {index + 1:,} / {len(manifest):,}")

        metadata["shard_file"] = f"chips_{shard_id:05d}.npz"
        metadata["shard_index"] = len(xs)
        metadata_rows.append(metadata)

        for row in footprint_metadata:
            row["shard_file"] = metadata["shard_file"]
            row["shard_index"] = metadata["shard_index"]
            footprint_metadata_rows.append(row)

        xs.append(sample["X"])
        weights.append(sample["aggregate_weight_map"])
        targets.append(sample["y_aggregate"])
        footprint_counts.append(sample["n_footprints"])
        aggregation_ids.append(str(metadata["aggregation_id"]))
        if len(xs) == SHARD_SIZE:
            flush_shard()

    payloads = (
        make_payload(
            index,
            manifest_row,
            assignment_groups[str(manifest_row["aggregation_id"])],
        )
        for index, manifest_row in manifest.iterrows()
    )

    if N_WORKERS == 1:
        for payload in payloads:
            handle_result(build_sample_from_payload(payload))
    else:
        with ProcessPoolExecutor(max_workers=N_WORKERS) as executor:
            for result in executor.map(
                build_sample_from_payload,
                payloads,
                chunksize=1,
            ):
                handle_result(result)

    flush_shard()

    metadata = pd.DataFrame(metadata_rows)
    footprint_metadata = pd.DataFrame(footprint_metadata_rows)
    metadata.to_csv(OUTPUT_DIR / "chip_metadata.csv", index=False)
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
            "target_name": [AGGREGATE_TARGET_COLUMN],
            "assignment_target_name": [TARGET_COLUMN],
            "aggregation": ["equal arithmetic mean across assigned footprints"],
        }
    ).to_csv(OUTPUT_DIR / "target_names.csv", index=False)

    config = {
        "input_manifest": str(MANIFEST_PATH),
        "input_assignments": str(ASSIGNMENTS_PATH),
        "output_sample_count": int(len(metadata)),
        "chip_size_m": CHIP_SIZE_M,
        "chip_resolution_m": CHIP_RES_M,
        "chip_rows": CHIP_SIZE,
        "chip_cols": CHIP_SIZE,
        "minimum_footprints": MIN_FOOTPRINTS,
        "minimum_assigned_footprint_inside_fraction": (
            MIN_ASSIGNED_FOOTPRINT_INSIDE_FRACTION
        ),
        "minimum_sentinel_valid_fraction_before_fill": (
            MIN_SENTINEL_VALID_FRACTION
        ),
        "sentinel_reflectance_bands": list(SENTINEL_REFLECTANCE_BANDS),
        "sentinel_quality_bands_included": False,
        "sentinel_reflectance_units": "unitless BOA reflectance",
        "channels": CHANNEL_NAMES,
        "target": AGGREGATE_TARGET_COLUMN,
        "assignment_target": TARGET_COLUMN,
        "target_range_filter": "none; retain the complete accepted SIF range",
        "target_aggregation": "equal arithmetic mean across assigned footprints",
        "supervision_weight_formula": (
            "mean_i(mask_i / sum_pixels(mask_i)); model scalar prediction is "
            "sum_pixels(weight_map * predicted_sif_map)"
        ),
        "mask_oversample": MASK_OVERSAMPLE,
        "fapar_modis_tiles": list(FAPAR_MODIS_TILES),
        "fapar_temporal_matching": "containing 8-day composite interval",
        "par_temporal_matching": "exact SIF Delta_Date",
        "par_units": "W/m2",
        "par_accepted_qa_codes": list(predictors.PAR_ACCEPTED_QA_CODES),
        "par_resampling": "bilinear after native-grid QA masking",
        "predictor_storage_dtype": "float16",
        "weight_map_storage_dtype": "float16",
        "target_storage_dtype": "float32",
        "predictor_nan_policy": (
            "preserve; after splitting, compute nan-aware normalization from "
            "training data only, normalize, then replace normalized NaNs with zero"
        ),
        "dataset_split": "not assigned by this script",
        "location_channels_or_lat_lon_metadata": False,
        "ancillary_spatial_check": (
            "each generated FAPAR, PAR and crop channel must contain at least "
            "one finite pixel; no separate coverage scan"
        ),
        "channel_reader_module": (
            "model_data_prep/prepare_sentinel2_multisif_cnn_chips.py"
        ),
    }
    with (OUTPUT_DIR / "dataset_config.json").open(
        "w", encoding="ascii"
    ) as file:
        json.dump(config, file, indent=2)

    print(f"Done. Wrote {len(metadata):,} density-window chips to {OUTPUT_DIR}")
    print(
        "Training reminder: create train/validation/test partitions first. "
        "Calculate channel statistics from training data only, normalize X, "
        "replace normalized NaNs with zero, and renormalize each weight map in "
        "float32 before calculating the weighted SIF prediction."
    )


if __name__ == "__main__":
    prepare_chips()
