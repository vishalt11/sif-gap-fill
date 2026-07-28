"""Prepare one 2 km x 2 km pixel table for tabular SIF map inference.

The XGBoost and Random Forest models in
``pmeans-tabular-xgboost-randomforest.ipynb`` were fitted to polygon-mean
predictors. This script creates the analogous predictors on a 100 x 100 grid
at 20 m resolution around one observed OCO-2 footprint. The resulting table is
intended for an illustrative, pixel-wise application of those fitted models.

Candidate footprints are ranked by polygon-level active-crop coverage and
distance from the MGRS tile centre. A bounded candidate scan then selects the
2 km window containing the most model-ready crop pixels.

Outputs:
  - pixel_features_2km_20m.csv: 10,000 rows with model features and coordinates
  - pixel_grid_2km_20m.npz: the same features as 100 x 100 arrays
  - selected_sif_row.csv: metadata and original polygon-mean values
  - training_feature_ranges.csv: model-data ranges for extrapolation checks
  - sample_metadata.json: grid, CRS, date and source information

Non-crop pixels are retained for mapping context but ``predictable_pixel`` is
False for them. The tabular models should normally be applied only where
``predictable_pixel`` is True because ``active_crop_pct`` was defined over
crop pixels in the training data.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
import re
from typing import Iterable

import numpy as np
import pandas as pd
import rasterio
from rasterio.enums import Resampling
from rasterio.transform import Affine, from_origin
from rasterio.warp import (
    reproject,
    transform as warp_coordinates,
    transform_bounds,
)
from shapely.geometry import box

import prepare_sentinel2_multisif_cnn_chips as sentinel


# ---------------------------------------------------------------------------
# Configuration

INPUT_CSV = Path(
    "data/pmeans_model_data/"
    "model_data_pmeans_6tiles_with_modis_sif_fapar_par_cleaned.csv"
)
MGRS_REFERENCE_DIR = Path("data/temp_data/mgrs_tifs")
OUTPUT_DIR = Path("data/pmeans_model_data/tabular_2km_pixel_sample")

TARGET_COLUMN = "target_modis_sif"
FINAL_CHECK_COLUMN = "final_check_modis_sif"
TARGET_TILE = "T32UNA"

WINDOW_SIZE_M = 2000.0
PIXEL_SIZE_M = 20.0
GRID_SIZE = 100
TEN_METRE_GRID_SIZE = 200
MAX_CANDIDATES_TO_EVALUATE = 50
FULL_COVERAGE_PIXEL_COUNT = GRID_SIZE * GRID_SIZE

NON_CROP_LOW_NDVI_THRESHOLD = 0.45

BASE_MODEL_FEATURES = [
    "mean_ndvi",
    "mean_nirv",
    "mean_fapar",
    "mean_par",
    "apar",
    "active_crop_pct",
]

MODEL_FEATURES = [
    *BASE_MODEL_FEATURES,
    "month_sin",
    "month_cos",
]

CORNER_COLUMNS = [
    "Lat_corner1",
    "Lat_corner2",
    "Lat_corner3",
    "Lat_corner4",
    "Lon_corner1",
    "Lon_corner2",
    "Lon_corner3",
    "Lon_corner4",
]

SENTINEL_BAND_FLAGS = {
    "B4": "R1",
    "B8": "R1",
}

# Reuse the already-tested FAPAR, PAR and footprint helpers on this grid.
sentinel.CHIP_SIZE_M = WINDOW_SIZE_M
sentinel.CHIP_RES_M = PIXEL_SIZE_M
sentinel.CHIP_SIZE = GRID_SIZE
sentinel.MASK_OVERSAMPLE = 4


# ---------------------------------------------------------------------------
# Generic helpers

def require_columns(
    frame: pd.DataFrame,
    columns: Iterable[str],
    label: str,
) -> None:
    missing = [column for column in columns if column not in frame.columns]
    if missing:
        raise ValueError(f"{label} is missing required columns: {missing}")


def normalize_mgrs_tile(value: object) -> str:
    tile = str(value).strip().upper()
    if not tile.startswith("T"):
        tile = f"T{tile}"
    if re.fullmatch(r"T\d{2}[A-Z]{3}", tile) is None:
        raise ValueError(f"Invalid MGRS tile: {value}")
    return tile


def resolve_product_path(value: object) -> Path:
    path = Path(str(value).strip())
    candidates = [path]
    if not path.is_absolute():
        candidates.append(Path.cwd() / path)

    # Some downloaded products previously contained one duplicate inner
    # directory. Supporting it here makes row selection robust to either form.
    for candidate in list(candidates):
        candidates.append(candidate / candidate.name)

    for candidate in candidates:
        if candidate.is_dir():
            return candidate
    raise FileNotFoundError(f"Sentinel product directory not found: {value}")


def reference_tif_for_tile(
    tile: str,
    product_paths: Iterable[object] | None = None,
) -> Path:
    matches = sorted(
        MGRS_REFERENCE_DIR.rglob(f"*_T{tile[1:]}_*_FRC_B5.tif")
    )
    if matches:
        return matches[0]

    # The convenience reference directory does not necessarily contain every
    # tile represented by the polygon-means table. Any actual product B5 for
    # the same tile has the same grid extent and is therefore an equivalent
    # reference for ranking footprints by distance from the tile centre.
    if product_paths is not None:
        unique_paths = pd.Series(product_paths).dropna().astype(str).unique()
        for value in unique_paths:
            try:
                product_path = resolve_product_path(value)
                if sentinel.sentinel_product_tile(product_path) != tile:
                    continue
                return sentinel.sentinel_band_path(product_path, "B5")
            except (FileNotFoundError, ValueError):
                continue

    raise FileNotFoundError(
        f"No B5 reference or usable Sentinel product found for {tile}"
    )


def grid_bounds(
    transform: Affine,
    height: int,
    width: int,
) -> tuple[float, float, float, float]:
    xmin = transform.c
    ymax = transform.f
    xmax = xmin + width * transform.a
    ymin = ymax + height * transform.e
    return float(xmin), float(ymin), float(xmax), float(ymax)


def normalized_difference(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    denominator = a + b
    output = np.full(a.shape, np.nan, dtype=np.float32)
    valid = (
        np.isfinite(a)
        & np.isfinite(b)
        & (np.abs(denominator) > 1e-6)
    )
    output[valid] = (a[valid] - b[valid]) / denominator[valid]
    return output


def reproject_float_array(
    source: np.ndarray,
    source_transform: Affine,
    source_crs,
    destination_transform: Affine,
    destination_crs,
    destination_shape: tuple[int, int],
    resampling: Resampling,
) -> np.ndarray:
    source_values = np.where(
        np.isfinite(source),
        source,
        sentinel.WARP_NODATA,
    ).astype(np.float32)
    destination = np.full(destination_shape, np.nan, dtype=np.float32)
    reproject(
        source=source_values,
        destination=destination,
        src_transform=source_transform,
        src_crs=source_crs,
        src_nodata=sentinel.WARP_NODATA,
        dst_transform=destination_transform,
        dst_crs=destination_crs,
        dst_nodata=np.nan,
        resampling=resampling,
    )
    return destination


def valid_fraction(array: np.ndarray) -> float:
    return float(np.isfinite(array).mean())


# ---------------------------------------------------------------------------
# Row selection and grid construction

def load_candidate_rows() -> pd.DataFrame:
    frame = pd.read_csv(INPUT_CSV, low_memory=False)
    frame["source_csv_row"] = np.arange(len(frame), dtype=np.int64)

    require_columns(
        frame,
        [
            "Delta_Date",
            "mgrs_tile",
            "product_path",
            TARGET_COLUMN,
            *CORNER_COLUMNS,
            *BASE_MODEL_FEATURES,
        ],
        "Polygon-means input",
    )

    frame["Delta_Date"] = pd.to_datetime(
        frame["Delta_Date"],
        errors="coerce",
    )
    frame["mgrs_tile_t"] = frame["mgrs_tile"].map(normalize_mgrs_tile)
    frame = frame.loc[frame["mgrs_tile_t"].eq(TARGET_TILE)].copy()
    if frame.empty:
        raise ValueError(f"No input rows belong to TARGET_TILE={TARGET_TILE}")

    month_angle = 2.0 * np.pi * frame["Delta_Date"].dt.month / 12.0
    frame["month_sin"] = np.sin(month_angle)
    frame["month_cos"] = np.cos(month_angle)

    numeric_columns = [
        TARGET_COLUMN,
        *CORNER_COLUMNS,
        *BASE_MODEL_FEATURES,
    ]
    for column in numeric_columns:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    frame[numeric_columns] = frame[numeric_columns].replace(
        [np.inf, -np.inf],
        np.nan,
    )

    keep = frame["Delta_Date"].notna()
    keep &= frame[numeric_columns].notna().all(axis=1)
    if FINAL_CHECK_COLUMN in frame.columns:
        keep &= (
            frame[FINAL_CHECK_COLUMN]
            .astype(str)
            .str.strip()
            .str.lower()
            .eq("accept")
        )
    frame = frame.loc[keep].copy()

    # The mean of the four corners is the centroid for these parallelogram
    # footprints and is sufficient for ranking candidates by tile-centre
    # distance. The exact projected polygon centroid is used for the final grid.
    frame["candidate_lon"] = frame[
        [f"Lon_corner{i}" for i in range(1, 5)]
    ].mean(axis=1)
    frame["candidate_lat"] = frame[
        [f"Lat_corner{i}" for i in range(1, 5)]
    ].mean(axis=1)

    distance_parts = []
    for tile, group in frame.groupby("mgrs_tile_t", sort=False):
        reference_path = reference_tif_for_tile(
            tile,
            group["product_path"],
        )
        with rasterio.open(reference_path) as reference:
            if reference.crs is None:
                raise ValueError(f"Reference raster has no CRS: {reference_path}")
            centre_x = (reference.bounds.left + reference.bounds.right) / 2.0
            centre_y = (reference.bounds.bottom + reference.bounds.top) / 2.0
            projected_x, projected_y = warp_coordinates(
                "EPSG:4326",
                reference.crs,
                group["candidate_lon"].to_numpy(dtype=float).tolist(),
                group["candidate_lat"].to_numpy(dtype=float).tolist(),
            )

        ranked = group.copy()
        ranked["candidate_center_distance_m"] = np.hypot(
            np.asarray(projected_x) - centre_x,
            np.asarray(projected_y) - centre_y,
        )
        distance_parts.append(ranked)

    candidates = pd.concat(distance_parts, ignore_index=True)
    candidates = candidates.sort_values(
        [
            "active_crop_pct",
            "candidate_center_distance_m",
            "source_csv_row",
        ],
        ascending=[False, True, True],
        kind="stable",
    ).reset_index(drop=True)

    if candidates.empty:
        raise ValueError("No valid candidate rows remain after filtering")
    return candidates


def centred_grid_for_row(
    row: pd.Series,
    product_path: Path,
) -> dict:
    b5_path = sentinel.sentinel_band_path(product_path, "B5")
    with rasterio.open(b5_path) as reference:
        if reference.crs is None:
            raise ValueError(f"Sentinel B5 raster has no CRS: {b5_path}")
        if not (
            math.isclose(abs(reference.transform.a), PIXEL_SIZE_M, abs_tol=1e-6)
            and math.isclose(abs(reference.transform.e), PIXEL_SIZE_M, abs_tol=1e-6)
        ):
            raise ValueError(
                f"Expected 20 m B5 pixels, found {reference.transform.a}, "
                f"{reference.transform.e}: {b5_path}"
            )

        footprint = sentinel.build_sif_polygon_projected(row, reference.crs)
        if footprint.is_empty or footprint.area <= 0:
            raise ValueError("Selected footprint is empty after projection")
        centroid = footprint.centroid

        centre_col = (centroid.x - reference.transform.c) / reference.transform.a
        centre_row = (
            reference.transform.f - centroid.y
        ) / abs(reference.transform.e)
        col_offset = int(round(centre_col - GRID_SIZE / 2.0))
        row_offset = int(round(centre_row - GRID_SIZE / 2.0))

        if (
            col_offset < 0
            or row_offset < 0
            or col_offset + GRID_SIZE > reference.width
            or row_offset + GRID_SIZE > reference.height
        ):
            raise ValueError("Centred 2 km window falls outside Sentinel product")

        transform = reference.transform * Affine.translation(
            col_offset,
            row_offset,
        )
        bounds = grid_bounds(transform, GRID_SIZE, GRID_SIZE)

    footprint_mask = sentinel.fractional_footprint_mask(
        footprint,
        transform,
    )
    footprint_inside_fraction = float(
        footprint.intersection(box(*bounds)).area
        / footprint.area
    )

    return {
        "transform": transform,
        "bounds": bounds,
        "crs": reference.crs,
        "footprint": footprint,
        "footprint_mask": footprint_mask,
        "footprint_inside_fraction": footprint_inside_fraction,
        "centroid_x": float(centroid.x),
        "centroid_y": float(centroid.y),
        "col_offset": col_offset,
        "row_offset": row_offset,
    }


def validate_source_files(
    row: pd.Series,
    product_path: Path,
) -> None:
    for band_id, flag_resolution in SENTINEL_BAND_FLAGS.items():
        sentinel.sentinel_band_path(product_path, band_id)
        sentinel.sentinel_flag_path(product_path, flag_resolution)

    date = pd.Timestamp(row["Delta_Date"])
    composite_doy = sentinel.fapar_composite_doy(date.dayofyear)
    for tile in sentinel.MODIS_TILES:
        sentinel.fapar_path(date.year, composite_doy, tile)
    sentinel.par_paths(date.strftime("%Y-%m-%d"))
    sentinel.crop_path_for_year(date.year)


# ---------------------------------------------------------------------------
# Sentinel-2 and crop features

def read_scaled_land_band_to_grid(
    product_path: Path,
    band_id: str,
    bounds: tuple[float, float, float, float],
    destination_transform: Affine,
    destination_crs,
    destination_shape: tuple[int, int],
) -> np.ndarray:
    flag_resolution = SENTINEL_BAND_FLAGS[band_id]
    band_path = sentinel.sentinel_band_path(product_path, band_id)
    flag_path = sentinel.sentinel_flag_path(product_path, flag_resolution)

    with rasterio.open(band_path) as band_src, rasterio.open(flag_path) as flag_src:
        if band_src.crs is None:
            raise ValueError(f"Band raster has no CRS: {band_path}")
        if (
            band_src.crs != flag_src.crs
            or band_src.transform != flag_src.transform
            or band_src.width != flag_src.width
            or band_src.height != flag_src.height
        ):
            raise ValueError(
                f"Band and flag grids differ: {band_path.name}, {flag_path.name}"
            )

        source_bounds = transform_bounds(
            destination_crs,
            band_src.crs,
            *bounds,
            densify_pts=21,
        )
        window = sentinel.expanded_window_for_bounds(
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
        flag = flag_src.read(
            1,
            window=window,
            boundless=True,
            fill_value=0,
        )
        source_transform = band_src.window_transform(window)

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

        destination = np.full(destination_shape, np.nan, dtype=np.float32)
        reproject(
            source=source,
            destination=destination,
            src_transform=source_transform,
            src_crs=band_src.crs,
            src_nodata=sentinel.WARP_NODATA,
            dst_transform=destination_transform,
            dst_crs=destination_crs,
            dst_nodata=np.nan,
            resampling=Resampling.average,
        )
    return destination


def read_crop_source(
    year: int,
    bounds: tuple[float, float, float, float],
    destination_crs,
) -> tuple[np.ndarray, Affine, object]:
    crop_path = sentinel.crop_path_for_year(year)
    with rasterio.open(crop_path) as src:
        if src.crs is None:
            raise ValueError(f"Crop raster has no CRS: {crop_path}")
        source_bounds = transform_bounds(
            destination_crs,
            src.crs,
            *bounds,
            densify_pts=21,
        )
        window = sentinel.expanded_window_for_bounds(
            source_bounds,
            src.transform,
            pad_pixels=2,
        )
        values = src.read(
            1,
            window=window,
            boundless=True,
            fill_value=sentinel.NON_CROP_CODE,
        ).astype(np.int16)
        if src.nodata is not None:
            values[values == src.nodata] = sentinel.NON_CROP_CODE
        source_transform = src.window_transform(window)
        source_crs = src.crs
    return values, source_transform, source_crs


def crop_codes_on_grid(
    crop_values: np.ndarray,
    crop_transform: Affine,
    crop_crs,
    destination_transform: Affine,
    destination_crs,
    destination_shape: tuple[int, int],
) -> np.ndarray:
    destination = np.full(
        destination_shape,
        sentinel.NON_CROP_CODE,
        dtype=np.int16,
    )
    reproject(
        source=crop_values,
        destination=destination,
        src_transform=crop_transform,
        src_crs=crop_crs,
        src_nodata=None,
        dst_transform=destination_transform,
        dst_crs=destination_crs,
        dst_nodata=sentinel.NON_CROP_CODE,
        resampling=Resampling.nearest,
    )
    return destination


def crop_indicator_fraction(
    crop_values: np.ndarray,
    crop_transform: Affine,
    crop_crs,
    codes: Iterable[int],
    destination_transform: Affine,
    destination_crs,
) -> np.ndarray:
    source = np.isin(crop_values, list(codes)).astype(np.float32)
    destination = np.zeros((GRID_SIZE, GRID_SIZE), dtype=np.float32)
    reproject(
        source=source,
        destination=destination,
        src_transform=crop_transform,
        src_crs=crop_crs,
        dst_transform=destination_transform,
        dst_crs=destination_crs,
        dst_nodata=0.0,
        resampling=Resampling.average,
    )
    return np.clip(destination, 0.0, 1.0)


def build_local_features(
    row: pd.Series,
    product_path: Path,
    grid: dict,
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    transform_20m = grid["transform"]
    bounds = grid["bounds"]
    crs = grid["crs"]
    transform_10m = from_origin(
        bounds[0],
        bounds[3],
        10.0,
        10.0,
    )

    b4_10m = read_scaled_land_band_to_grid(
        product_path,
        "B4",
        bounds,
        transform_10m,
        crs,
        (TEN_METRE_GRID_SIZE, TEN_METRE_GRID_SIZE),
    )
    b8_10m = read_scaled_land_band_to_grid(
        product_path,
        "B8",
        bounds,
        transform_10m,
        crs,
        (TEN_METRE_GRID_SIZE, TEN_METRE_GRID_SIZE),
    )

    ndvi_10m = normalized_difference(b8_10m, b4_10m)
    nirv_10m = b8_10m * ndvi_10m

    date = pd.Timestamp(row["Delta_Date"])
    crop_values, crop_transform, crop_crs = read_crop_source(
        date.year,
        bounds,
        crs,
    )
    crop_10m = crop_codes_on_grid(
        crop_values,
        crop_transform,
        crop_crs,
        transform_10m,
        crs,
        (TEN_METRE_GRID_SIZE, TEN_METRE_GRID_SIZE),
    )

    non_crop_10m = ~np.isin(crop_10m, sentinel.KNOWN_CROP_CODES)
    contamination_10m = (
        non_crop_10m
        & np.isfinite(ndvi_10m)
        & (ndvi_10m <= NON_CROP_LOW_NDVI_THRESHOLD)
    )

    raw_ndvi_20m = reproject_float_array(
        ndvi_10m,
        transform_10m,
        crs,
        transform_20m,
        crs,
        (GRID_SIZE, GRID_SIZE),
        Resampling.average,
    )

    ndvi_masked_10m = ndvi_10m.copy()
    nirv_masked_10m = nirv_10m.copy()
    ndvi_masked_10m[contamination_10m] = np.nan
    nirv_masked_10m[contamination_10m] = np.nan

    ndvi_20m = reproject_float_array(
        ndvi_masked_10m,
        transform_10m,
        crs,
        transform_20m,
        crs,
        (GRID_SIZE, GRID_SIZE),
        Resampling.average,
    )
    nirv_20m = reproject_float_array(
        nirv_masked_10m,
        transform_10m,
        crs,
        transform_20m,
        crs,
        (GRID_SIZE, GRID_SIZE),
        Resampling.average,
    )

    contamination_fraction_20m = reproject_float_array(
        contamination_10m.astype(np.float32),
        transform_10m,
        crs,
        transform_20m,
        crs,
        (GRID_SIZE, GRID_SIZE),
        Resampling.average,
    )

    composite_doy = sentinel.fapar_composite_doy(date.dayofyear)
    fapar_20m = sentinel.fapar_chip(
        date.year,
        composite_doy,
        bounds,
        transform_20m,
        crs,
    )
    par_20m = sentinel.par_chip(
        date,
        bounds,
        transform_20m,
        crs,
    )
    apar_20m = fapar_20m * par_20m

    crop_fraction = crop_indicator_fraction(
        crop_values,
        crop_transform,
        crop_crs,
        sentinel.KNOWN_CROP_CODES,
        transform_20m,
        crs,
    )
    active_codes = sentinel.active_codes_for_month(date.month)
    active_crop_area_fraction = crop_indicator_fraction(
        crop_values,
        crop_transform,
        crop_crs,
        active_codes,
        transform_20m,
        crs,
    )
    active_crop_pct = np.full(
        (GRID_SIZE, GRID_SIZE),
        np.nan,
        dtype=np.float32,
    )
    np.divide(
        active_crop_area_fraction,
        crop_fraction,
        out=active_crop_pct,
        where=crop_fraction > 0,
    )
    active_crop_pct = np.clip(active_crop_pct, 0.0, 1.0)

    month_angle = 2.0 * np.pi * date.month / 12.0
    month_sin = np.full(
        (GRID_SIZE, GRID_SIZE),
        np.sin(month_angle),
        dtype=np.float32,
    )
    month_cos = np.full(
        (GRID_SIZE, GRID_SIZE),
        np.cos(month_angle),
        dtype=np.float32,
    )

    model_features = {
        "mean_ndvi": ndvi_20m.astype(np.float32),
        "mean_nirv": nirv_20m.astype(np.float32),
        "mean_fapar": fapar_20m.astype(np.float32),
        "mean_par": par_20m.astype(np.float32),
        "apar": apar_20m.astype(np.float32),
        "active_crop_pct": active_crop_pct.astype(np.float32),
        "month_sin": month_sin,
        "month_cos": month_cos,
    }
    auxiliary = {
        "raw_ndvi_20m": raw_ndvi_20m.astype(np.float32),
        "crop_fraction": crop_fraction.astype(np.float32),
        "active_crop_area_fraction": active_crop_area_fraction.astype(
            np.float32
        ),
        "contamination_fraction": contamination_fraction_20m.astype(
            np.float32
        ),
    }
    return model_features, auxiliary


# ---------------------------------------------------------------------------
# Output

def feature_ranges(
    training_frame: pd.DataFrame,
) -> pd.DataFrame:
    records = []
    for feature in MODEL_FEATURES:
        values = pd.to_numeric(
            training_frame[feature],
            errors="coerce",
        ).to_numpy(dtype=np.float64)
        values = values[np.isfinite(values)]
        records.append(
            {
                "feature": feature,
                "training_min": float(np.min(values)),
                "training_q01": float(np.quantile(values, 0.01)),
                "training_mean": float(np.mean(values)),
                "training_q99": float(np.quantile(values, 0.99)),
                "training_max": float(np.max(values)),
            }
        )
    return pd.DataFrame(records)


def pixel_coordinate_table(
    transform: Affine,
    crs,
) -> pd.DataFrame:
    rows, columns = np.indices((GRID_SIZE, GRID_SIZE))
    x = (
        transform.c
        + (columns.ravel() + 0.5) * transform.a
        + (rows.ravel() + 0.5) * transform.b
    )
    y = (
        transform.f
        + (columns.ravel() + 0.5) * transform.d
        + (rows.ravel() + 0.5) * transform.e
    )
    longitude, latitude = warp_coordinates(
        crs,
        "EPSG:4326",
        x.tolist(),
        y.tolist(),
    )
    return pd.DataFrame(
        {
            "pixel_row": rows.ravel().astype(np.int16),
            "pixel_col": columns.ravel().astype(np.int16),
            "x_center": x,
            "y_center": y,
            "longitude": longitude,
            "latitude": latitude,
        }
    )


def selected_row_output(
    row: pd.Series,
    product_path: Path,
    grid: dict,
) -> pd.DataFrame:
    fields = [
        "source_csv_row",
        "Delta_Date",
        "mgrs_tile",
        "mgrs_tile_t",
        "Latitude",
        "Longitude",
        TARGET_COLUMN,
        FINAL_CHECK_COLUMN,
        *CORNER_COLUMNS,
        *BASE_MODEL_FEATURES,
    ]
    output = {
        field: row[field]
        for field in fields
        if field in row.index
    }
    output.update(
        {
            "resolved_product_path": str(product_path),
            "tile_center_distance_m": float(
                row["candidate_center_distance_m"]
            ),
            "footprint_inside_2km_fraction": grid[
                "footprint_inside_fraction"
            ],
        }
    )
    return pd.DataFrame([output])


def prepare_best_coverage_candidate(
    candidates: pd.DataFrame,
) -> tuple[
    pd.Series,
    Path,
    dict,
    dict[str, np.ndarray],
    dict[str, np.ndarray],
    pd.DataFrame,
]:
    skipped = []
    diagnostics = []
    best_payload = None
    best_quality = None
    scan_limit = min(MAX_CANDIDATES_TO_EVALUATE, len(candidates))

    for scan_rank, (_, row) in enumerate(
        candidates.head(scan_limit).iterrows(),
        start=1,
    ):
        try:
            product_path = resolve_product_path(row["product_path"])
            product_tile = sentinel.sentinel_product_tile(product_path)
            if product_tile != row["mgrs_tile_t"]:
                raise ValueError(
                    f"Product tile {product_tile} differs from "
                    f"{row['mgrs_tile_t']}"
                )
            validate_source_files(row, product_path)
            grid = centred_grid_for_row(row, product_path)
            features, auxiliary = build_local_features(
                row,
                product_path,
                grid,
            )

            finite_crop_pixels = (
                auxiliary["crop_fraction"] > 0
            )
            for feature in MODEL_FEATURES:
                finite_crop_pixels &= np.isfinite(features[feature])
            n_predictable = int(finite_crop_pixels.sum())
            if n_predictable == 0:
                raise ValueError("No model-ready crop pixels in the 2 km grid")

            diagnostics.append(
                {
                    "scan_rank": scan_rank,
                    "source_csv_row": int(row["source_csv_row"]),
                    "status": "usable",
                    "n_predictable_pixels": n_predictable,
                    "predictable_fraction": (
                        n_predictable / FULL_COVERAGE_PIXEL_COUNT
                    ),
                    "active_crop_pct_polygon": float(
                        row["active_crop_pct"]
                    ),
                    "footprint_inside_2km_fraction": float(
                        grid["footprint_inside_fraction"]
                    ),
                    "tile_center_distance_m": float(
                        row["candidate_center_distance_m"]
                    ),
                    "reason": "",
                }
            )

            quality = (
                n_predictable,
                float(grid["footprint_inside_fraction"]),
                -float(row["candidate_center_distance_m"]),
            )
            if best_quality is None or quality > best_quality:
                best_quality = quality
                best_payload = (
                    row,
                    product_path,
                    grid,
                    features,
                    auxiliary,
                )
                print(
                    f"New best candidate {scan_rank}/{scan_limit}: "
                    f"source row {int(row['source_csv_row'])}, "
                    f"{n_predictable:,} / "
                    f"{FULL_COVERAGE_PIXEL_COUNT:,} predictable pixels"
                )

            if n_predictable == FULL_COVERAGE_PIXEL_COUNT:
                print("Found a full-coverage 10,000-pixel window.")
                break
        except (FileNotFoundError, ValueError, rasterio.errors.RasterioError) as exc:
            skipped.append(
                {
                    "source_csv_row": int(row["source_csv_row"]),
                    "reason": str(exc),
                }
            )
            diagnostics.append(
                {
                    "scan_rank": scan_rank,
                    "source_csv_row": int(row["source_csv_row"]),
                    "status": "skipped",
                    "n_predictable_pixels": 0,
                    "predictable_fraction": 0.0,
                    "active_crop_pct_polygon": float(
                        row["active_crop_pct"]
                    ),
                    "footprint_inside_2km_fraction": np.nan,
                    "tile_center_distance_m": float(
                        row["candidate_center_distance_m"]
                    ),
                    "reason": str(exc),
                }
            )
            if len(skipped) <= 5:
                print(
                    f"Skipping source row {row['source_csv_row']}: {exc}"
                )

    if best_payload is None:
        raise RuntimeError(
            "No candidate in the bounded scan could be prepared. "
            f"First skipped rows: {skipped[:5]}"
        )

    diagnostics_frame = pd.DataFrame(diagnostics).sort_values(
        [
            "n_predictable_pixels",
            "footprint_inside_2km_fraction",
            "tile_center_distance_m",
        ],
        ascending=[False, False, True],
        kind="stable",
        na_position="last",
    )
    diagnostics_frame["selected"] = False
    selected_source_row = int(best_payload[0]["source_csv_row"])
    diagnostics_frame.loc[
        diagnostics_frame["source_csv_row"].eq(selected_source_row),
        "selected",
    ] = True

    best_count = int(best_quality[0])
    if best_count < FULL_COVERAGE_PIXEL_COUNT:
        print(
            "No complete 10,000-pixel crop window was found among "
            f"{scan_limit} candidates; using the best window with "
            f"{best_count:,} predictable pixels."
        )

    return (*best_payload, diagnostics_frame)


def main() -> None:
    if not INPUT_CSV.exists():
        raise FileNotFoundError(INPUT_CSV)
    if not MGRS_REFERENCE_DIR.exists():
        raise FileNotFoundError(MGRS_REFERENCE_DIR)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    candidates = load_candidate_rows()
    training_ranges = feature_ranges(candidates)

    row, product_path, grid, features, auxiliary, candidate_diagnostics = (
        prepare_best_coverage_candidate(candidates)
    )

    pixel_table = pixel_coordinate_table(grid["transform"], grid["crs"])
    for feature in MODEL_FEATURES:
        pixel_table[feature] = features[feature].ravel()
    for name, values in auxiliary.items():
        pixel_table[name] = values.ravel()
    pixel_table["footprint_fraction"] = grid["footprint_mask"].ravel()

    feature_matrix = pixel_table[MODEL_FEATURES].to_numpy(dtype=np.float64)
    finite_features = np.isfinite(feature_matrix).all(axis=1)
    is_crop_pixel = pixel_table["crop_fraction"].to_numpy() > 0
    pixel_table["is_crop_pixel"] = is_crop_pixel
    pixel_table["predictable_pixel"] = finite_features & is_crop_pixel

    within_training_range = np.ones(len(pixel_table), dtype=bool)
    range_lookup = training_ranges.set_index("feature")
    for feature in MODEL_FEATURES:
        values = pixel_table[feature].to_numpy(dtype=np.float64)
        within_training_range &= (
            np.isfinite(values)
            & (values >= range_lookup.loc[feature, "training_min"])
            & (values <= range_lookup.loc[feature, "training_max"])
        )
    pixel_table["within_training_minmax"] = within_training_range
    pixel_table["recommended_prediction_pixel"] = (
        pixel_table["predictable_pixel"]
        & pixel_table["within_training_minmax"]
    )

    selected_row = selected_row_output(
        row,
        product_path,
        grid,
    )

    pixel_csv = OUTPUT_DIR / "pixel_features_2km_20m.csv"
    selected_csv = OUTPUT_DIR / "selected_sif_row.csv"
    ranges_csv = OUTPUT_DIR / "training_feature_ranges.csv"
    diagnostics_csv = OUTPUT_DIR / "candidate_scan_diagnostics.csv"
    npz_path = OUTPUT_DIR / "pixel_grid_2km_20m.npz"
    metadata_path = OUTPUT_DIR / "sample_metadata.json"

    pixel_table.to_csv(pixel_csv, index=False)
    selected_row.to_csv(selected_csv, index=False)
    training_ranges.to_csv(ranges_csv, index=False)
    candidate_diagnostics.to_csv(diagnostics_csv, index=False)

    feature_cube = np.stack(
        [features[name] for name in MODEL_FEATURES],
        axis=0,
    ).astype(np.float32)
    np.savez_compressed(
        npz_path,
        feature_cube=feature_cube,
        feature_names=np.asarray(MODEL_FEATURES),
        raw_ndvi_20m=auxiliary["raw_ndvi_20m"],
        crop_fraction=auxiliary["crop_fraction"],
        active_crop_area_fraction=auxiliary[
            "active_crop_area_fraction"
        ],
        contamination_fraction=auxiliary["contamination_fraction"],
        footprint_mask=grid["footprint_mask"].astype(np.float32),
        predictable_mask=(
            pixel_table["predictable_pixel"]
            .to_numpy()
            .reshape(GRID_SIZE, GRID_SIZE)
        ),
        within_training_minmax_mask=(
            pixel_table["within_training_minmax"]
            .to_numpy()
            .reshape(GRID_SIZE, GRID_SIZE)
        ),
        recommended_prediction_mask=(
            pixel_table["recommended_prediction_pixel"]
            .to_numpy()
            .reshape(GRID_SIZE, GRID_SIZE)
        ),
    )

    date = pd.Timestamp(row["Delta_Date"])
    metadata = {
        "source_csv_row": int(row["source_csv_row"]),
        "observed_target_modis_sif": float(row[TARGET_COLUMN]),
        "delta_date": date.strftime("%Y-%m-%d"),
        "month": int(date.month),
        "mgrs_tile_t": row["mgrs_tile_t"],
        "product_path": str(product_path),
        "grid_width": GRID_SIZE,
        "grid_height": GRID_SIZE,
        "pixel_size_m": PIXEL_SIZE_M,
        "window_size_m": WINDOW_SIZE_M,
        "bounds": list(grid["bounds"]),
        "affine_transform": list(grid["transform"])[:6],
        "crs_wkt": grid["crs"].to_wkt(),
        "footprint_inside_2km_fraction": grid[
            "footprint_inside_fraction"
        ],
        "n_crop_pixels": int(pixel_table["is_crop_pixel"].sum()),
        "n_predictable_pixels": int(
            pixel_table["predictable_pixel"].sum()
        ),
        "n_recommended_prediction_pixels": int(
            pixel_table["recommended_prediction_pixel"].sum()
        ),
        "candidate_scan_limit": MAX_CANDIDATES_TO_EVALUATE,
        "n_candidates_evaluated": int(len(candidate_diagnostics)),
        "feature_valid_fractions": {
            name: valid_fraction(values)
            for name, values in features.items()
        },
        "model_feature_order": MODEL_FEATURES,
        "notes": [
            "The tabular models were trained on footprint means; pixel-wise "
            "predictions are an illustrative spatial application.",
            "Non-crop pixels are not recommended for model inference.",
            "FAPAR and PAR are resampled to 20 m but retain their coarser "
            "native information content.",
        ],
    }
    metadata_path.write_text(
        json.dumps(metadata, indent=2),
        encoding="utf-8",
    )

    print(f"Selected source CSV row: {row['source_csv_row']}")
    print(f"Selected tile/date: {row['mgrs_tile_t']} / {date.date()}")
    print(
        "Distance from tile centre: "
        f"{row['candidate_center_distance_m'] / 1000.0:.2f} km"
    )
    print(
        "Footprint inside 2 km window: "
        f"{grid['footprint_inside_fraction']:.1%}"
    )
    print(
        "Predictable crop pixels: "
        f"{pixel_table['predictable_pixel'].sum():,} / {len(pixel_table):,}"
    )
    print(f"Saved: {pixel_csv}")
    print(f"Saved: {npz_path}")
    print(f"Saved: {selected_csv}")
    print(f"Saved: {ranges_csv}")
    print(f"Saved: {diagnostics_csv}")
    print(f"Saved: {metadata_path}")


if __name__ == "__main__":
    main()
