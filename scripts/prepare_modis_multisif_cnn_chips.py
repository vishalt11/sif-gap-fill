"""
Prepare multi-footprint MODIS-scale CNN chips from density-cluster SIF groups.

Inputs are produced by diagnose_multisif_density_cluster_chips.R:
  - one row per 24x24 chip in multi_sif_24px_5to10_chip_manifest.csv
  - one row per assigned SIF footprint in multi_sif_24px_5to10_chip_assignments.csv

Each output sample is one chip with multiple SIF supervision points:
  X:               [channels, 24, 24]
  footprint_masks: [max_footprints, 24, 24]
  y_targets:       [max_footprints, 3]
  target_accept:   [max_footprints, 3]
  footprint_valid: [max_footprints]

Training should average loss per chip, then average over the batch, so chips with
more footprints do not dominate chips with fewer footprints.
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from functools import lru_cache
import glob
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import rasterio
from rasterio.enums import Resampling
from rasterio.features import rasterize
from rasterio.transform import Affine, from_origin
from rasterio.warp import reproject, transform_geom
from rasterio.windows import Window, from_bounds
from shapely.geometry import MultiPoint, Polygon, box, mapping, shape


# ---------------------------------------------------------------------------
# Config

CHIP_MANIFEST_PATH = Path(
    "data/cnn_multisif_chip_diagnostics/multi_sif_24px_5to10_chip_manifest.csv"
)
CHIP_ASSIGNMENTS_PATH = Path(
    "data/cnn_multisif_chip_diagnostics/multi_sif_24px_5to10_chip_assignments.csv"
)

GEOTIFF_ROOT = Path("data/glass_geotiff")
FAPAR_DIR = GEOTIFF_ROOT / "fapar"
EVI_DIR = GEOTIFF_ROOT / "evi"
NDVI_DIR = GEOTIFF_ROOT / "ndvi"
PAR_DIR = GEOTIFF_ROOT / "par"
CROP_DIR = Path("data/crop_type_tif")

OUTPUT_DIR = Path("data/cnn_modis_chips/multisif_24x24_5to10_active_crop")

MODIS_TILES = ("h18v03", "h18v04")

TARGET_COLUMNS = ["Daily_SIF_757nm", "Daily_SIF_771nm", "target_modis_sif"]
FINAL_CHECK_COLUMNS = ["final_check_757", "final_check_771", "final_check_modis_sif"]

CHIP_SIZE = 24
CHIP_RES_M = 250.0
MASK_OVERSAMPLE = 8
MAX_FOOTPRINTS = 10
MIN_VALID_FOOTPRINTS = 5
MIN_FOOTPRINT_INSIDE_FRACTION = 0.90
SHARD_SIZE = 512

# Set to a small integer for a smoke test.
MAX_CHIPS: int | None = None
#MAX_CHIPS = 100
#N_WORKERS = 1
# Use 1 for debugging. 2-4 is usually a reasonable Windows starting point.
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
    "fapar",
    "evi",
    "ndvi",
    "par",
    "apar",
    *CROP_GROUPS.keys(),
    "active_crop_fraction",
    "non_crop_fraction",
    "month_sin",
    "month_cos",
]


# ---------------------------------------------------------------------------
# Raster containers and file lookup


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


def find_one_file(pattern: str) -> Path:
    matches = sorted(Path(match) for match in glob.glob(pattern))
    if len(matches) != 1:
        raise FileNotFoundError(f"Expected 1 file for pattern {pattern}, found {len(matches)}")
    return matches[0]


@lru_cache(maxsize=512)
def load_modis_tile(product: str, year: int, doy: int, tile: str) -> RasterSource:
    cfg = PRODUCTS[product]
    pattern = str(cfg["directory"] / tile / str(year) / f"*A{year}{doy:03d}.{tile}.*.tif")
    return read_geotiff_source(find_one_file(pattern))


@lru_cache(maxsize=256)
def load_par(year: int, doy: int) -> RasterSource:
    pattern = str(PAR_DIR / str(year) / f"*A{year}{doy:03d}.*.tif")
    return read_geotiff_source(find_one_file(pattern))


# ---------------------------------------------------------------------------
# Manifest loading


def require_columns(df: pd.DataFrame, columns: list[str], label: str) -> None:
    missing = [column for column in columns if column not in df.columns]
    if missing:
        raise ValueError(f"{label} is missing required columns: {missing}")


def load_chip_tables() -> tuple[pd.DataFrame, dict[str, pd.DataFrame]]:
    manifest = pd.read_csv(CHIP_MANIFEST_PATH)
    assignments = pd.read_csv(CHIP_ASSIGNMENTS_PATH)

    require_columns(
        manifest,
        [
            "chip_id",
            "Delta_Date",
            "sif_year",
            "sif_doy",
            "composite_doy",
            "par_doy",
            "chip_pixels",
            "chip_xmin",
            "chip_ymin",
            "chip_xmax",
            "chip_ymax",
            "n_sif",
        ],
        "chip manifest",
    )
    require_columns(
        assignments,
        [
            "chip_id",
            "sif_row_id",
            "Delta_Date",
            "sif_year",
            "sif_doy",
            "composite_doy",
            "par_doy",
            "Lat_corner1",
            "Lat_corner2",
            "Lat_corner3",
            "Lat_corner4",
            "Lon_corner1",
            "Lon_corner2",
            "Lon_corner3",
            "Lon_corner4",
            *TARGET_COLUMNS,
            *FINAL_CHECK_COLUMNS,
        ],
        "chip assignments",
    )

    manifest = manifest.copy()
    assignments = assignments.copy()

    manifest["Delta_Date"] = pd.to_datetime(manifest["Delta_Date"], errors="coerce").dt.date
    assignments["Delta_Date"] = pd.to_datetime(assignments["Delta_Date"], errors="coerce").dt.date

    numeric_manifest_cols = [
        "sif_year",
        "sif_doy",
        "composite_doy",
        "par_doy",
        "chip_pixels",
        "chip_xmin",
        "chip_ymin",
        "chip_xmax",
        "chip_ymax",
        "n_sif",
    ]
    for column in numeric_manifest_cols:
        manifest[column] = pd.to_numeric(manifest[column], errors="coerce")

    numeric_assignment_cols = [
        "sif_row_id",
        "sif_year",
        "sif_doy",
        "composite_doy",
        "par_doy",
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
        *TARGET_COLUMNS,
    ]
    for column in numeric_assignment_cols:
        if column in assignments.columns:
            assignments[column] = pd.to_numeric(assignments[column], errors="coerce")

    manifest = manifest[manifest["chip_pixels"].astype(int) == CHIP_SIZE].copy()
    manifest = manifest.sort_values(["sif_year", "composite_doy", "Delta_Date", "chip_id"])

    if MAX_CHIPS is not None:
        manifest = manifest.head(MAX_CHIPS).copy()

    assignments = assignments[assignments["chip_id"].isin(manifest["chip_id"])].copy()

    if "track_score" in assignments.columns:
        assignments = assignments.sort_values(["chip_id", "track_score", "sif_row_id"])
    else:
        assignments = assignments.sort_values(["chip_id", "sif_row_id"])

    assignment_groups = {
        chip_id: group.reset_index(drop=True)
        for chip_id, group in assignments.groupby("chip_id", sort=False)
    }

    print(
        f"Loaded {len(manifest):,} candidate chips and "
        f"{len(assignments):,} assigned SIF footprints"
    )
    return manifest.reset_index(drop=True), assignment_groups


# ---------------------------------------------------------------------------
# Chip generation


def clean_product_array(arr: np.ndarray, product: str) -> np.ndarray:
    cfg = PRODUCTS[product]
    out = arr.astype(np.float32, copy=True)
    out[(out < cfg["valid_min"]) | (out > cfg["valid_max"])] = np.nan
    return out


def clean_par_array(arr: np.ndarray) -> np.ndarray:
    out = arr.astype(np.float32, copy=True)
    out[out < 0] = np.nan
    return out


def make_chip_transform(chip_row: pd.Series) -> Affine:
    return from_origin(
        float(chip_row["chip_xmin"]),
        float(chip_row["chip_ymax"]),
        CHIP_RES_M,
        CHIP_RES_M,
    )


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


def build_sif_polygon_modis(row: pd.Series):
    polygon = build_sif_polygon_wgs84(row)
    geom_json = transform_geom(
        "EPSG:4326",
        MODIS_SINUSOIDAL_CRS,
        mapping(polygon),
        precision=6,
    )
    return shape(geom_json)


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


def mask_inside_fraction(mask: np.ndarray, geom_modis) -> float:
    expected_mask_sum = float(geom_modis.area) / (CHIP_RES_M * CHIP_RES_M)
    if expected_mask_sum <= 0:
        return 0.0
    return float(min(mask.sum() / expected_mask_sum, 1.0))


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
    group_codes = list(group_codes)
    if len(group_codes) == 0:
        return np.zeros((CHIP_SIZE, CHIP_SIZE), dtype=np.float32)

    group_mask = np.isin(crop_values, group_codes).astype(np.float32)
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


def active_codes_for_month(month: int) -> list[int]:
    return [
        code
        for code, active_months in ACTIVE_GROWTH_MONTHS_BY_CODE.items()
        if month in active_months
    ]


def crop_fractions_chip(year: int, month: int, dst_transform: Affine) -> list[np.ndarray]:
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

        active_crop = crop_group_fraction_chip(
            crop_values,
            crop_transform,
            src.crs,
            active_codes_for_month(month),
            dst_transform,
        )

        non_crop_values = (~np.isin(crop_values, KNOWN_CROP_CODES)).astype(np.float32)
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

    crop_channels.append(np.clip(active_crop, 0.0, 1.0).astype(np.float32))
    crop_channels.append(np.clip(non_crop, 0.0, 1.0).astype(np.float32))
    return crop_channels


def constant_channel(value: float) -> np.ndarray:
    return np.full((CHIP_SIZE, CHIP_SIZE), value, dtype=np.float32)


def to_float(value) -> float:
    if pd.isna(value):
        return np.nan
    return float(value)


def accept_to_bool(value) -> bool:
    return str(value).strip().lower() == "accept"


def build_predictor_chip(chip_row: pd.Series, dst_transform: Affine) -> tuple[np.ndarray, dict]:
    year = int(chip_row["sif_year"])
    doy = int(chip_row["composite_doy"])
    par_doy = int(chip_row["par_doy"])
    month = pd.to_datetime(chip_row["Delta_Date"]).month

    fapar = mosaic_modis_product_chip("fapar", year, doy, dst_transform)
    evi = mosaic_modis_product_chip("evi", year, doy, dst_transform)
    ndvi = mosaic_modis_product_chip("ndvi", year, doy, dst_transform)
    par = par_chip(year, par_doy, dst_transform)
    apar = fapar * par

    crop_channels = crop_fractions_chip(year, month, dst_transform)

    month_angle = 2.0 * np.pi * month / 12.0
    month_sin = constant_channel(np.sin(month_angle))
    month_cos = constant_channel(np.cos(month_angle))

    channels = [
        fapar,
        evi,
        ndvi,
        par,
        apar,
        *crop_channels,
        month_sin,
        month_cos,
    ]

    x = np.stack(channels, axis=0).astype(np.float32)

    valid_fraction = {
        f"{name}_valid_fraction": float(np.isfinite(channel).mean())
        for name, channel in zip(CHANNEL_NAMES, channels)
        if name not in {"month_sin", "month_cos"}
    }

    x = np.nan_to_num(x, nan=0.0, posinf=0.0, neginf=0.0)
    return x, valid_fraction


def build_multisif_masks(footprints: pd.DataFrame, dst_transform: Affine) -> tuple[dict, list[dict]]:
    masks = np.zeros((MAX_FOOTPRINTS, CHIP_SIZE, CHIP_SIZE), dtype=np.float32)
    y_targets = np.full((MAX_FOOTPRINTS, len(TARGET_COLUMNS)), np.nan, dtype=np.float32)
    target_accept = np.zeros((MAX_FOOTPRINTS, len(TARGET_COLUMNS)), dtype=np.uint8)
    footprint_valid = np.zeros(MAX_FOOTPRINTS, dtype=np.uint8)
    sif_row_ids = np.full(MAX_FOOTPRINTS, -1, dtype=np.int64)
    mask_sums = np.zeros(MAX_FOOTPRINTS, dtype=np.float32)
    inside_fractions = np.zeros(MAX_FOOTPRINTS, dtype=np.float32)

    kept_rows: list[dict] = []
    skipped_rows: list[dict] = []

    if footprints.empty:
        payload = {
            "footprint_masks": masks,
            "y_targets": y_targets,
            "target_accept": target_accept,
            "footprint_valid": footprint_valid,
            "sif_row_ids": sif_row_ids,
            "mask_sums": mask_sums,
            "mask_inside_fractions": inside_fractions,
            "n_valid_footprints": 0,
            "n_skipped_footprints": 0,
        }
        return payload, []

    if "track_score" in footprints.columns:
        footprints = footprints.sort_values(["track_score", "sif_row_id"])
    else:
        footprints = footprints.sort_values(["sif_row_id"])

    for _, footprint in footprints.iterrows():
        geom_modis = build_sif_polygon_modis(footprint)
        mask = fractional_footprint_mask(geom_modis, dst_transform)
        mask_sum = float(mask.sum())
        inside_fraction = mask_inside_fraction(mask, geom_modis)

        skip_reason = None
        if mask_sum <= 0:
            skip_reason = "empty_mask"
        elif inside_fraction < MIN_FOOTPRINT_INSIDE_FRACTION:
            skip_reason = "clipped_footprint"

        footprint_meta = {
            "sif_row_id": int(footprint["sif_row_id"]),
            "mask_sum": mask_sum,
            "mask_inside_fraction": inside_fraction,
            "Daily_SIF_757nm": to_float(footprint.get("Daily_SIF_757nm", np.nan)),
            "Daily_SIF_771nm": to_float(footprint.get("Daily_SIF_771nm", np.nan)),
            "target_modis_sif": to_float(footprint.get("target_modis_sif", np.nan)),
            "final_check_757": footprint.get("final_check_757", ""),
            "final_check_771": footprint.get("final_check_771", ""),
            "final_check_modis_sif": footprint.get("final_check_modis_sif", ""),
            "state": footprint.get("state", ""),
            "hzs": footprint.get("hzs", ""),
        }

        if skip_reason is not None:
            skipped_rows.append({**footprint_meta, "slot": -1, "skip_reason": skip_reason})
            continue

        slot = len(kept_rows)
        masks[slot] = mask
        y_targets[slot] = [to_float(footprint[column]) for column in TARGET_COLUMNS]
        target_accept[slot] = [accept_to_bool(footprint[column]) for column in FINAL_CHECK_COLUMNS]
        footprint_valid[slot] = 1
        sif_row_ids[slot] = int(footprint["sif_row_id"])
        mask_sums[slot] = mask_sum
        inside_fractions[slot] = inside_fraction
        kept_rows.append({**footprint_meta, "slot": slot, "skip_reason": ""})

        if len(kept_rows) == MAX_FOOTPRINTS:
            break

    payload = {
        "footprint_masks": masks,
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


def build_chip(chip_row: pd.Series, footprints: pd.DataFrame) -> tuple[dict, dict, list[dict]]:
    dst_transform = make_chip_transform(chip_row)
    x, valid_fraction = build_predictor_chip(chip_row, dst_transform)
    mask_payload, footprint_metadata = build_multisif_masks(footprints, dst_transform)

    chip_metadata = {
        "chip_id": chip_row["chip_id"],
        "Delta_Date": chip_row["Delta_Date"],
        "sif_year": int(chip_row["sif_year"]),
        "sif_doy": int(chip_row["sif_doy"]),
        "composite_doy": int(chip_row["composite_doy"]),
        "par_doy": int(chip_row["par_doy"]),
        "chip_pixels": int(chip_row["chip_pixels"]),
        "chip_size_m": float(chip_row.get("chip_size_m", CHIP_SIZE * CHIP_RES_M)),
        "chip_xmin": float(chip_row["chip_xmin"]),
        "chip_ymin": float(chip_row["chip_ymin"]),
        "chip_xmax": float(chip_row["chip_xmax"]),
        "chip_ymax": float(chip_row["chip_ymax"]),
        "n_sif_assigned": int(len(footprints)),
        "n_valid_footprints": mask_payload["n_valid_footprints"],
        "n_skipped_footprints": mask_payload["n_skipped_footprints"],
        "states": chip_row.get("states", ""),
        "hzs_values": chip_row.get("hzs_values", ""),
        "mean_target_modis_sif": to_float(chip_row.get("mean_target_modis_sif", np.nan)),
        "min_target_modis_sif": to_float(chip_row.get("min_target_modis_sif", np.nan)),
        "max_target_modis_sif": to_float(chip_row.get("max_target_modis_sif", np.nan)),
        **valid_fraction,
    }

    for target_i, target_name in enumerate(TARGET_COLUMNS):
        chip_metadata[f"{target_name}_accepted_footprints"] = int(
            (mask_payload["target_accept"][:, target_i] * mask_payload["footprint_valid"]).sum()
        )

    chip_payload = {
        "X": x,
        **mask_payload,
    }
    return chip_payload, chip_metadata, footprint_metadata


# ---------------------------------------------------------------------------
# Optional multiprocessing helpers


def make_payload(idx: int, chip_row: pd.Series, footprints: pd.DataFrame) -> tuple[int, dict, list[dict]]:
    return idx, chip_row.to_dict(), footprints.to_dict(orient="records")


def build_chip_from_payload(payload: tuple[int, dict, list[dict]]) -> tuple[int, dict, dict, list[dict]]:
    idx, chip_dict, footprint_dicts = payload
    chip_row = pd.Series(chip_dict)
    footprints = pd.DataFrame(footprint_dicts)
    chip_payload, chip_metadata, footprint_metadata = build_chip(chip_row, footprints)
    return idx, chip_payload, chip_metadata, footprint_metadata


# ---------------------------------------------------------------------------
# Writing shards


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
        X=np.stack(xs, axis=0).astype(np.float32),
        footprint_masks=np.stack(masks, axis=0).astype(np.float32),
        y_targets=np.stack(y_targets, axis=0).astype(np.float32),
        target_accept=np.stack(target_accepts, axis=0).astype(np.uint8),
        footprint_valid=np.stack(footprint_valids, axis=0).astype(np.uint8),
        sif_row_ids=np.stack(sif_row_ids, axis=0).astype(np.int64),
        mask_sums=np.stack(mask_sums, axis=0).astype(np.float32),
        mask_inside_fractions=np.stack(mask_inside_fractions, axis=0).astype(np.float32),
        chip_id=np.asarray(chip_ids),
        target_names=np.asarray(TARGET_COLUMNS),
        final_check_names=np.asarray(FINAL_CHECK_COLUMNS),
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
        nonlocal shard_id, xs, masks, y_targets, target_accepts
        nonlocal footprint_valids, sif_row_ids, mask_sums, mask_inside_fractions, chip_ids
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
        idx, chip_payload, chip_metadata, footprint_metadata = result
        if idx % 100 == 0:
            print(f"Prepared chip {idx + 1:,} / {len(manifest):,}")

        chip_id = chip_metadata["chip_id"]
        for row in footprint_metadata:
            footprint_metadata_rows.append({**row, "chip_id": chip_id})

        if chip_metadata["n_valid_footprints"] < MIN_VALID_FOOTPRINTS:
            print(
                f"Skipping chip_id={chip_id} because only "
                f"{chip_metadata['n_valid_footprints']} valid footprints remain"
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
        make_payload(idx, chip_row, assignment_groups.get(chip_row["chip_id"], pd.DataFrame()))
        for idx, chip_row in manifest.iterrows()
    )

    if N_WORKERS == 1:
        for payload in payloads:
            handle_result(build_chip_from_payload(payload))
    else:
        with ProcessPoolExecutor(max_workers=N_WORKERS) as executor:
            for result in executor.map(build_chip_from_payload, payloads, chunksize=4):
                handle_result(result)

    flush_shard()

    metadata = pd.DataFrame(metadata_rows)
    metadata.to_csv(OUTPUT_DIR / "chip_metadata.csv", index=False)

    footprint_metadata = pd.DataFrame(footprint_metadata_rows)
    footprint_metadata.to_csv(OUTPUT_DIR / "footprint_metadata.csv", index=False)

    channel_table = pd.DataFrame(
        {"channel_index": np.arange(len(CHANNEL_NAMES)), "channel_name": CHANNEL_NAMES}
    )
    channel_table.to_csv(OUTPUT_DIR / "channel_names.csv", index=False)

    target_table = pd.DataFrame(
        {
            "target_index": np.arange(len(TARGET_COLUMNS)),
            "target_name": TARGET_COLUMNS,
            "final_check_column": FINAL_CHECK_COLUMNS,
        }
    )
    target_table.to_csv(OUTPUT_DIR / "target_names.csv", index=False)

    print(f"Done. Wrote {len(metadata):,} chips to {OUTPUT_DIR}")
    print(
        "Reminder for training: compute footprint losses inside each chip, "
        "average within chip first, then average across batch."
    )


if __name__ == "__main__":
    prepare_chips()
