"""Prepare centred 4 km chips for external native-footprint SIF validation.

The two external MGRS tiles were not used for model fitting. Each output sample
contains one accepted OCO-2 footprint centred within a 200 x 200 Sentinel-2
predictor chip at 20 m resolution:

  X:                 [19, 200, 200], stored as float16
  footprint_mask:    [200, 200], fractional coverage stored as float16
  observed_sif:      scalar target_modis_sif, stored as float32

Only footprints whose centroids fall within the inner 99.8 km square of their
109.8 km Sentinel product tile are eligible. This removes a 5 km border on all
sides, ensuring that the centred 4 km chip remains inside the source raster.

Predictor NaNs are intentionally preserved. Inference must load channel means,
standard deviations, target statistics, channel order and model weights from
the saved training checkpoint. Normalize X first and then replace normalized
NaNs with zero.
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
from shapely.geometry import box

import prepare_sentinel2_multisif_cnn_chips as sentinel


# ---------------------------------------------------------------------------
# Config

SIF_CSV_PATH = Path(
    "data/main_sif_data/2tiles_2_7_M01_QF01_inoutrange_PARrm.csv"
)
MGRS_REFERENCE_DIR = Path("data/temp_data/mgrs_tifs")

OUTPUT_DIR = Path(
    "data/cnn_sentinel2_chips/"
    "external_2tiles_native_footprints_4km_20m_1p"
)

TARGET_COLUMN = "target_modis_sif"
FINAL_CHECK_COLUMN = "final_check_modis_sif"
TARGET_TILES = ("T32UQV", "T32UMC")

CHIP_SIZE_M = 4000.0
CHIP_RES_M = 20.0
CHIP_SIZE = 200
TILE_BORDER_EXCLUSION_M = 5000.0
MIN_FOOTPRINT_INSIDE_FRACTION = 0.99
# Supersampled rasterization approximates polygon area on 5 m subpixels. The
# exact vector intersection above enforces the 99% containment requirement;
# this looser ratio only catches unexpectedly malformed rasterized masks.
MIN_RASTERIZED_MASK_AREA_RATIO = 0.95
MASK_OVERSAMPLE = 4
# Five 20 m pixels form one 100 m pixel. Aligning chip origins to this lattice
# preserves compatibility with a later fixed-grid 100 m aggregation experiment.
WINDOW_ALIGNMENT_PIXELS = 5

SAMPLES_PER_TILE = 1000
STRATIFY_COLUMNS = ("sif_month", "measurement_mode", "date_align")
RANDOM_SEED = 42

# Sixteen 19-channel chips form a moderate compressed shard.
SHARD_SIZE = 16

# Set to a small positive integer for a smoke test. The limit is applied after
# the full stratified sample has been selected, so sampling summaries remain
# representative of the planned full run.
MAX_CHIPS: int | None = None

# Four workers matched the completed 4 km training-chip preparation. Reduce to
# 1 or 2 if the large Sentinel products cause storage contention.
N_WORKERS = 4

# Prevent stale shards from different runs being mixed together.
FAIL_IF_OUTPUT_EXISTS = True

CHANNEL_NAMES = sentinel.CHANNEL_NAMES

# Shared channel and mask readers use these module globals at runtime.
sentinel.CHIP_SIZE_M = CHIP_SIZE_M
sentinel.CHIP_RES_M = CHIP_RES_M
sentinel.CHIP_SIZE = CHIP_SIZE
sentinel.MASK_OVERSAMPLE = MASK_OVERSAMPLE


# ---------------------------------------------------------------------------
# Generic helpers

def require_columns(df: pd.DataFrame, columns: Iterable[str], label: str) -> None:
    missing = [column for column in columns if column not in df.columns]
    if missing:
        raise ValueError(f"{label} is missing required columns: {missing}")


def normalize_mgrs_tile(value: object) -> str:
    tile = str(value).strip().upper()
    if not tile.startswith("T"):
        tile = f"T{tile}"
    if re.fullmatch(r"T\d{2}[A-Z]{3}", tile) is None:
        raise ValueError(f"Invalid MGRS tile value: {value}")
    return tile


def reference_tif_for_tile(tile: str) -> Path:
    matches = sorted(MGRS_REFERENCE_DIR.glob(f"*_T{tile[1:]}_*_FRC_B5.tif"))
    if len(matches) != 1:
        raise FileNotFoundError(
            f"Expected one B5 reference TIFF for {tile} in {MGRS_REFERENCE_DIR}; "
            f"found {len(matches)}"
        )
    return matches[0]


def reference_grid_for_tile(tile: str) -> dict:
    path = reference_tif_for_tile(tile)
    with rasterio.open(path) as src:
        if src.crs is None:
            raise ValueError(f"Reference raster has no CRS: {path}")
        if not (
            math.isclose(abs(src.transform.a), CHIP_RES_M, abs_tol=1e-6)
            and math.isclose(abs(src.transform.e), CHIP_RES_M, abs_tol=1e-6)
        ):
            raise ValueError(
                f"Expected a {CHIP_RES_M:g} m reference grid for {tile}, found "
                f"{src.transform.a}, {src.transform.e}: {path}"
            )

        bounds = src.bounds
        inner_bounds = (
            bounds.left + TILE_BORDER_EXCLUSION_M,
            bounds.bottom + TILE_BORDER_EXCLUSION_M,
            bounds.right - TILE_BORDER_EXCLUSION_M,
            bounds.top - TILE_BORDER_EXCLUSION_M,
        )
        inner_width = inner_bounds[2] - inner_bounds[0]
        inner_height = inner_bounds[3] - inner_bounds[1]
        if inner_width <= 0 or inner_height <= 0:
            raise ValueError(f"Border exclusion removes the entire tile: {tile}")

        return {
            "tile": tile,
            "path": path,
            "crs": src.crs,
            "transform": src.transform,
            "width": int(src.width),
            "height": int(src.height),
            "bounds": bounds,
            "inner_bounds": inner_bounds,
            "inner_width_m": float(inner_width),
            "inner_height_m": float(inner_height),
        }


def centred_chip_geometry(row: pd.Series, reference: dict) -> dict:
    geometry = sentinel.build_sif_polygon_projected(row, reference["crs"])
    if geometry.is_empty or not np.isfinite(geometry.area) or geometry.area <= 0:
        return {"eligible": False, "reason": "invalid_projected_footprint"}

    centroid = geometry.centroid
    inner_left, inner_bottom, inner_right, inner_top = reference["inner_bounds"]
    centroid_in_inner = (
        inner_left <= centroid.x <= inner_right
        and inner_bottom <= centroid.y <= inner_top
    )
    if not centroid_in_inner:
        return {"eligible": False, "reason": "centroid_outside_inner_tile"}

    transform = reference["transform"]
    centre_col = (centroid.x - transform.c) / transform.a
    centre_row = (transform.f - centroid.y) / abs(transform.e)
    raw_col_off = centre_col - CHIP_SIZE / 2.0
    raw_row_off = centre_row - CHIP_SIZE / 2.0
    col_off = (
        int(round(raw_col_off / WINDOW_ALIGNMENT_PIXELS))
        * WINDOW_ALIGNMENT_PIXELS
    )
    row_off = (
        int(round(raw_row_off / WINDOW_ALIGNMENT_PIXELS))
        * WINDOW_ALIGNMENT_PIXELS
    )

    inside_raster = (
        col_off >= 0
        and row_off >= 0
        and col_off + CHIP_SIZE <= reference["width"]
        and row_off + CHIP_SIZE <= reference["height"]
    )
    if not inside_raster:
        return {"eligible": False, "reason": "centred_chip_outside_raster"}

    chip_xmin = transform.c + col_off * transform.a
    chip_ymax = transform.f + row_off * transform.e
    chip_xmax = chip_xmin + CHIP_SIZE * CHIP_RES_M
    chip_ymin = chip_ymax - CHIP_SIZE * CHIP_RES_M
    chip_polygon = box(chip_xmin, chip_ymin, chip_xmax, chip_ymax)
    geometry_inside_fraction = float(
        geometry.intersection(chip_polygon).area / geometry.area
    )
    if geometry_inside_fraction < MIN_FOOTPRINT_INSIDE_FRACTION:
        return {"eligible": False, "reason": "footprint_not_inside_centred_chip"}

    chip_center_x = chip_xmin + CHIP_SIZE_M / 2.0
    chip_center_y = chip_ymin + CHIP_SIZE_M / 2.0
    return {
        "eligible": True,
        "reason": "eligible",
        "centroid_x": float(centroid.x),
        "centroid_y": float(centroid.y),
        "chip_center_x": float(chip_center_x),
        "chip_center_y": float(chip_center_y),
        "centering_error_m": float(
            math.hypot(centroid.x - chip_center_x, centroid.y - chip_center_y)
        ),
        "chip_col_off": col_off,
        "chip_row_off": row_off,
        "chip_xmin": float(chip_xmin),
        "chip_ymin": float(chip_ymin),
        "chip_xmax": float(chip_xmax),
        "chip_ymax": float(chip_ymax),
        "geometry_inside_fraction": geometry_inside_fraction,
        "projected_footprint_area_km2": float(geometry.area / 1_000_000.0),
    }


# ---------------------------------------------------------------------------
# Loading, eligibility and sampling

def load_candidates() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    df = pd.read_csv(SIF_CSV_PATH)
    require_columns(
        df,
        [
            TARGET_COLUMN,
            FINAL_CHECK_COLUMN,
            "Delta_Date",
            "sif_doy",
            "date_align",
            "Metadata.MeasurementMode",
            "Quality_Flag",
            "mgrs_tile",
            "product_path",
            "Lat_corner1",
            "Lat_corner2",
            "Lat_corner3",
            "Lat_corner4",
            "Lon_corner1",
            "Lon_corner2",
            "Lon_corner3",
            "Lon_corner4",
        ],
        "external SIF CSV",
    )

    df = df.copy()
    df["source_csv_row"] = np.arange(1, len(df) + 1, dtype=np.int64)
    df["sif_row_id"] = df["source_csv_row"]
    df["Delta_Date"] = pd.to_datetime(df["Delta_Date"], errors="raise").dt.date
    df["sif_year"] = pd.to_datetime(df["Delta_Date"]).dt.year.astype(int)
    df["sif_month"] = pd.to_datetime(df["Delta_Date"]).dt.month.astype(int)
    df["sif_doy"] = pd.to_numeric(df["sif_doy"], errors="raise").astype(int)
    df["measurement_mode"] = pd.to_numeric(
        df["Metadata.MeasurementMode"], errors="raise"
    ).astype(int)
    df["Quality_Flag"] = pd.to_numeric(
        df["Quality_Flag"], errors="raise"
    ).astype(int)
    df[TARGET_COLUMN] = pd.to_numeric(df[TARGET_COLUMN], errors="coerce")
    df["date_align"] = df["date_align"].astype(str).str.strip().str.lower()
    df["mgrs_tile_t"] = df["mgrs_tile"].map(normalize_mgrs_tile)

    corner_columns = [
        *[f"Lat_corner{i}" for i in range(1, 5)],
        *[f"Lon_corner{i}" for i in range(1, 5)],
    ]
    for column in corner_columns:
        df[column] = pd.to_numeric(df[column], errors="coerce")

    invalid_modes = sorted(set(df["measurement_mode"]) - {0, 1})
    if invalid_modes:
        raise ValueError(f"Unexpected measurement modes: {invalid_modes}")
    invalid_quality = sorted(set(df["Quality_Flag"]) - {0, 1})
    if invalid_quality:
        raise ValueError(f"Unexpected Quality_Flag values: {invalid_quality}")
    invalid_alignment = sorted(set(df["date_align"]) - {"inrange", "outrange"})
    if invalid_alignment:
        raise ValueError(f"Unexpected date_align values: {invalid_alignment}")

    filter_rows: list[dict] = [{"stage": "input", "n_rows": int(len(df))}]

    accepted = df[FINAL_CHECK_COLUMN].astype(str).str.strip().str.lower().eq("accept")
    finite_target = np.isfinite(df[TARGET_COLUMN])
    target_tile = df["mgrs_tile_t"].isin(TARGET_TILES)
    finite_corners = np.isfinite(df[corner_columns].to_numpy()).all(axis=1)
    df = df[accepted & finite_target & target_tile & finite_corners].copy()
    filter_rows.append(
        {"stage": "accepted_finite_target_tiles_corners", "n_rows": int(len(df))}
    )

    tile_references = {tile: reference_grid_for_tile(tile) for tile in TARGET_TILES}
    geometry_rows: list[dict] = []
    for _, row in df.iterrows():
        tile = row["mgrs_tile_t"]
        product_tile = sentinel.sentinel_product_tile(Path(str(row["product_path"])))
        if product_tile != tile:
            geometry_rows.append(
                {"eligible": False, "reason": "product_tile_mismatch"}
            )
            continue
        geometry_rows.append(centred_chip_geometry(row, tile_references[tile]))

    geometry_table = pd.DataFrame(geometry_rows, index=df.index)
    df = pd.concat([df, geometry_table], axis=1)

    eligibility_summary = (
        df.groupby(["mgrs_tile_t", "reason"], dropna=False)
        .size()
        .rename("n_rows")
        .reset_index()
    )
    df = df[df["eligible"]].copy()
    filter_rows.append(
        {"stage": "centroid_and_centred_chip_eligible", "n_rows": int(len(df))}
    )

    for tile in TARGET_TILES:
        if (df["mgrs_tile_t"] == tile).sum() < SAMPLES_PER_TILE:
            print(
                f"Warning: {tile} has fewer than {SAMPLES_PER_TILE:,} eligible "
                "footprints; all eligible rows will be used."
            )

    return df.reset_index(drop=True), pd.DataFrame(filter_rows), eligibility_summary


def allocate_stratum_samples(counts: pd.Series, target_n: int) -> pd.Series:
    counts = counts.astype(int)
    if target_n >= int(counts.sum()):
        return counts.copy()

    quota = counts.astype(float) * float(target_n) / float(counts.sum())
    allocation = np.floor(quota).astype(int)

    # Preserve every observed month/mode/alignment combination when possible.
    if target_n >= len(counts):
        allocation[(counts > 0) & (allocation == 0)] = 1

    while int(allocation.sum()) > target_n:
        removable = allocation > 1
        if not removable.any():
            removable = allocation > 0
        excess = allocation.astype(float) - quota
        key = excess.where(removable, -np.inf).idxmax()
        allocation.loc[key] -= 1

    while int(allocation.sum()) < target_n:
        available = allocation < counts
        if not available.any():
            break
        deficit = quota - allocation.astype(float)
        key = deficit.where(available, -np.inf).idxmax()
        allocation.loc[key] += 1

    if int(allocation.sum()) != target_n:
        raise RuntimeError(
            f"Could not allocate {target_n} samples across strata; allocated "
            f"{int(allocation.sum())}"
        )
    if (allocation > counts).any() or (allocation < 0).any():
        raise RuntimeError("Invalid proportional stratum allocation")
    return allocation.astype(int)


def stratified_sample(candidates: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    selected_parts: list[pd.DataFrame] = []
    summary_rows: list[dict] = []

    for tile_index, tile in enumerate(TARGET_TILES):
        tile_df = candidates[candidates["mgrs_tile_t"] == tile].copy()
        target_n = min(SAMPLES_PER_TILE, len(tile_df))
        counts = tile_df.groupby(list(STRATIFY_COLUMNS), dropna=False).size()
        allocation = allocate_stratum_samples(counts, target_n)

        grouped = tile_df.groupby(list(STRATIFY_COLUMNS), dropna=False, sort=True)
        for stratum_index, (stratum_key, group) in enumerate(grouped):
            if not isinstance(stratum_key, tuple):
                stratum_key = (stratum_key,)
            n_select = int(allocation.loc[stratum_key])
            rng = np.random.default_rng(
                RANDOM_SEED + tile_index * 10000 + stratum_index
            )
            chosen = rng.choice(group.index.to_numpy(), size=n_select, replace=False)
            selected_parts.append(tile_df.loc[chosen].copy())

            summary = {
                "mgrs_tile_t": tile,
                **dict(zip(STRATIFY_COLUMNS, stratum_key)),
                "available_n": int(len(group)),
                "selected_n": n_select,
                "available_fraction_within_tile": float(len(group) / len(tile_df)),
                "selected_fraction_within_tile": float(n_select / target_n),
            }
            summary_rows.append(summary)

    selected = pd.concat(selected_parts, ignore_index=True)
    selected = selected.sort_values(
        ["mgrs_tile_t", "sif_year", "sif_doy", "source_csv_row"]
    ).reset_index(drop=True)
    selected["sample_id"] = selected.apply(
        lambda row: (
            f"external_{row['mgrs_tile_t']}_{row['sif_year']}"
            f"{int(row['sif_doy']):03d}_row{int(row['source_csv_row'])}"
        ),
        axis=1,
    )

    if selected["sample_id"].duplicated().any():
        raise ValueError("Duplicate external inference sample IDs were generated")
    for tile in TARGET_TILES:
        expected = min(
            SAMPLES_PER_TILE,
            int((candidates["mgrs_tile_t"] == tile).sum()),
        )
        actual = int((selected["mgrs_tile_t"] == tile).sum())
        if actual != expected:
            raise RuntimeError(
                f"Sampling count mismatch for {tile}: expected {expected}, got {actual}"
            )

    return selected, pd.DataFrame(summary_rows)


# ---------------------------------------------------------------------------
# Chip assembly

def chip_row_for_predictors(row: pd.Series) -> pd.Series:
    chip_row = row.copy()
    chip_row["chip_id"] = row["sample_id"]
    chip_row["chip_rows"] = CHIP_SIZE
    chip_row["chip_cols"] = CHIP_SIZE
    chip_row["chip_size_m"] = CHIP_SIZE_M
    chip_row["chip_status"] = "centred_external_footprint"
    return chip_row


def build_sample(row: pd.Series) -> tuple[dict, dict]:
    chip_row = chip_row_for_predictors(row)
    dst_transform = sentinel.chip_transform(chip_row)
    footprint_df = pd.DataFrame([row])

    x, valid_fractions, dst_crs, product_path, fapar_doy = (
        sentinel.build_predictor_chip(chip_row, footprint_df, dst_transform)
    )

    geometry = sentinel.build_sif_polygon_projected(row, dst_crs)
    mask = sentinel.fractional_footprint_mask(geometry, dst_transform)
    mask_sum = float(mask.sum(dtype=np.float64))
    mask_inside_fraction = sentinel.mask_inside_fraction(mask, geometry)
    if not np.isfinite(mask_sum) or mask_sum <= 0:
        raise ValueError(f"Empty footprint mask for sample_id={row['sample_id']}")
    if mask_inside_fraction < MIN_RASTERIZED_MASK_AREA_RATIO:
        raise ValueError(
            f"Rasterized mask area ratio is {mask_inside_fraction:.6f} for "
            f"sample_id={row['sample_id']}"
        )

    metadata = {
        "sample_id": row["sample_id"],
        "sif_row_id": int(row["sif_row_id"]),
        "source_csv_row": int(row["source_csv_row"]),
        "Delta_Date": row["Delta_Date"],
        "sif_year": int(row["sif_year"]),
        "sif_month": int(row["sif_month"]),
        "sif_doy": int(row["sif_doy"]),
        "measurement_mode": int(row["measurement_mode"]),
        "Quality_Flag": int(row["Quality_Flag"]),
        "date_align": str(row["date_align"]),
        "mgrs_tile_t": row["mgrs_tile_t"],
        "product_path": str(product_path),
        "fapar_composite_doy": int(fapar_doy),
        "par_date": row["Delta_Date"],
        "observed_sif": float(row[TARGET_COLUMN]),
        "final_check": str(row[FINAL_CHECK_COLUMN]),
        "state": row.get("state", ""),
        "hzs": row.get("hzs", ""),
        "Latitude": sentinel.to_float(row.get("Latitude", np.nan)),
        "Longitude": sentinel.to_float(row.get("Longitude", np.nan)),
        "centroid_x": float(row["centroid_x"]),
        "centroid_y": float(row["centroid_y"]),
        "chip_center_x": float(row["chip_center_x"]),
        "chip_center_y": float(row["chip_center_y"]),
        "centering_error_m": float(row["centering_error_m"]),
        "chip_col_off": int(row["chip_col_off"]),
        "chip_row_off": int(row["chip_row_off"]),
        "chip_xmin": float(row["chip_xmin"]),
        "chip_ymin": float(row["chip_ymin"]),
        "chip_xmax": float(row["chip_xmax"]),
        "chip_ymax": float(row["chip_ymax"]),
        "projected_footprint_area_km2": float(row["projected_footprint_area_km2"]),
        "geometry_inside_fraction": float(row["geometry_inside_fraction"]),
        "mask_sum_pixels": mask_sum,
        "mask_area_km2": mask_sum * CHIP_RES_M**2 / 1_000_000.0,
        "mask_inside_fraction": float(mask_inside_fraction),
        **valid_fractions,
    }

    payload = {
        "X": x,
        "footprint_mask": mask.astype(np.float16),
        "observed_sif": np.float32(row[TARGET_COLUMN]),
        "sif_row_id": np.int64(row["sif_row_id"]),
    }
    return payload, metadata


def make_payload(index: int, row: pd.Series) -> tuple[int, dict]:
    return index, row.to_dict()


def build_sample_from_payload(payload: tuple[int, dict]) -> tuple[int, dict, dict]:
    index, row_dict = payload
    sample, metadata = build_sample(pd.Series(row_dict))
    return index, sample, metadata


# ---------------------------------------------------------------------------
# Shards and outputs

def check_output_directory() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    generated_names = [
        "chip_metadata.csv",
        "sampling_summary.csv",
        "filter_summary.csv",
        "eligibility_summary.csv",
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
    masks: list[np.ndarray],
    targets: list[np.float32],
    sif_row_ids: list[np.int64],
    sample_ids: list[str],
) -> Path:
    path = OUTPUT_DIR / f"chips_{shard_id:05d}.npz"
    np.savez_compressed(
        path,
        X=np.stack(xs, axis=0).astype(np.float16),
        footprint_mask=np.stack(masks, axis=0).astype(np.float16),
        observed_sif=np.asarray(targets, dtype=np.float32),
        sif_row_id=np.asarray(sif_row_ids, dtype=np.int64),
        sample_id=np.asarray(sample_ids),
        channel_names=np.asarray(CHANNEL_NAMES),
        target_name=np.asarray([TARGET_COLUMN]),
    )
    return path


def prepare_chips() -> None:
    check_output_directory()
    candidates, filter_summary, eligibility_summary = load_candidates()
    selected, sampling_summary = stratified_sample(candidates)
    planned_sample_count = len(selected)

    if MAX_CHIPS is not None:
        selected = selected.head(MAX_CHIPS).copy().reset_index(drop=True)
        print(
            f"Smoke-test limit active: preparing {len(selected):,} of "
            f"{planned_sample_count:,} selected footprints"
        )

    print(
        f"Selected {planned_sample_count:,} external footprints from "
        f"{len(candidates):,} eligible candidates"
    )
    print(
        selected.groupby("mgrs_tile_t").size().rename("n_selected").to_string()
    )

    metadata_rows: list[dict] = []
    xs: list[np.ndarray] = []
    masks: list[np.ndarray] = []
    targets: list[np.float32] = []
    sif_row_ids: list[np.int64] = []
    sample_ids: list[str] = []
    shard_id = 0

    def flush_shard() -> None:
        nonlocal shard_id, xs, masks, targets, sif_row_ids, sample_ids
        if not xs:
            return
        written = write_shard(
            shard_id,
            xs,
            masks,
            targets,
            sif_row_ids,
            sample_ids,
        )
        print(f"Wrote {written}")
        shard_id += 1
        xs = []
        masks = []
        targets = []
        sif_row_ids = []
        sample_ids = []

    def handle_result(result: tuple[int, dict, dict]) -> None:
        index, sample, metadata = result
        if index % 10 == 0:
            print(f"Prepared footprint {index + 1:,} / {len(selected):,}")

        metadata["shard_file"] = f"chips_{shard_id:05d}.npz"
        metadata["shard_index"] = len(xs)
        metadata_rows.append(metadata)

        xs.append(sample["X"])
        masks.append(sample["footprint_mask"])
        targets.append(sample["observed_sif"])
        sif_row_ids.append(sample["sif_row_id"])
        sample_ids.append(str(metadata["sample_id"]))
        if len(xs) == SHARD_SIZE:
            flush_shard()

    payloads = (
        make_payload(index, row)
        for index, (_, row) in enumerate(selected.iterrows())
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
    metadata.to_csv(OUTPUT_DIR / "chip_metadata.csv", index=False)
    sampling_summary.to_csv(OUTPUT_DIR / "sampling_summary.csv", index=False)
    filter_summary.to_csv(OUTPUT_DIR / "filter_summary.csv", index=False)
    eligibility_summary.to_csv(
        OUTPUT_DIR / "eligibility_summary.csv", index=False
    )

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
            "sample_support": ["one native OCO-2 footprint"],
        }
    ).to_csv(OUTPUT_DIR / "target_names.csv", index=False)

    tile_config = {}
    for tile in TARGET_TILES:
        reference = reference_grid_for_tile(tile)
        tile_config[tile] = {
            "reference_tif": str(reference["path"]),
            "raster_width_pixels": reference["width"],
            "raster_height_pixels": reference["height"],
            "inner_width_m": reference["inner_width_m"],
            "inner_height_m": reference["inner_height_m"],
        }

    config = {
        "input_sif_csv": str(SIF_CSV_PATH),
        "output_sample_count": int(len(metadata)),
        "planned_full_sample_count": int(planned_sample_count),
        "target_tiles": list(TARGET_TILES),
        "samples_per_tile": SAMPLES_PER_TILE,
        "sampling_seed": RANDOM_SEED,
        "sampling_strata": list(STRATIFY_COLUMNS),
        "chip_size_m": CHIP_SIZE_M,
        "chip_resolution_m": CHIP_RES_M,
        "chip_rows": CHIP_SIZE,
        "chip_cols": CHIP_SIZE,
        "window_alignment_pixels": WINDOW_ALIGNMENT_PIXELS,
        "window_alignment_m": WINDOW_ALIGNMENT_PIXELS * CHIP_RES_M,
        "tile_border_exclusion_m_per_side": TILE_BORDER_EXCLUSION_M,
        "minimum_footprint_inside_fraction": MIN_FOOTPRINT_INSIDE_FRACTION,
        "minimum_rasterized_mask_area_ratio": MIN_RASTERIZED_MASK_AREA_RATIO,
        "mask_oversample": MASK_OVERSAMPLE,
        "channels": CHANNEL_NAMES,
        "target": TARGET_COLUMN,
        "final_check": FINAL_CHECK_COLUMN,
        "predictor_storage_dtype": "float16",
        "mask_storage_dtype": "float16",
        "target_storage_dtype": "float32",
        "predictor_nan_policy": (
            "preserve; use checkpoint channel statistics, normalize, then fill "
            "normalized NaNs with zero"
        ),
        "normalization_source": "saved model checkpoint only",
        "footprint_prediction": (
            "normalize the float32 fractional mask to sum to one, then compute "
            "sum(mask * denormalized predicted 20 m SIF map)"
        ),
        "evi_valid_range": [sentinel.EVI_VALID_MIN, sentinel.EVI_VALID_MAX],
        "evi_invalid_policy": "mask as NaN; do not clip",
        "par_units": "W/m2",
        "par_accepted_qa_codes": list(sentinel.PAR_ACCEPTED_QA_CODES),
        "par_resampling": "bilinear after native-grid QA masking",
        "fapar_temporal_matching": "containing 8-day composite interval",
        "tile_references": tile_config,
        "channel_reader_module": "prepare_sentinel2_multisif_cnn_chips.py",
    }
    with (OUTPUT_DIR / "dataset_config.json").open(
        "w", encoding="ascii"
    ) as file:
        json.dump(config, file, indent=2)

    print(f"Done. Wrote {len(metadata):,} footprint chips to {OUTPUT_DIR}")
    print(
        "Inference reminder: load channel order, channel mean/std, target "
        "mean/std and model weights from the saved .pt checkpoint. Cast X and "
        "masks to float32; normalize X before replacing NaNs with zero; "
        "renormalize each footprint mask to sum to one."
    )


if __name__ == "__main__":
    prepare_chips()
