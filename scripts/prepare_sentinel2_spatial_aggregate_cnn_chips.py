"""Prepare 4 km Sentinel-2 CNN chips with spatially aggregated SIF labels.

Inputs are produced by spatial_grouping.R:
  - fixed_grid_4000m_aggregate_manifest.csv: one row per tile/date/mode/cell
  - fixed_grid_4000m_sif_assignments.csv: one row per contributing footprint

Each output sample contains:
  X:                    [19, 200, 200], float16 predictor channels
  aggregate_weight_map: [200, 200], float16 supervision weights
  y_aggregate:          scalar equal-mean footprint SIF target, float32

For footprint masks m_i, the aggregate supervision map is

    w(p) = mean_i(m_i(p) / sum_p(m_i(p)))

and therefore the model's scalar prediction is

    y_hat = sum_p(w(p) * predicted_sif_map(p)).

This exactly gives every represented footprint equal weight. Predictor NaNs
are intentionally preserved. Training statistics should ignore NaNs; after
normalization, replace remaining NaNs with zero (the normalized channel mean).
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

import prepare_sentinel2_multisif_cnn_chips as sentinel


# ---------------------------------------------------------------------------
# Config

AGGREGATION_DIR = Path("data/sentinel2_spatial_aggregation_4000m_original_tiles_inout")
MANIFEST_PATH = AGGREGATION_DIR / "fixed_grid_4000m_aggregate_manifest.csv"
ASSIGNMENTS_PATH = AGGREGATION_DIR / "fixed_grid_4000m_sif_assignments.csv"

OUTPUT_DIR = Path(
    "data/cnn_sentinel2_chips/"
    "spatial_aggregate_4km_20m_indices_updated_inout_evi_masked"
)

TARGET_COLUMN = "target_modis_sif"
AGGREGATE_TARGET_COLUMN = "aggregated_target_modis_sif"
FINAL_CHECK_COLUMN = "final_check_modis_sif"

CHIP_SIZE_M = 4000.0
CHIP_RES_M = 20.0
CHIP_SIZE = 200
MIN_FOOTPRINTS = 5
MASK_OVERSAMPLE = 4

# Sixteen 19-channel 200 x 200 samples form a reasonably sized compressed shard.
SHARD_SIZE = 16

# Set to a small integer for a smoke test, otherwise leave as None.
MAX_CHIPS: int | None = None

# Keep at 1 initially. Increasing this can cause heavy contention on the large
# Sentinel-2 products even when more CPU cores are available.
N_WORKERS = 4

# Fail instead of mixing shards from separate runs in the same directory.
FAIL_IF_OUTPUT_EXISTS = True

CHANNEL_NAMES = sentinel.CHANNEL_NAMES

# The shared channel readers use these module globals to define output arrays.
# They are changed once, before any chips are built.
sentinel.CHIP_SIZE_M = CHIP_SIZE_M
sentinel.CHIP_RES_M = CHIP_RES_M
sentinel.CHIP_SIZE = CHIP_SIZE
sentinel.MASK_OVERSAMPLE = MASK_OVERSAMPLE


# ---------------------------------------------------------------------------
# Input tables


def require_columns(df: pd.DataFrame, columns: Iterable[str], label: str) -> None:
    missing = [column for column in columns if column not in df.columns]
    if missing:
        raise ValueError(f"{label} is missing required columns: {missing}")


def parse_bool(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False)
    return series.astype(str).str.strip().str.lower().isin({"true", "t", "1", "yes"})


def load_aggregation_tables() -> tuple[pd.DataFrame, dict[str, pd.DataFrame]]:
    manifest = pd.read_csv(MANIFEST_PATH)
    assignments = pd.read_csv(ASSIGNMENTS_PATH)

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
            "Delta_Date",
            "sif_year",
            "sif_month",
            "sif_doy",
            "measurement_mode",
            "n_footprints",
            AGGREGATE_TARGET_COLUMN,
            "n_product_paths",
            "eligible_n5",
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
            "Quality_Flag",
            "date_align",
            "mgrs_tile_t",
            "product_path",
            TARGET_COLUMN,
            FINAL_CHECK_COLUMN,
            "Lat_corner1",
            "Lat_corner2",
            "Lat_corner3",
            "Lat_corner4",
            "Lon_corner1",
            "Lon_corner2",
            "Lon_corner3",
            "Lon_corner4",
        ],
        "SIF assignments",
    )

    manifest = manifest.copy()
    assignments = assignments.copy()
    manifest["Delta_Date"] = pd.to_datetime(manifest["Delta_Date"], errors="raise").dt.date
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
        "n_product_paths",
    ]
    assignment_numeric = [
        "sif_row_id",
        "source_csv_row",
        "sif_year",
        "sif_month",
        "sif_doy",
        "measurement_mode",
        "Quality_Flag",
        TARGET_COLUMN,
        "Lat_corner1",
        "Lat_corner2",
        "Lat_corner3",
        "Lat_corner4",
        "Lon_corner1",
        "Lon_corner2",
        "Lon_corner3",
        "Lon_corner4",
    ]
    for column in manifest_numeric:
        manifest[column] = pd.to_numeric(manifest[column], errors="raise")
    for column in assignment_numeric:
        assignments[column] = pd.to_numeric(assignments[column], errors="raise")

    assignments["date_align"] = (
        assignments["date_align"].astype(str).str.strip().str.lower()
    )
    invalid_quality_flags = sorted(
        set(assignments["Quality_Flag"].astype(int)) - {0, 1}
    )
    if invalid_quality_flags:
        raise ValueError(
            "Unexpected Quality_Flag values in SIF assignments: "
            f"{invalid_quality_flags}"
        )

    invalid_date_align = sorted(
        set(assignments["date_align"]) - {"inrange", "outrange"}
    )
    if invalid_date_align:
        raise ValueError(
            "Unexpected date_align values in SIF assignments: "
            f"{invalid_date_align}"
        )

    manifest = manifest[
        parse_bool(manifest["eligible_n5"])
        & (manifest["n_footprints"] >= MIN_FOOTPRINTS)
        & np.isclose(manifest["aggregation_size_m"], CHIP_SIZE_M)
        & (manifest["cell_pixels"].astype(int) == CHIP_SIZE)
        & (manifest["n_product_paths"].astype(int) == 1)
    ].copy()

    if "all_source_tiles_match_input" in manifest.columns:
        manifest = manifest[parse_bool(manifest["all_source_tiles_match_input"])].copy()

    manifest = manifest.sort_values(
        ["sif_year", "sif_doy", "mgrs_tile_t", "cell_id", "aggregation_id"]
    ).reset_index(drop=True)
    if manifest["aggregation_id"].duplicated().any():
        duplicate = manifest.loc[
            manifest["aggregation_id"].duplicated(), "aggregation_id"
        ].iloc[0]
        raise ValueError(f"Duplicate aggregation_id in manifest: {duplicate}")

    if MAX_CHIPS is not None:
        manifest = manifest.head(MAX_CHIPS).copy()

    assignments = assignments[
        assignments["aggregation_id"].isin(manifest["aggregation_id"])
    ].copy()
    assignments = assignments.sort_values(
        ["aggregation_id", "sif_row_id"]
    ).reset_index(drop=True)
    if assignments["sif_row_id"].duplicated().any():
        duplicate = assignments.loc[
            assignments["sif_row_id"].duplicated(), "sif_row_id"
        ].iloc[0]
        raise ValueError(
            f"sif_row_id={int(duplicate)} is assigned to more than one 4 km group"
        )

    groups: dict[str, pd.DataFrame] = {}
    manifest_lookup = manifest.set_index("aggregation_id", drop=False)

    for aggregation_id, group in assignments.groupby("aggregation_id", sort=False):
        group = group.reset_index(drop=True)
        manifest_row = manifest_lookup.loc[aggregation_id]

        homogeneous_columns = [
            "Delta_Date",
            "mgrs_tile_t",
            "measurement_mode",
            "product_path",
            "date_align",
        ]
        invalid = {
            column: int(group[column].nunique(dropna=False))
            for column in homogeneous_columns
            if group[column].nunique(dropna=False) != 1
        }
        if invalid:
            raise ValueError(
                f"aggregation_id={aggregation_id} is not homogeneous: {invalid}"
            )

        if len(group) != int(manifest_row["n_footprints"]):
            raise ValueError(
                f"Footprint count mismatch for {aggregation_id}: manifest="
                f"{int(manifest_row['n_footprints'])}, assignments={len(group)}"
            )
        if group.iloc[0]["Delta_Date"] != manifest_row["Delta_Date"]:
            raise ValueError(f"Date mismatch for {aggregation_id}")
        if str(group.iloc[0]["mgrs_tile_t"]) != str(manifest_row["mgrs_tile_t"]):
            raise ValueError(f"MGRS tile mismatch for {aggregation_id}")
        if float(group.iloc[0]["measurement_mode"]) != float(
            manifest_row["measurement_mode"]
        ):
            raise ValueError(f"Measurement mode mismatch for {aggregation_id}")

        final_checks = group[FINAL_CHECK_COLUMN].astype(str).str.strip().str.lower()
        if not final_checks.eq("accept").all():
            raise ValueError(f"Non-accepted SIF target found in {aggregation_id}")
        if not np.isfinite(group[TARGET_COLUMN].to_numpy(dtype=np.float64)).all():
            raise ValueError(f"Non-finite SIF target found in {aggregation_id}")

        assignment_mean = float(group[TARGET_COLUMN].mean())
        manifest_mean = float(manifest_row[AGGREGATE_TARGET_COLUMN])
        if not np.isclose(assignment_mean, manifest_mean, rtol=1e-6, atol=1e-7):
            raise ValueError(
                f"Aggregate target mismatch for {aggregation_id}: "
                f"manifest={manifest_mean}, assignment mean={assignment_mean}"
            )

        groups[str(aggregation_id)] = group

    missing = sorted(set(manifest["aggregation_id"].astype(str)) - set(groups))
    if missing:
        raise ValueError(f"Manifest aggregations have no assignments: {missing[:5]}")

    print(
        f"Loaded {len(manifest):,} eligible 4 km groups and "
        f"{len(assignments):,} contributing SIF footprints"
    )
    return manifest, groups


# ---------------------------------------------------------------------------
# Aggregate supervision


def chip_row_for_predictors(manifest_row: pd.Series) -> pd.Series:
    row = manifest_row.copy()
    row["chip_id"] = row["aggregation_id"]
    row["chip_xmin"] = float(row["cell_xmin"])
    row["chip_ymin"] = float(row["cell_ymin"])
    row["chip_xmax"] = float(row["cell_xmax"])
    row["chip_ymax"] = float(row["cell_ymax"])
    row["chip_rows"] = CHIP_SIZE
    row["chip_cols"] = CHIP_SIZE
    row["chip_size_m"] = CHIP_SIZE_M
    row["chip_status"] = "inside"
    return row


def build_aggregate_weight_map(
    footprints: pd.DataFrame,
    dst_transform,
    dst_crs,
) -> tuple[np.ndarray, list[dict], dict]:
    normalized_masks: list[np.ndarray] = []
    footprint_rows: list[dict] = []

    for slot, (_, footprint) in enumerate(footprints.iterrows()):
        geometry = sentinel.build_sif_polygon_projected(footprint, dst_crs)
        mask = sentinel.fractional_footprint_mask(geometry, dst_transform)
        mask_sum = float(mask.sum())
        inside_fraction = sentinel.mask_inside_fraction(mask, geometry)

        # Centroid assignment should guarantee intersection. Failing loudly is
        # preferable to changing the aggregate target by silently dropping rows.
        if not np.isfinite(mask_sum) or mask_sum <= 0:
            raise ValueError(
                f"Empty footprint mask for sif_row_id={int(footprint['sif_row_id'])}"
            )

        normalized_masks.append(mask / mask_sum)
        footprint_rows.append(
            {
                "slot": slot,
                "sif_row_id": int(footprint["sif_row_id"]),
                "source_csv_row": int(footprint["source_csv_row"]),
                "observed_sif": float(footprint[TARGET_COLUMN]),
                "mask_sum_pixels": mask_sum,
                "mask_area_km2_on_chip": mask_sum * CHIP_RES_M**2 / 1_000_000.0,
                "mask_inside_fraction": inside_fraction,
                "Delta_Date": footprint["Delta_Date"],
                "measurement_mode": int(footprint["measurement_mode"]),
                "quality_flag": int(footprint["Quality_Flag"]),
                "date_align": str(footprint["date_align"]),
                "mgrs_tile_t": footprint["mgrs_tile_t"],
                "state": footprint.get("state", ""),
                "hzs": footprint.get("hzs", ""),
                "Latitude": sentinel.to_float(footprint.get("Latitude", np.nan)),
                "Longitude": sentinel.to_float(footprint.get("Longitude", np.nan)),
                "Lat_corner1": float(footprint["Lat_corner1"]),
                "Lat_corner2": float(footprint["Lat_corner2"]),
                "Lat_corner3": float(footprint["Lat_corner3"]),
                "Lat_corner4": float(footprint["Lat_corner4"]),
                "Lon_corner1": float(footprint["Lon_corner1"]),
                "Lon_corner2": float(footprint["Lon_corner2"]),
                "Lon_corner3": float(footprint["Lon_corner3"]),
                "Lon_corner4": float(footprint["Lon_corner4"]),
            }
        )

    # Each normalized footprint mask sums to one. Their arithmetic mean gives
    # each footprint equal influence regardless of area or footprint overlap.
    aggregate_weight = np.mean(np.stack(normalized_masks, axis=0), axis=0).astype(
        np.float32
    )
    weight_sum = float(aggregate_weight.sum(dtype=np.float64))
    if not np.isclose(weight_sum, 1.0, rtol=1e-6, atol=1e-6):
        raise ValueError(f"Aggregate weight map sums to {weight_sum}, expected 1")

    # Renormalize once to remove tiny accumulation error before float16 storage.
    aggregate_weight /= aggregate_weight.sum(dtype=np.float64)
    positive = aggregate_weight > 0
    diagnostics = {
        "aggregate_weight_sum_float32": float(
            aggregate_weight.sum(dtype=np.float64)
        ),
        "aggregate_support_pixels": int(positive.sum()),
        "aggregate_support_fraction": float(positive.mean()),
        "effective_weighted_pixels": float(
            1.0 / np.square(aggregate_weight.astype(np.float64)).sum()
        ),
        "mean_mask_inside_fraction": float(
            np.mean([row["mask_inside_fraction"] for row in footprint_rows])
        ),
        "min_mask_inside_fraction": float(
            np.min([row["mask_inside_fraction"] for row in footprint_rows])
        ),
        "max_mask_inside_fraction": float(
            np.max([row["mask_inside_fraction"] for row in footprint_rows])
        ),
    }
    return aggregate_weight.astype(np.float16), footprint_rows, diagnostics


# ---------------------------------------------------------------------------
# Complete sample assembly


def build_sample(
    manifest_row: pd.Series,
    footprints: pd.DataFrame,
) -> tuple[dict, dict, list[dict]]:
    chip_row = chip_row_for_predictors(manifest_row)
    dst_transform = sentinel.chip_transform(chip_row)

    x, valid_fractions, dst_crs, product_path, fapar_doy = (
        sentinel.build_predictor_chip(chip_row, footprints, dst_transform)
    )
    weight_map, footprint_metadata, weight_diagnostics = build_aggregate_weight_map(
        footprints,
        dst_transform,
        dst_crs,
    )

    assignment_target = float(footprints[TARGET_COLUMN].mean())
    manifest_target = float(manifest_row[AGGREGATE_TARGET_COLUMN])
    if not np.isclose(assignment_target, manifest_target, rtol=1e-6, atol=1e-7):
        raise ValueError(
            f"Target changed while building {manifest_row['aggregation_id']}"
        )

    quality_flags = footprints["Quality_Flag"].astype(int)
    quality_flag_values = ";".join(
        str(value) for value in sorted(quality_flags.unique().tolist())
    )
    n_quality_flag_0 = int((quality_flags == 0).sum())
    n_quality_flag_1 = int((quality_flags == 1).sum())
    date_align = str(footprints.iloc[0]["date_align"])

    metadata = {
        "aggregation_id": manifest_row["aggregation_id"],
        "cell_id": manifest_row["cell_id"],
        "Delta_Date": manifest_row["Delta_Date"],
        "sif_year": int(manifest_row["sif_year"]),
        "sif_month": int(manifest_row["sif_month"]),
        "sif_doy": int(manifest_row["sif_doy"]),
        "measurement_mode": int(manifest_row["measurement_mode"]),
        "date_align": date_align,
        "quality_flag_values": quality_flag_values,
        "n_quality_flag_0": n_quality_flag_0,
        "n_quality_flag_1": n_quality_flag_1,
        "quality_flag_1_fraction": n_quality_flag_1 / len(footprints),
        "mgrs_tile_t": manifest_row["mgrs_tile_t"],
        "product_path": str(product_path),
        "fapar_composite_doy": int(fapar_doy),
        "par_date": manifest_row["Delta_Date"],
        "n_footprints": int(len(footprints)),
        "y_aggregate": assignment_target,
        "manifest_aggregated_target_modis_sif": manifest_target,
        "target_sd": sentinel.to_float(
            manifest_row.get("sd_target_modis_sif", np.nan)
        ),
        "target_se": sentinel.to_float(
            manifest_row.get("se_target_modis_sif", np.nan)
        ),
        "states": manifest_row.get("states", ""),
        "hzs_values": manifest_row.get("hzs_values", ""),
        "cell_xmin": float(manifest_row["cell_xmin"]),
        "cell_ymin": float(manifest_row["cell_ymin"]),
        "cell_xmax": float(manifest_row["cell_xmax"]),
        "cell_ymax": float(manifest_row["cell_ymax"]),
        **weight_diagnostics,
        **valid_fractions,
    }

    for row in footprint_metadata:
        row["aggregation_id"] = manifest_row["aggregation_id"]

    payload = {
        "X": x,
        "aggregate_weight_map": weight_map,
        "y_aggregate": np.float32(assignment_target),
        "n_footprints": np.int16(len(footprints)),
    }
    return payload, metadata, footprint_metadata


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
# Shards and metadata


def check_output_directory() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    existing = list(OUTPUT_DIR.glob("chips_*.npz"))
    existing.extend(
        path
        for path in [
            OUTPUT_DIR / "chip_metadata.csv",
            OUTPUT_DIR / "footprint_metadata.csv",
            OUTPUT_DIR / "channel_names.csv",
            OUTPUT_DIR / "target_names.csv",
            OUTPUT_DIR / "dataset_config.json",
        ]
        if path.exists()
    )
    if FAIL_IF_OUTPUT_EXISTS and existing:
        names = ", ".join(path.name for path in existing[:5])
        raise FileExistsError(
            f"Output directory already contains generated files ({names}). "
            "Use an empty output directory before starting a new run."
        )


def write_shard(
    shard_id: int,
    xs: list[np.ndarray],
    weights: list[np.ndarray],
    targets: list[np.float32],
    n_footprints: list[np.int16],
    aggregation_ids: list[str],
) -> Path:
    path = OUTPUT_DIR / f"chips_{shard_id:05d}.npz"
    np.savez_compressed(
        path,
        X=np.stack(xs, axis=0).astype(np.float16),
        aggregate_weight_map=np.stack(weights, axis=0).astype(np.float16),
        y_aggregate=np.asarray(targets, dtype=np.float32),
        n_footprints=np.asarray(n_footprints, dtype=np.int16),
        aggregation_id=np.asarray(aggregation_ids),
        channel_names=np.asarray(CHANNEL_NAMES),
        target_name=np.asarray([AGGREGATE_TARGET_COLUMN]),
    )
    return path


def prepare_chips() -> None:
    check_output_directory()
    manifest, assignment_groups = load_aggregation_tables()

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
            print(f"Prepared group {index + 1:,} / {len(manifest):,}")

        current_shard_id = shard_id
        current_shard_index = len(xs)
        metadata["shard_file"] = f"chips_{current_shard_id:05d}.npz"
        metadata["shard_index"] = current_shard_index
        metadata_rows.append(metadata)

        for row in footprint_metadata:
            row["shard_file"] = metadata["shard_file"]
            row["shard_index"] = current_shard_index
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
        "channels": CHANNEL_NAMES,
        "target": AGGREGATE_TARGET_COLUMN,
        "assignment_target": TARGET_COLUMN,
        "target_aggregation": "equal arithmetic mean across assigned footprints",
        "supervision_weight_formula": (
            "mean_i(mask_i / sum_pixels(mask_i)); model scalar prediction is "
            "sum_pixels(weight_map * predicted_sif_map)"
        ),
        "mask_oversample": MASK_OVERSAMPLE,
        "partial_footprint_policy": (
            "retain centroid-assigned footprints; masks are clipped to the 4 km "
            "cell and mask_inside_fraction is recorded per footprint"
        ),
        "predictor_storage_dtype": "float16",
        "weight_map_storage_dtype": "float16",
        "target_storage_dtype": "float32",
        "predictor_nan_policy": (
            "preserve; compute nan-aware training statistics, normalize, then fill "
            "normalized NaNs with zero"
        ),
        "evi_valid_range": [sentinel.EVI_VALID_MIN, sentinel.EVI_VALID_MAX],
        "evi_invalid_policy": "mask as NaN; do not clip",
        "par_units": "W/m2",
        "par_accepted_qa_codes": list(sentinel.PAR_ACCEPTED_QA_CODES),
        "par_resampling": "bilinear after native-grid QA masking",
        "fapar_temporal_matching": "containing 8-day composite interval",
        "split_requirement": (
            "split by whole Delta_Date before training; reserve untouched dates "
            "for native-footprint evaluation"
        ),
        "metadata_fields": {
            "footprint": ["quality_flag", "date_align"],
            "chip": [
                "date_align",
                "quality_flag_values",
                "n_quality_flag_0",
                "n_quality_flag_1",
                "quality_flag_1_fraction",
            ],
        },
        "channel_reader_module": "prepare_sentinel2_multisif_cnn_chips.py",
    }
    with (OUTPUT_DIR / "dataset_config.json").open("w", encoding="ascii") as file:
        json.dump(config, file, indent=2)

    print(f"Done. Wrote {len(metadata):,} aggregate chips to {OUTPUT_DIR}")
    print(
        "Training reminder: cast X/aggregate_weight_map to float32, calculate "
        "nan-aware normalization statistics from training dates only, normalize "
        "X, replace normalized NaNs with zero, and renormalize each float32 "
        "weight map to sum to one before computing the weighted prediction."
    )


if __name__ == "__main__":
    prepare_chips()
