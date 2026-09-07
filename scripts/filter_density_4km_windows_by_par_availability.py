"""Remove density-aggregated 4 km windows with confirmed missing PAR data.

The source valid95 manifest and assignments are left unchanged. Dates marked
``no_granules`` in the VIIRS PAR download checklist are excluded, and a new
PAR-available manifest/assignment folder is written for CNN chip preparation.

Checklist dates with ``error`` or missing records are not silently removed.
They raise an error because they may represent an incomplete download rather
than a date on which the VIIRS product genuinely had no granules.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Configuration

PROJECT_ROOT = Path(__file__).resolve().parent
DENSITY_ROOT = PROJECT_ROOT / "data" / "density_aggregation"

INPUT_DIR = (
    DENSITY_ROOT
    / "sentinel2_spatial_aggregation_density_4000m_landcover_combined_"
    "mask60_min4_s2valid95"
)
OUTPUT_DIR = (
    DENSITY_ROOT
    / "sentinel2_spatial_aggregation_density_4000m_landcover_combined_"
    "mask60_min4_s2valid95_par_available"
)

MANIFEST_FILENAME = "density_cluster_4000m_aggregate_manifest.csv"
ASSIGNMENTS_FILENAME = "density_cluster_4000m_sif_assignments.csv"

INPUT_MANIFEST_PATH = INPUT_DIR / MANIFEST_FILENAME
INPUT_ASSIGNMENTS_PATH = INPUT_DIR / ASSIGNMENTS_FILENAME
PAR_CHECKLIST_PATH = (
    PROJECT_ROOT
    / "data"
    / "viirs_vnp18a2_daily_mean_par_germany_native"
    / "vnp18a2_daily_mean_par_download_checklist.csv"
)

OUTPUT_MANIFEST_PATH = OUTPUT_DIR / MANIFEST_FILENAME
OUTPUT_ASSIGNMENTS_PATH = OUTPUT_DIR / ASSIGNMENTS_FILENAME
EXCLUDED_WINDOWS_PATH = OUTPUT_DIR / "par_no_granules_excluded_windows.csv"
FILTER_SUMMARY_PATH = OUTPUT_DIR / "par_availability_filter_summary.csv"

TARGET_COLUMN = "target_modis_sif"
AGGREGATE_TARGET_COLUMN = "aggregated_target_modis_sif"
AVAILABLE_STATUS = "completed"
NO_DATA_STATUS = "no_granules"


# ---------------------------------------------------------------------------
# Helpers

def require_columns(
    table: pd.DataFrame,
    columns: Iterable[str],
    table_name: str,
) -> None:
    missing = [column for column in columns if column not in table.columns]
    if missing:
        raise ValueError(f"{table_name} is missing required columns: {missing}")


def atomic_write_csv(table: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = path.with_suffix(path.suffix + ".partial")
    try:
        table.to_csv(temporary_path, index=False)
        os.replace(temporary_path, path)
    finally:
        if temporary_path.exists():
            temporary_path.unlink()


def validate_retained_tables(
    manifest: pd.DataFrame,
    assignments: pd.DataFrame,
) -> None:
    manifest_ids = set(manifest["aggregation_id"])
    assignment_ids = set(assignments["aggregation_id"])
    if manifest_ids != assignment_ids:
        raise ValueError("Retained manifest and assignment aggregation IDs differ")

    expected_counts = manifest.set_index("aggregation_id")["n_footprints"].astype(int)
    observed_counts = assignments.groupby("aggregation_id").size()
    if not expected_counts.sort_index().equals(observed_counts.sort_index()):
        raise ValueError("Retained manifest and assignment footprint counts differ")

    assignment_means = assignments.groupby("aggregation_id")[TARGET_COLUMN].mean()
    manifest_means = manifest.set_index("aggregation_id")[AGGREGATE_TARGET_COLUMN]
    if not np.allclose(
        assignment_means.sort_index().to_numpy(dtype=float),
        manifest_means.sort_index().to_numpy(dtype=float),
        rtol=1e-6,
        atol=1e-7,
    ):
        raise ValueError("Retained aggregate targets do not match assignment means")


# ---------------------------------------------------------------------------
# Filtering

def load_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    manifest = pd.read_csv(INPUT_MANIFEST_PATH, low_memory=False)
    assignments = pd.read_csv(INPUT_ASSIGNMENTS_PATH, low_memory=False)
    checklist = pd.read_csv(PAR_CHECKLIST_PATH, low_memory=False)

    require_columns(
        manifest,
        [
            "aggregation_id",
            "Delta_Date",
            "mgrs_tile_t",
            "measurement_mode",
            "n_footprints",
            AGGREGATE_TARGET_COLUMN,
        ],
        "density manifest",
    )
    require_columns(
        assignments,
        ["aggregation_id", "sif_row_id", "Delta_Date", TARGET_COLUMN],
        "density assignments",
    )
    require_columns(
        checklist,
        [
            "date",
            "status",
            "granules_found",
            "par_output",
            "quality_output",
            "valid_par_pixels",
            "error",
        ],
        "PAR download checklist",
    )

    manifest = manifest.copy()
    assignments = assignments.copy()
    checklist = checklist.copy()

    manifest["aggregation_id"] = manifest["aggregation_id"].astype(str)
    assignments["aggregation_id"] = assignments["aggregation_id"].astype(str)
    manifest["Delta_Date"] = pd.to_datetime(
        manifest["Delta_Date"], errors="raise"
    ).dt.date
    assignments["Delta_Date"] = pd.to_datetime(
        assignments["Delta_Date"], errors="raise"
    ).dt.date
    checklist["date"] = pd.to_datetime(checklist["date"], errors="raise").dt.date
    checklist["status"] = checklist["status"].astype(str).str.strip().str.lower()

    if manifest["aggregation_id"].duplicated().any():
        raise ValueError("The source manifest contains duplicate aggregation IDs")
    if checklist["date"].duplicated().any():
        duplicate_dates = checklist.loc[
            checklist["date"].duplicated(keep=False), "date"
        ].drop_duplicates()
        raise ValueError(
            "The PAR checklist contains duplicate dates: "
            + ", ".join(str(value) for value in duplicate_dates.head(5))
        )

    return manifest, assignments, checklist


def filter_by_par_availability(
    manifest: pd.DataFrame,
    assignments: pd.DataFrame,
    checklist: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    checklist_fields = checklist[
        [
            "date",
            "status",
            "granules_found",
            "par_output",
            "quality_output",
            "valid_par_pixels",
            "error",
        ]
    ].rename(
        columns={
            "date": "par_checklist_date",
            "status": "par_checklist_status",
            "granules_found": "par_granules_found",
            "par_output": "par_output_path",
            "quality_output": "par_quality_output_path",
            "valid_par_pixels": "par_valid_pixels",
            "error": "par_checklist_error",
        }
    )

    joined = manifest.merge(
        checklist_fields,
        left_on="Delta_Date",
        right_on="par_checklist_date",
        how="left",
        validate="many_to_one",
    )

    unresolved = joined["par_checklist_status"].isna()
    failed = joined["par_checklist_status"].notna() & ~joined[
        "par_checklist_status"
    ].isin({AVAILABLE_STATUS, NO_DATA_STATUS})
    if unresolved.any() or failed.any():
        problem_rows = joined.loc[
            unresolved | failed,
            [
                "aggregation_id",
                "Delta_Date",
                "par_checklist_status",
                "par_checklist_error",
            ],
        ].head()
        raise ValueError(
            "Some window dates have a missing or unresolved PAR checklist "
            f"record. First rows:\n{problem_rows}"
        )

    no_granules = joined["par_checklist_status"].eq(NO_DATA_STATUS)
    retained_manifest = joined.loc[~no_granules].copy()
    excluded_windows = joined.loc[no_granules].copy()

    retained_ids = set(retained_manifest["aggregation_id"])
    retained_assignments = assignments[
        assignments["aggregation_id"].isin(retained_ids)
    ].copy()
    validate_retained_tables(retained_manifest, retained_assignments)

    retained_manifest = retained_manifest.sort_values(
        ["Delta_Date", "mgrs_tile_t", "measurement_mode", "aggregation_id"]
    ).reset_index(drop=True)
    retained_assignments = retained_assignments.sort_values(
        ["aggregation_id", "sif_row_id"]
    ).reset_index(drop=True)
    excluded_windows = excluded_windows.sort_values(
        ["Delta_Date", "mgrs_tile_t", "measurement_mode", "aggregation_id"]
    ).reset_index(drop=True)

    summary = pd.DataFrame(
        [
            {
                "input_windows": len(manifest),
                "input_unique_dates": manifest["Delta_Date"].nunique(),
                "retained_windows": len(retained_manifest),
                "retained_unique_dates": retained_manifest["Delta_Date"].nunique(),
                "excluded_no_granules_windows": len(excluded_windows),
                "excluded_no_granules_unique_dates": excluded_windows[
                    "Delta_Date"
                ].nunique(),
                "retained_assignments": len(retained_assignments),
                "retained_window_fraction": len(retained_manifest) / len(manifest),
            }
        ]
    )
    return retained_manifest, retained_assignments, excluded_windows, summary


def main() -> None:
    manifest, assignments, checklist = load_inputs()
    retained_manifest, retained_assignments, excluded_windows, summary = (
        filter_by_par_availability(manifest, assignments, checklist)
    )

    atomic_write_csv(retained_manifest, OUTPUT_MANIFEST_PATH)
    atomic_write_csv(retained_assignments, OUTPUT_ASSIGNMENTS_PATH)
    atomic_write_csv(excluded_windows, EXCLUDED_WINDOWS_PATH)
    atomic_write_csv(summary, FILTER_SUMMARY_PATH)

    record = summary.iloc[0]
    print(f"Input windows: {int(record['input_windows']):,}")
    print(
        "Excluded windows on confirmed no-granules dates: "
        f"{int(record['excluded_no_granules_windows']):,} across "
        f"{int(record['excluded_no_granules_unique_dates']):,} dates"
    )
    print(
        f"Retained PAR-available windows: {int(record['retained_windows']):,} "
        f"({float(record['retained_window_fraction']):.1%})"
    )
    print(f"Retained footprint assignments: {len(retained_assignments):,}")
    print(f"Filtered manifest: {OUTPUT_MANIFEST_PATH}")
    print(f"Filtered assignments: {OUTPUT_ASSIGNMENTS_PATH}")


if __name__ == "__main__":
    main()
