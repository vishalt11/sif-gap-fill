"""Combine density-based 4 km datasets and retain clear Sentinel-2 windows.

The four source aggregation folders and downloaded GeoTIFFs are never changed.
This script creates one derived manifest and assignment table for windows whose
Sentinel-2 L2A raster had at least 95% valid pixels before neighbour filling.
"""

from __future__ import annotations

import os
from pathlib import Path
import re
from typing import Iterable

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Configuration

PROJECT_ROOT = Path(__file__).resolve().parent
DENSITY_ROOT = PROJECT_ROOT / "data" / "density_aggregation"

SOURCE_DIRS = {
    "original_9tiles": DENSITY_ROOT
    / "sentinel2_spatial_aggregation_density_4000m_landcover_mask60_min4",
    "redtiles": DENSITY_ROOT
    / "sentinel2_spatial_aggregation_density_4000m_landcover_redtiles_mask60_min4",
    "yellowtiles_32": DENSITY_ROOT
    / "sentinel2_spatial_aggregation_density_4000m_landcover_yellowtiles_32_mask60_min4",
    "yellowtiles_33": DENSITY_ROOT
    / "sentinel2_spatial_aggregation_density_4000m_landcover_yellowtiles_33_mask60_min4",
}

DOWNLOAD_MANIFEST_PATH = (
    PROJECT_ROOT
    / "data"
    / "sentinel2_l2a"
    / "sentinel2_l2a_download_manifest.csv"
)

OUTPUT_DIR = (
    DENSITY_ROOT
    / "sentinel2_spatial_aggregation_density_4000m_landcover_combined_"
    "mask60_min4_s2valid95"
)
OUTPUT_MANIFEST_PATH = (
    OUTPUT_DIR / "density_cluster_4000m_aggregate_manifest.csv"
)
OUTPUT_ASSIGNMENTS_PATH = (
    OUTPUT_DIR / "density_cluster_4000m_sif_assignments.csv"
)
EXCLUDED_WINDOWS_PATH = (
    OUTPUT_DIR / "density_cluster_4000m_sentinel2_excluded_windows.csv"
)
FILTER_SUMMARY_PATH = OUTPUT_DIR / "sentinel2_quality_filter_summary.csv"
EXCLUSION_SUMMARY_PATH = OUTPUT_DIR / "sentinel2_exclusion_reason_summary.csv"

MIN_VALID_FRACTION_BEFORE_FILL = 0.95
WINDOW_SIZE_M = 4000.0
WINDOW_PIXELS = 200
MIN_FOOTPRINTS = 4
SUCCESS_STATUSES = {"completed", "skipped_existing"}

MANIFEST_FILENAME = "density_cluster_4000m_aggregate_manifest.csv"
ASSIGNMENTS_FILENAME = "density_cluster_4000m_sif_assignments.csv"

MANIFEST_REQUIRED_COLUMNS = [
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
    "measurement_mode",
    "n_footprints",
    "aggregated_target_modis_sif",
]

ASSIGNMENT_REQUIRED_COLUMNS = [
    "aggregation_id",
    "sif_row_id",
    "source_csv_row",
    "Delta_Date",
    "mgrs_tile_t",
    "measurement_mode",
    "target_modis_sif",
    "Lat_corner1",
    "Lat_corner2",
    "Lat_corner3",
    "Lat_corner4",
    "Lon_corner1",
    "Lon_corner2",
    "Lon_corner3",
    "Lon_corner4",
]

DOWNLOAD_REQUIRED_COLUMNS = [
    "aggregation_id",
    "target_date",
    "mgrs_tile",
    "output_path",
    "status",
    "error",
    "valid_fraction_before_fill",
    "invalid_pixels_before_fill",
    "neighbour_fill_iterations",
    "neighbour_fill_complete",
    "source_date_counts",
]


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


def parse_bool(series: pd.Series) -> pd.Series:
    return (
        series.astype("string")
        .str.strip()
        .str.lower()
        .isin({"true", "t", "1", "yes"})
    )


def normalize_mgrs_tile(series: pd.Series) -> pd.Series:
    tile = (
        series.astype("string")
        .str.strip()
        .str.upper()
        .str.replace(r"^T", "", regex=True)
    )
    valid = tile.str.fullmatch(r"[0-9]{2}[C-X][A-Z]{2}", na=False)
    return ("T" + tile).where(valid)


def crs_from_mgrs_tile(tile: str) -> str:
    match = re.fullmatch(r"T?([0-9]{2})([C-X])[A-Z]{2}", str(tile).upper())
    if match is None:
        raise ValueError(f"Invalid MGRS tile: {tile}")

    zone = int(match.group(1))
    latitude_band = match.group(2)
    epsg = (32600 if latitude_band >= "N" else 32700) + zone
    return f"EPSG:{epsg}"


def project_relative_path(value: object) -> tuple[str, bool]:
    if pd.isna(value) or not str(value).strip():
        return "", False

    path = Path(str(value).strip())
    if not path.is_absolute():
        path = PROJECT_ROOT / path
    resolved = path.resolve()

    try:
        stored_path = resolved.relative_to(PROJECT_ROOT).as_posix()
    except ValueError:
        stored_path = str(resolved)
    return stored_path, resolved.is_file()


def atomic_write_csv(table: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = path.with_suffix(path.suffix + ".partial")
    try:
        table.to_csv(temporary_path, index=False)
        os.replace(temporary_path, path)
    finally:
        if temporary_path.exists():
            temporary_path.unlink()


def add_source_metadata(
    table: pd.DataFrame,
    source_name: str,
    source_dir: Path,
) -> pd.DataFrame:
    table = table.copy()
    table["source_dataset"] = source_name
    table["source_aggregation_dir"] = source_dir.relative_to(
        PROJECT_ROOT
    ).as_posix()
    leading = ["source_dataset", "source_aggregation_dir"]
    return table[leading + [column for column in table.columns if column not in leading]]


def add_and_validate_window_crs(table: pd.DataFrame, table_name: str) -> pd.DataFrame:
    table = table.copy()
    normalized_tiles = normalize_mgrs_tile(table["mgrs_tile_t"])
    if normalized_tiles.isna().any():
        bad = table.loc[normalized_tiles.isna(), "mgrs_tile_t"].head().tolist()
        raise ValueError(f"{table_name} has invalid MGRS tile IDs: {bad}")
    table["mgrs_tile_t"] = normalized_tiles

    expected_crs = table["mgrs_tile_t"].map(crs_from_mgrs_tile)
    if "window_crs" not in table.columns:
        table["window_crs"] = expected_crs
        return table

    existing_crs = (
        table["window_crs"]
        .astype("string")
        .str.strip()
        .str.upper()
        .replace("", pd.NA)
    )
    mismatch = existing_crs.notna() & (existing_crs != expected_crs)
    if mismatch.any():
        columns = ["aggregation_id", "mgrs_tile_t", "window_crs"]
        raise ValueError(
            f"{table_name} has a window_crs/MGRS mismatch. First rows:\n"
            f"{table.loc[mismatch, columns].head()}"
        )
    table["window_crs"] = existing_crs.fillna(expected_crs)
    return table


# ---------------------------------------------------------------------------
# Source aggregation tables


def read_source_tables(
    source_name: str,
    source_dir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    manifest_path = source_dir / MANIFEST_FILENAME
    assignments_path = source_dir / ASSIGNMENTS_FILENAME
    if not manifest_path.is_file():
        raise FileNotFoundError(manifest_path)
    if not assignments_path.is_file():
        raise FileNotFoundError(assignments_path)

    manifest = pd.read_csv(manifest_path, low_memory=False)
    assignments = pd.read_csv(assignments_path, low_memory=False)
    require_columns(manifest, MANIFEST_REQUIRED_COLUMNS, f"{source_name} manifest")
    require_columns(
        assignments,
        ASSIGNMENT_REQUIRED_COLUMNS,
        f"{source_name} assignments",
    )

    manifest["aggregation_id"] = manifest["aggregation_id"].astype(str)
    assignments["aggregation_id"] = assignments["aggregation_id"].astype(str)
    if manifest["aggregation_id"].duplicated().any():
        duplicate = manifest.loc[
            manifest["aggregation_id"].duplicated(), "aggregation_id"
        ].iloc[0]
        raise ValueError(f"Duplicate aggregation_id in {source_name}: {duplicate}")

    for column in ["aggregation_size_m", "cell_pixels", "n_footprints"]:
        manifest[column] = pd.to_numeric(manifest[column], errors="raise")
    manifest["aggregated_target_modis_sif"] = pd.to_numeric(
        manifest["aggregated_target_modis_sif"], errors="raise"
    )
    assignments["target_modis_sif"] = pd.to_numeric(
        assignments["target_modis_sif"], errors="raise"
    )

    invalid_windows = (
        ~np.isclose(manifest["aggregation_size_m"], WINDOW_SIZE_M)
        | (manifest["cell_pixels"].astype(int) != WINDOW_PIXELS)
        | (manifest["n_footprints"].astype(int) < MIN_FOOTPRINTS)
        | ~np.isfinite(manifest["aggregated_target_modis_sif"])
    )
    if invalid_windows.any():
        raise ValueError(
            f"{source_name} contains {int(invalid_windows.sum())} invalid 4 km windows"
        )

    unknown_assignment_ids = set(assignments["aggregation_id"]) - set(
        manifest["aggregation_id"]
    )
    if unknown_assignment_ids:
        raise ValueError(
            f"{source_name} assignments contain IDs absent from its manifest"
        )

    expected_counts = manifest.set_index("aggregation_id")["n_footprints"].astype(int)
    found_counts = assignments.groupby("aggregation_id").size()
    count_check = expected_counts.to_frame("expected").join(
        found_counts.rename("found")
    )
    mismatched = count_check[
        count_check["found"].isna()
        | (count_check["expected"] != count_check["found"])
    ]
    if not mismatched.empty:
        raise ValueError(
            f"{source_name} manifest/assignment counts differ. First rows:\n"
            f"{mismatched.head()}"
        )

    if "passes_window_coverage_filter" in assignments.columns:
        passing = parse_bool(assignments["passes_window_coverage_filter"])
        if not passing.all():
            raise ValueError(
                f"{source_name} contains assignments that failed the 60% filter"
            )

    manifest = add_and_validate_window_crs(manifest, f"{source_name} manifest")
    assignments = add_and_validate_window_crs(
        assignments,
        f"{source_name} assignments",
    )
    manifest = add_source_metadata(manifest, source_name, source_dir)
    assignments = add_source_metadata(assignments, source_name, source_dir)
    return manifest, assignments


def combine_source_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    manifests: list[pd.DataFrame] = []
    assignments: list[pd.DataFrame] = []
    for source_name, source_dir in SOURCE_DIRS.items():
        source_manifest, source_assignments = read_source_tables(
            source_name,
            source_dir,
        )
        manifests.append(source_manifest)
        assignments.append(source_assignments)

    combined_manifest = pd.concat(manifests, ignore_index=True, sort=False)
    combined_assignments = pd.concat(assignments, ignore_index=True, sort=False)

    duplicated = combined_manifest["aggregation_id"].duplicated(keep=False)
    if duplicated.any():
        columns = ["aggregation_id", "source_dataset", "mgrs_tile_t", "Delta_Date"]
        raise ValueError(
            "Aggregation IDs overlap across source datasets. First rows:\n"
            f"{combined_manifest.loc[duplicated, columns].head()}"
        )

    assignment_keys = ["aggregation_id", "sif_row_id"]
    if combined_assignments.duplicated(assignment_keys).any():
        raise ValueError("Duplicate aggregation_id/sif_row_id assignment keys found")

    return combined_manifest, combined_assignments


# ---------------------------------------------------------------------------
# Sentinel-2 quality join and filtering


def read_download_manifest() -> pd.DataFrame:
    if not DOWNLOAD_MANIFEST_PATH.is_file():
        raise FileNotFoundError(DOWNLOAD_MANIFEST_PATH)

    download = pd.read_csv(DOWNLOAD_MANIFEST_PATH, low_memory=False)
    require_columns(download, DOWNLOAD_REQUIRED_COLUMNS, "Sentinel-2 download manifest")
    download["aggregation_id"] = download["aggregation_id"].astype(str)
    if download["aggregation_id"].duplicated().any():
        duplicate = download.loc[
            download["aggregation_id"].duplicated(), "aggregation_id"
        ].iloc[0]
        raise ValueError(f"Duplicate aggregation_id in download manifest: {duplicate}")

    rename = {
        column: f"sentinel2_{column}"
        for column in download.columns
        if column != "aggregation_id"
    }
    return download.rename(columns=rename)


def add_exclusion_reasons(table: pd.DataFrame) -> pd.DataFrame:
    table = table.copy()
    matched = table["download_record_found"]
    valid_fraction = table["sentinel2_valid_fraction_before_fill"]

    reason_masks = {
        "missing_download_record": ~matched,
        "unsuccessful_download_status": matched & ~table["successful_status"],
        "download_error_recorded": matched & ~table["download_error_is_empty"],
        "invalid_valid_fraction": matched
        & (~np.isfinite(valid_fraction) | ~valid_fraction.between(0.0, 1.0)),
        "valid_fraction_below_0.95": matched
        & valid_fraction.between(0.0, 1.0)
        & (valid_fraction < MIN_VALID_FRACTION_BEFORE_FILL),
        "neighbour_fill_incomplete": matched
        & ~table["sentinel2_neighbour_fill_complete_bool"],
        "sentinel2_output_file_missing": matched & ~table["sentinel2_file_exists"],
        "target_date_mismatch": matched & ~table["sentinel2_target_date_matches"],
        "mgrs_tile_mismatch": matched & ~table["sentinel2_mgrs_tile_matches"],
    }

    reasons: list[str] = []
    for row_index in table.index:
        reasons.append(
            ";".join(
                reason
                for reason, mask in reason_masks.items()
                if bool(mask.loc[row_index])
            )
        )
    table["sentinel2_exclusion_reason"] = reasons
    return table


def join_and_filter_sentinel2(
    manifest: pd.DataFrame,
    download: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    joined = manifest.merge(
        download,
        on="aggregation_id",
        how="left",
        validate="one_to_one",
        indicator=True,
    )
    joined["download_record_found"] = joined["_merge"].eq("both")
    joined = joined.drop(columns="_merge")

    joined["sentinel2_status"] = (
        joined["sentinel2_status"].astype("string").str.strip().str.lower()
    )
    joined["successful_status"] = joined["sentinel2_status"].isin(
        SUCCESS_STATUSES
    )
    joined["download_error_is_empty"] = (
        joined["sentinel2_error"].fillna("").astype(str).str.strip().eq("")
    )
    joined["sentinel2_valid_fraction_before_fill"] = pd.to_numeric(
        joined["sentinel2_valid_fraction_before_fill"], errors="coerce"
    )
    joined["sentinel2_neighbour_fill_complete_bool"] = parse_bool(
        joined["sentinel2_neighbour_fill_complete"]
    )

    resolved_paths = joined["sentinel2_output_path"].map(project_relative_path)
    joined["sentinel2_raster_path"] = resolved_paths.map(lambda value: value[0])
    joined["sentinel2_file_exists"] = resolved_paths.map(lambda value: value[1])

    aggregate_dates = pd.to_datetime(joined["Delta_Date"], errors="coerce").dt.date
    sentinel_dates = pd.to_datetime(
        joined["sentinel2_target_date"], errors="coerce"
    ).dt.date
    joined["sentinel2_target_date_matches"] = (
        aggregate_dates.eq(sentinel_dates).fillna(False).astype(bool)
    )

    aggregate_tiles = normalize_mgrs_tile(joined["mgrs_tile_t"])
    sentinel_tiles = normalize_mgrs_tile(joined["sentinel2_mgrs_tile"])
    joined["sentinel2_mgrs_tile_matches"] = (
        aggregate_tiles.eq(sentinel_tiles).fillna(False).astype(bool)
    )

    joined = add_exclusion_reasons(joined)
    joined["passes_sentinel2_quality_filter"] = joined[
        "sentinel2_exclusion_reason"
    ].eq("")

    retained = joined[joined["passes_sentinel2_quality_filter"]].copy()
    excluded = joined[~joined["passes_sentinel2_quality_filter"]].copy()
    return retained, excluded


# ---------------------------------------------------------------------------
# Output validation and summaries


def filter_and_validate_assignments(
    retained_manifest: pd.DataFrame,
    combined_assignments: pd.DataFrame,
) -> pd.DataFrame:
    retained_ids = set(retained_manifest["aggregation_id"])
    retained_assignments = combined_assignments[
        combined_assignments["aggregation_id"].isin(retained_ids)
    ].copy()

    found_ids = set(retained_assignments["aggregation_id"])
    if found_ids != retained_ids:
        raise ValueError("At least one retained window has no assignment rows")

    expected_counts = (
        retained_manifest.set_index("aggregation_id")["n_footprints"].astype(int)
    )
    found_counts = retained_assignments.groupby("aggregation_id").size()
    if not expected_counts.sort_index().equals(found_counts.sort_index()):
        raise ValueError("Retained manifest and assignment footprint counts differ")

    assignment_means = retained_assignments.groupby("aggregation_id")[
        "target_modis_sif"
    ].mean()
    manifest_means = retained_manifest.set_index("aggregation_id")[
        "aggregated_target_modis_sif"
    ]
    if not np.allclose(
        assignment_means.sort_index().to_numpy(dtype=float),
        manifest_means.sort_index().to_numpy(dtype=float),
        rtol=1e-6,
        atol=1e-7,
    ):
        raise ValueError("Retained aggregate targets do not match assignment means")

    return retained_assignments


def summarize_filter(joined: pd.DataFrame) -> pd.DataFrame:
    records: list[dict] = []

    def append_summary(grouping: str, value: object, group: pd.DataFrame) -> None:
        retained_n = int(group["passes_sentinel2_quality_filter"].sum())
        records.append(
            {
                "grouping": grouping,
                "group_value": str(value),
                "input_windows": len(group),
                "retained_windows": retained_n,
                "excluded_windows": len(group) - retained_n,
                "retained_fraction": retained_n / len(group),
            }
        )

    append_summary("overall", "all", joined)
    for column in [
        "source_dataset",
        "mgrs_tile_t",
        "window_crs",
        "sif_year",
        "sif_month",
        "majority_land_cover_class",
    ]:
        if column not in joined.columns:
            continue
        for value, group in joined.groupby(column, dropna=False, sort=True):
            append_summary(column, value, group)

    return pd.DataFrame(records)


def exclusion_reason_summary(excluded: pd.DataFrame) -> pd.DataFrame:
    if excluded.empty:
        return pd.DataFrame(columns=["exclusion_reason", "n_windows"])
    return (
        excluded["sentinel2_exclusion_reason"]
        .str.split(";")
        .explode()
        .value_counts()
        .rename_axis("exclusion_reason")
        .reset_index(name="n_windows")
    )


def sort_outputs(
    manifest: pd.DataFrame,
    assignments: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    manifest = manifest.sort_values(
        ["Delta_Date", "mgrs_tile_t", "measurement_mode", "aggregation_id"]
    ).reset_index(drop=True)
    assignments = assignments.sort_values(
        ["aggregation_id", "sif_row_id"]
    ).reset_index(drop=True)
    return manifest, assignments


def main() -> None:
    combined_manifest, combined_assignments = combine_source_tables()
    download_manifest = read_download_manifest()
    retained_manifest, excluded_windows = join_and_filter_sentinel2(
        combined_manifest,
        download_manifest,
    )
    retained_assignments = filter_and_validate_assignments(
        retained_manifest,
        combined_assignments,
    )

    all_joined = pd.concat(
        [retained_manifest, excluded_windows],
        ignore_index=True,
        sort=False,
    )
    filter_summary = summarize_filter(all_joined)
    reason_summary = exclusion_reason_summary(excluded_windows)
    retained_manifest, retained_assignments = sort_outputs(
        retained_manifest,
        retained_assignments,
    )

    atomic_write_csv(retained_manifest, OUTPUT_MANIFEST_PATH)
    atomic_write_csv(retained_assignments, OUTPUT_ASSIGNMENTS_PATH)
    atomic_write_csv(excluded_windows, EXCLUDED_WINDOWS_PATH)
    atomic_write_csv(filter_summary, FILTER_SUMMARY_PATH)
    atomic_write_csv(reason_summary, EXCLUSION_SUMMARY_PATH)

    input_windows = len(combined_manifest)
    retained_windows = len(retained_manifest)
    print(f"Combined source windows: {input_windows:,}")
    print(
        f"Retained Sentinel-2 windows (valid fraction >= "
        f"{MIN_VALID_FRACTION_BEFORE_FILL:.0%}): {retained_windows:,} "
        f"({retained_windows / input_windows:.1%})"
    )
    print(f"Excluded windows: {input_windows - retained_windows:,}")
    print(f"Retained footprint assignments: {len(retained_assignments):,}")
    print(f"Combined manifest: {OUTPUT_MANIFEST_PATH}")
    print(f"Combined assignments: {OUTPUT_ASSIGNMENTS_PATH}")


if __name__ == "__main__":
    main()
