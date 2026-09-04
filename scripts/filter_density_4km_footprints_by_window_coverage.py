"""Filter density-based 4 km SIF windows by footprint coverage.

This is a lightweight stage between density clustering and CNN chip creation.
It reads only the density-cluster manifest and footprint-assignment CSV files;
no Sentinel-2 or other predictor rasters are opened.

For every assigned OCO-2 footprint, the script reconstructs its polygon from
the four corner coordinates, projects it to EPSG:32633, and rasterizes it on
the same 200 x 200, 20 m grid used by the 4 km CNN chips. Footprints with less
than 60% of their area represented inside the assigned window are removed.
Windows are retained when at least four passing footprints remain.

All footprint-dependent manifest statistics and the aggregate SIF target are
recalculated from the retained assignments. The original density-clustering
folder is never modified. Filtered inputs and audit tables are written to a
new output folder.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from pyproj import Transformer
from rasterio.features import rasterize
from rasterio.transform import from_origin
from shapely.geometry import MultiPoint, Polygon, mapping
from shapely.ops import transform


# ---------------------------------------------------------------------------
# Configuration

INPUT_DIR = Path(
    "data/density_aggregation/"
    "sentinel2_spatial_aggregation_density_4000m_landcover_yellowtiles_33"
)
OUTPUT_DIR = Path(
    "data/density_aggregation/"
    "sentinel2_spatial_aggregation_density_4000m_landcover_yellowtiles_33_mask60_min4"
)

INPUT_MANIFEST_PATH = (
    INPUT_DIR / "density_cluster_4000m_aggregate_manifest.csv"
)
INPUT_ASSIGNMENTS_PATH = (
    INPUT_DIR / "density_cluster_4000m_sif_assignments.csv"
)

# These filenames match the original density output names so the filtered
# folder can be used as a direct input source by the later CNN preparation.
OUTPUT_MANIFEST_PATH = (
    OUTPUT_DIR / "density_cluster_4000m_aggregate_manifest.csv"
)
OUTPUT_ASSIGNMENTS_PATH = (
    OUTPUT_DIR / "density_cluster_4000m_sif_assignments.csv"
)
REMOVED_FOOTPRINTS_PATH = (
    OUTPUT_DIR / "density_cluster_4000m_removed_footprints.csv"
)
REMOVED_WINDOWS_PATH = (
    OUTPUT_DIR / "density_cluster_4000m_removed_windows.csv"
)
FILTER_SUMMARY_PATH = OUTPUT_DIR / "coverage_filter_summary.csv"

WINDOW_SIZE_M = 4000.0
PIXEL_SIZE_M = 20.0
WINDOW_PIXELS = 200
MASK_OVERSAMPLE = 4

MIN_MASK_INSIDE_FRACTION = 0.60
MIN_RETAINED_FOOTPRINTS = 4

WINDOW_CRS = "EPSG:32633"
WGS84_CRS = "EPSG:4326"
PROGRESS_EVERY_WINDOWS = 250


# ---------------------------------------------------------------------------
# Validation and geometry helpers


def require_columns(
    table: pd.DataFrame,
    columns: Iterable[str],
    table_name: str,
) -> None:
    missing = [column for column in columns if column not in table.columns]
    if missing:
        raise ValueError(f"{table_name} is missing required columns: {missing}")


def build_footprint_wgs84(row: pd.Series) -> Polygon:
    coordinates = [
        (float(row["Lon_corner1"]), float(row["Lat_corner1"])),
        (float(row["Lon_corner2"]), float(row["Lat_corner2"])),
        (float(row["Lon_corner3"]), float(row["Lat_corner3"])),
        (float(row["Lon_corner4"]), float(row["Lat_corner4"])),
    ]
    polygon = Polygon(coordinates)
    if not polygon.is_valid or polygon.area <= 0:
        polygon = MultiPoint(coordinates).convex_hull
    if polygon.is_empty or polygon.area <= 0:
        raise ValueError(
            "Could not construct footprint polygon for "
            f"sif_row_id={row['sif_row_id']}"
        )
    return polygon


def rasterized_inside_fraction(
    footprint_projected,
    cell_xmin: float,
    cell_ymax: float,
) -> tuple[float, float, float]:
    """Return window coverage fraction, inside area, and full polygon area."""
    high_pixel_size = PIXEL_SIZE_M / MASK_OVERSAMPLE
    high_pixels = WINDOW_PIXELS * MASK_OVERSAMPLE
    high_transform = from_origin(
        cell_xmin,
        cell_ymax,
        high_pixel_size,
        high_pixel_size,
    )
    high_mask = rasterize(
        [(mapping(footprint_projected), 1)],
        out_shape=(high_pixels, high_pixels),
        transform=high_transform,
        fill=0,
        dtype="uint8",
        all_touched=False,
    )

    full_area_m2 = float(footprint_projected.area)
    inside_area_m2 = float(high_mask.sum(dtype=np.float64)) * high_pixel_size**2
    if not np.isfinite(full_area_m2) or full_area_m2 <= 0:
        return 0.0, inside_area_m2, full_area_m2

    inside_fraction = float(np.clip(inside_area_m2 / full_area_m2, 0.0, 1.0))
    return inside_fraction, inside_area_m2, full_area_m2


def collapse_values(series: pd.Series) -> str:
    values = [str(value) for value in series.dropna().unique()]
    return ";".join(sorted(values))


def modal_value(series: pd.Series):
    values = series.dropna()
    if values.empty:
        return np.nan
    counts = values.value_counts()
    return counts.index[0]


def safe_numeric_summary(series: pd.Series, operation: str) -> float:
    values = pd.to_numeric(series, errors="coerce").dropna()
    if values.empty:
        return np.nan
    if operation == "mean":
        return float(values.mean())
    if operation == "min":
        return float(values.min())
    if operation == "max":
        return float(values.max())
    if operation == "sum":
        return float(values.sum())
    raise ValueError(f"Unsupported summary operation: {operation}")


# ---------------------------------------------------------------------------
# Input loading


def load_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    if not INPUT_MANIFEST_PATH.exists():
        raise FileNotFoundError(INPUT_MANIFEST_PATH)
    if not INPUT_ASSIGNMENTS_PATH.exists():
        raise FileNotFoundError(INPUT_ASSIGNMENTS_PATH)

    manifest = pd.read_csv(INPUT_MANIFEST_PATH)
    assignments = pd.read_csv(INPUT_ASSIGNMENTS_PATH)

    require_columns(
        manifest,
        [
            "aggregation_id",
            "aggregation_size_m",
            "cell_pixels",
            "cell_xmin",
            "cell_ymin",
            "cell_xmax",
            "cell_ymax",
            "n_footprints",
            "aggregated_target_modis_sif",
        ],
        "Density aggregate manifest",
    )
    require_columns(
        assignments,
        [
            "aggregation_id",
            "sif_row_id",
            "source_csv_row",
            "target_modis_sif",
            "Lat_corner1",
            "Lat_corner2",
            "Lat_corner3",
            "Lat_corner4",
            "Lon_corner1",
            "Lon_corner2",
            "Lon_corner3",
            "Lon_corner4",
        ],
        "Density SIF assignments",
    )

    manifest_numeric = [
        "aggregation_size_m",
        "cell_pixels",
        "cell_xmin",
        "cell_ymin",
        "cell_xmax",
        "cell_ymax",
        "n_footprints",
        "aggregated_target_modis_sif",
    ]
    assignment_numeric = [
        "sif_row_id",
        "source_csv_row",
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
    for column in manifest_numeric:
        manifest[column] = pd.to_numeric(manifest[column], errors="raise")
    for column in assignment_numeric:
        assignments[column] = pd.to_numeric(assignments[column], errors="raise")

    manifest["aggregation_id"] = manifest["aggregation_id"].astype(str)
    assignments["aggregation_id"] = assignments["aggregation_id"].astype(str)

    manifest = manifest[
        np.isclose(manifest["aggregation_size_m"], WINDOW_SIZE_M)
        & (manifest["cell_pixels"].astype(int) == WINDOW_PIXELS)
        & (manifest["n_footprints"].astype(int) >= MIN_RETAINED_FOOTPRINTS)
    ].copy()
    if manifest.empty:
        raise ValueError("No eligible 4 km density windows were found.")
    if manifest["aggregation_id"].duplicated().any():
        duplicate_id = manifest.loc[
            manifest["aggregation_id"].duplicated(), "aggregation_id"
        ].iloc[0]
        raise ValueError(f"Duplicate aggregation_id in manifest: {duplicate_id}")

    assignments = assignments[
        assignments["aggregation_id"].isin(manifest["aggregation_id"])
    ].copy()

    expected_counts = manifest.set_index("aggregation_id")["n_footprints"].astype(int)
    found_counts = assignments.groupby("aggregation_id").size()
    count_check = expected_counts.to_frame("expected").join(
        found_counts.rename("found")
    )
    mismatched = count_check[count_check["expected"] != count_check["found"]]
    if not mismatched.empty:
        raise ValueError(
            "Manifest/assignment footprint count mismatch. First rows:\n"
            f"{mismatched.head()}"
        )

    return (
        manifest.sort_values("aggregation_id").reset_index(drop=True),
        assignments.sort_values(["aggregation_id", "sif_row_id"]).reset_index(
            drop=True
        ),
    )


# ---------------------------------------------------------------------------
# Coverage filtering and manifest rebuilding


def add_coverage_measurements(
    manifest: pd.DataFrame,
    assignments: pd.DataFrame,
) -> pd.DataFrame:
    transformer = Transformer.from_crs(
        WGS84_CRS,
        WINDOW_CRS,
        always_xy=True,
    )
    window_bounds = manifest.set_index("aggregation_id")[[
        "cell_xmin",
        "cell_ymin",
        "cell_xmax",
        "cell_ymax",
    ]]

    output_rows: list[dict] = []
    total_windows = len(manifest)
    grouped = assignments.groupby("aggregation_id", sort=False)

    print(
        f"Filtering {total_windows:,} density-based 4 km windows and "
        f"{len(assignments):,} assigned footprints"
    )
    print(
        f"Coverage threshold: >= {MIN_MASK_INSIDE_FRACTION:.0%}; "
        f"minimum retained footprints per window: {MIN_RETAINED_FOOTPRINTS}"
    )

    for window_index, (aggregation_id, footprints) in enumerate(grouped, start=1):
        bounds = window_bounds.loc[aggregation_id]
        cell_xmin = float(bounds["cell_xmin"])
        cell_ymin = float(bounds["cell_ymin"])
        cell_xmax = float(bounds["cell_xmax"])
        cell_ymax = float(bounds["cell_ymax"])
        width = cell_xmax - cell_xmin
        height = cell_ymax - cell_ymin
        if not np.isclose(width, WINDOW_SIZE_M) or not np.isclose(
            height, WINDOW_SIZE_M
        ):
            raise ValueError(
                f"Window {aggregation_id} is {width} x {height} m; expected "
                f"{WINDOW_SIZE_M} x {WINDOW_SIZE_M} m."
            )

        for _, footprint in footprints.iterrows():
            footprint_wgs84 = build_footprint_wgs84(footprint)
            footprint_projected = transform(
                transformer.transform,
                footprint_wgs84,
            )
            inside_fraction, inside_area_m2, full_area_m2 = (
                rasterized_inside_fraction(
                    footprint_projected,
                    cell_xmin,
                    cell_ymax,
                )
            )
            row = footprint.to_dict()
            row.update(
                {
                    "window_mask_inside_fraction": inside_fraction,
                    "window_mask_inside_pct": inside_fraction * 100.0,
                    "window_rasterized_inside_area_km2": (
                        inside_area_m2 / 1_000_000.0
                    ),
                    "window_projected_full_area_km2": (
                        full_area_m2 / 1_000_000.0
                    ),
                    "passes_window_coverage_filter": bool(
                        inside_fraction >= MIN_MASK_INSIDE_FRACTION
                    ),
                }
            )
            output_rows.append(row)

        if window_index == 1 or window_index % PROGRESS_EVERY_WINDOWS == 0:
            print(f"Measured window {window_index:,} / {total_windows:,}")

    measured = pd.DataFrame(output_rows)
    return measured.sort_values(["aggregation_id", "sif_row_id"]).reset_index(
        drop=True
    )


def summarize_retained_assignments(
    group: pd.DataFrame,
    original_n: int,
) -> dict:
    target = pd.to_numeric(group["target_modis_sif"], errors="raise")
    n_footprints = len(group)
    result: dict = {
        "original_n_footprints": int(original_n),
        "n_footprints": int(n_footprints),
        "removed_by_window_coverage_n": int(original_n - n_footprints),
        "retained_footprint_fraction": float(n_footprints / original_n),
        "aggregated_target_modis_sif": float(target.mean()),
        "median_target_modis_sif": float(target.median()),
        "min_target_modis_sif": float(target.min()),
        "max_target_modis_sif": float(target.max()),
        "sd_target_modis_sif": (
            float(target.std(ddof=1)) if n_footprints > 1 else np.nan
        ),
        "se_target_modis_sif": (
            float(target.std(ddof=1) / np.sqrt(n_footprints))
            if n_footprints > 1
            else np.nan
        ),
        "mean_window_mask_inside_fraction": float(
            group["window_mask_inside_fraction"].mean()
        ),
        "min_window_mask_inside_fraction": float(
            group["window_mask_inside_fraction"].min()
        ),
        "max_window_mask_inside_fraction": float(
            group["window_mask_inside_fraction"].max()
        ),
        "sif_row_ids": ",".join(group["sif_row_id"].astype(int).astype(str)),
        "source_csv_rows": ",".join(
            group["source_csv_row"].astype(int).astype(str)
        ),
        "eligible_n3": bool(n_footprints >= 3),
        "eligible_n4": bool(n_footprints >= 4),
        "eligible_n5": bool(n_footprints >= 5),
        "eligible_n6": bool(n_footprints >= 6),
        "eligible_n8": bool(n_footprints >= 8),
        "eligible_n10": bool(n_footprints >= 10),
    }

    if "sif_area_km2_evi" in group.columns:
        result.update(
            {
                "mean_sif_area_km2": safe_numeric_summary(
                    group["sif_area_km2_evi"], "mean"
                ),
                "min_sif_area_km2": safe_numeric_summary(
                    group["sif_area_km2_evi"], "min"
                ),
                "max_sif_area_km2": safe_numeric_summary(
                    group["sif_area_km2_evi"], "max"
                ),
                "total_sif_area_km2": safe_numeric_summary(
                    group["sif_area_km2_evi"], "sum"
                ),
            }
        )
    if "tile_inside_fraction" in group.columns:
        result.update(
            {
                "mean_tile_inside_fraction": safe_numeric_summary(
                    group["tile_inside_fraction"], "mean"
                ),
                "min_tile_inside_fraction": safe_numeric_summary(
                    group["tile_inside_fraction"], "min"
                ),
            }
        )
    if "state" in group.columns:
        result["states"] = collapse_values(group["state"])
    if "hzs" in group.columns:
        result["hzs_values"] = collapse_values(group["hzs"])
    if "BKR10_ID" in group.columns:
        result.update(
            {
                "BKR10_ID_values": collapse_values(group["BKR10_ID"]),
                "majority_BKR10_ID": modal_value(group["BKR10_ID"]),
                "n_BKR10": int(group["BKR10_ID"].nunique(dropna=True)),
            }
        )
    if "BKR_NAME" in group.columns:
        result.update(
            {
                "BKR_NAME_values": collapse_values(group["BKR_NAME"]),
                "majority_BKR_NAME": modal_value(group["BKR_NAME"]),
            }
        )
    if "Quality_Flag" in group.columns:
        quality = pd.to_numeric(group["Quality_Flag"], errors="coerce")
        result.update(
            {
                "quality_flag_values": collapse_values(quality),
                "n_quality_flag_0": int((quality == 0).sum()),
                "n_quality_flag_1": int((quality == 1).sum()),
                "quality_flag_1_fraction": float((quality == 1).mean()),
            }
        )
    if "date_align" in group.columns:
        date_align = group["date_align"].astype("string")
        result.update(
            {
                "date_align_values": collapse_values(date_align),
                "n_date_inrange": int((date_align == "inrange").sum()),
                "n_date_outrange": int((date_align == "outrange").sum()),
                "date_outrange_fraction": float(
                    (date_align == "outrange").mean()
                ),
            }
        )
    source_tile_column = (
        "sentinel_source_tile_t"
        if "sentinel_source_tile_t" in group.columns
        else None
    )
    if source_tile_column:
        result.update(
            {
                "n_sentinel_source_tiles": int(
                    group[source_tile_column].nunique(dropna=True)
                ),
                "sentinel_source_tiles": collapse_values(
                    group[source_tile_column]
                ),
            }
        )
    if "product_path" in group.columns:
        result.update(
            {
                "n_product_paths": int(group["product_path"].nunique(dropna=True)),
                "product_paths": collapse_values(group["product_path"]),
            }
        )
    if "source_matches_input_tile" in group.columns:
        matches = group["source_matches_input_tile"].astype("string").str.lower()
        result["all_source_tiles_match_input"] = bool(
            matches.isin(["true", "1"]).all()
        )

    return result


def rebuild_manifest(
    original_manifest: pd.DataFrame,
    measured_assignments: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    passing = measured_assignments[
        measured_assignments["passes_window_coverage_filter"]
    ].copy()
    retained_counts = passing.groupby("aggregation_id").size()
    retained_ids = set(
        retained_counts[
            retained_counts >= MIN_RETAINED_FOOTPRINTS
        ].index.astype(str)
    )

    retained_assignments = passing[
        passing["aggregation_id"].isin(retained_ids)
    ].copy()
    removed_footprints = measured_assignments[
        ~measured_assignments["passes_window_coverage_filter"]
        | ~measured_assignments["aggregation_id"].isin(retained_ids)
    ].copy()
    removed_footprints["removal_reason"] = np.where(
        ~removed_footprints["passes_window_coverage_filter"],
        "footprint_below_60pct_window_coverage",
        "window_has_fewer_than_4_passing_footprints",
    )

    manifest_by_id = original_manifest.set_index("aggregation_id", drop=False)
    rebuilt_rows: list[dict] = []
    for aggregation_id, group in retained_assignments.groupby(
        "aggregation_id", sort=False
    ):
        original_row = manifest_by_id.loc[aggregation_id].to_dict()
        original_n = int(original_row["n_footprints"])
        original_target = float(original_row["aggregated_target_modis_sif"])
        original_row.update(
            summarize_retained_assignments(group, original_n=original_n)
        )
        original_row["original_aggregated_target_modis_sif"] = original_target
        original_row["target_mean_change_after_coverage_filter"] = (
            original_row["aggregated_target_modis_sif"] - original_target
        )
        original_row["coverage_filter_threshold"] = MIN_MASK_INSIDE_FRACTION
        original_row["coverage_filter_min_footprints"] = (
            MIN_RETAINED_FOOTPRINTS
        )
        rebuilt_rows.append(original_row)

    filtered_manifest = pd.DataFrame(rebuilt_rows)
    original_order = list(original_manifest.columns)
    added_columns = [
        column for column in filtered_manifest.columns
        if column not in original_order
    ]
    filtered_manifest = filtered_manifest[
        [column for column in original_order if column in filtered_manifest.columns]
        + added_columns
    ]
    filtered_manifest = filtered_manifest.sort_values(
        [
            column for column in [
                "mgrs_tile_t",
                "Delta_Date",
                "measurement_mode",
                "n_footprints",
                "aggregation_id",
            ]
            if column in filtered_manifest.columns
        ],
        ascending=[True, True, True, False, True],
    ).reset_index(drop=True)
    retained_assignments = retained_assignments.sort_values(
        ["aggregation_id", "sif_row_id"]
    ).reset_index(drop=True)

    removed_window_ids = set(original_manifest["aggregation_id"]) - retained_ids
    removed_windows = original_manifest[
        original_manifest["aggregation_id"].isin(removed_window_ids)
    ].copy()
    removed_windows["passing_n_footprints"] = (
        removed_windows["aggregation_id"].map(retained_counts).fillna(0).astype(int)
    )
    removed_windows["removed_n_footprints"] = (
        removed_windows["n_footprints"].astype(int)
        - removed_windows["passing_n_footprints"]
    )
    removed_windows["removal_reason"] = (
        "fewer_than_4_footprints_after_60pct_coverage_filter"
    )

    return (
        filtered_manifest,
        retained_assignments,
        removed_footprints,
        removed_windows,
    )


def validate_filtered_outputs(
    manifest: pd.DataFrame,
    assignments: pd.DataFrame,
) -> None:
    if manifest.empty or assignments.empty:
        raise ValueError("Coverage filtering produced empty outputs.")
    if manifest["aggregation_id"].duplicated().any():
        raise ValueError("Filtered manifest contains duplicate aggregation IDs.")
    if not assignments["aggregation_id"].isin(manifest["aggregation_id"]).all():
        raise ValueError("Filtered assignments contain an unknown aggregation ID.")
    if (
        assignments["window_mask_inside_fraction"]
        < MIN_MASK_INSIDE_FRACTION
    ).any():
        raise ValueError("A retained footprint is below the coverage threshold.")

    found_counts = assignments.groupby("aggregation_id").size().sort_index()
    expected_counts = (
        manifest.set_index("aggregation_id")["n_footprints"]
        .astype(int)
        .sort_index()
    )
    if not found_counts.equals(expected_counts):
        raise ValueError("Filtered manifest/assignment counts do not match.")
    if (expected_counts < MIN_RETAINED_FOOTPRINTS).any():
        raise ValueError("A retained window has fewer than four footprints.")

    assignment_targets = assignments.groupby("aggregation_id")[
        "target_modis_sif"
    ].mean().sort_index()
    manifest_targets = (
        manifest.set_index("aggregation_id")["aggregated_target_modis_sif"]
        .astype(float)
        .sort_index()
    )
    if not np.allclose(
        assignment_targets.to_numpy(),
        manifest_targets.to_numpy(),
        rtol=0,
        atol=1e-10,
    ):
        raise ValueError("Recalculated manifest targets do not match assignments.")


def write_outputs(
    original_manifest: pd.DataFrame,
    original_assignments: pd.DataFrame,
    filtered_manifest: pd.DataFrame,
    filtered_assignments: pd.DataFrame,
    removed_footprints: pd.DataFrame,
    removed_windows: pd.DataFrame,
) -> None:
    if OUTPUT_DIR.resolve() == INPUT_DIR.resolve():
        raise ValueError("OUTPUT_DIR must differ from INPUT_DIR.")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    existing_primary_outputs = [
        path for path in [OUTPUT_MANIFEST_PATH, OUTPUT_ASSIGNMENTS_PATH]
        if path.exists()
    ]
    if existing_primary_outputs:
        raise FileExistsError(
            "Filtered primary output already exists. Remove or rename the output "
            f"folder before rerunning: {existing_primary_outputs[0]}"
        )

    original_windows = len(original_manifest)
    retained_windows = len(filtered_manifest)
    original_footprints = len(original_assignments)
    retained_footprints = len(filtered_assignments)
    summary = pd.DataFrame(
        [
            {
                "input_manifest": str(INPUT_MANIFEST_PATH),
                "input_assignments": str(INPUT_ASSIGNMENTS_PATH),
                "mask_inside_fraction_threshold": MIN_MASK_INSIDE_FRACTION,
                "minimum_retained_footprints": MIN_RETAINED_FOOTPRINTS,
                "original_windows": original_windows,
                "retained_windows": retained_windows,
                "lost_windows": original_windows - retained_windows,
                "retained_window_fraction": retained_windows / original_windows,
                "original_footprints": original_footprints,
                "retained_footprints": retained_footprints,
                "removed_footprints": original_footprints - retained_footprints,
                "retained_footprint_fraction": (
                    retained_footprints / original_footprints
                ),
            }
        ]
    )

    filtered_manifest.to_csv(OUTPUT_MANIFEST_PATH, index=False)
    filtered_assignments.to_csv(OUTPUT_ASSIGNMENTS_PATH, index=False)
    removed_footprints.to_csv(REMOVED_FOOTPRINTS_PATH, index=False)
    removed_windows.to_csv(REMOVED_WINDOWS_PATH, index=False)
    summary.to_csv(FILTER_SUMMARY_PATH, index=False)

    print("\nCoverage filtering complete")
    print(f"Original windows: {original_windows:,}")
    print(
        f"Retained windows: {retained_windows:,} "
        f"({retained_windows / original_windows:.1%})"
    )
    print(f"Lost windows: {original_windows - retained_windows:,}")
    print(f"Original assignments: {original_footprints:,}")
    print(
        f"Retained assignments: {retained_footprints:,} "
        f"({retained_footprints / original_footprints:.1%})"
    )
    print(f"Removed assignments: {original_footprints - retained_footprints:,}")
    print(f"Filtered manifest: {OUTPUT_MANIFEST_PATH}")
    print(f"Filtered assignments: {OUTPUT_ASSIGNMENTS_PATH}")
    print(f"Removed-footprint audit: {REMOVED_FOOTPRINTS_PATH}")
    print(f"Removed-window audit: {REMOVED_WINDOWS_PATH}")
    print(f"Filter summary: {FILTER_SUMMARY_PATH}")


def main() -> None:
    manifest, assignments = load_tables()
    measured_assignments = add_coverage_measurements(manifest, assignments)
    (
        filtered_manifest,
        filtered_assignments,
        removed_footprints,
        removed_windows,
    ) = rebuild_manifest(manifest, measured_assignments)
    validate_filtered_outputs(filtered_manifest, filtered_assignments)
    write_outputs(
        original_manifest=manifest,
        original_assignments=assignments,
        filtered_manifest=filtered_manifest,
        filtered_assignments=filtered_assignments,
        removed_footprints=removed_footprints,
        removed_windows=removed_windows,
    )


if __name__ == "__main__":
    main()
