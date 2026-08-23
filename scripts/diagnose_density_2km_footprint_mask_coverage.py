"""Diagnose footprint coverage inside density-based 2 km SIF windows.

This script reads only the density-clustering manifest and SIF assignments. It
does not read Sentinel-2 or any other predictor raster and does not create CNN
chips. Each OCO-2 footprint is reconstructed from its corner coordinates,
projected to the 2 km window CRS, and rasterized on the same 20 m grid that
would later be used for model preparation.

A footprint passes when at least 80% of its rasterized area is inside its
assigned window. A window is retained after removing failed footprints when at
least three passing footprints remain. The script reports both this primary
retention count and the stricter number of windows in which every assigned
footprint passes the threshold. Results are printed to the terminal only.
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

AGGREGATION_DIR = Path(
    "data/density_aggregation/"
    "sentinel2_spatial_aggregation_density_2000m_landcover"
)
MANIFEST_PATH = (
    AGGREGATION_DIR / "density_cluster_2000m_aggregate_manifest.csv"
)
ASSIGNMENTS_PATH = (
    AGGREGATION_DIR / "density_cluster_2000m_sif_assignments.csv"
)

WINDOW_SIZE_M = 2000.0
PIXEL_SIZE_M = 20.0
WINDOW_PIXELS = 100
MASK_OVERSAMPLE = 4

MIN_MASK_INSIDE_FRACTION = 0.80
MIN_RETAINED_FOOTPRINTS = 3

# All study tiles used here are Sentinel-2 MGRS zone 32N products. No raster is
# opened merely to retrieve this CRS.
WINDOW_CRS = "EPSG:32632"
WGS84_CRS = "EPSG:4326"

PROGRESS_EVERY_WINDOWS = 250


# ---------------------------------------------------------------------------
# Input validation and geometry helpers


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
            f"Could not construct footprint polygon for "
            f"sif_row_id={row['sif_row_id']}"
        )
    return polygon


def rasterized_inside_fraction(
    footprint_projected,
    cell_xmin: float,
    cell_ymax: float,
) -> tuple[float, float, float]:
    """Return inside fraction, rasterized inside area, and full polygon area."""
    high_pixel_size = PIXEL_SIZE_M / MASK_OVERSAMPLE
    high_pixels = WINDOW_PIXELS * MASK_OVERSAMPLE
    high_transform = from_origin(
        cell_xmin,
        cell_ymax,
        high_pixel_size,
        high_pixel_size,
    )

    # Rasterization is automatically clipped to the 2 km output extent.
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


# ---------------------------------------------------------------------------
# Diagnostic


def load_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    if not MANIFEST_PATH.exists():
        raise FileNotFoundError(MANIFEST_PATH)
    if not ASSIGNMENTS_PATH.exists():
        raise FileNotFoundError(ASSIGNMENTS_PATH)

    manifest = pd.read_csv(MANIFEST_PATH)
    assignments = pd.read_csv(ASSIGNMENTS_PATH)

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

    numeric_manifest_columns = [
        "aggregation_size_m",
        "cell_pixels",
        "cell_xmin",
        "cell_ymin",
        "cell_xmax",
        "cell_ymax",
        "n_footprints",
    ]
    numeric_assignment_columns = [
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
    for column in numeric_manifest_columns:
        manifest[column] = pd.to_numeric(manifest[column], errors="raise")
    for column in numeric_assignment_columns:
        assignments[column] = pd.to_numeric(assignments[column], errors="raise")

    manifest["aggregation_id"] = manifest["aggregation_id"].astype(str)
    assignments["aggregation_id"] = assignments["aggregation_id"].astype(str)

    manifest = manifest[
        np.isclose(manifest["aggregation_size_m"], WINDOW_SIZE_M)
        & (manifest["cell_pixels"].astype(int) == WINDOW_PIXELS)
        & (manifest["n_footprints"].astype(int) >= MIN_RETAINED_FOOTPRINTS)
    ].copy()

    if manifest["aggregation_id"].duplicated().any():
        duplicate_id = manifest.loc[
            manifest["aggregation_id"].duplicated(), "aggregation_id"
        ].iloc[0]
        raise ValueError(f"Duplicate aggregation_id in manifest: {duplicate_id}")

    assignments = assignments[
        assignments["aggregation_id"].isin(manifest["aggregation_id"])
    ].copy()

    manifest_ids = set(manifest["aggregation_id"])
    assignment_ids = set(assignments["aggregation_id"])
    missing_assignment_ids = sorted(manifest_ids - assignment_ids)
    if missing_assignment_ids:
        raise ValueError(
            "Manifest windows without footprint assignments: "
            f"{missing_assignment_ids[:5]}"
        )

    assignment_counts = assignments.groupby("aggregation_id").size()
    expected_counts = manifest.set_index("aggregation_id")["n_footprints"].astype(int)
    count_comparison = expected_counts.to_frame("expected").join(
        assignment_counts.rename("found")
    )
    mismatched = count_comparison[
        count_comparison["expected"] != count_comparison["found"]
    ]
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


def diagnose_coverage() -> None:
    manifest, assignments = load_tables()
    assignment_groups = {
        aggregation_id: group.reset_index(drop=True)
        for aggregation_id, group in assignments.groupby(
            "aggregation_id", sort=False
        )
    }

    transformer = Transformer.from_crs(
        WGS84_CRS,
        WINDOW_CRS,
        always_xy=True,
    )

    footprint_rows: list[dict] = []
    window_rows: list[dict] = []
    total_windows = len(manifest)

    print(
        f"Diagnosing {total_windows:,} density-based 2 km windows and "
        f"{len(assignments):,} assigned footprints"
    )
    print(
        f"Footprint threshold: >= {MIN_MASK_INSIDE_FRACTION:.0%} inside window; "
        f"minimum retained footprints: {MIN_RETAINED_FOOTPRINTS}"
    )

    for window_index, manifest_row in manifest.iterrows():
        aggregation_id = str(manifest_row["aggregation_id"])
        footprints = assignment_groups[aggregation_id]
        cell_xmin = float(manifest_row["cell_xmin"])
        cell_ymin = float(manifest_row["cell_ymin"])
        cell_xmax = float(manifest_row["cell_xmax"])
        cell_ymax = float(manifest_row["cell_ymax"])

        width = cell_xmax - cell_xmin
        height = cell_ymax - cell_ymin
        if not np.isclose(width, WINDOW_SIZE_M) or not np.isclose(
            height, WINDOW_SIZE_M
        ):
            raise ValueError(
                f"Window {aggregation_id} is {width} x {height} m, expected "
                f"{WINDOW_SIZE_M} x {WINDOW_SIZE_M} m"
            )

        window_footprint_rows: list[dict] = []
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
            passes = bool(inside_fraction >= MIN_MASK_INSIDE_FRACTION)

            diagnostic_row = {
                "aggregation_id": aggregation_id,
                "sif_row_id": int(footprint["sif_row_id"]),
                "source_csv_row": int(footprint["source_csv_row"]),
                "target_modis_sif": float(footprint["target_modis_sif"]),
                "mask_inside_fraction": inside_fraction,
                "mask_inside_pct": inside_fraction * 100.0,
                "rasterized_inside_area_km2": inside_area_m2 / 1_000_000.0,
                "projected_full_area_km2": full_area_m2 / 1_000_000.0,
                "passes_80pct_mask_inside": passes,
            }
            window_footprint_rows.append(diagnostic_row)
            footprint_rows.append(diagnostic_row)

        window_coverage = pd.DataFrame(window_footprint_rows)
        passing = window_coverage["passes_80pct_mask_inside"]
        retained_targets = window_coverage.loc[passing, "target_modis_sif"]
        original_targets = window_coverage["target_modis_sif"]

        original_n = len(window_coverage)
        retained_n = int(passing.sum())
        failed_n = original_n - retained_n
        retained_window = retained_n >= MIN_RETAINED_FOOTPRINTS
        strict_all_pass = failed_n == 0

        output_row = manifest_row.to_dict()
        output_row.update(
            {
                "original_n_footprints": original_n,
                "retained_n_footprints": retained_n,
                "removed_n_footprints": failed_n,
                "retained_footprint_fraction": retained_n / original_n,
                "mean_mask_inside_fraction": float(
                    window_coverage["mask_inside_fraction"].mean()
                ),
                "min_mask_inside_fraction": float(
                    window_coverage["mask_inside_fraction"].min()
                ),
                "max_mask_inside_fraction": float(
                    window_coverage["mask_inside_fraction"].max()
                ),
                "all_original_footprints_pass_80pct": strict_all_pass,
                "window_retained_after_80pct_filter": retained_window,
                "window_lost_after_80pct_filter": not retained_window,
                "original_target_modis_sif_mean": float(
                    original_targets.mean()
                ),
                "retained_target_modis_sif_mean": (
                    float(retained_targets.mean())
                    if retained_n > 0
                    else np.nan
                ),
            }
        )
        if retained_n > 0:
            output_row["target_mean_change_after_filter"] = (
                output_row["retained_target_modis_sif_mean"]
                - output_row["original_target_modis_sif_mean"]
            )
        else:
            output_row["target_mean_change_after_filter"] = np.nan

        window_rows.append(output_row)

        completed = window_index + 1
        if completed == 1 or completed % PROGRESS_EVERY_WINDOWS == 0:
            print(f"Processed window {completed:,} / {total_windows:,}")

    footprint_diagnostics = pd.DataFrame(footprint_rows)
    window_diagnostics = pd.DataFrame(window_rows)

    original_windows = len(window_diagnostics)
    retained_windows = int(
        window_diagnostics["window_retained_after_80pct_filter"].sum()
    )
    lost_windows = original_windows - retained_windows
    strict_windows = int(
        window_diagnostics["all_original_footprints_pass_80pct"].sum()
    )
    original_footprints = len(footprint_diagnostics)
    retained_footprints = int(
        footprint_diagnostics["passes_80pct_mask_inside"].sum()
    )
    removed_footprints = original_footprints - retained_footprints

    retained_count_distribution = (
        window_diagnostics.groupby(
            [
                "retained_n_footprints",
                "window_retained_after_80pct_filter",
            ],
            dropna=False,
        )
        .size()
        .reset_index(name="n_windows")
        .sort_values("retained_n_footprints")
    )

    print("\nCoverage diagnostic complete")
    print(f"Original windows: {original_windows:,}")
    print(
        "Windows retained after removing <80% footprints and requiring "
        f">={MIN_RETAINED_FOOTPRINTS}: {retained_windows:,} "
        f"({retained_windows / original_windows:.1%})"
    )
    print(
        f"Windows lost: {lost_windows:,} "
        f"({lost_windows / original_windows:.1%})"
    )
    print(
        "Strict windows where every assigned footprint passes 80%: "
        f"{strict_windows:,} ({strict_windows / original_windows:.1%})"
    )
    print(
        f"Footprints passing: {retained_footprints:,} / "
        f"{original_footprints:,} ({retained_footprints / original_footprints:.1%})"
    )
    print(f"Footprints failing: {removed_footprints:,}")
    print("\nRetained-footprint count distribution:")
    print(retained_count_distribution.to_string(index=False))


if __name__ == "__main__":
    diagnose_coverage()
