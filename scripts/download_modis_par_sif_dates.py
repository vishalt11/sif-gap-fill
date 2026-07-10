"""Download MODIS MCD18A2 PAR at 12:00 UTC for the dates in sif_dates.csv.

The script searches the active MCD18A2 Version 6.2 Earthdata collection for
each unique ``Delta_Date``, downloads only the MODIS tiles intersecting the
Germany bounding box, and writes two compact, date-specific GeoTIFFs:

* total incident PAR from ``GMT_1200_PAR`` (W m-2), and
* the corresponding ``PAR_Quality`` categorical layer.

Outputs use a fixed 1,000 m ETRS89 / LAEA Europe grid (EPSG:3035), so every
date is pixel-aligned. Raw HDF-EOS2 granules are temporary by default because
they are much larger than the Germany subset. Pass ``--keep-hdf`` to retain
them.

Authentication is handled by ``earthaccess``. It checks an Earthdata
``.netrc`` file or the EARTHDATA_USERNAME/EARTHDATA_PASSWORD environment
variables and can prompt interactively when necessary. Do not put credentials
in this file.

Required packages: earthaccess, numpy, rasterio. The GDAL installation used by
rasterio must include HDF4/HDF-EOS support for MCD18A2 files.
"""

from __future__ import annotations

import argparse
import csv
import math
import tempfile
from dataclasses import dataclass
from datetime import date, datetime, time, timezone
from pathlib import Path

import earthaccess
import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.errors import RasterioIOError
from rasterio.transform import Affine, from_origin
from rasterio.vrt import WarpedVRT
from rasterio.warp import transform_bounds


PROJECT_ROOT = Path(__file__).resolve().parent
DEFAULT_DATE_CSV = PROJECT_ROOT / "data" / "sif_dates.csv"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "data" / "modis_par_mcd18a2_germany_1200utc"

DATE_COLUMN = "Delta_Date"
COLLECTION_CONCEPT_ID = "C2486282714-LPCLOUD"
PRODUCT = "MCD18A2.062"
PAR_LAYER = "GMT_1200_PAR"
QUALITY_LAYER = "PAR_Quality"

# Padded national extent, including Germany's offshore islands.
# Order: west, south, east, north in WGS84 longitude/latitude.
DEFAULT_GERMANY_BOUNDS = (5.5, 47.0, 15.5, 55.2)

OUTPUT_CRS = "EPSG:3035"
OUTPUT_RESOLUTION_METERS = 1000.0
PAR_NODATA = -9999.0
QUALITY_NODATA = 255

CHECKLIST_NAME = "mcd18a2_1200utc_download_checklist.csv"
CHECKLIST_FIELDS = [
    "date",
    "status",
    "granules_found",
    "par_output",
    "quality_output",
    "valid_par_pixels",
    "error",
    "updated_at_utc",
]


class FatalProcessingError(RuntimeError):
    """A local configuration problem that would make every date fail."""


@dataclass(frozen=True)
class TargetGrid:
    crs: str
    transform: Affine
    width: int
    height: int
    bounds: tuple[float, float, float, float]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Download MCD18A2.062 total PAR at 12:00 UTC for unique "
            "Delta_Date values and subset the tiles to Germany."
        )
    )
    parser.add_argument(
        "--dates-csv",
        type=Path,
        default=DEFAULT_DATE_CSV,
        help=f"CSV containing {DATE_COLUMN} (default: {DEFAULT_DATE_CSV})",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Output directory (default: {DEFAULT_OUTPUT_DIR})",
    )
    parser.add_argument(
        "--bounds",
        type=float,
        nargs=4,
        metavar=("WEST", "SOUTH", "EAST", "NORTH"),
        default=DEFAULT_GERMANY_BOUNDS,
        help="WGS84 search/subset bounds (default: padded Germany extent)",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Parallel Earthdata downloads within one date (default: 4)",
    )
    parser.add_argument(
        "--keep-hdf",
        action="store_true",
        help="Retain the large source HDF files after creating GeoTIFFs",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Recreate existing completed GeoTIFFs",
    )
    return parser.parse_args()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def load_unique_dates(csv_path: Path) -> list[date]:
    if not csv_path.is_file():
        raise FileNotFoundError(f"Date CSV does not exist: {csv_path}")

    parsed_dates: set[date] = set()
    with csv_path.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or DATE_COLUMN not in reader.fieldnames:
            raise ValueError(
                f"{csv_path} must contain a {DATE_COLUMN!r} column; "
                f"found {reader.fieldnames}"
            )

        for row_number, row in enumerate(reader, start=2):
            raw_value = (row.get(DATE_COLUMN) or "").strip()
            if not raw_value:
                continue
            try:
                parsed_dates.add(date.fromisoformat(raw_value[:10]))
            except ValueError as error:
                raise ValueError(
                    f"Invalid {DATE_COLUMN} value on CSV row {row_number}: "
                    f"{raw_value!r}"
                ) from error

    if not parsed_dates:
        raise ValueError(f"No valid dates found in {csv_path}")

    return sorted(parsed_dates)


def validate_bounds(bounds: tuple[float, float, float, float]) -> None:
    west, south, east, north = bounds
    if not (-180.0 <= west < east <= 180.0):
        raise ValueError(f"Invalid west/east bounds: {west}, {east}")
    if not (-90.0 <= south < north <= 90.0):
        raise ValueError(f"Invalid south/north bounds: {south}, {north}")


def make_target_grid(
    bounds_wgs84: tuple[float, float, float, float],
) -> TargetGrid:
    projected = transform_bounds(
        "EPSG:4326",
        OUTPUT_CRS,
        *bounds_wgs84,
        densify_pts=41,
    )
    left = math.floor(projected[0] / OUTPUT_RESOLUTION_METERS) * OUTPUT_RESOLUTION_METERS
    bottom = math.floor(projected[1] / OUTPUT_RESOLUTION_METERS) * OUTPUT_RESOLUTION_METERS
    right = math.ceil(projected[2] / OUTPUT_RESOLUTION_METERS) * OUTPUT_RESOLUTION_METERS
    top = math.ceil(projected[3] / OUTPUT_RESOLUTION_METERS) * OUTPUT_RESOLUTION_METERS

    width = int(round((right - left) / OUTPUT_RESOLUTION_METERS))
    height = int(round((top - bottom) / OUTPUT_RESOLUTION_METERS))
    transform = from_origin(
        left,
        top,
        OUTPUT_RESOLUTION_METERS,
        OUTPUT_RESOLUTION_METERS,
    )
    return TargetGrid(
        crs=OUTPUT_CRS,
        transform=transform,
        width=width,
        height=height,
        bounds=(left, bottom, right, top),
    )


def load_checklist(path: Path) -> dict[str, dict[str, str]]:
    if not path.is_file():
        return {}

    with path.open("r", newline="", encoding="utf-8") as handle:
        return {
            row["date"]: row
            for row in csv.DictReader(handle)
            if row.get("date")
        }


def save_checklist(path: Path, rows: dict[str, dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = path.with_name(path.name + ".part")
    with temporary_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CHECKLIST_FIELDS)
        writer.writeheader()
        for date_string in sorted(rows):
            writer.writerow(rows[date_string])
    temporary_path.replace(path)


def checklist_row(
    target_date: date,
    status: str,
    *,
    granules_found: int | None = None,
    par_output: Path | None = None,
    quality_output: Path | None = None,
    valid_par_pixels: int | None = None,
    error: str | None = None,
) -> dict[str, str]:
    return {
        "date": target_date.isoformat(),
        "status": status,
        "granules_found": "" if granules_found is None else str(granules_found),
        "par_output": "" if par_output is None else str(par_output),
        "quality_output": "" if quality_output is None else str(quality_output),
        "valid_par_pixels": (
            "" if valid_par_pixels is None else str(valid_par_pixels)
        ),
        "error": error or "",
        "updated_at_utc": utc_now(),
    }


def output_paths(output_dir: Path, target_date: date) -> tuple[Path, Path]:
    year_dir = output_dir / str(target_date.year)
    prefix = f"{PRODUCT}_{target_date.isoformat()}"
    par_path = year_dir / f"{prefix}_{PAR_LAYER}_EPSG3035_1km.tif"
    quality_path = year_dir / f"{prefix}_{QUALITY_LAYER}_EPSG3035_1km.tif"
    return par_path, quality_path


def file_is_complete(path: Path) -> bool:
    return path.is_file() and path.stat().st_size > 0


def search_granules(
    target_date: date,
    bounds_wgs84: tuple[float, float, float, float],
):
    start = datetime.combine(target_date, time.min, tzinfo=timezone.utc)
    end = datetime.combine(target_date, time.max, tzinfo=timezone.utc)
    return earthaccess.search_data(
        concept_id=COLLECTION_CONCEPT_ID,
        temporal=(start, end),
        bounding_box=bounds_wgs84,
        count=20,
    )


def normalize_layer_name(subdataset: str) -> str:
    return subdataset.rsplit(":", maxsplit=1)[-1].strip('"').replace(" ", "")


def subdataset_for(hdf_path: Path, layer_name: str) -> str:
    try:
        with rasterio.open(hdf_path) as container:
            subdatasets = list(container.subdatasets)
    except RasterioIOError as error:
        raise FatalProcessingError(
            "Could not open an MCD18A2 HDF-EOS2 file. Ensure the GDAL build "
            "used by rasterio includes HDF4/HDF-EOS support. "
            f"File: {hdf_path}; error: {error}"
        ) from error

    if not subdatasets:
        raise FatalProcessingError(
            f"No HDF subdatasets were exposed by GDAL for {hdf_path}. "
            "HDF4/HDF-EOS support is probably unavailable."
        )

    normalized_target = layer_name.replace(" ", "")
    matches = [
        item
        for item in subdatasets
        if normalize_layer_name(item) == normalized_target
    ]
    if len(matches) != 1:
        available = sorted(normalize_layer_name(item) for item in subdatasets)
        raise FatalProcessingError(
            f"Expected exactly one {layer_name!r} subdataset in {hdf_path}, "
            f"found {len(matches)}. Available layers: {available}"
        )
    return matches[0]


def mosaic_par(hdf_paths: list[Path], grid: TargetGrid) -> np.ndarray:
    mosaic = np.full(
        (grid.height, grid.width),
        PAR_NODATA,
        dtype=np.float32,
    )

    for hdf_path in hdf_paths:
        subdataset = subdataset_for(hdf_path, PAR_LAYER)
        with rasterio.open(subdataset) as source:
            if source.crs is None:
                raise FatalProcessingError(
                    f"The {PAR_LAYER} subdataset has no CRS: {hdf_path}"
                )
            source_nodata = source.nodata if source.nodata is not None else -1.0
            with WarpedVRT(
                source,
                crs=grid.crs,
                transform=grid.transform,
                width=grid.width,
                height=grid.height,
                src_nodata=source_nodata,
                nodata=PAR_NODATA,
                dtype="float32",
                resampling=Resampling.bilinear,
            ) as warped:
                values = warped.read(1, masked=True)

        data = np.asarray(values.data, dtype=np.float32)
        valid = ~np.ma.getmaskarray(values)
        valid &= np.isfinite(data)
        valid &= (data >= 0.0) & (data <= 700.0)
        mosaic[valid] = data[valid]

    return mosaic


def mosaic_quality(hdf_paths: list[Path], grid: TargetGrid) -> np.ndarray:
    mosaic = np.full(
        (grid.height, grid.width),
        QUALITY_NODATA,
        dtype=np.uint8,
    )
    allowed_values = np.array([0, 1, 2, 4], dtype=np.uint8)

    for hdf_path in hdf_paths:
        subdataset = subdataset_for(hdf_path, QUALITY_LAYER)
        with rasterio.open(subdataset) as source:
            if source.crs is None:
                raise FatalProcessingError(
                    f"The {QUALITY_LAYER} subdataset has no CRS: {hdf_path}"
                )
            with WarpedVRT(
                source,
                crs=grid.crs,
                transform=grid.transform,
                width=grid.width,
                height=grid.height,
                nodata=QUALITY_NODATA,
                dtype="uint8",
                resampling=Resampling.nearest,
            ) as warped:
                values = warped.read(1, masked=True)

        data = np.asarray(values.data, dtype=np.uint8)
        valid = ~np.ma.getmaskarray(values)
        valid &= np.isin(data, allowed_values)
        mosaic[valid] = data[valid]

    return mosaic


def write_geotiff(
    path: Path,
    array: np.ndarray,
    grid: TargetGrid,
    *,
    nodata: float | int,
    layer_name: str,
    target_date: date,
    bounds_wgs84: tuple[float, float, float, float],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = path.with_name(path.name + ".part")
    predictor = 3 if np.issubdtype(array.dtype, np.floating) else 2

    with rasterio.open(
        temporary_path,
        "w",
        driver="GTiff",
        width=grid.width,
        height=grid.height,
        count=1,
        dtype=array.dtype.name,
        crs=grid.crs,
        transform=grid.transform,
        nodata=nodata,
        compress="DEFLATE",
        predictor=predictor,
        tiled=True,
        blockxsize=256,
        blockysize=256,
        BIGTIFF="IF_SAFER",
    ) as destination:
        destination.write(array, 1)
        destination.set_band_description(1, layer_name)
        destination.update_tags(
            product=PRODUCT,
            collection_concept_id=COLLECTION_CONCEPT_ID,
            layer=layer_name,
            date_utc=target_date.isoformat(),
            observation_time_utc="12:00:00" if layer_name == PAR_LAYER else "",
            units="W m-2" if layer_name == PAR_LAYER else "category",
            germany_bounds_wgs84=",".join(str(value) for value in bounds_wgs84),
            output_crs=grid.crs,
            output_resolution_m=str(int(OUTPUT_RESOLUTION_METERS)),
            quality_meaning=(
                "0=no valid surface reflectance; 1=MCD43 surface reflectance; "
                "2=climatology surface reflectance; 4=non-land"
                if layer_name == QUALITY_LAYER
                else ""
            ),
        )

    temporary_path.replace(path)


def downloaded_hdf_paths(download_results) -> list[Path]:
    paths = [Path(item) for item in download_results]
    hdf_paths = [path for path in paths if path.suffix.lower() == ".hdf"]
    if not hdf_paths:
        raise RuntimeError(
            f"Earthdata returned no HDF files. Download results: {paths}"
        )
    return sorted(hdf_paths)


def process_date(
    target_date: date,
    granules,
    output_dir: Path,
    grid: TargetGrid,
    bounds_wgs84: tuple[float, float, float, float],
    *,
    threads: int,
    keep_hdf: bool,
) -> tuple[Path, Path, int]:
    par_path, quality_path = output_paths(output_dir, target_date)
    temporary_download: tempfile.TemporaryDirectory[str] | None = None

    if keep_hdf:
        download_dir = (
            output_dir
            / "raw_hdf"
            / str(target_date.year)
            / target_date.isoformat()
        )
        download_dir.mkdir(parents=True, exist_ok=True)
    else:
        temporary_root = output_dir / "temporary_hdf"
        temporary_root.mkdir(parents=True, exist_ok=True)
        temporary_download = tempfile.TemporaryDirectory(
            prefix=f"{target_date.isoformat()}_",
            dir=temporary_root,
        )
        download_dir = Path(temporary_download.name)

    try:
        downloaded = earthaccess.download(
            granules,
            local_path=str(download_dir),
            threads=threads,
            show_progress=True,
        )
        hdf_paths = downloaded_hdf_paths(downloaded)

        par = mosaic_par(hdf_paths, grid)
        quality = mosaic_quality(hdf_paths, grid)
        valid_par_pixels = int(np.count_nonzero(par != PAR_NODATA))
        if valid_par_pixels == 0:
            raise RuntimeError(
                f"No valid {PAR_LAYER} pixels were found for {target_date}"
            )

        write_geotiff(
            par_path,
            par,
            grid,
            nodata=PAR_NODATA,
            layer_name=PAR_LAYER,
            target_date=target_date,
            bounds_wgs84=bounds_wgs84,
        )
        write_geotiff(
            quality_path,
            quality,
            grid,
            nodata=QUALITY_NODATA,
            layer_name=QUALITY_LAYER,
            target_date=target_date,
            bounds_wgs84=bounds_wgs84,
        )
        return par_path, quality_path, valid_par_pixels
    finally:
        if temporary_download is not None:
            temporary_download.cleanup()


def main() -> None:
    args = parse_args()
    if args.threads < 1:
        raise ValueError("--threads must be at least 1")

    bounds_wgs84 = tuple(float(value) for value in args.bounds)
    validate_bounds(bounds_wgs84)
    dates = load_unique_dates(args.dates_csv.resolve())

    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    checklist_path = output_dir / CHECKLIST_NAME
    checklist = load_checklist(checklist_path)
    grid = make_target_grid(bounds_wgs84)

    print(f"Unique SIF dates: {len(dates)}")
    print(f"First/last date: {dates[0]} / {dates[-1]}")
    print(f"Earthdata product/layer: {PRODUCT} / {PAR_LAYER}")
    print(f"UTC time: 12:00")
    print(f"Germany WGS84 bounds: {bounds_wgs84}")
    print(
        f"Output grid: {grid.width} x {grid.height} pixels, "
        f"{int(OUTPUT_RESOLUTION_METERS)} m, {grid.crs}"
    )
    print(f"Output directory: {output_dir}")
    print("Authenticating with NASA Earthdata...")
    earthaccess.login()

    for index, target_date in enumerate(dates, start=1):
        par_path, quality_path = output_paths(output_dir, target_date)
        print(f"[{index}/{len(dates)}] {target_date.isoformat()}")

        if (
            not args.overwrite
            and file_is_complete(par_path)
            and file_is_complete(quality_path)
        ):
            previous = checklist.get(target_date.isoformat(), {})
            checklist[target_date.isoformat()] = checklist_row(
                target_date,
                "completed",
                granules_found=(
                    int(previous["granules_found"])
                    if previous.get("granules_found")
                    else None
                ),
                par_output=par_path,
                quality_output=quality_path,
                valid_par_pixels=(
                    int(previous["valid_par_pixels"])
                    if previous.get("valid_par_pixels")
                    else None
                ),
            )
            save_checklist(checklist_path, checklist)
            print("  Existing outputs are complete; skipped.")
            continue

        granules = search_granules(target_date, bounds_wgs84)
        granule_count = len(granules)
        print(f"  Matching Germany tiles: {granule_count}")
        if granule_count == 0:
            checklist[target_date.isoformat()] = checklist_row(
                target_date,
                "no_granules",
                granules_found=0,
                par_output=par_path,
                quality_output=quality_path,
                error="No MCD18A2.062 granules intersected the date and bounds",
            )
            save_checklist(checklist_path, checklist)
            continue

        try:
            par_path, quality_path, valid_pixels = process_date(
                target_date,
                granules,
                output_dir,
                grid,
                bounds_wgs84,
                threads=args.threads,
                keep_hdf=args.keep_hdf,
            )
            checklist[target_date.isoformat()] = checklist_row(
                target_date,
                "completed",
                granules_found=granule_count,
                par_output=par_path,
                quality_output=quality_path,
                valid_par_pixels=valid_pixels,
            )
            print(f"  Wrote: {par_path.name}")
            print(f"  Wrote: {quality_path.name}")
        except FatalProcessingError as error:
            checklist[target_date.isoformat()] = checklist_row(
                target_date,
                "fatal_processing_error",
                granules_found=granule_count,
                par_output=par_path,
                quality_output=quality_path,
                error=str(error),
            )
            save_checklist(checklist_path, checklist)
            raise
        except Exception as error:
            checklist[target_date.isoformat()] = checklist_row(
                target_date,
                "error",
                granules_found=granule_count,
                par_output=par_path,
                quality_output=quality_path,
                error=repr(error),
            )
            print(f"  ERROR: {error!r}")

        save_checklist(checklist_path, checklist)

    completed = sum(
        row.get("status") == "completed" for row in checklist.values()
    )
    print(f"Finished. Completed dates in checklist: {completed}/{len(dates)}")
    print(f"Checklist: {checklist_path}")


if __name__ == "__main__":
    main()
