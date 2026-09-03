"""Download all 4 km Sentinel-2 L2A nearest-clear-pixel composites.

The script reads 4 km windows from the density aggregation manifest and
requests Sentinel-2 L2A observations around each target date. The configured
window is target_date - 8 days through target_date + 8 days. Each request is
restricted to the MGRS tile assigned to that window in the manifest.

The Process API returns the six required reflectance bands as scaled UINT16
values to reduce processing-unit cost. They are decoded locally and the final
GeoTIFF stores unitless BOA reflectance as FLOAT32, preserving the convenient
physical-value representation used by downstream analysis.

Pixel selection order is:

1. A good-quality pixel on the target date.
2. A good-quality pixel on the nearest available date.
3. If dates are equally distant, prefer the date after the target date.
4. If multiple scenes remain tied, prefer the lower cloud probability.

A pixel is considered good when it has data, its Sentinel-2 Scene
Classification Layer (SCL) class is vegetation, bare soil, or water, and its
cloud probability is at most CLOUD_PROBABILITY_MAX_PERCENT. Cloud shadow,
unclassified/low-probability cloud, medium/high-probability cloud, cirrus,
snow/ice, saturated/defective, and no-data pixels are rejected.

Authentication
--------------
Create an OAuth client in the Copernicus Data Space Ecosystem dashboard, then
provide credentials through either:

* environment variables CDSE_SH_CLIENT_ID and CDSE_SH_CLIENT_SECRET; or
* an existing sentinelhub-py profile named by SH_PROFILE.

By default all manifest rows are processed. Set TEST_AGGREGATION_ID below, pass
--aggregation-id, or pass --limit for a smaller run. Completed, structurally
valid outputs are skipped so an interrupted run can be resumed safely.
"""

from __future__ import annotations

import argparse
import json
import os
import random
import time
from datetime import timedelta
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from rasterio.transform import from_bounds
from scipy.ndimage import convolve
from sentinelhub import (
    BBox,
    DataCollection,
    MimeType,
    MosaickingOrder,
    SentinelHubRequest,
    SHConfig,
    bbox_to_dimensions,
)


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent

MANIFEST_PATH = (
    PROJECT_ROOT
    / "data"
    / "density_aggregation"
    / "sentinel2_spatial_aggregation_density_4000m_landcover_redtiles_mask60_min4"
    / "density_cluster_4000m_aggregate_manifest.csv"
)
OUTPUT_DIR = PROJECT_ROOT / "data" / "sentinel2_l2a"

# Leave as None to process the complete manifest. This can also be overridden
# with --aggregation-id or --limit on the command line.
TEST_AGGREGATION_ID: str | None = None
TEST_WINDOW_COUNT: int | None = None

DAYS_EITHER_SIDE = 8
RESOLUTION_M = 20
CLOUD_PROBABILITY_MAX_PERCENT = 40
PREFER_FUTURE_ON_EQUAL_DISTANCE = True
LOW_VALID_FRACTION_THRESHOLD = 0.98
MAX_REQUEST_ATTEMPTS = 4
RETRY_BASE_DELAY_SECONDS = 5
PROGRESS_SAVE_EVERY = 10

# Conservative clear classes from the Sentinel-2 L2A Scene Classification
# Layer: 4 = vegetation, 5 = bare soil, 6 = water.
GOOD_SCL_CLASSES = (4, 5, 6)

REFLECTANCE_BANDS = (
    "B02",
    "B04",
    "B05",
    "B08",
    "B8A",
    "B11",
)
LEGACY_REFLECTANCE_BANDS = (
    "B02",
    "B03",
    "B04",
    "B05",
    "B06",
    "B07",
    "B08",
    "B8A",
    "B11",
    "B12",
)
QUALITY_BANDS = (
    "source_day_offset_from_target",
    "selected_SCL",
    "selected_cloud_probability_percent",
    "valid_pixel",
)
REFLECTANCE_SCALE = 10000
SOURCE_DAY_OFFSET_BIAS = 128
UINT16_NODATA = 65535

CDSE_BASE_URL = "https://sh.dataspace.copernicus.eu"
CDSE_TOKEN_URL = (
    "https://identity.dataspace.copernicus.eu/auth/realms/CDSE/"
    "protocol/openid-connect/token"
)

# sentinelhub-py 3.11.5 hard-codes the commercial Sentinel Hub deployment in
# DataCollection.SENTINEL2_L2A. A derived collection is required so that the
# Process request and its OAuth token both use the CDSE deployment.
CDSE_SENTINEL2_L2A = DataCollection.SENTINEL2_L2A.define_from(
    "SENTINEL2_L2A_CDSE",
    service_url=CDSE_BASE_URL,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--aggregation-id",
        default=TEST_AGGREGATION_ID,
        help="Process only this aggregation_id; default: all manifest rows.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=TEST_WINDOW_COUNT,
        help="Process only the first N manifest rows; default: all rows.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=OUTPUT_DIR,
        help=f"Output directory; default: {OUTPUT_DIR}",
    )
    return parser.parse_args()


def load_test_windows(aggregation_id: str | None, limit: int | None) -> pd.DataFrame:
    required_columns = [
        "aggregation_id",
        "cell_xmin",
        "cell_ymin",
        "cell_xmax",
        "cell_ymax",
        "mgrs_tile_t",
        "Delta_Date",
    ]
    manifest = pd.read_csv(MANIFEST_PATH, usecols=required_columns)

    if manifest.empty:
        raise RuntimeError(f"Manifest is empty: {MANIFEST_PATH}")

    if aggregation_id is None:
        if limit is not None and limit < 1:
            raise ValueError(f"--limit must be at least 1, got {limit}")
        rows = manifest.copy() if limit is None else manifest.head(limit).copy()
    else:
        matches = manifest.loc[manifest["aggregation_id"] == aggregation_id]
        if matches.empty:
            raise KeyError(f"aggregation_id not found in manifest: {aggregation_id}")
        if len(matches) != 1:
            raise RuntimeError(
                f"Expected one row for {aggregation_id}, found {len(matches)}"
            )
        rows = matches.copy()

    rows["Delta_Date"] = pd.to_datetime(rows["Delta_Date"]).dt.normalize()

    widths_m = rows["cell_xmax"].astype(float) - rows["cell_xmin"].astype(float)
    heights_m = rows["cell_ymax"].astype(float) - rows["cell_ymin"].astype(float)
    invalid = ~(np.isclose(widths_m, 4000.0) & np.isclose(heights_m, 4000.0))
    if invalid.any():
        bad_ids = rows.loc[invalid, "aggregation_id"].tolist()
        raise ValueError(
            "Selected manifest rows contain non-4 km windows: "
            + ", ".join(str(value) for value in bad_ids)
        )

    return rows.reset_index(drop=True)


def epsg_from_mgrs_tile(mgrs_tile: str) -> int:
    """Return the UTM EPSG code encoded by a tile such as T32ULA."""
    tile = str(mgrs_tile).upper()
    if tile.startswith("T"):
        tile = tile[1:]

    if len(tile) != 5 or not tile[:2].isdigit():
        raise ValueError(f"Unexpected MGRS tile value: {mgrs_tile}")

    zone = int(tile[:2])
    latitude_band = tile[2]
    if not 1 <= zone <= 60 or not "C" <= latitude_band <= "X":
        raise ValueError(f"Unexpected MGRS tile value: {mgrs_tile}")

    # MGRS latitude bands N through X lie in the northern hemisphere.
    return (32600 if latitude_band >= "N" else 32700) + zone


def build_cdse_config() -> SHConfig:
    profile_name = os.environ.get("SH_PROFILE")
    config = SHConfig(profile_name) if profile_name else SHConfig()

    client_id = os.environ.get("CDSE_SH_CLIENT_ID")
    client_secret = os.environ.get("CDSE_SH_CLIENT_SECRET")
    if client_id:
        config.sh_client_id = client_id
    if client_secret:
        config.sh_client_secret = client_secret

    config.sh_base_url = CDSE_BASE_URL
    config.sh_token_url = CDSE_TOKEN_URL

    if not config.sh_client_id or not config.sh_client_secret:
        raise RuntimeError(
            "Missing Copernicus Data Space Sentinel Hub OAuth credentials. "
            "Set CDSE_SH_CLIENT_ID and CDSE_SH_CLIENT_SECRET, or set "
            "SH_PROFILE to a configured sentinelhub-py profile."
        )

    return config


def make_evalscript(target_date: pd.Timestamp, mgrs_tile: str) -> str:
    target_iso = target_date.strftime("%Y-%m-%dT00:00:00Z")
    good_scl = ", ".join(str(value) for value in GOOD_SCL_CLASSES)
    prefer_future = str(PREFER_FUTURE_ON_EQUAL_DISTANCE).lower()
    output_band_count = len(REFLECTANCE_BANDS) + 4
    target_tile = str(mgrs_tile).upper()
    if not target_tile.startswith("T"):
        target_tile = f"T{target_tile}"
    tile_token = f"_{target_tile}_"

    return f"""//VERSION=3
const TARGET_TIME = Date.parse("{target_iso}");
const MILLISECONDS_PER_DAY = 24 * 60 * 60 * 1000;
const GOOD_SCL = [{good_scl}];
const MAX_CLOUD_PROBABILITY = {CLOUD_PROBABILITY_MAX_PERCENT};
const PREFER_FUTURE_ON_TIE = {prefer_future};
const TARGET_TILE_TOKEN = "{tile_token}";
const REFLECTANCE_SCALE = {REFLECTANCE_SCALE};
const SOURCE_DAY_OFFSET_BIAS = {SOURCE_DAY_OFFSET_BIAS};
const UINT16_NODATA = {UINT16_NODATA};

function setup() {{
  return {{
    input: [{{
      bands: [
        "B02", "B04", "B05", "B08", "B8A", "B11",
        "SCL", "CLD", "dataMask"
      ],
      units: [
        "REFLECTANCE", "REFLECTANCE", "REFLECTANCE", "REFLECTANCE",
        "REFLECTANCE", "REFLECTANCE", "DN", "PERCENT", "DN"
      ]
    }}],
    output: {{
      id: "default",
      bands: {output_band_count},
      sampleType: "UINT16"
    }},
    mosaicking: "TILE"
  }};
}}

function preProcessScenes(collections) {{
  collections.scenes.tiles = collections.scenes.tiles.filter(function(tile) {{
    const productId = tile.productId
      || tile.sentinel2ProductId
      || tile.tileOriginalId
      || "";
    return productId.includes(TARGET_TILE_TOKEN);
  }});
  return collections;
}}

function calendarDayOffset(sceneDate) {{
  const dateOnly = sceneDate.split("T")[0] + "T00:00:00Z";
  return Math.round((Date.parse(dateOnly) - TARGET_TIME) / MILLISECONDS_PER_DAY);
}}

function isGoodPixel(sample) {{
  return sample.dataMask === 1
    && GOOD_SCL.includes(sample.SCL)
    && sample.CLD <= MAX_CLOUD_PROBABILITY;
}}

function encodeReflectance(value) {{
  return Math.max(0, Math.min(UINT16_NODATA - 1, Math.round(value * REFLECTANCE_SCALE)));
}}

function compareCandidates(a, b) {{
  const distanceDifference = Math.abs(a.dayOffset) - Math.abs(b.dayOffset);
  if (distanceDifference !== 0) {{
    return distanceDifference;
  }}

  if (a.dayOffset !== b.dayOffset) {{
    return PREFER_FUTURE_ON_TIE
      ? b.dayOffset - a.dayOffset
      : a.dayOffset - b.dayOffset;
  }}

  const cloudDifference = a.sample.CLD - b.sample.CLD;
  return cloudDifference !== 0 ? cloudDifference : a.index - b.index;
}}

function evaluatePixel(samples, scenes) {{
  const candidates = [];

  for (let index = 0; index < samples.length; index++) {{
    if (isGoodPixel(samples[index])) {{
      candidates.push({{
        index: index,
        sample: samples[index],
        dayOffset: calendarDayOffset(scenes.tiles[index].date)
      }});
    }}
  }}

  if (candidates.length === 0) {{
    return [
      UINT16_NODATA, UINT16_NODATA, UINT16_NODATA,
      UINT16_NODATA, UINT16_NODATA, UINT16_NODATA,
      UINT16_NODATA, 0, UINT16_NODATA, 0
    ];
  }}

  candidates.sort(compareCandidates);
  const selected = candidates[0];
  const sample = selected.sample;

  return [
    encodeReflectance(sample.B02), encodeReflectance(sample.B04),
    encodeReflectance(sample.B05), encodeReflectance(sample.B08),
    encodeReflectance(sample.B8A), encodeReflectance(sample.B11),
    selected.dayOffset + SOURCE_DAY_OFFSET_BIAS,
    sample.SCL, Math.round(sample.CLD), 1
  ];
}}
"""


def request_composite(
    config: SHConfig,
    bbox: BBox,
    size: tuple[int, int],
    target_date: pd.Timestamp,
    mgrs_tile: str,
    start_date: pd.Timestamp,
    end_exclusive: pd.Timestamp,
) -> np.ndarray:
    request = SentinelHubRequest(
        evalscript=make_evalscript(target_date, mgrs_tile),
        input_data=[
            SentinelHubRequest.input_data(
                data_collection=CDSE_SENTINEL2_L2A,
                time_interval=(
                    start_date.to_pydatetime(),
                    end_exclusive.to_pydatetime(),
                ),
                maxcc=1.0,
                mosaicking_order=MosaickingOrder.LEAST_CC,
                other_args={"processing": {"harmonizeValues": True}},
            )
        ],
        responses=[SentinelHubRequest.output_response("default", MimeType.TIFF)],
        bbox=bbox,
        size=size,
        config=config,
    )

    request_url = request.download_list[0].url
    expected_url_prefix = f"{CDSE_BASE_URL}/"
    if not request_url.startswith(expected_url_prefix):
        raise RuntimeError(
            "Refusing to send the request to a non-CDSE endpoint: "
            f"{request_url}"
        )

    print(f"Process API endpoint: {request_url}")
    last_error: Exception | None = None
    for attempt in range(1, MAX_REQUEST_ATTEMPTS + 1):
        try:
            result = request.get_data()[0]
            break
        except Exception as exc:
            last_error = exc
            if attempt == MAX_REQUEST_ATTEMPTS:
                raise

            delay = RETRY_BASE_DELAY_SECONDS * (2 ** (attempt - 1))
            delay += random.uniform(0, RETRY_BASE_DELAY_SECONDS)
            print(
                f"Process request attempt {attempt}/{MAX_REQUEST_ATTEMPTS} failed: "
                f"{type(exc).__name__}: {exc}"
            )
            print(f"Retrying in {delay:.1f} seconds")
            time.sleep(delay)
    else:
        raise RuntimeError("Process request failed without an exception") from last_error

    expected_bands = len(REFLECTANCE_BANDS) + 4
    if result.ndim != 3 or result.shape[2] != expected_bands:
        raise RuntimeError(
            "Unexpected Process API array shape: "
            f"{result.shape}; expected (height, width, {expected_bands})"
        )

    encoded = result.astype(np.int32, copy=False)
    decoded = np.empty(encoded.shape, dtype=np.float32)
    reflectance_count = len(REFLECTANCE_BANDS)
    valid = encoded[:, :, -1] == 1

    decoded[:, :, :reflectance_count] = (
        encoded[:, :, :reflectance_count] / REFLECTANCE_SCALE
    )
    decoded[~valid, :reflectance_count] = np.nan
    decoded[:, :, reflectance_count] = np.where(
        valid,
        encoded[:, :, reflectance_count] - SOURCE_DAY_OFFSET_BIAS,
        255,
    )
    decoded[:, :, reflectance_count + 1] = np.where(
        valid, encoded[:, :, reflectance_count + 1], 0
    )
    decoded[:, :, reflectance_count + 2] = np.where(
        valid, encoded[:, :, reflectance_count + 2], 255
    )
    decoded[:, :, reflectance_count + 3] = valid.astype(np.float32)
    return decoded


def fill_missing_reflectance(
    data: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    """Fill missing reflectance pixels with iterative 3x3 neighbour means.

    Only the six reflectance bands are filled. Existing diagnostic bands are
    left unchanged. The original validity array is returned for CSV reporting.
    """
    reflectance = data[:, :, : len(REFLECTANCE_BANDS)]
    original_valid = (data[:, :, -1] == 1) & np.all(
        np.isfinite(reflectance), axis=2
    )
    filled_valid = original_valid.copy()
    filled_data = data.copy()

    if filled_valid.all() or not filled_valid.any():
        return filled_data, original_valid, filled_valid, 0

    kernel = np.ones((3, 3), dtype=np.float32)
    kernel[1, 1] = 0
    iterations = 0

    while not filled_valid.all():
        neighbour_count = convolve(
            filled_valid.astype(np.float32),
            kernel,
            mode="constant",
            cval=0.0,
        )
        fillable = (~filled_valid) & (neighbour_count > 0)
        if not fillable.any():
            break

        for band_index in range(len(REFLECTANCE_BANDS)):
            band = filled_data[:, :, band_index]
            neighbour_sum = convolve(
                np.where(filled_valid, band, 0.0).astype(np.float32),
                kernel,
                mode="constant",
                cval=0.0,
            )
            band[fillable] = neighbour_sum[fillable] / neighbour_count[fillable]

        filled_valid[fillable] = True
        iterations += 1

    return filled_data, original_valid, filled_valid, iterations


def write_output(
    output_path: Path,
    data: np.ndarray,
    dataset_valid_mask: np.ndarray,
    bounds: tuple[float, float, float, float],
    epsg: int,
    row: pd.Series,
    target_date: pd.Timestamp,
    start_date: pd.Timestamp,
    end_inclusive: pd.Timestamp,
    valid_fraction_before_fill: float,
    neighbour_fill_iterations: int,
    neighbour_fill_complete: bool,
) -> None:
    height, width, band_count = data.shape
    transform = from_bounds(*bounds, width=width, height=height)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(
        output_path,
        "w",
        driver="GTiff",
        width=width,
        height=height,
        count=band_count,
        dtype="float32",
        crs=f"EPSG:{epsg}",
        transform=transform,
        compress="deflate",
        predictor=3,
        tiled=True,
        blockxsize=128,
        blockysize=128,
    ) as dst:
        dst.write(np.moveaxis(data.astype(np.float32, copy=False), -1, 0))
        dst.write_mask(dataset_valid_mask.astype(np.uint8) * 255)

        descriptions = list(REFLECTANCE_BANDS) + list(QUALITY_BANDS)
        for band_index, description in enumerate(descriptions, start=1):
            dst.set_band_description(band_index, description)

        dst.update_tags(
            aggregation_id=str(row["aggregation_id"]),
            mgrs_tile=str(row["mgrs_tile_t"]),
            target_date=target_date.date().isoformat(),
            search_start_date=start_date.date().isoformat(),
            search_end_date=end_inclusive.date().isoformat(),
            selection_rule=(
                "nearest good pixel by absolute day distance; future preferred "
                "on equal distance; lower CLD preferred within the same day"
            ),
            good_scl_classes=",".join(str(value) for value in GOOD_SCL_CLASSES),
            max_cloud_probability_percent=str(CLOUD_PROBABILITY_MAX_PERCENT),
            reflectance_units="unitless BOA reflectance",
            valid_fraction_before_neighbour_fill=f"{valid_fraction_before_fill:.8f}",
            neighbour_fill_iterations=str(neighbour_fill_iterations),
            neighbour_fill_complete=str(neighbour_fill_complete),
        )


def source_date_counts(
    target_date: pd.Timestamp,
    source_offsets: np.ndarray,
    original_valid: np.ndarray,
) -> dict[str, int]:
    selected_offsets = source_offsets[original_valid]
    if selected_offsets.size == 0:
        return {}

    offsets, counts = np.unique(selected_offsets.astype(np.int16), return_counts=True)
    return {
        (target_date + timedelta(days=int(offset))).date().isoformat(): int(count)
        for offset, count in zip(offsets, counts)
    }


def convert_legacy_output_if_needed(
    output_path: Path,
    row: pd.Series,
    target_date: pd.Timestamp,
) -> bool:
    """Reduce an existing ten-reflectance-band TIFF without another API call."""
    expected_descriptions = list(REFLECTANCE_BANDS) + list(QUALITY_BANDS)
    try:
        with rasterio.open(output_path) as src:
            descriptions = list(src.descriptions)
            if descriptions == expected_descriptions:
                return False

            legacy_descriptions = list(LEGACY_REFLECTANCE_BANDS) + list(QUALITY_BANDS)
            if descriptions != legacy_descriptions:
                return False

            selected_indexes = [
                descriptions.index(description) + 1
                for description in expected_descriptions
            ]
            selected_data = np.moveaxis(src.read(selected_indexes), 0, -1)
            dataset_valid_mask = src.dataset_mask() > 0
            tags = src.tags()
    except (OSError, ValueError, rasterio.errors.RasterioError):
        return False

    original_valid = selected_data[:, :, -1] == 1
    valid_fraction = float(original_valid.mean())
    fill_iterations = int(tags.get("neighbour_fill_iterations", "0"))
    start_date = target_date - timedelta(days=DAYS_EITHER_SIDE)
    end_inclusive = target_date + timedelta(days=DAYS_EITHER_SIDE)
    bounds = (
        float(row["cell_xmin"]),
        float(row["cell_ymin"]),
        float(row["cell_xmax"]),
        float(row["cell_ymax"]),
    )

    temp_path = output_path.with_suffix(".legacy_conversion.partial.tif")
    try:
        write_output(
            output_path=temp_path,
            data=selected_data,
            dataset_valid_mask=dataset_valid_mask,
            bounds=bounds,
            epsg=epsg_from_mgrs_tile(row["mgrs_tile_t"]),
            row=row,
            target_date=target_date,
            start_date=start_date,
            end_inclusive=end_inclusive,
            valid_fraction_before_fill=valid_fraction,
            neighbour_fill_iterations=fill_iterations,
            neighbour_fill_complete=bool(dataset_valid_mask.all()),
        )
        os.replace(temp_path, output_path)
    finally:
        if temp_path.exists():
            temp_path.unlink()

    print(f"Converted existing legacy TIFF locally: {output_path}")
    return True


def inspect_existing_output(
    output_path: Path,
    row: pd.Series,
    target_date: pd.Timestamp,
) -> dict | None:
    """Validate a completed output and recover its summary for resume mode."""
    try:
        with rasterio.open(output_path) as src:
            expected_band_count = len(REFLECTANCE_BANDS) + 4
            if src.width != 200 or src.height != 200 or src.count != expected_band_count:
                return None
            expected_descriptions = tuple(REFLECTANCE_BANDS) + tuple(QUALITY_BANDS)
            if src.descriptions != expected_descriptions:
                return None
            if src.crs is None or src.crs.to_epsg() != epsg_from_mgrs_tile(row["mgrs_tile_t"]):
                return None

            tags = src.tags()
            if tags.get("aggregation_id") != str(row["aggregation_id"]):
                return None
            if tags.get("target_date") != target_date.date().isoformat():
                return None

            source_offsets = src.read(len(REFLECTANCE_BANDS) + 1)
            original_valid = src.read(len(REFLECTANCE_BANDS) + 4) == 1
            final_valid = src.dataset_mask() > 0
            reflectance = src.read(
                list(range(1, len(REFLECTANCE_BANDS) + 1))
            )
            if np.any(~np.isfinite(reflectance[:, final_valid])):
                return None
            counts = source_date_counts(target_date, source_offsets, original_valid)

            valid_fraction = float(original_valid.mean())
            return {
                "aggregation_id": str(row["aggregation_id"]),
                "target_date": target_date.date().isoformat(),
                "mgrs_tile": str(row["mgrs_tile_t"]),
                "output_path": str(output_path),
                "status": "skipped_existing",
                "error": "",
                "valid_fraction_before_fill": valid_fraction,
                "below_valid_fraction_0_98": valid_fraction < LOW_VALID_FRACTION_THRESHOLD,
                "invalid_pixels_before_fill": int((~original_valid).sum()),
                "neighbour_fill_iterations": int(
                    tags.get("neighbour_fill_iterations", "0")
                ),
                "neighbour_fill_complete": bool(final_valid.all()),
                "source_date_counts": counts,
            }
    except (OSError, ValueError, rasterio.errors.RasterioError):
        return None


def process_window(
    row: pd.Series,
    config: SHConfig,
    output_dir: Path,
    window_number: int,
    window_count: int,
) -> dict:
    target_date = row["Delta_Date"]
    start_date = target_date - timedelta(days=DAYS_EITHER_SIDE)
    end_inclusive = target_date + timedelta(days=DAYS_EITHER_SIDE)

    # The API time interval ends at midnight. Advancing by one day makes the
    # requested target + 8 calendar day fully inclusive. The CDSE service can
    # include an adjacent calendar day at the endpoint, which is acceptable
    # for this workflow.
    end_exclusive = end_inclusive + timedelta(days=1)

    epsg = epsg_from_mgrs_tile(row["mgrs_tile_t"])
    bounds = (
        float(row["cell_xmin"]),
        float(row["cell_ymin"]),
        float(row["cell_xmax"]),
        float(row["cell_ymax"]),
    )
    bbox = BBox(bounds, crs=f"EPSG:{epsg}")
    size = bbox_to_dimensions(bbox, resolution=RESOLUTION_M)

    expected_size = (
        round((bounds[2] - bounds[0]) / RESOLUTION_M),
        round((bounds[3] - bounds[1]) / RESOLUTION_M),
    )
    if size != expected_size:
        raise RuntimeError(
            f"Unexpected output size {size}; expected {expected_size} for the "
            f"manifest bounds at {RESOLUTION_M} m"
        )

    print("\n" + "=" * 80)
    print(f"WINDOW {window_number}/{window_count}")
    print("=" * 80)
    print(f"Aggregation ID: {row['aggregation_id']}")
    print(f"MGRS tile / CRS: {row['mgrs_tile_t']} / EPSG:{epsg}")
    print(f"Bounds: {bounds}")
    print(f"Output size: {size[0]} x {size[1]} pixels at {RESOLUTION_M} m")
    print(
        "Date search: "
        f"{start_date.date().isoformat()} through "
        f"{end_inclusive.date().isoformat()} (target {target_date.date().isoformat()})"
    )

    safe_id = str(row["aggregation_id"])
    stem = f"{safe_id}_target_{target_date.strftime('%Y%m%d')}_pm8d_nearest_clear"
    output_path = output_dir / f"{stem}.tif"
    existing_summary = None
    if output_path.exists():
        convert_legacy_output_if_needed(output_path, row, target_date)
        existing_summary = inspect_existing_output(output_path, row, target_date)
    if existing_summary is not None:
        print(f"Verified existing output; skipping download: {output_path}")
        return existing_summary

    data = request_composite(
        config=config,
        bbox=bbox,
        size=size,
        target_date=target_date,
        mgrs_tile=str(row["mgrs_tile_t"]),
        start_date=start_date,
        end_exclusive=end_exclusive,
    )
    filled_data, original_valid, final_valid, fill_iterations = fill_missing_reflectance(
        data
    )
    valid_fraction = float(original_valid.mean())
    counts = source_date_counts(target_date, data[:, :, -4], original_valid)

    temp_output_path = output_path.with_suffix(".partial.tif")
    try:
        write_output(
            output_path=temp_output_path,
            data=filled_data,
            dataset_valid_mask=final_valid,
            bounds=bounds,
            epsg=epsg,
            row=row,
            target_date=target_date,
            start_date=start_date,
            end_inclusive=end_inclusive,
            valid_fraction_before_fill=valid_fraction,
            neighbour_fill_iterations=fill_iterations,
            neighbour_fill_complete=bool(final_valid.all()),
        )
        os.replace(temp_output_path, output_path)
    finally:
        if temp_output_path.exists():
            temp_output_path.unlink()

    print(f"Wrote composite: {output_path}")
    print(f"Valid clear-pixel fraction before neighbour fill: {valid_fraction:.4f}")
    print(
        "Neighbour fill: "
        f"{int((~original_valid).sum())} input gaps, "
        f"{fill_iterations} iteration(s), complete={bool(final_valid.all())}"
    )
    if counts:
        print("Observed-pixel source dates:")
        for source_date, count in counts.items():
            fraction = count / original_valid.sum()
            print(f"  {source_date}: {count} pixels ({fraction:.2%})")
    else:
        print("No good-quality observations were found in the temporal window.")

    return {
        "aggregation_id": str(row["aggregation_id"]),
        "target_date": target_date.date().isoformat(),
        "mgrs_tile": str(row["mgrs_tile_t"]),
        "output_path": str(output_path),
        "status": "completed",
        "error": "",
        "valid_fraction_before_fill": valid_fraction,
        "below_valid_fraction_0_98": valid_fraction < LOW_VALID_FRACTION_THRESHOLD,
        "invalid_pixels_before_fill": int((~original_valid).sum()),
        "neighbour_fill_iterations": fill_iterations,
        "neighbour_fill_complete": bool(final_valid.all()),
        "source_date_counts": counts,
    }


def load_progress_records(progress_path: Path) -> dict[str, dict]:
    if not progress_path.exists():
        return {}

    try:
        progress = pd.read_csv(progress_path, dtype=str, keep_default_na=False)
    except (OSError, pd.errors.ParserError):
        return {}

    if "aggregation_id" not in progress.columns:
        return {}
    return {
        str(record["aggregation_id"]): record
        for record in progress.to_dict(orient="records")
    }


def write_progress_records(progress_path: Path, records: dict[str, dict]) -> None:
    columns = [
        "aggregation_id",
        "target_date",
        "mgrs_tile",
        "output_path",
        "status",
        "error",
        "valid_fraction_before_fill",
        "below_valid_fraction_0_98",
        "invalid_pixels_before_fill",
        "neighbour_fill_iterations",
        "neighbour_fill_complete",
        "source_date_counts",
    ]
    csv_records = []
    for record in records.values():
        csv_record = dict(record)
        counts = csv_record.get("source_date_counts", {})
        if isinstance(counts, dict):
            csv_record["source_date_counts"] = json.dumps(counts, sort_keys=True)
        csv_records.append(csv_record)

    progress_path.parent.mkdir(parents=True, exist_ok=True)
    temp_path = progress_path.with_suffix(".partial.csv")
    try:
        pd.DataFrame(csv_records, columns=columns).to_csv(temp_path, index=False)
        os.replace(temp_path, progress_path)
    finally:
        if temp_path.exists():
            temp_path.unlink()


def main() -> None:
    args = parse_args()
    rows = load_test_windows(args.aggregation_id, args.limit)
    config = build_cdse_config()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    progress_path = output_dir / "sentinel2_l2a_download_manifest.csv"
    progress_records = load_progress_records(progress_path)

    print(f"Processing {len(rows)} window(s)")
    print(f"Output directory: {output_dir}")
    print(f"Progress CSV: {progress_path}")
    print(f"Reflectance bands: {', '.join(REFLECTANCE_BANDS)}")
    print("Process response: scaled UINT16; local GeoTIFF reflectance: FLOAT32")
    print("Scene filter: use each manifest row's assigned MGRS tile")

    for row_index, (_, row) in enumerate(rows.iterrows(), start=1):
        try:
            summary = process_window(
                row=row,
                config=config,
                output_dir=output_dir,
                window_number=row_index,
                window_count=len(rows),
            )
        except Exception as exc:
            summary = {
                "aggregation_id": str(row["aggregation_id"]),
                "target_date": row["Delta_Date"].date().isoformat(),
                "mgrs_tile": str(row["mgrs_tile_t"]),
                "output_path": "",
                "status": "failed",
                "error": f"{type(exc).__name__}: {exc}",
                "valid_fraction_before_fill": "",
                "below_valid_fraction_0_98": "",
                "invalid_pixels_before_fill": "",
                "neighbour_fill_iterations": "",
                "neighbour_fill_complete": "",
                "source_date_counts": {},
            }
            print(
                f"FAILED {row['aggregation_id']}: "
                f"{type(exc).__name__}: {exc}"
            )

        progress_records[str(row["aggregation_id"])] = summary
        if row_index % PROGRESS_SAVE_EVERY == 0 or summary["status"] == "failed":
            write_progress_records(progress_path, progress_records)

    write_progress_records(progress_path, progress_records)

    processed_ids = set(rows["aggregation_id"].astype(str))
    run_records = [
        record
        for aggregation_id, record in progress_records.items()
        if aggregation_id in processed_ids
    ]
    completed = sum(record.get("status") == "completed" for record in run_records)
    skipped = sum(record.get("status") == "skipped_existing" for record in run_records)
    failed = sum(record.get("status") == "failed" for record in run_records)
    low_valid = sum(
        str(record.get("below_valid_fraction_0_98", "")).lower() == "true"
        for record in run_records
    )

    print("\n" + "=" * 80)
    print("DOWNLOAD SUMMARY")
    print("=" * 80)
    print(f"Downloaded this run: {completed}")
    print(f"Verified and skipped: {skipped}")
    print(f"Failed: {failed}")
    print(f"Pre-fill valid fraction below 0.98: {low_valid}")
    print(f"Run records: {len(run_records)} of {len(rows)}")
    print(f"Manifest: {progress_path}")


if __name__ == "__main__":
    main()
