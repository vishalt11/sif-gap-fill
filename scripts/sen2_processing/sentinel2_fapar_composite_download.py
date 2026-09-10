"""Download one coherent nearest-clear Sentinel-2 composite per window for SNAP.

Download only: no SNAP executable is required or launched. Uses the existing
manifest/authentication/request helpers; leave those scripts alongside this one.
Default manifest: data/density_aggregation/sen2_spataggr_fapar/
density_cluster_4000m_aggregate_manifest.csv.

Server-side selection: nearest clear UTC calendar date within +/-8 days,
future on an equal-distance tie, then lower cloud probability on the same day.
Eight reflectances, four angles and three QA values come from ONE acquisition
per pixel. No spatial gap filling and no SNAP-flag-based temporal selection.
The final sensor_id band identifies the network needed later: 1=S2A, 2=S2B.

One 16-band UINT16 TIFF per window, 200 x 200 pixels at 20 m in its manifest
UTM CRS. Mask 65535 before decoding: reflectances /10000, angles /100, QA and
sensor_id unscaled. A mixed-sensor TIFF must NOT be processed with one sensor's
SNAP network across all pixels; apply each network to its corresponding pixels.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from rasterio.transform import from_bounds
from sentinelhub import BBox

import sentinel2_fapar_snap_window_download as previous


MANIFEST_PATH = previous.MANIFEST_PATH
OUTPUT_DIR = previous.PROJECT_ROOT / "data" / "sentinel2_fapar_composite_inputs"
INPUT_BANDS = previous.INPUT_BANDS
BANDS = INPUT_BANDS + ("sensor_id",)
SIZE = previous.SIZE
NODATA = previous.ENCODED_NODATA
SCALES = (1 / previous.REFLECTANCE_SCALE,) * 8 + (1 / previous.ANGLE_SCALE,) * 4 + (1.0,) * 4
SENSOR_CODES = {"S2A": 1, "S2B": 2}
VERSION = "s2-fapar-coherent-nearest-clear-composite-v1"
SETTINGS = {
    "version": VERSION, "bands": BANDS, "scales": SCALES, "nodata": NODATA,
    "sensor_codes": SENSOR_CODES, "size": SIZE,
    "good_scl": previous.GOOD_SCL, "max_cld_percent": previous.MAX_CLD,
    "days_either_side": previous.DAYS_EITHER_SIDE,
    "selection": "nearest_UTC_calendar_date_future_tie_lower_CLD_then_scene_order",
    "scene_order": "acquisition_timestamp_then_scene_id",
    "tile_selection": "assigned_MGRS_only_latest_reprocessing_per_acquisition",
    "harmonizeValues": True, "resampling": "NEAREST", "gap_fill": "none",
    "snap_executed": False,
}


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=MANIFEST_PATH)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    parser.add_argument("--aggregation-id", help="Download exactly one selected window.")
    parser.add_argument("--limit", type=int, help="Download the first N selected windows.")
    parser.add_argument("--overwrite", action="store_true",
                        help="Re-query and redownload even verified composites; consumes PU.")
    return parser.parse_args()


def window_definition(row, bounds, crs):
    return {
        "aggregation_id": row["aggregation_id"],
        "target_date": row["Delta_Date"].date().isoformat(),
        "mgrs_tile": row["mgrs_tile_t"], "bounds": bounds,
        "crs": crs.to_string(), "settings": SETTINGS,
    }


def signature(definition):
    return hashlib.sha256(json.dumps(definition, sort_keys=True).encode("utf-8")).hexdigest()


def make_evalscript(scenes):
    # Only known, deduplicated acquisitions are admitted. The mapping is keyed
    # by shId rather than array position, which the API may reorder.
    selected = {
        str(scene["scene_id"]): {
            "offset": int(scene["day_offset"]),
            "sensor": SENSOR_CODES[scene["sensor"]], "rank": rank,
        }
        for rank, scene in enumerate(scenes)
    }
    units = ["REFLECTANCE"] * 8 + ["DEGREES"] * 4 + ["DN", "PERCENT", "DN"]
    return f"""//VERSION=3
const allowed = {json.dumps(selected)};
const bands = {json.dumps(list(INPUT_BANDS))};
const missing = {NODATA};
function setup() {{
  return {{input: [{{bands: bands, units: {json.dumps(units)}}}],
    output: {{id: "default", bands: {len(BANDS)}, sampleType: "UINT16", nodataValue: missing}},
    mosaicking: "TILE"}};
}}
function preProcessScenes(collections) {{
  collections.scenes.tiles = collections.scenes.tiles.filter(
    t => Object.prototype.hasOwnProperty.call(allowed, String(t.shId)));
  return collections;
}}
function encode(value, scale) {{
  if (!Number.isFinite(value) || value < 0) return missing;
  const encoded = Math.round(value * scale);
  return encoded >= missing ? missing : encoded;
}}
function clearAndComplete(s) {{
  if (s.dataMask !== 1 || {json.dumps(list(previous.GOOD_SCL))}.indexOf(s.SCL) < 0
      || !Number.isFinite(s.CLD) || s.CLD < 0 || s.CLD > {previous.MAX_CLD}) return false;
  // These are physical-input checks, NOT SNAP's training-domain validation.
  if (!(s.sunZenithAngles >= 0 && s.sunZenithAngles < 90
        && s.viewZenithMean >= 0 && s.viewZenithMean < 90
        && s.sunAzimuthAngles >= 0 && s.sunAzimuthAngles <= 360
        && s.viewAzimuthMean >= 0 && s.viewAzimuthMean <= 360)) return false;
  for (let j = 0; j < 12; j++) {{
    if (encode(s[bands[j]], j < 8 ? {previous.REFLECTANCE_SCALE} : {previous.ANGLE_SCALE}) === missing)
      return false;
  }}
  // Reject values that would round to an invalid 90-degree zenith after packing.
  return Math.round(s.sunZenithAngles * {previous.ANGLE_SCALE}) < 90 * {previous.ANGLE_SCALE}
      && Math.round(s.viewZenithMean * {previous.ANGLE_SCALE}) < 90 * {previous.ANGLE_SCALE};
}}
function evaluatePixel(samples, scenes) {{
  let best = -1, bestInfo = null, bestCloud = Infinity;
  for (let i = 0; i < samples.length; i++) {{
    const info = allowed[String(scenes.tiles[i].shId)];
    const s = samples[i];
    if (!info || !clearAndComplete(s)) continue;
    const distance = Math.abs(info.offset);
    const oldDistance = bestInfo ? Math.abs(bestInfo.offset) : Infinity;
    const better = best < 0 || distance < oldDistance
      || (distance === oldDistance && info.offset > bestInfo.offset)
      || (distance === oldDistance && info.offset === bestInfo.offset && s.CLD < bestCloud)
      || (distance === oldDistance && info.offset === bestInfo.offset
          && s.CLD === bestCloud && info.rank < bestInfo.rank);
    if (better) {{ best = i; bestInfo = info; bestCloud = s.CLD; }}
  }}
  if (best < 0) return new Array({len(BANDS)}).fill(missing);
  // Copy every channel from this ONE sample; never select bands independently.
  const s = samples[best];
  const result = bands.map((b, j) => encode(s[b],
    j < 8 ? {previous.REFLECTANCE_SCALE} : (j < 12 ? {previous.ANGLE_SCALE} : 1)));
  result.push(bestInfo.sensor);
  return result;
}}
function updateOutputMetadata(scenes, inputMetadata, outputMetadata) {{
  outputMetadata.userData = {{tiles: scenes.tiles}};
}}
"""


def validate_pixels(data):
    if data.shape != (len(BANDS), SIZE[1], SIZE[0]) or data.dtype != np.uint16:
        raise ValueError(f"Unexpected composite array: {data.shape}, {data.dtype}")
    missing = data == NODATA
    # Selected pixels have every channel present; unselected ones have none.
    if not np.array_equal(missing.any(axis=0), missing.all(axis=0)):
        raise ValueError("Composite contains partially missing input vectors")
    valid = ~missing[0]
    if (not np.isin(data[BANDS.index("sensor_id")][valid], tuple(SENSOR_CODES.values())).all()
            or not (data[BANDS.index("dataMask")][valid] == 1).all()
            or not np.isin(data[BANDS.index("SCL")][valid], previous.GOOD_SCL).all()
            or np.any(data[BANDS.index("CLD")][valid] > previous.MAX_CLD)):
        raise ValueError("Composite violates sensor/cloud/land selection rules")
    for name, limit, inclusive in (
        ("sunZenithAngles", 90, False), ("viewZenithMean", 90, False),
        ("sunAzimuthAngles", 360, True), ("viewAzimuthMean", 360, True),
    ):
        values = data[BANDS.index(name)][valid]
        maximum = limit * previous.ANGLE_SCALE
        if np.any(values > maximum if inclusive else values >= maximum):
            raise ValueError(f"Invalid packed acquisition angles: {name}")
    return {
        "clear_pixels": int(valid.sum()), "clear_fraction": float(valid.mean()),
        "nodata_pixels": int((~valid).sum()),
        "s2a_pixels": int((data[-1] == SENSOR_CODES["S2A"]).sum()),
        "s2b_pixels": int((data[-1] == SENSOR_CODES["S2B"]).sum()),
    }


def inspect_composite(path, definition, bounds, crs):
    if not path.is_file():
        return None
    try:
        with rasterio.open(path) as src:
            if (src.count != len(BANDS) or src.shape != (SIZE[1], SIZE[0])
                    or src.dtypes != ("uint16",) * len(BANDS) or src.nodata != NODATA
                    or src.descriptions != BANDS or src.crs != crs
                    or not src.transform.almost_equals(from_bounds(*bounds, *SIZE))
                    or not np.allclose(src.scales, SCALES, rtol=0, atol=1e-12)
                    or src.offsets != (0.0,) * len(BANDS)
                    or src.tags().get("signature") != signature(definition)):
                return None
            tags = src.tags()
            sources = json.loads(tags["source_acquisitions_json"])
            if not isinstance(sources, list) or len(sources) != int(tags["acquisitions"]):
                return None
            summary = validate_pixels(src.read())  # Also detects corrupt TIFF blocks.
            source_sensors = {scene["sensor"] for scene in sources}
            if ((summary["s2a_pixels"] and "S2A" not in source_sensors)
                    or (summary["s2b_pixels"] and "S2B" not in source_sensors)):
                return None
            return {**summary, "api_candidates": int(tags["api_candidates"]),
                    "acquisitions": len(sources), "file_bytes": path.stat().st_size}
    except (OSError, ValueError, KeyError, TypeError, rasterio.errors.RasterioError):
        return None


def write_composite(path, data, scenes, candidates, definition, bounds, crs):
    validate_pixels(data)
    # Retain the source inventory in TIFF metadata, without extra per-window
    # files. This is NOT pixel-level date provenance; only sensor_id is spatial.
    source_keys = ("scene_id", "sensor", "date", "day_offset", "production_timestamp", "acquisition_key")
    sources = [{key: scene[key] for key in source_keys} for scene in scenes]
    temporary = path.with_suffix(".partial.tif")
    try:
        with rasterio.open(
            temporary, "w", driver="GTiff", width=SIZE[0], height=SIZE[1],
            count=len(BANDS), dtype="uint16", nodata=NODATA, crs=crs,
            transform=from_bounds(*bounds, *SIZE), compress="deflate", predictor=2, tiled=True,
        ) as dst:
            dst.write(data)
            dst.scales = SCALES
            dst.offsets = (0.0,) * len(BANDS)
            for index, (name, scale) in enumerate(zip(BANDS, SCALES), 1):
                units = "reflectance" if index <= 8 else "degrees" if index <= 12 else "percent" if name == "CLD" else "code"
                dst.set_band_description(index, name)
                dst.update_tags(index, units=units, scale=str(scale),
                                decoding="physical_value = stored_value * scale",
                                nodata_before_scaling=str(NODATA))
            dst.update_tags(len(BANDS), codes="1=S2A;2=S2B;65535=NoData")
            dst.update_tags(
                signature=signature(definition), workflow_version=VERSION,
                aggregation_id=definition["aggregation_id"], target_date=definition["target_date"],
                mgrs_tile=definition["mgrs_tile"], selection=SETTINGS["selection"],
                window_definition_json=json.dumps(definition, sort_keys=True),
                source_acquisitions_json=json.dumps(sources, sort_keys=True),
                api_candidates=str(candidates), acquisitions=str(len(scenes)),
                sensor_codes_json=json.dumps(SENSOR_CODES), gap_fill="none",
                snap_executed="false", angle_source_resolution_m="5000",
                requires_sensor_specific_snap="true", pixel_source_date_retained="false",
            )
        if inspect_composite(temporary, definition, bounds, crs) is None:
            raise RuntimeError("Composite failed TIFF read-back validation")
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def download_composite(row, path, definition, bounds, crs, config):
    target = row["Delta_Date"]
    start = target - pd.Timedelta(days=previous.DAYS_EITHER_SIDE)
    end = target + pd.Timedelta(days=previous.DAYS_EITHER_SIDE + 1)
    bbox = BBox(bounds, crs=crs.to_string())
    # A tiny 4x4 dataMask request discovers scene metadata at supported 1 km
    # resolution. Python then applies the same strict acquisition filter as the
    # previous downloader before the SINGLE full-resolution composite request.
    _, metadata = previous.snap.process_request(
        config, bbox, (4, 4), start, end, previous.make_evalscript()
    )
    scenes, _ = previous.snap.select_acquisitions(metadata["tiles"], target, row["mgrs_tile_t"])
    candidates = len(metadata["tiles"])
    print(f"  API candidates: {candidates}; unique assigned-tile acquisitions: {len(scenes)}", flush=True)
    unsupported = sorted({scene["sensor"] for scene in scenes} - SENSOR_CODES.keys())
    if unsupported:
        raise ValueError(f"Unsupported SNAP sensor(s): {unsupported}; refusing to substitute S2A/S2B")
    if len({str(scene["scene_id"]) for scene in scenes}) != len(scenes):
        raise ValueError("Duplicate scene IDs after acquisition selection")
    if scenes:
        raw, returned = previous.snap.process_request(config, bbox, SIZE, start, end, make_evalscript(scenes))
        expected = {str(scene["scene_id"]): pd.to_datetime(scene["date"], utc=True) for scene in scenes}
        actual = returned["tiles"]
        # Do not silently write a composite built from an incomplete/different
        # inventory if sources changed between discovery and download.
        if len(actual) != len(expected) or {str(scene.get("shId")) for scene in actual} != set(expected):
            raise RuntimeError("Composite response does not contain exactly the selected acquisitions")
        if any(pd.to_datetime(scene["date"], utc=True) != expected[str(scene["shId"])] for scene in actual):
            raise RuntimeError("Composite response acquisition timestamps differ from discovery")
        if raw.shape != (SIZE[1], SIZE[0], len(BANDS)) or raw.dtype != np.uint16:
            raise RuntimeError(f"Unexpected UINT16 response: {raw.shape}, {raw.dtype}")
        data = np.moveaxis(raw, -1, 0)
    else:
        data = np.full((len(BANDS), SIZE[1], SIZE[0]), NODATA, dtype="uint16")
    write_composite(path, data, scenes, candidates, definition, bounds, crs)
    summary = inspect_composite(path, definition, bounds, crs)
    if summary is None:
        raise RuntimeError("Saved composite failed validation")
    return summary


def main():
    args = parse_args()
    rows = previous.load_rows(args)
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    checklist = output / "composite_download_checklist.csv"
    records = previous.load_checklist(checklist)
    config = None
    failures = 0
    print(f"DOWNLOAD ONLY: {len(rows)} windows\nManifest: {args.manifest.resolve()}\nOutput: {output}", flush=True)
    print("One 16-band UINT16 TIFF per window; nearest-clear server composite; no SNAP or gap filling", flush=True)
    for number, (_, row) in enumerate(rows.iterrows(), 1):
        window_id = row["aggregation_id"]
        path = output / f"{window_id}_snap_inputs.tif"
        bounds, crs = previous.geometry(row)
        definition = window_definition(row, bounds, crs)
        print(f"\nWINDOW {number}/{len(rows)}: {window_id}", flush=True)
        record = {"aggregation_id": window_id, "target_date": definition["target_date"],
                  "mgrs_tile": definition["mgrs_tile"], "composite_path": str(path), "error": ""}
        summary = None if args.overwrite else inspect_composite(path, definition, bounds, crs)
        if summary is not None:
            record.update(summary, status="skipped_existing")
            print("  Verified existing composite; skipped", flush=True)
        else:
            records[window_id] = {**record, "status": "downloading"}
            previous.save_checklist(checklist, records)
            try:
                if config is None:
                    config = previous.snap.base.build_cdse_config()
                summary = download_composite(row, path, definition, bounds, crs, config)
                status = "completed" if summary["clear_pixels"] else "completed_empty"
                record.update(summary, status=status)
            except KeyboardInterrupt:
                records[window_id] = {**record, "status": "interrupted", "error": "Interrupted; rerun to resume"}
                previous.save_checklist(checklist, records)
                raise
            except Exception as exc:
                record.update(status="failed", error=f"{type(exc).__name__}: {exc}")
                failures += 1
                print(f"  FAILED: {record['error']}", flush=True)
        if summary is not None:
            print(f"  Clear input coverage: {summary['clear_fraction']:.2%}; "
                  f"S2A pixels: {summary['s2a_pixels']}; S2B pixels: {summary['s2b_pixels']}", flush=True)
        records[window_id] = record
        previous.save_checklist(checklist, records)
    print(f"\nFinished {len(rows)} windows; failures: {failures}\nChecklist: {checklist}", flush=True)
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
