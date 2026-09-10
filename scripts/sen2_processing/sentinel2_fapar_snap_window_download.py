"""Produce two 20 m FAPAR/flag TIFFs per manifest window using local ESA SNAP.

Downloads eight L2A reflectances, four angles and SCL/CLD/dataMask as UINT16.
Runs SNAP separately per acquisition in the inclusive +/-8 UTC calendar days.
Flag 0 has priority over all fallback flags (1, 2, 3), regardless of date.
Within each priority group: nearest date, future on a day-distance tie, then
lower CLD on the same day. No spatial gap filling or source-date bands.

Default: process every row of data/density_aggregation/sen2_spataggr_fapar/
density_cluster_4000m_aggregate_manifest.csv. Verified output pairs are skipped.
Requires the two existing Sentinel-2 Python scripts in this same directory,
the existing CDSE credentials, and ESA SNAP's Optical Toolbox/gpt executable.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import tempfile
from uuid import uuid4

import numpy as np
import pandas as pd
import rasterio
from rasterio.transform import from_bounds
from sentinelhub import BBox

import sentinel2_fapar_snap_one_window_test as snap


PROJECT_ROOT = snap.base.PROJECT_ROOT
MANIFEST_PATH = (PROJECT_ROOT / "data" / "density_aggregation" / "sen2_spataggr_fapar"
                 / "density_cluster_4000m_aggregate_manifest.csv")
OUTPUT_DIR = PROJECT_ROOT / "data" / "sentinel2_fapar_snap"
REFLECTANCE_BANDS = ("B03", "B04", "B05", "B06", "B07", "B8A", "B11", "B12")
ANGLE_BANDS = snap.ANGLE_BANDS
INPUT_BANDS = REFLECTANCE_BANDS + ANGLE_BANDS + ("SCL", "CLD", "dataMask")
REFLECTANCE_SCALE = 10000
ANGLE_SCALE = 100  # 0.01 degree precision; 360 degrees fits in UINT16.
ENCODED_NODATA = 65535
FAPAR_NODATA = -9999.0
FLAG_NODATA = 255
GOOD_SCL = (4, 5)
MAX_CLD = 40
ACCEPTED_FLAGS = (0, 1, 2, 3)
DAYS_EITHER_SIDE = 8
SIZE = (200, 200)
WORKFLOW_VERSION = "fapar-production-uint16-flag0-first-v1"
CHECKLIST_NAME = "fapar_download_checklist.csv"
SETTINGS = {
    "workflow_version": WORKFLOW_VERSION, "input_bands": INPUT_BANDS,
    "reflectance_scale": REFLECTANCE_SCALE, "angle_scale": ANGLE_SCALE,
    "good_scl": GOOD_SCL, "max_cld": MAX_CLD, "accepted_flags": ACCEPTED_FLAGS,
    "days_either_side": DAYS_EITHER_SIDE, "size": SIZE,
    "selection": "flag0_first_then_nearest_date_future_tie_lower_cld",
    "harmonizeValues": True, "resampling": "NEAREST", "gap_fill": "none",
}


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=MANIFEST_PATH)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    parser.add_argument("--gpt", default=os.environ.get("SNAP_GPT"))
    parser.add_argument("--aggregation-id", help="Optional single-window selection.")
    parser.add_argument("--limit", type=int, help="Optional first N selected rows.")
    parser.add_argument("--overwrite", action="store_true", help="Regenerate even verified output pairs.")
    return parser.parse_args()


def load_rows(args):
    required = ["aggregation_id", "Delta_Date", "mgrs_tile_t", "window_crs",
                "cell_xmin", "cell_ymin", "cell_xmax", "cell_ymax"]
    rows = pd.read_csv(args.manifest, usecols=required, dtype={"aggregation_id": str})
    if rows.empty or rows[required].isna().any().any():
        raise ValueError("Manifest is empty or contains missing required values")
    if rows.aggregation_id.duplicated().any():
        raise ValueError("Manifest aggregation_id values must be unique")
    # Use IDs verbatim as filenames, without silent sanitization/collisions.
    if not rows.aggregation_id.str.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*").all():
        raise ValueError("aggregation_id contains characters unsafe for output filenames")
    if args.aggregation_id:
        rows = rows.loc[rows.aggregation_id == args.aggregation_id]
        if rows.empty:
            raise ValueError(f"Unknown aggregation_id: {args.aggregation_id}")
    if args.limit is not None:
        if args.limit < 1:
            raise ValueError("--limit must be positive")
        rows = rows.head(args.limit)
    rows = rows.copy()
    rows["Delta_Date"] = pd.to_datetime(rows.Delta_Date, utc=True, errors="raise").dt.normalize()
    for name in required[4:]:
        rows[name] = pd.to_numeric(rows[name], errors="raise")
    # Validate all chosen bounds and CRSs before making any download request.
    for _, row in rows.iterrows():
        geometry(row)
    return rows.reset_index(drop=True)


def geometry(row):
    bounds = tuple(float(row[key]) for key in ("cell_xmin", "cell_ymin", "cell_xmax", "cell_ymax"))
    if not np.isfinite(bounds).all() or not np.allclose(
        [bounds[2] - bounds[0], bounds[3] - bounds[1]], [4000, 4000], rtol=0, atol=0.01
    ):
        raise ValueError(f"Invalid 4 km bounds for {row['aggregation_id']}: {bounds}")
    crs = rasterio.crs.CRS.from_user_input(row["window_crs"])
    if crs.to_epsg() != snap.base.epsg_from_mgrs_tile(row["mgrs_tile_t"]):
        raise ValueError(f"CRS/MGRS disagreement for {row['aggregation_id']}")
    return bounds, crs


def signature(row, bounds, crs):
    config = {**SETTINGS, "id": row["aggregation_id"], "date": row["Delta_Date"].isoformat(),
              "tile": row["mgrs_tile_t"], "bounds": bounds, "crs": crs.to_string()}
    return hashlib.sha256(json.dumps(config, sort_keys=True).encode()).hexdigest()


def make_evalscript(scene_id=None):
    discovery = scene_id is None
    bands = ("dataMask",) if discovery else INPUT_BANDS
    units = ["DN"] if discovery else (["REFLECTANCE"] * 8 + ["DEGREES"] * 4 + ["DN", "PERCENT", "DN"])
    filter_js = "" if discovery else f"""
function preProcessScenes(collections) {{
  collections.scenes.tiles = collections.scenes.tiles.filter(t => String(t.shId) === {json.dumps(str(scene_id))});
  return collections;
}}
"""
    pixels = "return [0];" if discovery else f"""
if (samples.length !== 1 || samples[0].dataMask !== 1)
  return new Array({len(INPUT_BANDS)}).fill({ENCODED_NODATA});
const s = samples[0];
return {json.dumps(list(INPUT_BANDS))}.map((b, i) =>
  encode(s[b], i < 8 ? {REFLECTANCE_SCALE} : (i < 12 ? {ANGLE_SCALE} : 1)));
"""
    return f"""//VERSION=3
function setup() {{
  return {{input: [{{bands: {json.dumps(list(bands))}, units: {json.dumps(units)}}}],
    output: {{id: "default", bands: {len(bands)}, sampleType: "UINT16", noDataValue: {ENCODED_NODATA}}},
    mosaicking: "TILE"}};
}}
function encode(value, scale) {{
  if (!Number.isFinite(value) || value < 0) return {ENCODED_NODATA};
  const encoded = Math.round(value * scale);
  // Mark overflow missing instead of silently clipping it to a plausible input.
  return encoded >= {ENCODED_NODATA} ? {ENCODED_NODATA} : encoded;
}}
{filter_js}
function evaluatePixel(samples) {{ {pixels} }}
function updateOutputMetadata(scenes, inputMetadata, outputMetadata) {{
  outputMetadata.userData = {{tiles: scenes.tiles}};
}}
"""


def download_scene(config, bbox, scene):
    day = pd.to_datetime(scene["date"], utc=True).normalize()
    raw, metadata = snap.process_request(
        config, bbox, SIZE, day, day + pd.Timedelta(days=1), make_evalscript(scene["scene_id"])
    )
    returned = metadata["tiles"]
    if len(returned) != 1 or str(returned[0].get("shId")) != scene["scene_id"]:
        raise RuntimeError("Download did not return exactly the selected scene")
    if pd.to_datetime(returned[0]["date"], utc=True).normalize() != day:
        raise RuntimeError("Returned acquisition date does not match selected date")
    if raw.shape != (SIZE[1], SIZE[0], len(INPUT_BANDS)) or raw.dtype != np.uint16:
        raise RuntimeError(f"Unexpected UINT16 response: {raw.shape}, {raw.dtype}")
    encoded = np.moveaxis(raw, -1, 0)
    missing = encoded == ENCODED_NODATA
    data = encoded.astype("float32")
    data[:8] /= REFLECTANCE_SCALE
    data[8:12] /= ANGLE_SCALE
    data[missing] = snap.NODATA
    return data


def pixel_counts(flags):
    counts = {f"flag_{value}_pixels": int((flags == value).sum()) for value in ACCEPTED_FLAGS}
    counts["valid_pixels"] = int((flags != FLAG_NODATA).sum())
    counts["nodata_pixels"] = int((flags == FLAG_NODATA).sum())
    counts["valid_fraction"] = counts["valid_pixels"] / flags.size
    # Denominator is the whole window, so this can be used for window filtering.
    counts["flag_0_fraction"] = counts["flag_0_pixels"] / flags.size
    return counts


def output_paths(output, window_id):
    return output / f"{window_id}_fapar.tif", output / f"{window_id}_snap_flag.tif"


def inspect_pair(paths, row, bounds, crs, fingerprint):
    if not all(path.is_file() for path in paths):
        return None
    try:
        with rasterio.open(paths[0]) as f, rasterio.open(paths[1]) as q:
            expected_transform = from_bounds(*bounds, *SIZE)
            for src, description, dtype, nodata in (
                (f, "FAPAR", "float32", FAPAR_NODATA),
                (q, "SNAP_flag", "uint8", FLAG_NODATA),
            ):
                if (src.count != 1 or src.shape != (SIZE[1], SIZE[0]) or src.crs != crs
                        or not src.transform.almost_equals(expected_transform)
                        or src.descriptions != (description,) or src.dtypes != (dtype,)
                        or src.nodata != nodata or src.tags().get("signature") != fingerprint
                        or src.tags().get("aggregation_id") != row["aggregation_id"]):
                    return None
            # A crash between the two replacements must not yield a valid pair.
            if not f.tags().get("pair_id") or f.tags()["pair_id"] != q.tags().get("pair_id"):
                return None
            values, flags = f.read(1), q.read(1)
            good = flags != FLAG_NODATA
            if (not np.isin(flags, (*ACCEPTED_FLAGS, FLAG_NODATA)).all()
                    or not np.array_equal(values != FAPAR_NODATA, good)
                    or not np.isfinite(values[good]).all()
                    or np.any((values[good] < 0) | (values[good] > 1))):
                return None
            return {**pixel_counts(flags), "acquisitions": int(f.tags()["acquisitions"]),
                    "api_candidates": int(f.tags()["api_candidates"])}
    except (OSError, ValueError, KeyError, rasterio.errors.RasterioError):
        return None


def write_pair(paths, values, flags, row, bounds, crs, fingerprint, candidates, acquisitions):
    tags = {"signature": fingerprint, "pair_id": uuid4().hex,
            "aggregation_id": row["aggregation_id"], "target_date": row["Delta_Date"].date().isoformat(),
            "workflow_version": WORKFLOW_VERSION, "selection": SETTINGS["selection"],
            "accepted_flags": "0,1,2,3", "good_scl": "4,5", "max_cld_percent": MAX_CLD,
            "api_candidates": candidates, "acquisitions": acquisitions,
            "gap_fill": "none", "days_either_side": DAYS_EITHER_SIDE}
    temporary = [path.with_suffix(".partial.tif") for path in paths]
    try:
        for path, array, dtype, nodata, name in (
            (temporary[0], values, "float32", FAPAR_NODATA, "FAPAR"),
            (temporary[1], flags, "uint8", FLAG_NODATA, "SNAP_flag"),
        ):
            with rasterio.open(
                path, "w", driver="GTiff", width=SIZE[0], height=SIZE[1], count=1,
                dtype=dtype, crs=crs, transform=from_bounds(*bounds, *SIZE),
                nodata=nodata, compress="deflate", predictor=3 if dtype == "float32" else 2,
                tiled=True,
            ) as dst:
                dst.write(array.astype(dtype, copy=False), 1)
                dst.set_band_description(1, name)
                dst.update_tags(**{key: str(value) for key, value in tags.items()})
                dst.update_tags(1, units="fraction" if name == "FAPAR" else "bit_flags")
        for src, dst in zip(temporary, paths):
            os.replace(src, dst)
    finally:
        for path in temporary:
            path.unlink(missing_ok=True)


def process_window(row, config, gpt, output):
    bounds, crs = geometry(row)
    fingerprint = signature(row, bounds, crs)
    target = row["Delta_Date"]
    bbox = BBox(bounds, crs=crs.to_string())
    start = target - pd.Timedelta(days=DAYS_EITHER_SIDE)
    end = target + pd.Timedelta(days=DAYS_EITHER_SIDE + 1)
    _, metadata = snap.process_request(config, bbox, (4, 4), start, end, make_evalscript())
    scenes, _ = snap.select_acquisitions(metadata["tiles"], target, row["mgrs_tile_t"])
    print(f"  API candidates: {len(metadata['tiles'])}; unique acquisitions: {len(scenes)}", flush=True)
    if any(scene["sensor"] not in ("S2A", "S2B") for scene in scenes):
        raise RuntimeError("Only S2A/S2B SNAP models are supported; no automatic sensor substitution")
    result = np.full((SIZE[1], SIZE[0]), FAPAR_NODATA, dtype="float32")
    selected_flags = np.full(result.shape, FLAG_NODATA, dtype="uint8")
    best_priority = np.full(result.shape, 2, dtype="uint8")
    best_distance = np.full(result.shape, 999, dtype="int16")
    best_offset = np.full(result.shape, -999, dtype="int16")
    best_cloud = np.full(result.shape, np.inf, dtype="float32")
    for index, scene in enumerate(scenes, 1):
        print(f"  [{index}/{len(scenes)}] {scene['date']} {scene['sensor']} offset={scene['day_offset']:+d}", flush=True)
        data = download_scene(config, bbox, scene)
        cloud = data[INPUT_BANDS.index("CLD")]
        clear = (snap.numeric_validity(data, INPUT_BANDS)
                 & np.isin(data[INPUT_BANDS.index("SCL")], GOOD_SCL)
                 & np.isfinite(cloud) & (cloud >= 0) & (cloud <= MAX_CLD))
        if not clear.any():
            print("    No clear land pixels; skipping local SNAP calculation", flush=True)
            continue
        # Temporary files are outside the final output directory. The context
        # removes only the directory it created, including after a SNAP error.
        with tempfile.TemporaryDirectory(prefix="s2_fapar_snap_") as work:
            fapar, flags, _ = snap.run_snap(gpt, Path(work), data, bounds, crs,
                                          scene["sensor"], input_bands=INPUT_BANDS)
        usable = clear & np.isfinite(fapar) & (fapar >= 0) & (fapar <= 1) & np.isin(flags, ACCEPTED_FLAGS)
        priority = np.where(flags == 0, 0, 1).astype("uint8")
        offset = scene["day_offset"]
        distance = abs(offset)
        nearer = ((distance < best_distance)
                  | ((distance == best_distance) & (offset > best_offset))
                  | ((distance == best_distance) & (offset == best_offset) & (cloud < best_cloud)))
        better = usable & ((priority < best_priority) | ((priority == best_priority) & nearer))
        # Any eligible flag-0 candidate replaces a fallback; a fallback can
        # NEVER replace a flag-0 candidate, even if it is closer to the target.
        result[better] = fapar[better]
        selected_flags[better] = flags[better].astype("uint8")
        best_priority[better] = priority[better]
        best_distance[better] = distance
        best_offset[better] = offset
        best_cloud[better] = cloud[better]
        print(f"    Clear land: {clear.mean():.2%}; accepted FAPAR: {usable.mean():.2%}; "
              f"flag 0: {(usable & (flags == 0)).mean():.2%}", flush=True)
    paths = output_paths(output, row["aggregation_id"])
    write_pair(paths, result, selected_flags, row, bounds, crs, fingerprint,
               len(metadata["tiles"]), len(scenes))
    summary = inspect_pair(paths, row, bounds, crs, fingerprint)
    if summary is None:
        raise RuntimeError("Final TIFF pair failed read-back validation")
    print(f"  Final valid: {summary['valid_fraction']:.2%}; flag 0: {summary['flag_0_fraction']:.2%}", flush=True)
    return summary


def load_checklist(path):
    if not path.exists():
        return {}
    # Do not silently discard a damaged checklist; output files remain safe.
    frame = pd.read_csv(path, dtype=str, keep_default_na=False)
    if "aggregation_id" not in frame or frame.aggregation_id.duplicated().any():
        raise ValueError(f"Invalid progress checklist: {path}")
    return {item["aggregation_id"]: item for item in frame.to_dict("records")}


def save_checklist(path, records):
    temporary = path.with_suffix(".partial.csv")
    pd.DataFrame(records.values()).to_csv(temporary, index=False)
    os.replace(temporary, path)


def main():
    args = parse_args()
    rows = load_rows(args)
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    checklist = output / CHECKLIST_NAME
    records = load_checklist(checklist)
    config = gpt = None  # Initialize only if a window actually needs processing.
    failures = 0
    print(f"Manifest: {args.manifest.resolve()}\nWindows: {len(rows)}\nOutput: {output}", flush=True)
    print("UINT16 download; flag 0 first, then flags 1/2/3; no spatial gap filling", flush=True)
    for number, (_, row) in enumerate(rows.iterrows(), 1):
        window_id = row["aggregation_id"]
        paths = output_paths(output, window_id)
        bounds, crs = geometry(row)
        print(f"\nWINDOW {number}/{len(rows)}: {window_id}", flush=True)
        record = {"aggregation_id": window_id, "target_date": row["Delta_Date"].date().isoformat(),
                  "mgrs_tile": row["mgrs_tile_t"], "fapar_path": str(paths[0]),
                  "snap_flag_path": str(paths[1]), "status": "processing", "error": ""}
        existing = None if args.overwrite else inspect_pair(paths, row, bounds, crs, signature(row, bounds, crs))
        if existing is not None:
            record.update(existing, status="skipped_existing")
            print("  Verified output pair; skipped", flush=True)
        else:
            if config is None:
                # Fail before downloading if local SNAP or credentials are unavailable.
                gpt = snap.find_gpt(args.gpt)
                with tempfile.TemporaryDirectory(prefix="s2_fapar_help_") as work:
                    help_text = snap.run_command([gpt, "BiophysicalOp", "-h"], Path(work) / "help.txt", 120)
                if any(name not in help_text for name in ("computeFapar", "sensor")):
                    raise RuntimeError("Installed SNAP BiophysicalOp lacks required parameters")
                config = snap.base.build_cdse_config()
            records[window_id] = record.copy()
            save_checklist(checklist, records)
            try:
                summary = process_window(row, config, gpt, output)
                record.update(summary, status="completed" if summary["valid_pixels"] else "completed_empty")
            except KeyboardInterrupt:
                record.update(status="interrupted", error="Interrupted by user; incomplete window will be retried")
                records[window_id] = record
                save_checklist(checklist, records)
                raise
            except Exception as exc:
                record.update(status="failed", error=f"{type(exc).__name__}: {exc}")
                failures += 1
                print(f"  FAILED: {record['error']}", flush=True)
        records[window_id] = record
        save_checklist(checklist, records)
    print(f"\nFinished {len(rows)} windows; failures: {failures}\nChecklist: {checklist}", flush=True)
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
