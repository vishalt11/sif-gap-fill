"""Download and retain Sentinel-2 acquisition chips for later SNAP processing.

No SNAP executable is needed or launched. Uses the production downloader's
manifest reader, CDSE request builder and scene-selection helpers without
changing either previous downloader. Keep the existing scripts alongside this.

Default: all rows in data/density_aggregation/sen2_spataggr_fapar/
density_cluster_4000m_aggregate_manifest.csv, inclusive +/-8 UTC calendar days,
assigned MGRS tile only, latest reprocessing per acquisition.

Each acquisition is stored as one 15-band UINT16 GeoTIFF: eight reflectances,
four angles, SCL, CLD and dataMask. Reflectance = stored / 10000; angles =
stored / 100; QA bands are unscaled. Mask 65535 BEFORE applying these scales.
Cloudy acquisitions are retained. FAPAR and temporal selection are deferred.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re

import numpy as np
import pandas as pd
import rasterio
from rasterio.transform import from_bounds
from sentinelhub import BBox

import sentinel2_fapar_snap_window_download as previous


MANIFEST_PATH = previous.MANIFEST_PATH
OUTPUT_DIR = previous.PROJECT_ROOT / "data" / "sentinel2_fapar_acquisitions"
BANDS = previous.INPUT_BANDS
SIZE = previous.SIZE
NODATA = previous.ENCODED_NODATA
SCALES = (1 / previous.REFLECTANCE_SCALE,) * 8 + (1 / previous.ANGLE_SCALE,) * 4 + (1.0,) * 3
VERSION = "s2-fapar-acquisition-cache-v1"
SETTINGS = {
    "version": VERSION, "bands": BANDS, "scale_factors": SCALES,
    "nodata": NODATA, "size": SIZE, "days_either_side": previous.DAYS_EITHER_SIDE,
    "harmonizeValues": True, "upsampling": "NEAREST", "downsampling": "NEAREST",
}


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=MANIFEST_PATH)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    parser.add_argument("--aggregation-id", help="Download one selected window.")
    parser.add_argument("--limit", type=int, help="Download the first N selected windows.")
    parser.add_argument("--refresh-scenes", action="store_true",
                        help="Query current scenes again; still reuse valid acquisition TIFFs.")
    parser.add_argument("--overwrite", action="store_true",
                        help="Redownload acquisition TIFFs even when valid; consumes download PU.")
    return parser.parse_args()


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode("utf-8")).hexdigest()


def save_json(path, value):
    temporary = path.with_suffix(".partial.json")
    try:
        temporary.write_text(json.dumps(value, indent=2), encoding="utf-8")
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def window_definition(row, bounds, crs):
    target = row["Delta_Date"]
    days = previous.DAYS_EITHER_SIDE
    return {
        "aggregation_id": row["aggregation_id"], "target_date": target.date().isoformat(),
        "mgrs_tile": row["mgrs_tile_t"], "bounds": bounds, "crs": crs.to_string(),
        "start_date_inclusive": (target - pd.Timedelta(days=days)).date().isoformat(),
        "end_date_inclusive": (target + pd.Timedelta(days=days)).date().isoformat(),
        "settings": SETTINGS,
    }


def get_inventory(path, definition, bbox, config, refresh):
    signature = digest(definition)
    if path.is_file() and not refresh:
        try:
            inventory = json.loads(path.read_text(encoding="utf-8"))
            if (inventory.get("signature") == signature
                    and isinstance(inventory.get("selected_acquisitions"), list)):
                print("  Reusing saved acquisition inventory", flush=True)
                return inventory
        except (ValueError, OSError):
            pass
    start = pd.Timestamp(definition["start_date_inclusive"], tz="UTC")
    end = pd.Timestamp(definition["end_date_inclusive"], tz="UTC") + pd.Timedelta(days=1)
    _, metadata = previous.snap.process_request(
        config(), bbox, (4, 4), start, end, previous.make_evalscript()
    )
    scenes, audit = previous.snap.select_acquisitions(
        metadata["tiles"], pd.Timestamp(definition["target_date"], tz="UTC"), definition["mgrs_tile"]
    )
    inventory = {
        "signature": signature, "window": definition,
        "api_candidates": len(metadata["tiles"]), "selected_acquisitions": scenes,
        "scene_selection_audit": audit, "download_status": "pending",
    }
    # Persist discovery before downloads, so interruptions do not repeat this request.
    save_json(path, inventory)
    return inventory


def acquisition_path(output, window_id, scene):
    scene_id = str(scene["scene_id"])
    sensor = str(scene["sensor"])
    if not re.fullmatch(r"\d+", scene_id) or not re.fullmatch(r"S2[A-Z]", sensor):
        raise ValueError("Unexpected scene ID/sensor in acquisition inventory")
    stamp = pd.to_datetime(scene["date"], utc=True).strftime("%Y%m%dT%H%M%S")
    return output / "acquisitions" / f"{window_id}_{stamp}_{sensor}_{scene_id}.tif"


def inspect_acquisition(path, fingerprint, bounds, crs):
    if not path.is_file():
        return None
    try:
        with rasterio.open(path) as src:
            if (src.shape != (SIZE[1], SIZE[0]) or src.count != len(BANDS)
                    or src.descriptions != BANDS or src.dtypes != ("uint16",) * len(BANDS)
                    or src.nodata != NODATA or src.crs != crs
                    or not src.transform.almost_equals(from_bounds(*bounds, *SIZE))
                    or not np.allclose(src.scales, SCALES, rtol=0, atol=1e-12)
                    or src.tags().get("cache_signature") != fingerprint):
                return None
            # Read every band to catch truncated/corrupt compressed image blocks.
            data = src.read()
            mask = data[BANDS.index("dataMask")]
            if not np.isin(mask, (0, 1, NODATA)).all():
                return None
            if not src.tags().get("source_scene_json"):
                return None
            return {"file_bytes": path.stat().st_size, "data_pixels": int((mask == 1).sum())}
    except (OSError, ValueError, rasterio.errors.RasterioError):
        return None


def request_acquisition(config, bbox, scene):
    day = pd.to_datetime(scene["date"], utc=True).normalize()
    raw, metadata = previous.snap.process_request(
        config(), bbox, SIZE, day, day + pd.Timedelta(days=1),
        previous.make_evalscript(scene["scene_id"]),
    )
    returned = metadata["tiles"]
    if len(returned) != 1 or str(returned[0].get("shId")) != str(scene["scene_id"]):
        raise RuntimeError("Response does not contain exactly the requested acquisition")
    if pd.to_datetime(returned[0]["date"], utc=True) != pd.to_datetime(scene["date"], utc=True):
        raise RuntimeError("Returned acquisition timestamp differs from saved inventory; try --refresh-scenes")
    if raw.shape != (SIZE[1], SIZE[0], len(BANDS)) or raw.dtype != np.uint16:
        raise RuntimeError(f"Unexpected UINT16 response: {raw.shape}, {raw.dtype}")
    return np.moveaxis(raw, -1, 0), returned[0]


def write_acquisition(path, data, returned, scene, definition, fingerprint, bounds, crs):
    temporary = path.with_suffix(".partial.tif")
    try:
        with rasterio.open(
            temporary, "w", driver="GTiff", width=SIZE[0], height=SIZE[1],
            count=len(BANDS), dtype="uint16", nodata=NODATA, crs=crs,
            transform=from_bounds(*bounds, *SIZE), compress="deflate", predictor=2,
            tiled=True,
        ) as dst:
            dst.write(data)
            dst.scales = SCALES
            dst.offsets = (0.0,) * len(BANDS)
            for index, (name, scale) in enumerate(zip(BANDS, SCALES), 1):
                dst.set_band_description(index, name)
                units = "reflectance" if index <= 8 else "degrees" if index <= 12 else "percent" if name == "CLD" else "code"
                dst.update_tags(index, units=units, decoding="physical_value = stored_value * scale",
                                scale=str(scale), nodata_before_scaling=str(NODATA))
            dst.update_tags(
                cache_signature=fingerprint, cache_version=VERSION,
                aggregation_id=definition["aggregation_id"], target_date=definition["target_date"],
                mgrs_tile=definition["mgrs_tile"], sensor=scene["sensor"],
                scene_id=str(scene["scene_id"]), acquisition_time_utc=scene["date"],
                day_offset_from_target=str(scene["day_offset"]),
                source_scene_json=json.dumps(returned, sort_keys=True),
                harmonizeValues="true", resampling="NEAREST", cloud_filter_applied="false",
                angle_source_resolution_m="5000",
            )
        if inspect_acquisition(temporary, fingerprint, bounds, crs) is None:
            raise RuntimeError("Downloaded acquisition failed TIFF read-back validation")
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def process_window(row, output, config, args):
    bounds, crs = previous.geometry(row)
    definition = window_definition(row, bounds, crs)
    inventory_path = output / "window_metadata" / f"{row['aggregation_id']}.json"
    bbox = BBox(bounds, crs=crs.to_string())
    inventory = get_inventory(inventory_path, definition, bbox, config, args.refresh_scenes)
    scenes = inventory["selected_acquisitions"]
    print(f"  API candidates: {inventory['api_candidates']}; selected acquisitions: {len(scenes)}", flush=True)
    inventory["download_status"] = "downloading"
    save_json(inventory_path, inventory)
    downloaded = reused = size_bytes = 0
    files = []
    for index, scene in enumerate(scenes, 1):
        path = acquisition_path(output, row["aggregation_id"], scene)
        # Ignore discovery ordering (__idx) and other incidental metadata when
        # refreshing inventories: the same source/grid can reuse the same TIFF.
        fingerprint = digest({"window_signature": inventory["signature"],
                              "scene_id": str(scene["scene_id"]), "sensor": scene["sensor"],
                              "acquisition_time_utc": scene["date"]})
        info = None if args.overwrite else inspect_acquisition(path, fingerprint, bounds, crs)
        if info is None:
            print(f"  [{index}/{len(scenes)}] Download {scene['date']} {scene['sensor']}", flush=True)
            data, returned = request_acquisition(config, bbox, scene)
            write_acquisition(path, data, returned, scene, definition, fingerprint, bounds, crs)
            info = inspect_acquisition(path, fingerprint, bounds, crs)
            if info is None:
                raise RuntimeError(f"Saved acquisition failed validation: {path}")
            downloaded += 1
        else:
            print(f"  [{index}/{len(scenes)}] Verified; reuse {scene['date']} {scene['sensor']}", flush=True)
            reused += 1
        size_bytes += info["file_bytes"]
        files.append({"path": path.relative_to(output).as_posix(), "scene_id": scene["scene_id"],
                      "sensor": scene["sensor"], "acquisition_time_utc": scene["date"],
                      "day_offset_from_target": scene["day_offset"], **info})
    # Completion means every scene in this inventory has a validated local TIFF.
    inventory.update(download_status="completed" if scenes else "no_acquisitions", files=files,
                     total_file_bytes=size_bytes)
    save_json(inventory_path, inventory)
    print(f"  Downloaded: {downloaded}; reused: {reused}; stored: {size_bytes / 1024**2:.2f} MiB", flush=True)
    return {
        "status": "completed" if scenes else "no_acquisitions",
        "api_candidates": inventory["api_candidates"], "acquisitions": len(scenes),
        "downloaded_this_run": downloaded, "reused_this_run": reused,
        "stored_bytes": size_bytes, "inventory_path": str(inventory_path),
    }


def main():
    args = parse_args()
    rows = previous.load_rows(args)
    output = args.output_dir.resolve()
    (output / "acquisitions").mkdir(parents=True, exist_ok=True)
    (output / "window_metadata").mkdir(exist_ok=True)
    checklist = output / "acquisition_download_checklist.csv"
    records = previous.load_checklist(checklist)
    session = None

    def config():
        nonlocal session
        if session is None:
            session = previous.snap.base.build_cdse_config()
        return session

    print(f"DOWNLOAD ONLY: {len(rows)} windows\nManifest: {args.manifest.resolve()}\nOutput: {output}", flush=True)
    print("15-band UINT16 chips at 20 m; no cloud filtering, SNAP or temporal compositing", flush=True)
    failures = 0
    downloaded = reused = size_bytes = 0
    for number, (_, row) in enumerate(rows.iterrows(), 1):
        window_id = row["aggregation_id"]
        print(f"\nWINDOW {number}/{len(rows)}: {window_id}", flush=True)
        record = {"aggregation_id": window_id, "target_date": row["Delta_Date"].date().isoformat(),
                  "mgrs_tile": row["mgrs_tile_t"], "status": "downloading", "error": ""}
        records[window_id] = record.copy()
        previous.save_checklist(checklist, records)
        try:
            summary = process_window(row, output, config, args)
            record.update(summary)
            downloaded += summary["downloaded_this_run"]
            reused += summary["reused_this_run"]
            size_bytes += summary["stored_bytes"]
        except KeyboardInterrupt:
            record.update(status="interrupted", error="User interrupted; saved acquisitions can be reused")
            records[window_id] = record
            previous.save_checklist(checklist, records)
            raise
        except Exception as exc:
            record.update(status="failed", error=f"{type(exc).__name__}: {exc}")
            failures += 1
            print(f"  FAILED: {record['error']}", flush=True)
        records[window_id] = record
        previous.save_checklist(checklist, records)
    print(f"\nCompleted-window acquisition downloads: {downloaded}; reused: {reused}; failed windows: {failures}", flush=True)
    print(f"Stored TIFF size for successful selected windows: {size_bytes / 1024**3:.3f} GiB\nChecklist: {checklist}", flush=True)
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
