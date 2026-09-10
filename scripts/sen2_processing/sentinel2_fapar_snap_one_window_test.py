"""Download individual L2A acquisitions, run SNAP FAPAR, then select nearest clear.

Exactly ONE manifest window is processed (first row unless --aggregation-id is
given). Search bounds are target - 8 through target + 8 UTC calendar days.
Uses the same CDSE credentials as sentinel2_l2a_nearest_clear_window_download.py.

Requires ESA SNAP with Optical Toolbox and its gpt executable, in addition to
the existing Python environment. Set SNAP_GPT or pass --gpt if not on PATH.
This uses the actual installed BiophysicalOp, not a Python approximation of it.
An ENVI intermediate preserves the band names required by the SNAP operator.

Example, from the scripts directory:
    python sen2_processing/sentinel2_fapar_snap_one_window_test.py
    python sen2_processing/sentinel2_fapar_snap_one_window_test.py --gpt "C:/Program Files/esa-snap/bin/gpt.exe"

See sentinel2_fapar_snap_one_window_test.md for outputs and scientific caveats.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import time
from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd
import rasterio
from rasterio.transform import from_bounds
from sentinelhub import BBox, MimeType, MosaickingOrder, SentinelHubRequest

import sentinel2_l2a_nearest_clear_window_download as base


MANIFEST_PATH = (
    base.PROJECT_ROOT / "data" / "density_aggregation"
    / "sentinel2_spatial_aggregation_density_4000m_landcover_combined_mask60_min4_s2valid95_par_available"
    / "density_cluster_4000m_aggregate_manifest.csv"
)
OUTPUT_DIR = base.PROJECT_ROOT / "data" / "sentinel2_fapar_snap_one_window_test"
DAYS_EITHER_SIDE = 8
RESOLUTION_M = 20
# All L2A spectral bands available from this API. B10 is not an L2A band.
REFLECTANCE_BANDS = (
    "B01", "B02", "B03", "B04", "B05", "B06", "B07", "B08",
    "B8A", "B09", "B11", "B12",
)
ANGLE_BANDS = (
    "sunAzimuthAngles", "sunZenithAngles", "viewAzimuthMean", "viewZenithMean",
)
QA_BANDS = ("SCL", "CLD", "dataMask")
ALL_BANDS = REFLECTANCE_BANDS + ANGLE_BANDS + QA_BANDS
SNAP_INPUTS = {
    "B3": "B03", "B4": "B04", "B5": "B05", "B6": "B06",
    "B7": "B07", "B8A": "B8A", "B11": "B11", "B12": "B12",
    "view_zenith_mean": "viewZenithMean", "sun_zenith": "sunZenithAngles",
    "sun_azimuth": "sunAzimuthAngles", "view_azimuth_mean": "viewAzimuthMean",
}
# Conservative land FAPAR selection. Water (SCL 6) is not a canopy retrieval.
GOOD_SCL_CLASSES = (4, 5)
CLOUD_PROBABILITY_MAX_PERCENT = 40
# Strict smoke-test policy: reject all SNAP warning/error flags, including
# tolerance clipping flags. Raw outputs and flags are still retained.
ACCEPTED_FAPAR_FLAGS = (0,)
NODATA = -9999.0
WORKFLOW_VERSION = "1"


def save_json(path: Path, value) -> None:
    temporary = path.with_suffix(path.suffix + ".partial")
    temporary.write_text(json.dumps(value, indent=2, default=str), encoding="utf-8")
    os.replace(temporary, path)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--aggregation-id", help="Choose exactly one aggregation_id.")
    parser.add_argument("--manifest", type=Path, default=MANIFEST_PATH)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    parser.add_argument("--gpt", default=os.environ.get("SNAP_GPT"))
    parser.add_argument("--download-only", action="store_true",
                        help="Save all acquisition inputs; defer SNAP and compositing.")
    return parser.parse_args()


def load_window(path: Path, aggregation_id: str | None):
    columns = ["aggregation_id", "Delta_Date", "mgrs_tile_t", "cell_xmin",
               "cell_ymin", "cell_xmax", "cell_ymax", "window_crs"]
    frame = pd.read_csv(path, usecols=columns)
    if aggregation_id:
        frame = frame.loc[frame.aggregation_id == aggregation_id]
        if len(frame) != 1:
            raise ValueError(f"Expected one row for {aggregation_id}; found {len(frame)}")
    if frame.empty:
        raise ValueError(f"Empty manifest: {path}")
    row = frame.iloc[0].copy()
    row["Delta_Date"] = pd.to_datetime(row["Delta_Date"], utc=True).normalize()
    bounds = tuple(float(row[name]) for name in columns[3:7])
    if not np.isfinite(bounds).all() or not np.allclose(
        [bounds[2] - bounds[0], bounds[3] - bounds[1]], [4000, 4000]
    ):
        raise ValueError(f"Expected a finite 4 km window; got {bounds}")
    crs = rasterio.crs.CRS.from_user_input(row["window_crs"])
    expected = base.epsg_from_mgrs_tile(row["mgrs_tile_t"])
    if crs.to_epsg() != expected:
        raise ValueError(f"Manifest CRS {crs} conflicts with MGRS EPSG:{expected}")
    return row, bounds, crs


def find_gpt(explicit: str | None) -> str:
    candidates = [explicit] if explicit else [
        shutil.which("gpt"),
        "C:/Program Files/esa-snap/bin/gpt.exe",
        "C:/Program Files/snap/bin/gpt.exe",
    ]
    for candidate in candidates:
        if candidate:
            resolved = shutil.which(candidate) or candidate
            if Path(resolved).is_file():
                return str(Path(resolved).resolve())
    raise FileNotFoundError(
        "SNAP gpt was not found. Pass --gpt with your SNAP bin/gpt.exe path, "
        "set SNAP_GPT, or use --download-only to defer the SNAP step."
    )


def run_command(command: list[str], log_path: Path, timeout: int) -> str:
    # No shell; keep Windows helpers hidden. Called only when the user runs this script.
    with log_path.open("w", encoding="utf-8") as log:
        result = subprocess.run(
            command, stdout=log, stderr=subprocess.STDOUT, check=False,
            timeout=timeout, creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0),
        )
    output = log_path.read_text(encoding="utf-8", errors="replace")
    if result.returncode:
        raise RuntimeError(f"SNAP exited with {result.returncode}. See {log_path}\n{output[-4000:]}")
    return output


def evalscript(scene_id: str | None = None) -> str:
    discovery = scene_id is None
    bands = ("dataMask",) if discovery else ALL_BANDS
    units = ["DN"] if discovery else (
        ["REFLECTANCE"] * len(REFLECTANCE_BANDS)
        + ["DEGREES"] * len(ANGLE_BANDS) + ["DN", "PERCENT", "DN"]
    )
    selection = "" if discovery else f"""
function preProcessScenes(collections) {{
  collections.scenes.tiles = collections.scenes.tiles.filter(
    t => String(t.shId) === {json.dumps(scene_id)});
  return collections;
}}
"""
    pixel = "return [0];" if discovery else f"""
if (samples.length !== 1 || samples[0].dataMask !== 1)
  return new Array({len(ALL_BANDS)}).fill({NODATA});
const s = samples[0];
return {json.dumps(list(ALL_BANDS))}.map(b => s[b]);
"""
    return f"""//VERSION=3
function setup() {{
  return {{input: [{{bands: {json.dumps(list(bands))}, units: {json.dumps(units)}}}],
    output: {{id: "default", bands: {len(bands)}, sampleType: "FLOAT32", noDataValue: {NODATA}}},
    mosaicking: "TILE"}};
}}
{selection}
function evaluatePixel(samples) {{ {pixel} }}
function updateOutputMetadata(scenes, inputMetadata, outputMetadata) {{
  outputMetadata.userData = {{tiles: scenes.tiles, normalizationFactor: inputMetadata.normalizationFactor}};
}}
"""


def process_request(config, bbox, size, start, end, script):
    request = SentinelHubRequest(
        evalscript=script,
        input_data=[SentinelHubRequest.input_data(
            data_collection=base.CDSE_SENTINEL2_L2A,
            time_interval=(start.isoformat(), end.isoformat()), maxcc=1.0,
            mosaicking_order=MosaickingOrder.MOST_RECENT,
            other_args={"processing": {
                "harmonizeValues": True, "upsampling": "NEAREST", "downsampling": "NEAREST",
            }},
        )],
        responses=[SentinelHubRequest.output_response("default", MimeType.TIFF),
                   SentinelHubRequest.output_response("userdata", MimeType.JSON)],
        bbox=bbox, size=size, config=config,
    )
    if not request.download_list[0].url.startswith(base.CDSE_BASE_URL + "/"):
        raise RuntimeError("Unexpected non-CDSE Process API endpoint")
    for attempt in range(1, 5):
        try:
            response = request.get_data()[0]
            break
        except Exception:
            if attempt == 4:
                raise
            print(f"  API request failed; retry {attempt}/3", flush=True)
            time.sleep(5 * 2 ** (attempt - 1))
    if not isinstance(response, dict) or not {"default.tif", "userdata.json"} <= response.keys():
        raise RuntimeError("Expected TIFF + userdata.json response from CDSE")
    metadata = response["userdata.json"]
    if not isinstance(metadata, dict) or not isinstance(metadata.get("tiles"), list):
        raise RuntimeError(f"Unexpected scene metadata: {metadata}")
    return response["default.tif"], metadata


def select_acquisitions(tiles, target, tile):
    """Keep separate acquisitions; remove older reprocessing of the same pass."""
    tile = str(tile).upper().removeprefix("T")
    groups = {}
    audit = []
    modern = re.compile(
        r"(S2[A-Z])_MSIL2A_(\d{8}T\d{6})_N(\d{4})_R(\d{3})_T(\d{2}[A-Z]{3})_(\d{8}T\d{6})"
    )
    legacy = re.compile(
        r"(S2[A-Z])_OPER_MSI_L2A_TL_.*?_(\d{8}T\d{6})_A(\d{6})_T(\d{2}[A-Z]{3})_N([\d.]+)"
    )
    for raw in tiles:
        item = dict(raw)
        text = " ".join(str(raw.get(key, "")) for key in
                        ("productId", "sentinel2ProductId", "tileOriginalId", "dataPath"))
        if not re.search(rf"_T{tile}(?:_|\b)", text):
            audit.append({**item, "selection": "other_tile"})
            continue
        stamp = pd.to_datetime(raw["date"], utc=True)
        offset = (stamp.normalize() - target).days
        if abs(offset) > DAYS_EITHER_SIDE:
            audit.append({**item, "selection": "outside_calendar_window"})
            continue
        if raw.get("shId") is None:
            raise ValueError(f"Missing shId for strict scene selection: {raw}")
        match = modern.search(text)
        if match:
            sensor, sensing, baseline, orbit, parsed_tile, production = match.groups()
            key = (sensor, sensing, orbit, parsed_tile)
            rank = (production, int(baseline), str(raw["shId"]))
        else:
            match = legacy.search(text)
            if not match:
                raise ValueError(f"Cannot identify acquisition/sensor from metadata: {raw}")
            sensor, production, orbit, parsed_tile, baseline = match.groups()
            key = (sensor, orbit, parsed_tile)
            rank = (production, int(baseline.replace(".", "")), str(raw["shId"]))
        item.update(sensor=sensor, day_offset=int(offset), scene_id=str(raw["shId"]),
                    acquisition_key=list(key), production_timestamp=production)
        groups.setdefault(key, []).append((rank, item))
    selected = []
    for entries in groups.values():
        entries.sort(key=lambda entry: entry[0], reverse=True)
        selected.append(entries[0][1])
        audit.append({**entries[0][1], "selection": "selected"})
        audit.extend({**entry[1], "selection": "older_reprocessing"} for entry in entries[1:])
    selected.sort(key=lambda item: (item["date"], item["scene_id"]))
    return selected, audit


def write_tif(path, values, names, bounds, crs, tags):
    values = np.asarray(values, dtype=np.float32)
    if values.ndim == 2:
        values = values[np.newaxis]
    count, height, width = values.shape
    if count != len(names):
        raise ValueError("Band count/name mismatch")
    temporary = path.with_suffix(".partial.tif")
    with rasterio.open(
        temporary, "w", driver="GTiff", width=width, height=height, count=count,
        dtype="float32", crs=crs, transform=from_bounds(*bounds, width, height),
        nodata=NODATA, compress="deflate", predictor=3, tiled=True,
    ) as dst:
        dst.write(np.where(np.isfinite(values), values, NODATA))
        for index, name in enumerate(names, 1):
            dst.set_band_description(index, name)
        dst.update_tags(**{key: str(value) for key, value in tags.items()})
    os.replace(temporary, path)


def band(data, name):
    return data[ALL_BANDS.index(name)]


def numeric_validity(data, input_bands=ALL_BANDS):
    get_band = lambda name: data[input_bands.index(name)]
    inputs = np.stack([get_band(name) for name in SNAP_INPUTS.values()])
    valid = (get_band("dataMask") == 1) & np.all(
        np.isfinite(inputs) & (inputs != NODATA), axis=0
    )
    for name in ("sunZenithAngles", "viewZenithMean"):
        valid &= (get_band(name) >= 0) & (get_band(name) < 90)
    for name in ("sunAzimuthAngles", "viewAzimuthMean"):
        valid &= (get_band(name) >= 0) & (get_band(name) <= 360)
    return valid


def run_snap(gpt, directory, data, bounds, crs, sensor, input_bands=ALL_BANDS):
    valid = numeric_validity(data, input_bands)
    height, width = valid.shape
    # ENVI's explicit band names avoid GeoTIFF-reader-dependent band_1 names.
    envi = directory / "snap_input.img"
    with rasterio.open(
        envi, "w", driver="ENVI", width=width, height=height,
        count=len(SNAP_INPUTS), dtype="float32", crs=crs,
        transform=from_bounds(*bounds, width, height), interleave="bsq",
    ) as dst:
        for index, (snap_name, api_name) in enumerate(SNAP_INPUTS.items(), 1):
            # Finite placeholders prevent NaNs entering SNAP's domain lookup.
            # These pixels are excluded from every usable output afterwards.
            values = data[input_bands.index(api_name)]
            dst.write(np.where(valid, values, 0).astype("float32"), index)
            dst.set_band_description(index, snap_name)
    header = envi.with_suffix(".hdr")
    if not header.is_file():
        raise RuntimeError(f"ENVI header was not created: {header}")
    header_text = header.read_text(encoding="utf-8", errors="replace")
    if not all(name in header_text for name in SNAP_INPUTS):
        raise RuntimeError("ENVI header does not preserve required SNAP band names")

    graph = ET.Element("graph", id="FAPAR_single_acquisition")
    ET.SubElement(graph, "version").text = "1.0"
    read = ET.SubElement(graph, "node", id="Read")
    ET.SubElement(read, "operator").text = "Read"
    params = ET.SubElement(read, "parameters")
    ET.SubElement(params, "file").text = str(header.resolve())
    bio = ET.SubElement(graph, "node", id="FAPAR")
    ET.SubElement(bio, "operator").text = "BiophysicalOp"
    sources = ET.SubElement(bio, "sources")
    ET.SubElement(sources, "source", refid="Read")
    params = ET.SubElement(bio, "parameters")
    for key, value in {
        "sensor": sensor, "computeLAI": "false", "computeFapar": "true",
        "computeFcover": "false", "computeCab": "false", "computeCw": "false",
    }.items():
        ET.SubElement(params, key).text = value
    write = ET.SubElement(graph, "node", id="Write")
    ET.SubElement(write, "operator").text = "Write"
    ET.SubElement(ET.SubElement(write, "sources"), "sourceProduct", refid="FAPAR")
    params = ET.SubElement(write, "parameters")
    ET.SubElement(params, "file").text = str((directory / "snap_fapar.dim").resolve())
    ET.SubElement(params, "formatName").text = "BEAM-DIMAP"
    graph_path = directory / "snap_graph.xml"
    ET.ElementTree(graph).write(graph_path, encoding="utf-8", xml_declaration=True)
    run_command([gpt, str(graph_path), "-c", "512M", "-q", "2"],
                directory / "snap_gpt.log", timeout=1800)
    arrays = []
    for name in ("fapar", "fapar_flags"):
        path = directory / "snap_fapar.data" / f"{name}.img"
        with rasterio.open(path) as src:
            if (src.height, src.width, src.count) != (height, width, 1):
                raise RuntimeError(f"Unexpected SNAP output dimensions: {path}")
            arrays.append(src.read(1))
    return arrays[0], arrays[1], valid


def download_acquisition(config, bbox, bounds, crs, scene, directory, tags):
    path = directory / "l2a_inputs_20m.tif"
    transform = from_bounds(*bounds, 200, 200)
    if path.is_file():
        with rasterio.open(path) as src:
            if (src.descriptions == ALL_BANDS and src.shape == (200, 200)
                    and src.crs == crs and src.transform.almost_equals(transform)
                    and src.tags().get("scene_id") == scene["scene_id"]
                    and src.tags().get("workflow_version") == WORKFLOW_VERSION):
                print("  Reusing verified acquisition input", flush=True)
                return src.read()
    day = pd.to_datetime(scene["date"], utc=True).normalize()
    script = evalscript(scene["scene_id"])
    (directory / "download_evalscript.js").write_text(script, encoding="utf-8")
    raw, metadata = process_request(config, bbox, (200, 200), day,
                                    day + pd.Timedelta(days=1), script)
    save_json(directory / "download_metadata.json", metadata)
    returned = metadata["tiles"]
    if len(returned) != 1 or str(returned[0].get("shId")) != scene["scene_id"]:
        raise RuntimeError("API did not return exactly the requested acquisition; see download_metadata.json")
    if pd.to_datetime(returned[0]["date"], utc=True).normalize() != day:
        raise RuntimeError("Returned acquisition is on the wrong UTC date")
    if raw.shape != (200, 200, len(ALL_BANDS)):
        raise RuntimeError(f"Unexpected acquisition array shape: {raw.shape}")
    data = np.moveaxis(raw.astype("float32"), -1, 0)
    write_tif(path, data, ALL_BANDS, bounds, crs, tags)
    return data


def main():
    args = parse_args()
    row, bounds, crs = load_window(args.manifest, args.aggregation_id)
    target = row["Delta_Date"]
    safe_id = re.sub(r"[^A-Za-z0-9_.-]", "_", str(row["aggregation_id"]))
    output = args.output_dir.resolve() / safe_id
    output.mkdir(parents=True, exist_ok=True)
    gpt = None if args.download_only else find_gpt(args.gpt)
    if gpt:
        help_text = run_command([gpt, "BiophysicalOp", "-h"], output / "snap_operator_help.txt", 120)
        for required in ("computeFapar", "sensor"):
            if required not in help_text:
                raise RuntimeError(f"SNAP BiophysicalOp help lacks {required}; see snap_operator_help.txt")
    config = base.build_cdse_config()
    start = target - pd.Timedelta(days=DAYS_EITHER_SIDE)
    end = target + pd.Timedelta(days=DAYS_EITHER_SIDE + 1)
    bbox = BBox(bounds, crs=crs.to_string())
    common_tags = dict(
        workflow_version=WORKFLOW_VERSION, aggregation_id=row["aggregation_id"],
        target_date=target.date().isoformat(), reflectance_units="unitless BOA reflectance",
        angle_units="degrees", angle_source_resolution_m=5000,
        upsampling="NEAREST", downsampling="NEAREST", harmonizeValues=True,
    )
    save_json(output / "run_configuration.json", {
        **common_tags, "manifest": str(args.manifest.resolve()), "bounds": bounds,
        "crs": crs.to_string(), "resolution_m": RESOLUTION_M, "gpt": gpt,
        "start_inclusive": start, "end_exclusive": end,
        "bands": ALL_BANDS, "snap_band_mapping": SNAP_INPUTS,
        "good_scl": GOOD_SCL_CLASSES, "max_cld": CLOUD_PROBABILITY_MAX_PERCENT,
        "accepted_fapar_flags": ACCEPTED_FAPAR_FLAGS,
    })
    print(f"One-window test: {row['aggregation_id']}", flush=True)
    print(f"Target: {target.date()}; inclusive search: {start.date()} to {(end - pd.Timedelta(days=1)).date()}")
    print(f"Output: {output}\nGrid: 200 x 200 at 20 m, {crs}", flush=True)
    discovery_script = evalscript()
    (output / "discovery_evalscript.js").write_text(discovery_script, encoding="utf-8")
    # S2L2A permits at most 1500 m/pixel. A 4 x 4 metadata preview over
    # this validated 4 km window is 1000 m/pixel; 1 x 1 would be rejected.
    # Acquisition downloads below still use the full 200 x 200, 20 m grid.
    _, metadata = process_request(config, bbox, (4, 4), start, end, discovery_script)
    save_json(output / "all_scene_metadata.json", metadata)
    scenes, audit = select_acquisitions(metadata["tiles"], target, row["mgrs_tile_t"])
    save_json(output / "scene_selection_audit.json", audit)
    print(f"API candidates: {len(metadata['tiles'])}; unique acquisitions on assigned tile: {len(scenes)}", flush=True)
    if not scenes:
        raise RuntimeError("No acquisitions in the requested tile/calendar window; see scene_selection_audit.json")
    unsupported = [scene for scene in scenes if scene["sensor"] not in ("S2A", "S2B")]
    if unsupported and not args.download_only:
        raise RuntimeError("This workflow supports SNAP's S2A/S2B models only; no S2C substitution is made.")
    records = []
    products = []
    # The presence of this file, not an old TIFF alone, marks a complete run.
    completion = output / "completed.json"
    if completion.exists():
        completion.unlink()
    for index, scene in enumerate(scenes, 1):
        stamp = pd.to_datetime(scene["date"], utc=True)
        directory = output / "acquisitions" / f"{stamp.strftime('%Y%m%dT%H%M%S')}_{scene['sensor']}_{scene['scene_id']}"
        directory.mkdir(parents=True, exist_ok=True)
        tags = {**common_tags, "scene_id": scene["scene_id"], "sensor": scene["sensor"],
                "acquisition_date_time": scene["date"], "day_offset": scene["day_offset"]}
        save_json(directory / "source_scene.json", scene)
        record = {**scene, "source_index": index, "directory": str(directory), "status": "started"}
        print(f"[{index}/{len(scenes)}] {scene['date']} {scene['sensor']} offset={scene['day_offset']:+d}", flush=True)
        try:
            data = download_acquisition(config, bbox, bounds, crs, scene, directory, tags)
            record["data_fraction"] = float((band(data, "dataMask") == 1).mean())
            if args.download_only:
                record["status"] = "downloaded"
            else:
                fapar, flags, numeric_valid = run_snap(gpt, directory, data, bounds, crs, scene["sensor"])
                clear = (numeric_valid & np.isin(band(data, "SCL"), GOOD_SCL_CLASSES)
                         & np.isfinite(band(data, "CLD")) & (band(data, "CLD") >= 0)
                         & (band(data, "CLD") <= CLOUD_PROBABILITY_MAX_PERCENT))
                usable = (clear & np.isfinite(fapar) & (fapar >= 0) & (fapar <= 1)
                          & np.isin(flags, ACCEPTED_FAPAR_FLAGS))
                write_tif(directory / "fapar_raw_20m.tif", np.where(numeric_valid, fapar, NODATA),
                          ["FAPAR"], bounds, crs, tags)
                write_tif(directory / "fapar_quality_20m.tif", np.stack([
                    np.where(numeric_valid, flags, NODATA), band(data, "SCL"), band(data, "CLD"),
                    clear.astype("float32"), usable.astype("float32"),
                ]), ["fapar_flags", "SCL", "CLD_percent", "clear_land_pixel", "usable_fapar"], bounds, crs, tags)
                write_tif(directory / "fapar_clear_20m.tif", np.where(usable, fapar, NODATA),
                          ["FAPAR"], bounds, crs, tags)
                values, counts = np.unique(flags[numeric_valid], return_counts=True)
                record.update(status="completed", clear_fraction=float(clear.mean()),
                              usable_fapar_fraction=float(usable.mean()),
                              flag_counts={str(int(k)): int(v) for k, v in zip(values, counts)})
                print(f"  Clear land: {clear.mean():.2%}; usable FAPAR: {usable.mean():.2%}; flags: {record['flag_counts']}", flush=True)
                products.append((scene, index, data, fapar, flags, usable))
        except Exception as exc:
            record.update(status="failed", error=f"{type(exc).__name__}: {exc}")
            records.append(record)
            save_json(output / "acquisition_summary.json", records)
            # Do not quietly composite only a subset of requested acquisitions.
            raise
        records.append(record)
        save_json(output / "acquisition_summary.json", records)
    if args.download_only:
        print("All inputs downloaded. Rerun without --download-only to run SNAP using cached inputs.")
        return

    result = np.full((200, 200), NODATA, dtype="float32")
    selection = np.full((5, 200, 200), NODATA, dtype="float32")
    best_distance = np.full((200, 200), np.inf)
    best_offset = np.full((200, 200), -999)
    best_cloud = np.full((200, 200), np.inf)
    for scene, index, data, fapar, flags, usable in products:
        offset = scene["day_offset"]
        distance = abs(offset)
        cloud = band(data, "CLD")
        better = usable & (
            (distance < best_distance)
            | ((distance == best_distance) & (offset > best_offset))
            | ((distance == best_distance) & (offset == best_offset) & (cloud < best_cloud))
        )
        result[better] = fapar[better]
        for destination, value in zip(selection, (index, offset, band(data, "SCL"), cloud, flags)):
            destination[better] = value[better] if isinstance(value, np.ndarray) else value
        best_distance[better], best_offset[better], best_cloud[better] = distance, offset, cloud[better]
    valid = result != NODATA
    tags = {**common_tags, "selection": "nearest usable FAPAR; future on day tie; lower CLD on same-day tie", "gap_fill": "none"}
    write_tif(output / "fapar_nearest_clear_20m.tif", result, ["FAPAR"], bounds, crs,
              {**tags, "units": "fraction 0 to 1"})
    write_tif(output / "fapar_nearest_clear_provenance_20m.tif",
              np.concatenate([selection, valid[np.newaxis].astype("float32")]),
              ["source_index", "source_day_offset_from_target", "selected_SCL", "selected_CLD_percent", "selected_fapar_flags", "valid_pixel"],
              bounds, crs, tags)
    # Carry exactly the SAME selected observation into predictor bands. This
    # prevents FAPAR and reflectance channels in a later CNN chip using different dates.
    selected_inputs = np.full((len(ALL_BANDS), 200, 200), NODATA, dtype="float32")
    for _, index, data, _, _, _ in products:
        chosen = selection[0] == index
        selected_inputs[:, chosen] = data[:, chosen]
    write_tif(output / "l2a_inputs_matching_fapar_selection_20m.tif", selected_inputs,
              ALL_BANDS, bounds, crs, tags)
    save_json(completion, {"status": "completed", "acquisitions": len(scenes),
                          "valid_fraction": float(valid.mean()), "output": str(output)})
    print(f"Final valid FAPAR: {valid.mean():.2%}. Remaining gaps are NoData.\nWrote: {output / 'fapar_nearest_clear_20m.tif'}", flush=True)


if __name__ == "__main__":
    main()
