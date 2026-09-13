"""Local sensor-aware SNAP FAPAR from packed coherent composites; no downloads.

Default smoke test: first five input TIFFs, one batch, one worker. Use --all
explicitly for the full directory. --compare-individual additionally checks
batch results against one-window graphs (up to 20 test windows).
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import threading
import time
from uuid import uuid4
import xml.etree.ElementTree as ET

import numpy as np
import rasterio


ROOT = Path(__file__).resolve().parents[1]
INPUT_DIR = ROOT / "data" / "sentinel2_fapar_composite_inputs"
OUTPUT_DIR = ROOT / "data" / "sentinel2_fapar_composite_snap"
REFLECTANCE = ("B03", "B04", "B05", "B06", "B07", "B8A", "B11", "B12")
ANGLES = ("sunAzimuthAngles", "sunZenithAngles", "viewAzimuthMean", "viewZenithMean")
BANDS = REFLECTANCE + ANGLES + ("SCL", "CLD", "dataMask", "sensor_id")
SCALES = (0.0001,) * 8 + (0.01,) * 4 + (1.0,) * 4
SNAP_NAMES = ("B3", "B4", "B5", "B6", "B7", "B8A", "B11", "B12",
              "sun_azimuth", "sun_zenith", "view_azimuth_mean", "view_zenith_mean")
ACCEPTED = (0, 1, 2, 3)
NODATA = -9999.0
FLAG_NODATA = 255
VERSION = "coherent-composite-snap-batch-v1"


def arguments():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input-dir", type=Path, default=INPUT_DIR)
    p.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    p.add_argument("--gpt", default=os.environ.get("SNAP_GPT"))
    selection = p.add_mutually_exclusive_group()
    selection.add_argument("--limit", type=int, help="First N windows; default is five.")
    selection.add_argument("--all", action="store_true", help="Explicitly process every input window.")
    p.add_argument("--aggregation-id", help="Select one WINDOW_ID, without _snap_inputs.tif.")
    p.add_argument("--batch-size", type=int, default=5)
    p.add_argument("--workers", type=int, default=1, help="Concurrent GPT processes; start with one.")
    p.add_argument("--snap-threads", type=int, default=2, help="GPT -q threads PER process.")
    p.add_argument("--cache", default="512M", help="GPT tile cache PER process, not total Java heap.")
    p.add_argument("--timeout", type=int, default=1800, help="Seconds allowed per GPT graph.")
    p.add_argument("--compare-individual", action="store_true")
    p.add_argument("--overwrite", action="store_true")
    p.add_argument("--keep-work", action="store_true", help="Keep intermediate ENVI and SNAP files.")
    args = p.parse_args()
    if any(x < 1 for x in (args.batch_size, args.workers, args.snap_threads, args.timeout)):
        p.error("Batch size, workers, threads and timeout must be positive")
    if args.limit is not None and args.limit < 1:
        p.error("--limit must be positive")
    if not re.fullmatch(r"[1-9][0-9]*[KMG]?", args.cache, re.IGNORECASE):
        p.error("--cache must be a positive size such as 512M or 1G")
    if args.compare_individual and args.all:
        p.error("Use --compare-individual on a small --limit test, not --all")
    return args


def find_gpt(explicit):
    candidates = [explicit] if explicit else [shutil.which("gpt"),
        "C:/Program Files/esa-snap/bin/gpt.exe", "C:/Program Files/snap/bin/gpt.exe"]
    for candidate in candidates:
        if candidate:
            path = Path(shutil.which(candidate) or candidate)
            if path.is_file():
                return path.resolve()
    raise FileNotFoundError("Pass --gpt with the installed SNAP bin/gpt.exe path")


def atomic_json(path, data):
    tmp = path.with_suffix(".partial.json")
    tmp.write_text(json.dumps(data, indent=2), encoding="utf-8")
    os.replace(tmp, path)


def file_hash(path):
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def run_gpt(command, log_path, timeout, stop):
    if stop.is_set():
        raise RuntimeError("Processing cancelled")
    began = time.perf_counter()
    with log_path.open("w", encoding="utf-8") as log:
        child = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT,
                                 creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
        try:
            while True:
                try:
                    code = child.wait(timeout=0.2)
                    break
                except subprocess.TimeoutExpired:
                    if stop.is_set():
                        raise RuntimeError("Processing cancelled")
                    if time.perf_counter() - began > timeout:
                        raise TimeoutError(f"GPT exceeded {timeout} seconds; see {log_path}")
        finally:
            if child.poll() is None:
                child.terminate()
                try:
                    child.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    child.kill()
                    child.wait()
    if code:
        tail = log_path.read_text(encoding="utf-8", errors="replace")[-3000:]
        raise RuntimeError(f"GPT exited {code}; see {log_path}\n{tail}")
    return time.perf_counter() - began


def load_input(path, runtime):
    window_id = path.name.removesuffix("_snap_inputs.tif")
    before = path.stat()
    checksum = file_hash(path)
    with rasterio.open(path) as src:
        if (src.descriptions != BANDS or src.dtypes != ("uint16",) * 16
                or src.nodata != 65535 or src.shape != (200, 200) or src.crs is None
                or not np.allclose(src.scales, SCALES, rtol=0, atol=1e-12)
                or src.offsets != (0.0,) * 16):
            raise ValueError("Expected the original 16-band packed UINT16 composite with its scales")
        transform = src.transform
        if not np.allclose((transform.a, transform.b, transform.d, transform.e), (20, 0, 0, -20)):
            raise ValueError("Expected a north-up 20 m input grid")
        tags = src.tags()
        if tags.get("aggregation_id") != window_id:
            raise ValueError("TIFF aggregation_id does not match its filename")
        # Rasterio read() returns stored integers: apply metadata scales ONCE.
        raw = src.read()
        missing = raw == 65535
        if not np.array_equal(missing.any(axis=0), missing.all(axis=0)):
            raise ValueError("Input contains partially missing observation vectors")
        data = raw.astype("float32") * np.asarray(src.scales, dtype="float32")[:, None, None]
        data[missing] = np.nan
        sensor = raw[15]
        present = ~missing[0]
        if not np.isin(sensor[present], (1, 2)).all():
            raise ValueError("Unsupported or missing sensor ID at non-NoData pixels")
        valid = (present & (raw[14] == 1) & np.isin(raw[12], (4, 5))
                 & (raw[13] <= 40) & np.isfinite(data[:12]).all(axis=0))
        for index in (9, 11):
            valid &= (data[index] >= 0) & (data[index] < 90)
        for index in (8, 10):
            valid &= (data[index] >= 0) & (data[index] <= 360)
        crs = src.crs
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise RuntimeError("Input changed while being read; rerun after downloads finish")
    definition = {"input_sha256": checksum, "runtime": runtime, "version": VERSION,
                  "accepted_flags": ACCEPTED, "good_scl": (4, 5), "max_cld": 40}
    fingerprint = hashlib.sha256(json.dumps(definition, sort_keys=True).encode()).hexdigest()
    return {"id": window_id, "path": path, "data": data, "valid": valid,
            "sensor": sensor, "sensors": [x for x in (1, 2) if (valid & (sensor == x)).any()],
            "transform": transform, "crs": crs, "signature": fingerprint,
            "input_sha256": checksum, "source_tags": tags}


def paths_for(output, item):
    return output / f"{item['id']}_fapar.tif", output / f"{item['id']}_snap_flag.tif"


def inspect_pair(output, item):
    paths = paths_for(output, item)
    if not all(p.is_file() for p in paths):
        return None
    try:
        with rasterio.open(paths[0]) as f, rasterio.open(paths[1]) as q:
            for src, name, dtype, nodata in ((f, "FAPAR", "float32", NODATA),
                                           (q, "SNAP_flag", "uint8", FLAG_NODATA)):
                if (src.shape != (200, 200) or src.count != 1 or src.descriptions != (name,)
                        or src.dtypes != (dtype,) or src.nodata != nodata
                        or src.crs != item["crs"] or src.transform != item["transform"]
                        or src.tags().get("signature") != item["signature"]):
                    return None
            if not f.tags().get("pair_id") or f.tags()["pair_id"] != q.tags().get("pair_id"):
                return None
            values, flags = f.read(1), q.read(1)
            good = flags != FLAG_NODATA
            if (not np.isin(flags, (*ACCEPTED, FLAG_NODATA)).all()
                    or not np.array_equal(values != NODATA, good)
                    or np.any(good & ~item["valid"])
                    or not np.isfinite(values[good]).all()
                    or np.any((values[good] < 0) | (values[good] > 1))):
                return None
            summary = json.loads(f.tags()["summary_json"])
            if not isinstance(summary, dict) or summary.get("valid_pixels") != int(good.sum()):
                return None
            return summary
    except (OSError, ValueError, KeyError, TypeError, rasterio.errors.RasterioError):
        return None


def prepare_envi(item, directory):
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / "inputs.img"
    with rasterio.open(path, "w", driver="ENVI", width=200, height=200, count=12,
                       dtype="float32", crs=item["crs"], transform=item["transform"],
                       interleave="bsq") as dst:
        for index, name in enumerate(SNAP_NAMES, 1):
            # Placeholder values outside the original valid mask never survive
            # final masking. Every model sees the same coherent input vectors.
            dst.write(np.where(item["valid"], item["data"][index - 1], 0).astype("float32"), index)
            dst.set_band_description(index, name)
    header = path.with_suffix(".hdr")
    if not header.is_file() or not all(name in header.read_text(errors="replace") for name in SNAP_NAMES):
        raise RuntimeError("ENVI input did not retain SNAP's required band names")
    item["header"] = header


def make_graph(items, directory, graph_path):
    directory.mkdir(parents=True, exist_ok=True)
    graph = ET.Element("graph", id="sensor_aware_FAPAR_batch")
    ET.SubElement(graph, "version").text = "1.0"
    destinations = {}
    for index, item in enumerate(items):
        if not item["sensors"]:
            continue
        read_id = f"Read_{index}"
        node = ET.SubElement(graph, "node", id=read_id)
        ET.SubElement(node, "operator").text = "Read"
        ET.SubElement(ET.SubElement(node, "parameters"), "file").text = str(item["header"].resolve())
        for sensor in item["sensors"]:
            bio_id = f"FAPAR_{index}_{sensor}"
            node = ET.SubElement(graph, "node", id=bio_id)
            ET.SubElement(node, "operator").text = "BiophysicalOp"
            ET.SubElement(ET.SubElement(node, "sources"), "source", refid=read_id)
            params = ET.SubElement(node, "parameters")
            for key, value in {"sensor": "S2A" if sensor == 1 else "S2B", "computeFapar": "true",
                               "computeLAI": "false", "computeFcover": "false",
                               "computeCab": "false", "computeCw": "false"}.items():
                ET.SubElement(params, key).text = value
            destination = directory / f"window_{index}_sensor_{sensor}.dim"
            destinations[item["id"], sensor] = destination
            node = ET.SubElement(graph, "node", id=f"Write_{index}_{sensor}")
            ET.SubElement(node, "operator").text = "Write"
            ET.SubElement(ET.SubElement(node, "sources"), "sourceProduct", refid=bio_id)
            params = ET.SubElement(node, "parameters")
            ET.SubElement(params, "file").text = str(destination.resolve())
            ET.SubElement(params, "formatName").text = "BEAM-DIMAP"
    ET.ElementTree(graph).write(graph_path, encoding="utf-8", xml_declaration=True)
    return destinations


def execute_graph(items, directory, graph_path, log_path, args, gpt, stop):
    destinations = make_graph(items, directory, graph_path)
    if not destinations:
        return destinations, 0.0
    seconds = run_gpt([str(gpt), str(graph_path), "-c", args.cache, "-q", str(args.snap_threads), "-e"],
                      log_path, args.timeout, stop)
    return destinations, seconds


def read_merged(item, destinations):
    values = np.full((200, 200), np.nan, dtype="float32")
    flags = np.full((200, 200), FLAG_NODATA, dtype="uint8")
    for sensor in item["sensors"]:
        mask = item["valid"] & (item["sensor"] == sensor)
        folder = destinations[item["id"], sensor].with_suffix(".data")
        arrays = []
        for name in ("fapar", "fapar_flags"):
            path = folder / f"{name}.img"
            with rasterio.open(path) as src:
                if (src.shape != (200, 200) or src.count != 1 or src.crs != item["crs"]
                        or not src.transform.almost_equals(item["transform"], precision=1e-6)):
                    raise RuntimeError(f"SNAP output grid differs from input: {path}")
                arrays.append(src.read(1))
        raw_flags = arrays[1][mask]
        if (not np.isfinite(raw_flags).all() or np.any(raw_flags != np.floor(raw_flags))
                or np.any((raw_flags < 0) | (raw_flags > 31))):
            raise RuntimeError("Unexpected SNAP flag values; inspect the retained intermediate files")
        values[mask] = arrays[0][mask]
        flags[mask] = raw_flags.astype("uint8")
    return values, flags


def finish_arrays(item, values, flags):
    valid = item["valid"]
    accepted = valid & np.isin(flags, ACCEPTED) & np.isfinite(values) & (values >= 0) & (values <= 1)
    unique, counts = np.unique(flags[valid], return_counts=True)
    summary = {"input_valid_pixels": int(valid.sum()), "valid_pixels": int(accepted.sum()),
               "valid_fraction": float(accepted.mean()), "nodata_pixels": int((~accepted).sum()),
               "raw_flag_counts": {str(int(k)): int(v) for k, v in zip(unique, counts)},
               "s2a_input_pixels": int((valid & (item["sensor"] == 1)).sum()),
               "s2b_input_pixels": int((valid & (item["sensor"] == 2)).sum())}
    for flag in ACCEPTED:
        summary[f"accepted_flag_{flag}_pixels"] = int((accepted & (flags == flag)).sum())
    summary["fapar_min"] = float(values[accepted].min()) if accepted.any() else None
    summary["fapar_max"] = float(values[accepted].max()) if accepted.any() else None
    return (np.where(accepted, values, NODATA).astype("float32"),
            np.where(accepted, flags, FLAG_NODATA).astype("uint8"), summary)


def write_pair(output, item, values, flags, summary):
    paths = paths_for(output, item)
    temporary = [p.with_suffix(".partial.tif") for p in paths]
    pair_id = uuid4().hex
    try:
        for path, data, name, dtype, nodata in (
            (temporary[0], values, "FAPAR", "float32", NODATA),
            (temporary[1], flags, "SNAP_flag", "uint8", FLAG_NODATA),
        ):
            with rasterio.open(path, "w", driver="GTiff", width=200, height=200, count=1,
                               dtype=dtype, nodata=nodata, transform=item["transform"], crs=item["crs"],
                               compress="deflate", predictor=3 if dtype == "float32" else 2, tiled=True) as dst:
                dst.write(data, 1)
                dst.set_band_description(1, name)
                dst.update_tags(1, units="fraction" if name == "FAPAR" else "bit_flags")
                dst.update_tags(signature=item["signature"], pair_id=pair_id, workflow_version=VERSION,
                                aggregation_id=item["id"], input_sha256=item["input_sha256"],
                                source_path=str(item["path"]), accepted_flags="0,1,2,3",
                                target_date=item["source_tags"].get("target_date", ""),
                                sensor_models="S2A for sensor_id=1; S2B for sensor_id=2",
                                summary_json=json.dumps(summary), gap_fill="none", reprojection="none")
        for source, destination in zip(temporary, paths):
            os.replace(source, destination)
    finally:
        for path in temporary:
            path.unlink(missing_ok=True)
    if inspect_pair(output, item) is None:
        raise RuntimeError("Final output pair failed read-back validation")


def compare(item, batched, reference):
    mask = item["valid"]
    a, fa = batched
    b, fb = reference
    equal_flags = np.array_equal(fa[mask], fb[mask])
    close = np.isclose(a[mask], b[mask], rtol=0, atol=1e-6, equal_nan=True)
    finite = mask & np.isfinite(a) & np.isfinite(b)
    maximum = float(np.max(np.abs(a[finite] - b[finite]))) if finite.any() else 0.0
    report = {"passed": bool(equal_flags and close.all()), "identical_flags": bool(equal_flags),
              "flag_different_pixels": int((fa[mask] != fb[mask]).sum()),
              "fapar_different_pixels": int((~close).sum()), "max_abs_difference": maximum,
              "absolute_tolerance": 1e-6}
    if not report["passed"]:
        raise RuntimeError(f"Batch/individual mismatch for {item['id']}: {report}")
    return report


def process_batch(number, paths, args, output, run_dir, runtime, gpt, stop):
    started = time.perf_counter()
    work_root = output / "_work"
    work_root.mkdir(exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix=f"batch_{number:05d}_", dir=work_root)).resolve()
    prefix = f"batch_{number:05d}"
    report = {"batch": number, "work_dir": str(work), "records": [], "gpt_seconds": 0.0,
              "individual_gpt_seconds": 0.0, "comparisons": {}}
    pending = []
    success = False
    try:
        for path in paths:
            if stop.is_set():
                raise RuntimeError("Processing cancelled")
            try:
                item = load_input(path, runtime)
                existing = None if args.overwrite or args.compare_individual else inspect_pair(output, item)
                if existing is not None:
                    report["records"].append({"aggregation_id": item["id"], "status": "skipped_existing", **existing})
                    continue
                prepare_envi(item, work / f"input_{len(pending)}")
                pending.append(item)
                print(f"  {item['id']}: input valid {item['valid'].mean():.2%}; "
                      f"models {[('S2A' if s == 1 else 'S2B') for s in item['sensors']]}", flush=True)
            except Exception as exc:
                report["records"].append({"aggregation_id": path.name.removesuffix("_snap_inputs.tif"),
                                          "status": "failed", "error": f"{type(exc).__name__}: {exc}"})
        destinations, seconds = execute_graph(pending, work / "batch_results", run_dir / f"{prefix}.xml",
                                              run_dir / f"{prefix}.log", args, gpt, stop)
        report["gpt_seconds"] = seconds
        merged = {item["id"]: read_merged(item, destinations) for item in pending}
        # Compare RAW results at every input-valid pixel, not only accepted
        # FAPAR pixels. Do not publish any new pair if a comparison fails.
        if args.compare_individual:
            for index, item in enumerate(pending):
                dest, seconds = execute_graph([item], work / f"individual_{index}",
                    run_dir / f"{prefix}_individual_{index}.xml", run_dir / f"{prefix}_individual_{index}.log",
                    args, gpt, stop)
                report["individual_gpt_seconds"] += seconds
                report["comparisons"][item["id"]] = compare(item, merged[item["id"]], read_merged(item, dest))
        for item in pending:
            values, flags, summary = finish_arrays(item, *merged[item["id"]])
            write_pair(output, item, values, flags, summary)
            record = {"aggregation_id": item["id"], "status": "completed" if summary["valid_pixels"] else "completed_empty",
                      "fapar_path": str(paths_for(output, item)[0]), "flag_path": str(paths_for(output, item)[1]), **summary}
            report["records"].append(record)
            print(f"  {item['id']}: accepted FAPAR {summary['valid_fraction']:.2%}; "
                  f"raw flags {summary['raw_flag_counts']}", flush=True)
        success = not any(r["status"] == "failed" for r in report["records"])
    except Exception as exc:
        done = {r["aggregation_id"] for r in report["records"]}
        for path in paths:
            window_id = path.name.removesuffix("_snap_inputs.tif")
            if window_id not in done:
                report["records"].append({"aggregation_id": window_id, "status": "failed",
                                          "error": f"{type(exc).__name__}: {exc}"})
        report["error"] = f"{type(exc).__name__}: {exc}"
    finally:
        report["elapsed_seconds"] = time.perf_counter() - started
        if success and not args.keep_work:
            # Delete only the temporary directory this batch created directly
            # under this script's _work root, never a user-supplied input path.
            if work.parent != work_root.resolve() or not work.name.startswith(f"batch_{number:05d}_"):
                raise RuntimeError("Refusing unsafe work-directory cleanup")
            try:
                shutil.rmtree(work)
                report["work_retained"] = False
            except OSError as exc:
                # Antivirus/indexing can briefly lock Windows intermediates;
                # that must not invalidate successfully written final TIFFs.
                report["work_retained"] = True
                report["cleanup_warning"] = str(exc)
        else:
            report["work_retained"] = True
        atomic_json(run_dir / f"{prefix}_report.json", report)
    return report


def update_checklist(path, records):
    columns = sorted({key for record in records.values() for key in record} - {"aggregation_id"})
    temporary = path.with_suffix(".partial.csv")
    with temporary.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=["aggregation_id", *columns])
        writer.writeheader()
        for record in records.values():
            writer.writerow({k: json.dumps(v) if isinstance(v, (dict, list)) else v for k, v in record.items()})
    os.replace(temporary, path)


def main():
    args = arguments()
    gpt = find_gpt(args.gpt)
    output = args.output_dir.resolve()
    if output == args.input_dir.resolve():
        raise ValueError("Use a different directory for outputs and downloaded inputs")
    files = sorted(args.input_dir.resolve().glob("*_snap_inputs.tif"))
    if args.aggregation_id:
        files = [p for p in files if p.name == f"{args.aggregation_id}_snap_inputs.tif"]
    if not args.all:
        files = files[:args.limit if args.limit is not None else 5]
    if not files:
        raise ValueError("No matching input composites found")
    if args.compare_individual and len(files) > 20:
        raise ValueError("Limit comparison tests to 20 windows or fewer")
    output.mkdir(parents=True, exist_ok=True)
    run_dir = output / "logs" / (time.strftime("%Y%m%d_%H%M%S") + "_" + uuid4().hex[:8])
    run_dir.mkdir(parents=True)
    stop = threading.Event()
    run_gpt([str(gpt), "BiophysicalOp", "-h"], run_dir / "snap_help.txt", 120, stop)
    help_text = (run_dir / "snap_help.txt").read_text(encoding="utf-8", errors="replace")
    if any(word not in help_text for word in ("sensor", "computeFapar", "S2A", "S2B")):
        raise RuntimeError(f"Required BiophysicalOp sensor/FAPAR support not found; see {run_dir / 'snap_help.txt'}")
    # Extract usage before startup/shutdown logs so timestamped logger output
    # cannot invalidate every existing output on a later run.
    usage = re.search(r"(?im)^\s*Usage:", help_text)
    if usage is None:
        raise RuntimeError("Could not identify stable GPT usage text; inspect snap_help.txt")
    parameter_help = help_text[usage.start():]
    end_help = re.search(r"(?m)^(?:INFO|WARNING|SEVERE|DEBUG):|^\d{4}-\d{2}-\d{2}[ T]", parameter_help)
    if end_help:
        parameter_help = parameter_help[:end_help.start()]
    stat = gpt.stat()
    runtime = {"gpt_path": str(gpt), "gpt_size": stat.st_size, "gpt_mtime_ns": stat.st_mtime_ns,
               "operator_help_sha256": hashlib.sha256(parameter_help.encode()).hexdigest()}
    atomic_json(run_dir / "run_configuration.json", {"runtime": runtime, "version": VERSION,
        "arguments": {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
        "files": [str(p) for p in files], "accepted_flags": ACCEPTED})
    checklist = output / "snap_processing_checklist.csv"
    records = {}
    if checklist.is_file():
        with checklist.open(newline="", encoding="utf-8") as stream:
            for record in csv.DictReader(stream):
                key = record["aggregation_id"]
                if key in records:
                    raise ValueError("Duplicate window in processing checklist")
                records[key] = record
    batches = [files[i:i + args.batch_size] for i in range(0, len(files), args.batch_size)]
    print(f"LOCAL SNAP: {len(files)} windows; {len(batches)} batches; {args.workers} workers\n"
          f"Output: {output}\nLogs: {run_dir}\nComparison: {args.compare_individual}", flush=True)
    reports = []
    began = time.perf_counter()
    executor = ThreadPoolExecutor(max_workers=args.workers)
    try:
        futures = [executor.submit(process_batch, n, paths, args, output, run_dir, runtime, gpt, stop)
                   for n, paths in enumerate(batches, 1)]
        for future in as_completed(futures):
            report = future.result()
            reports.append(report)
            for record in report["records"]:
                records[record["aggregation_id"]] = record
                if record["status"] == "failed":
                    print(f"FAILED {record['aggregation_id']}: {record.get('error', '')}", flush=True)
            update_checklist(checklist, records)
            print(f"Batch {report['batch']}: GPT {report['gpt_seconds']:.2f}s; "
                  f"individual GPT {report['individual_gpt_seconds']:.2f}s; "
                  f"total {report['elapsed_seconds']:.2f}s; comparisons {len(report['comparisons'])}", flush=True)
    except BaseException:
        stop.set()
        raise
    finally:
        executor.shutdown(wait=True, cancel_futures=True)
    failures = sum(r["status"] == "failed" for b in reports for r in b["records"])
    summary = {"windows": len(files), "failed_windows": failures, "elapsed_seconds": time.perf_counter() - began,
               "batch_gpt_seconds_sum": sum(b["gpt_seconds"] for b in reports),
               "individual_gpt_seconds_sum": sum(b["individual_gpt_seconds"] for b in reports),
               "comparison_windows_passed": sum(len(b["comparisons"]) for b in reports), "batches": reports}
    atomic_json(run_dir / "run_summary.json", summary)
    print(f"Finished: {failures} failed; elapsed {summary['elapsed_seconds']:.2f}s\n"
          f"Share this report: {run_dir / 'run_summary.json'}", flush=True)
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
