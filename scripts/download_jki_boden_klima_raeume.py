"""Download Germany's Boden-Klima-Raeume (BKR) polygons from the JKI WFS.

The Julius Kuehn Institute provides two services for this dataset:

* WFS: vector polygons and attributes for analysis and spatial joins.
* WMS: rendered map images intended mainly for visualisation.

This script uses WFS because the project needs the underlying BKR geometries.

Default output:
    data/boden_klima_raeume/jki_boden_klima_raeume.geojson

Example:
    python download_jki_boden_klima_raeume.py

To inspect the layers advertised by the service:
    python download_jki_boden_klima_raeume.py --list-layers
"""

from __future__ import annotations

import argparse
import json
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
from xml.etree import ElementTree

import requests


WFS_URL = "https://geoservices.julius-kuehn.de/geoserver/bkr/wfs"
WFS_VERSION = "1.1.0"
DEFAULT_OUTPUT = Path(
    "data/boden_klima_raeume/jki_boden_klima_raeume.geojson"
)
REQUEST_TIMEOUT_SECONDS = 120
USER_AGENT = "OCO2-SIF-thesis/1.0 (JKI BKR research-data download)"


def local_name(tag: str) -> str:
    """Return an XML tag without its namespace."""
    return tag.rsplit("}", 1)[-1]


def child_text(element: ElementTree.Element, wanted_name: str) -> str | None:
    """Read the first direct child whose namespace-free tag matches."""
    for child in element:
        if local_name(child.tag) == wanted_name and child.text:
            return child.text.strip()
    return None


def request_session() -> requests.Session:
    session = requests.Session()
    session.headers.update({"User-Agent": USER_AGENT})
    return session


def get_feature_types(session: requests.Session) -> list[dict[str, str]]:
    """Read the WFS capabilities document and return advertised layers."""
    response = session.get(
        WFS_URL,
        params={
            "service": "WFS",
            "version": WFS_VERSION,
            "request": "GetCapabilities",
        },
        timeout=REQUEST_TIMEOUT_SECONDS,
    )
    response.raise_for_status()

    try:
        root = ElementTree.fromstring(response.content)
    except ElementTree.ParseError as exc:
        raise RuntimeError(
            "The JKI service did not return a valid WFS capabilities document."
        ) from exc

    layers: list[dict[str, str]] = []
    for element in root.iter():
        if local_name(element.tag) != "FeatureType":
            continue

        name = child_text(element, "Name")
        if not name:
            continue

        layers.append(
            {
                "name": name,
                "title": child_text(element, "Title") or "",
                "crs": (
                    child_text(element, "DefaultSRS")
                    or child_text(element, "DefaultCRS")
                    or ""
                ),
            }
        )

    if not layers:
        raise RuntimeError("The JKI WFS did not advertise any feature layers.")
    return layers


def choose_bkr_layer(
    layers: list[dict[str, str]], requested_layer: str | None
) -> dict[str, str]:
    """Select an explicit layer or infer the most likely BKR layer."""
    if requested_layer:
        for layer in layers:
            if layer["name"] == requested_layer:
                return layer
        available = ", ".join(layer["name"] for layer in layers)
        raise RuntimeError(
            f"Layer {requested_layer!r} was not advertised by the WFS. "
            f"Available layers: {available}"
        )

    if len(layers) == 1:
        return layers[0]

    candidates = [
        layer
        for layer in layers
        if any(
            keyword in f"{layer['name']} {layer['title']}".casefold()
            for keyword in ("bkr", "boden", "klima")
        )
    ]
    if len(candidates) == 1:
        return candidates[0]

    available = "\n".join(
        f"  {layer['name']} | {layer['title']} | {layer['crs']}"
        for layer in layers
    )
    raise RuntimeError(
        "The BKR layer could not be selected unambiguously. Run with "
        "--list-layers, then pass the required name with --layer.\n"
        f"Advertised layers:\n{available}"
    )


def validate_geojson(path: Path) -> tuple[int, dict[str, Any]]:
    """Confirm that the downloaded file is a non-empty FeatureCollection."""
    try:
        with path.open("r", encoding="utf-8") as file:
            data = json.load(file)
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise RuntimeError(
            "The WFS response was not valid GeoJSON. The server may have "
            "returned an XML error message instead."
        ) from exc

    if data.get("type") != "FeatureCollection":
        raise RuntimeError("The downloaded GeoJSON is not a FeatureCollection.")

    features = data.get("features")
    if not isinstance(features, list) or not features:
        raise RuntimeError("The downloaded BKR FeatureCollection is empty.")

    return len(features), data


def download_bkr(
    session: requests.Session,
    layer: dict[str, str],
    output_path: Path,
    overwrite: bool,
) -> None:
    """Download the selected BKR layer as WGS84 GeoJSON."""
    if output_path.exists() and not overwrite:
        raise FileExistsError(
            f"Output already exists: {output_path}. Use --overwrite to replace it."
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = output_path.with_suffix(output_path.suffix + ".part")

    response = session.get(
        WFS_URL,
        params={
            "service": "WFS",
            "version": WFS_VERSION,
            "request": "GetFeature",
            "typeName": layer["name"],
            "outputFormat": "application/json",
            "srsName": "EPSG:4326",
            "maxFeatures": 10000,
        },
        timeout=REQUEST_TIMEOUT_SECONDS,
        stream=True,
    )
    response.raise_for_status()

    try:
        with temporary_path.open("wb") as file:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    file.write(chunk)

        feature_count, geojson = validate_geojson(temporary_path)
        os.replace(temporary_path, output_path)
    except Exception:
        temporary_path.unlink(missing_ok=True)
        raise

    metadata_path = output_path.with_suffix(".metadata.json")
    metadata = {
        "dataset": "Boden-Klima-Raeume (BKR)",
        "publisher": "Julius Kuehn Institute (JKI)",
        "service": WFS_URL,
        "wfs_version": WFS_VERSION,
        "layer": layer,
        "requested_crs": "EPSG:4326",
        "geojson_crs": geojson.get("crs"),
        "feature_count": feature_count,
        "downloaded_utc": datetime.now(timezone.utc).isoformat(),
        "request_url": response.url,
    }
    with metadata_path.open("w", encoding="utf-8") as file:
        json.dump(metadata, file, ensure_ascii=False, indent=2)
        file.write("\n")

    print(f"Layer: {layer['name']} ({layer['title']})")
    print(f"Features: {feature_count}")
    print(f"GeoJSON: {output_path}")
    print(f"Metadata: {metadata_path}")


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Download Germany's JKI Boden-Klima-Raeume polygons via WFS."
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"GeoJSON output path (default: {DEFAULT_OUTPUT})",
    )
    parser.add_argument(
        "--layer",
        help="Exact WFS layer name. Normally unnecessary because it is discovered.",
    )
    parser.add_argument(
        "--list-layers",
        action="store_true",
        help="Print advertised WFS layers without downloading data.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace an existing output and metadata file.",
    )
    return parser.parse_args()


def main() -> None:
    arguments = parse_arguments()

    with request_session() as session:
        layers = get_feature_types(session)

        if arguments.list_layers:
            print("Layers advertised by the JKI BKR WFS:")
            for layer in layers:
                print(
                    f"  {layer['name']} | title={layer['title']} | "
                    f"default_crs={layer['crs']}"
                )
            return

        selected_layer = choose_bkr_layer(layers, arguments.layer)
        download_bkr(
            session=session,
            layer=selected_layer,
            output_path=arguments.output,
            overwrite=arguments.overwrite,
        )


if __name__ == "__main__":
    main()
