import os
import random
import time
from datetime import datetime

import rasterio
from pystac_client import Client


STAC_URL = "https://geoservice.dlr.de/eoc/ogc/stac/v1"
COLLECTION = "S2_L3A_WASP"

# Test case: one exact monthly composite, one tile, one band.
YEAR = 2020
MONTH = 9
TILE = "32UNC"
BAND = "FRC_B8"

# Broad enough to intersect Niedersachsen / the target MGRS tiles.
NIEDERSACHSEN_BBOX = [6.3, 51.0, 11.8, 54.1]

OUT_DIR = os.path.join("data", "sentinel2_wasp_python_test")
MAX_RETRIES = 3


def target_date_token(year, month):
    return f"{year}{month:02d}15"


def month_window(year, month):
    start = f"{year}-{month:02d}-01"
    if month == 12:
        end = f"{year + 1}-01-01"
    else:
        end = f"{year}-{month + 1:02d}-01"
    return f"{start}/{end}"


def find_exact_wasp_item(year, month, tile):
    client = Client.open(STAC_URL)

    search = client.search(
        collections=[COLLECTION],
        bbox=NIEDERSACHSEN_BBOX,
        datetime=month_window(year, month),
        max_items=100,
    )

    date_token = target_date_token(year, month)
    tile_token = f"_T{tile}_"

    candidates = [
        item
        for item in search.items()
        if date_token in item.id and tile_token in item.id
    ]

    if not candidates:
        raise RuntimeError(
            f"No exact STAC item found for {year}-{month:02d}, tile {tile}. "
            f"Expected item id to contain {date_token} and {tile_token}."
        )

    if len(candidates) > 1:
        print("Multiple candidates found; using the first:")
        for item in candidates:
            print(f"  {item.id}")

    return candidates[0]


def output_path(item, band):
    os.makedirs(OUT_DIR, exist_ok=True)
    band_short = band.replace("FRC_", "")
    return os.path.join(OUT_DIR, f"{item.id}_{band_short}.tif")


def download_band_with_rasterio(item, band, max_retries=MAX_RETRIES):
    if band not in item.assets:
        raise RuntimeError(f"Asset {band} not found in item {item.id}")

    url = item.assets[band].href
    out_path = output_path(item, band)

    print(f"Item: {item.id}")
    print(f"Asset: {band}")
    print(f"URL: {url}")
    print(f"Output: {out_path}")

    if os.path.exists(out_path) and os.path.getsize(out_path) > 0:
        print("Output already exists; skipping.")
        return out_path

    for attempt in range(max_retries):
        try:
            time.sleep(random.uniform(0.5, 1.5))

            # This mirrors the old Colab method: GDAL/rasterio opens the
            # remote STAC asset URL, then we write a local compressed GeoTIFF.
            with rasterio.open(url) as src:
                data = src.read(1)
                profile = src.profile.copy()

            profile.update(
                driver="GTiff",
                compress="lzw",
                tiled=True,
                blockxsize=512,
                blockysize=512,
            )

            with rasterio.open(out_path, "w", **profile) as dst:
                dst.write(data, 1)

            print(f"Downloaded successfully at {datetime.now().isoformat(timespec='seconds')}")
            return out_path

        except Exception as exc:
            if attempt < max_retries - 1:
                wait_time = (2 ** attempt) + random.uniform(1, 3)
                print(
                    f"Retry {attempt + 1}/{max_retries} failed for {item.id} {band}: {exc}. "
                    f"Waiting {wait_time:.1f}s"
                )
                time.sleep(wait_time)
            else:
                raise


if __name__ == "__main__":
    item = find_exact_wasp_item(YEAR, MONTH, TILE)
    download_band_with_rasterio(item, BAND)
