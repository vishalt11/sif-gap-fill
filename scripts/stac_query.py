

from google.colab import drive
drive.mount('/content/drive')
     



import numpy as np
import rasterio
from pystac_client import Client
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
import os
import time
import json
import random
from datetime import datetime


CACHE_DIR = '/content/drive/MyDrive/Capstone Project/WASP_Cache'

YEARS = [2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024]
MONTHS = [3, 4, 5, 6]
MONTH_NAMES = {3: 'March', 4: 'April', 5: 'May', 6: 'June'}

BAVARIA_BBOX = [8.9, 47.2, 13.9, 50.6]
MAX_WORKERS = 8

os.makedirs(CACHE_DIR, exist_ok=True)

print(f"Cache directory: {CACHE_DIR}")
print(f"Years: {YEARS[0]} - {YEARS[-1]}")
print(f"Months: {[MONTH_NAMES[m] for m in MONTHS]}")
print(f"Parallel workers: {MAX_WORKERS}")

STAC_URL = "https://geoservice.dlr.de/eoc/ogc/stac/v1"

def find_wasp_tiles(year, month):
    """Find all WASP tiles for a given year and month."""
    client = Client.open(STAC_URL)
    start = f"{year}-{month:02d}-01"
    end = f"{year}-{month+1:02d}-01" if month < 12 else f"{year+1}-01-01"
    search = client.search(
        collections=["S2_L3A_WASP"],
        bbox=BAVARIA_BBOX,
        datetime=f"{start}/{end}"
    )
    return list(search.items())

print("Searching for tiles...\n")

all_tiles = {}
total_tiles = 0

for year in YEARS:
    all_tiles[year] = {}
    print(f"{year}:")
    for month in MONTHS:
        tiles = find_wasp_tiles(year, month)
        all_tiles[year][month] = tiles
        total_tiles += len(tiles)
        print(f"  {MONTH_NAMES[month]}: {len(tiles)} tiles")

print(f"\n{'='*50}")
print(f"Total tiles to download: {total_tiles}")
print(f"Total files (B4 + B8): {total_tiles * 2}")
print(f"Estimated size: ~{total_tiles * 2 * 150 / 1000:.0f} GB")
print(f"{'='*50}")


def get_tile_cache_path(year, month, tile_id, band):
    """Get the cache path for a tile band."""

    month_dir = os.path.join(CACHE_DIR, str(year), f"{month:02d}_{MONTH_NAMES[month]}")
    os.makedirs(month_dir, exist_ok=True)


    filename = f"{tile_id}_{band}.tif"
    return os.path.join(month_dir, filename)

def download_and_cache_tile(item, year, month, max_retries=3):
    """Download a tile's B4 and B8 bands with retry logic."""
    tile_id = item.id
    results = {'id': tile_id, 'year': year, 'month': month, 'status': 'success', 'skipped': []}

    for band in ['FRC_B4', 'FRC_B8']:
        band_name = band.split('_')[1]
        cache_path = get_tile_cache_path(year, month, tile_id, band_name)


        if os.path.exists(cache_path):
            results['skipped'].append(band_name)
            continue

        for attempt in range(max_retries):
            try:
                time.sleep(random.uniform(0.5, 1.5))

                url = item.assets[band].href
                with rasterio.open(url) as src:
                    data = src.read(1)
                    profile = src.profile.copy()

                profile.update(
                    driver='GTiff',
                    compress='lzw',
                    tiled=True,
                    blockxsize=512,
                    blockysize=512
                )

                with rasterio.open(cache_path, 'w', **profile) as dst:
                    dst.write(data, 1)

                break

            except Exception as e:
                if attempt < max_retries - 1:
                    wait_time = (2 ** attempt) + random.uniform(1, 3)
                    print(f"    ⚠️ Retry {attempt+1}/{max_retries} for {tile_id} {band_name} (waiting {wait_time:.1f}s)")
                    time.sleep(wait_time)
                else:
                    results['status'] = 'error'
                    results['error'] = f"{band_name}: {str(e)}"
                    return results

    return results


def download_month(year, month, tiles):
    """Download all tiles for a month with rate limiting."""
    if not tiles:
        return {'downloaded': 0, 'skipped': 0, 'errors': 0}

    stats = {'downloaded': 0, 'skipped': 0, 'errors': 0}

    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {}

        for i, tile in enumerate(tiles):
            if i > 0 and i % MAX_WORKERS == 0:
                time.sleep(1)

            future = executor.submit(download_and_cache_tile, tile, year, month)
            futures[future] = tile

        for future in tqdm(as_completed(futures), total=len(tiles),
                          desc=f"{year}-{MONTH_NAMES[month]}"):
            result = future.result()

            if result['status'] == 'error':
                stats['errors'] += 1
                print(f"\n Error: {result['id']} - {result.get('error', 'Unknown')}")
            elif len(result['skipped']) == 2:
                stats['skipped'] += 1
            else:
                stats['downloaded'] += 1

    return stats

print("Download functions ready!")

print("="*70)
print("STARTING TILE CACHE DOWNLOAD")
print(f"Started at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("="*70)

total_stats = {'downloaded': 0, 'skipped': 0, 'errors': 0}
start_time = time.time()

for year in YEARS:
    print(f"\n{'#'*70}")
    print(f"YEAR: {year}")
    print(f"{'#'*70}")

    for month in MONTHS:
        tiles = all_tiles[year][month]

        if not tiles:
            print(f"\n  {MONTH_NAMES[month]}: No tiles found")
            continue

        print(f"\n  {MONTH_NAMES[month]}: {len(tiles)} tiles")

        stats = download_month(year, month, tiles)

        total_stats['downloaded'] += stats['downloaded']
        total_stats['skipped'] += stats['skipped']
        total_stats['errors'] += stats['errors']

        print(f"    ✓ Downloaded: {stats['downloaded']} | Skipped: {stats['skipped']} | Errors: {stats['errors']}")

    progress_log = {
        'last_completed_year': year,
        'timestamp': datetime.now().isoformat(),
        'stats': total_stats
    }
    with open(os.path.join(CACHE_DIR, 'download_progress.json'), 'w') as f:
        json.dump(progress_log, f, indent=2)
    print(f"\n  💾 Progress saved for {year}")

elapsed = time.time() - start_time

print(f"\n{'='*70}")
print("DOWNLOAD COMPLETE")
print(f"{'='*70}")
print(f"\nTotal downloaded: {total_stats['downloaded']} tiles")
print(f"Total skipped (cached): {total_stats['skipped']} tiles")
print(f"Total errors: {total_stats['errors']} tiles")
print(f"\nTime elapsed: {elapsed/3600:.1f} hours")
print(f"\nCache location: {CACHE_DIR}")

print("Verifying cache...\n")

cache_stats = {'total_files': 0, 'total_size_gb': 0, 'by_year': {}}

for year in YEARS:
    year_dir = os.path.join(CACHE_DIR, str(year))
    if not os.path.exists(year_dir):
        print(f"{year}: Not found")
        continue

    year_files = 0
    year_size = 0

    for month in MONTHS:
        month_dir = os.path.join(year_dir, f"{month:02d}_{MONTH_NAMES[month]}")
        if os.path.exists(month_dir):
            files = [f for f in os.listdir(month_dir) if f.endswith('.tif')]
            size = sum(os.path.getsize(os.path.join(month_dir, f)) for f in files)
            year_files += len(files)
            year_size += size

    cache_stats['by_year'][year] = {'files': year_files, 'size_gb': year_size / 1e9}
    cache_stats['total_files'] += year_files
    cache_stats['total_size_gb'] += year_size / 1e9

    print(f"{year}: {year_files} files, {year_size/1e9:.1f} GB")

print(f"\n{'='*50}")
print(f"Total files: {cache_stats['total_files']}")
print(f"Total size: {cache_stats['total_size_gb']:.1f} GB")
print(f"{'='*50}")

expected_files = total_tiles * 2
if cache_stats['total_files'] == expected_files:
    print("\n✅ Cache complete!")
else:
    print(f"\n⚠️ Expected {expected_files} files, found {cache_stats['total_files']}")
    print("   Missing files will be downloaded on next run.")

print("Sample cached files:\n")

for year in YEARS[:2]:
    year_dir = os.path.join(CACHE_DIR, str(year))
    if os.path.exists(year_dir):
        print(f"{year}/")
        for month in MONTHS[:1]:
            month_dir = os.path.join(year_dir, f"{month:02d}_{MONTH_NAMES[month]}")
            if os.path.exists(month_dir):
                files = sorted(os.listdir(month_dir))[:4]
                print(f"  {month:02d}_{MONTH_NAMES[month]}/")
                for f in files:
                    size = os.path.getsize(os.path.join(month_dir, f)) / 1e6
                    print(f"    {f} ({size:.1f} MB)")
                if len(os.listdir(month_dir)) > 4:
                    print(f"    ... and {len(os.listdir(month_dir)) - 4} more files")