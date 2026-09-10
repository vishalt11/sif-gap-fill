# Coherent Sentinel-2 composites for later SNAP FAPAR

This **download-only** alternative leaves the existing per-acquisition and
SNAP-processing scripts unchanged. It requires the same Python packages and
CDSE credentials as those scripts, but does not launch SNAP or require `--gpt`.
Keep `sentinel2_fapar_snap_window_download.py`,
`sentinel2_fapar_snap_one_window_test.py` and
`sentinel2_l2a_nearest_clear_window_download.py` in this script directory: their
manifest, credential, request and acquisition-selection helpers are imported.

From `D:/UE/4_Semester/code/scripts`, first test one window:

```powershell
python sen2_processing/sentinel2_fapar_composite_download.py --limit 1
```

Then run all manifest windows; validated existing composites are skipped:

```powershell
python sen2_processing/sentinel2_fapar_composite_download.py
```

Options: `--aggregation-id WINDOW_ID`, `--limit N`, `--manifest PATH`,
`--output-dir PATH`, `--overwrite`. Overwrite re-queries and redownloads the
selected windows; it consumes PU. Run only one instance against a given output
directory/checklist. The older acquisition downloader uses a different default
directory and is not stopped or changed by this script. Running both at once
still shares your account's API quota and request limits.

## Inputs and outputs

Default manifest:
`data/density_aggregation/sen2_spataggr_fapar/density_cluster_4000m_aggregate_manifest.csv`.

Default output directory: `data/sentinel2_fapar_composite_inputs`.

Each window produces **one flat file**, `WINDOW_ID_snap_inputs.tif`, on the
manifest's UTM grid: 200 x 200 pixels, 20 m, 4 km x 4 km. There are no per-window
subdirectories or separate per-acquisition files. Progress and coverage counts
are recorded in one shared `composite_download_checklist.csv`.

The TIFF is UINT16 with NoData **65535 in every band** at unselected pixels.
Band descriptions, GDAL scales, sensor definitions, processing settings and the
candidate source-acquisition inventory are embedded in the TIFF. There is no
pixel-level source-date or source-scene-index band, as requested.

| Bands | Names, in file order | Decode valid stored values |
| --- | --- | --- |
| 1–8 | B03, B04, B05, B06, B07, B8A, B11, B12 | Divide by 10000 for surface reflectance |
| 9–12 | sunAzimuthAngles, sunZenithAngles, viewAzimuthMean, viewZenithMean | Divide by 100 for degrees |
| 13–15 | SCL, CLD, dataMask | Unscaled; CLD is percent |
| 16 | sensor_id | 1 = S2A; 2 = S2B; 65535 = NoData |

Mask NoData **before** scaling. Rasterio `read()` returns stored values, not
automatically scaled values; some other readers apply GDAL scales themselves.
Avoid double scaling. Angles are provided on the 20 m output grid but originate
from the API's much coarser angle layers, not independent 20 m measurements.

## Server-side selection

1. A small 4 x 4 dataMask request retrieves scene metadata over the inclusive
   target-date +/-8 UTC calendar-day interval.
2. Existing scene-selection logic keeps only the assigned MGRS tile and the
   latest reprocessing of each acquisition. Unsupported sensors fail the
   window explicitly; they are not silently treated as S2A/S2B.
3. One full-resolution Process API request admits only those scene IDs. At each
   pixel, eligibility requires dataMask=1, SCL in {4,5}, CLD in [0,40] percent,
   finite packable reflectances and physically valid angles. This is not SNAP
   input-domain validation. There is no canopy-model flag filtering yet.
4. Choose nearest UTC calendar day. On equal day distance, prefer the future
   day. Within the same day, prefer lower CLD; remaining ties use the earlier
   acquisition timestamp, then scene ID. Copy **all** bands and angles from the
   same chosen sample and attach its sensor ID. No bandwise averaging or
   synthetic spatial gap filling occurs.

This normally means two requests per window (one tiny discovery, one composite),
not one request per acquisition. It reduces local downloads and later SNAP
work but does not promise a particular speedup or PU reduction. No-acquisition
windows get an all-NoData TIFF without the full-resolution request.

## Subsequent SNAP processing

These are **inputs**, not FAPAR outputs. Do not feed raw packed integers directly
to the network. The later processing step must decode the bands, map their names
to SNAP's required input names, and respect sensor_id:

- Apply the S2A network to pixels with sensor_id=1.
- Apply the S2B network to pixels with sensor_id=2.
- Combine the corresponding FAPAR and SNAP-flag values and preserve NoData.

One mixed-sensor TIFF can therefore require two SNAP runs, followed by masking
and combination. Never choose one sensor for the entire mixed-sensor raster.
The old per-acquisition processing scripts do not automatically process this
new directory; the composite-aware local processing step is deferred.

The eventual accepted flags can still be 0,1,2,3, but **flag-0-first temporal
selection is no longer possible**: alternate acquisitions have already been
discarded. A selected input might yield an unacceptable FAPAR flag and leave a
gap despite another acquisition being suitable. Clear input coverage reported
by this downloader is therefore not validated FAPAR coverage.

## Resume and validation

The script reads and validates existing TIFFs rather than trusting the checklist
alone. It checks the grid, bands, scales, settings fingerprint, complete input
vectors, sensor values and screening rules. New TIFFs are written to temporary
files, read back, and atomically replaced only after validation. Failures are
recorded per window and retried on rerun. Empty completed composites are also
skipped; use `--aggregation-id WINDOW_ID --overwrite` to query one again.

No Python execution or live API test was performed during implementation;
verification was by source review, in accordance with the project instructions.

API reference for TILE sample/scene correspondence and output encoding:
[Sentinel Hub Evalscript V3](https://docs.sentinel-hub.com/api/latest/evalscript/v3/).
