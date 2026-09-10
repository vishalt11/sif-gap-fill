# Sentinel-2 acquisition downloads for later SNAP processing

Run from `D:/UE/4_Semester/code/scripts` in your existing Python environment:

```powershell
python sen2_processing/sentinel2_fapar_acquisitions_download.py
```

For a small download first, append `--limit 1`. Use `--aggregation-id ID` for a
particular window. There is no `--gpt` argument: this script never runs SNAP.
It reads the `sen2_spataggr_fapar` manifest used by the production script.

The previous scripts remain in place. Keep these alongside the new script,
because their existing request/metadata helpers are reused:

- `sentinel2_fapar_snap_window_download.py`
- `sentinel2_fapar_snap_one_window_test.py`
- `sentinel2_l2a_nearest_clear_window_download.py`

The same CDSE credentials/profile are used. No quota stop or parallel SNAP
processing is added in this download-only stage.

## Local files

The default output directory is `data/sentinel2_fapar_acquisitions`:

- `acquisitions/`: one 15-band GeoTIFF per window/acquisition, named
  `WINDOW_ID_TIMESTAMP_SENSOR_SCENEID.tif`.
- `window_metadata/`: one JSON inventory per window, recording bounds, CRS,
  target date, source acquisition metadata, deduplication decisions and the
  relative paths of the saved TIFFs. No separate directory per window.
- `acquisition_download_checklist.csv`: shared progress and error checklist.

Do not discard the inventories: they identify the complete set of inputs to
use for each window in the later SNAP/compositing stage. TIFFs also embed the
source scene, sensor, observation time, target date and encoding metadata.

## Band encoding

All TIFFs have a 200 x 200 grid at 20 m in the window's UTM CRS. They store
unsigned 16-bit integers, not already-decoded floating-point reflectances.

| Band positions | Names | Physical value |
| --- | --- | --- |
| 1-8 | B03, B04, B05, B06, B07, B8A, B11, B12 | stored value / 10000 |
| 9-12 | sunAzimuthAngles, sunZenithAngles, viewAzimuthMean, viewZenithMean | stored value / 100, in degrees |
| 13-15 | SCL, CLD, dataMask | unchanged; CLD is percent |

**65535 is NoData and must be masked before scaling.** GDAL band scale/offset
metadata is included. Some readers apply it automatically; others, including
`rasterio.read()`, return the stored integers. Apply scaling exactly once.
Angle quantization is 0.01 degrees; the underlying angle grids are much coarser
than 20 m. This is the same UINT16 request encoding as the production script.

## Resume and processing later

All selected acquisitions in the inclusive +/-8-day UTC window are downloaded,
including cloudy ones. Scene selection retains the assigned MGRS tile and latest
reprocessing of each acquisition. No FAPAR is calculated; no cloud/SNAP flags are
used to discard pixels or choose dates at this stage.

On a normal rerun, saved inventories avoid repeated discovery requests and
validated TIFFs avoid repeated downloads. An interrupted window resumes from
its remaining acquisitions. `--refresh-scenes` queries current scene availability
again but reuses matching TIFFs; `--overwrite` deliberately redownloads the chips.
An empty inventory can be refreshed with `--refresh-scenes` if necessary.

Local inputs are retained until you remove them. They are not automatically
deleted after the download run. Distinct windows get separate cropped chips,
even when the same satellite acquisition covers both windows.

This script does not import old final FAPAR/flag TIFFs as acquisition inputs.
The earlier production script did not retain its acquisition downloads, so those
cannot be recovered from its final FAPAR outputs.

The later processor can use these files without contacting CDSE, run the actual
SNAP model per acquisition, then choose flag 0 first and flags 1/2/3 as fallback.
It will need its own implementation; the existing production script does not yet
read this new cache.

Code reviewed statically, not executed by Codex, per project instructions.
