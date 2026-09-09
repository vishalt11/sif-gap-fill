# One-window Sentinel-2 SNAP FAPAR test

Run `sentinel2_fapar_snap_one_window_test.py` in your existing Python environment.
The original nearest-clear downloader is unchanged; this script imports only its
CDSE configuration and MGRS CRS helper. Keep both Python files in the same folder.

## Run

From `D:/UE/4_Semester/code/scripts`:

```powershell
python sen2_processing/sentinel2_fapar_snap_one_window_test.py
```

ESA SNAP, with its Optical Toolbox/BiophysicalOp installed, is required. The
script checks `SNAP_GPT`, PATH, and two common Windows installation locations.
If necessary, specify your actual executable:

```powershell
python sen2_processing/sentinel2_fapar_snap_one_window_test.py --gpt "C:/Program Files/esa-snap/bin/gpt.exe"
```

There is no dependency on Python's `esa_snappy` bridge. SNAP runs as a separate
command-line process. Its operator help is checked before imagery downloads.
The existing `CDSE_SH_CLIENT_ID` / `CDSE_SH_CLIENT_SECRET` or `SH_PROFILE`
authentication configuration is reused.

To download the inputs first and defer SNAP:

```powershell
python sen2_processing/sentinel2_fapar_snap_one_window_test.py --download-only
```

Rerun without that option to reuse verified input TIFFs and calculate FAPAR.
SNAP calculations are repeated on a rerun; valid acquisition inputs are cached.

## Window and scene selection

The actual manifest found in this project is:

`data/density_aggregation/sentinel2_spatial_aggregation_density_4000m_landcover_combined_mask60_min4_s2valid95_par_available/density_cluster_4000m_aggregate_manifest.csv`

The default is exactly the first row:

- Aggregation: `s2_density_4000m_T32UPV_20190206_m0_s00346_db001_w001_lc10`
- Target date: 2019-02-06
- Search: 2019-01-29 through 2019-02-14, inclusive UTC calendar dates
- Bounds: 695900, 5400440, 699900, 5404440 in EPSG:32632
- Output grid: 200 x 200 pixels, each 20 x 20 m

Use `--aggregation-id YOUR_ID` to select another single row, or `--manifest PATH`
to override the manifest. There is deliberately no all-windows option.

The discovery request records every returned scene. Only the assigned MGRS tile
and dates strictly within +/-8 calendar days are retained. Reprocessed versions
of an acquisition are grouped using product identity (sensor, sensing time,
orbit, tile), with the latest production timestamp retained. Legacy tile IDs use
sensor, absolute orbit and tile. This preserves distinct passes on the same day.
Each full-resolution request filters to one exact Sentinel Hub `shId` and checks
the returned metadata, so reflectance and angles cannot come from different scenes.
Unrecognized in-scope product identifiers cause an explicit error for inspection.

## Inputs and SNAP processing

Each acquisition gets a 19-channel float32 GeoTIFF:

| Channels | Meaning |
| --- | --- |
| B01, B02, B03, B04, B05, B06, B07, B08, B8A, B09, B11, B12 | Unitless L2A BOA reflectance |
| sunAzimuthAngles, sunZenithAngles, viewAzimuthMean, viewZenithMean | Acquisition angles in degrees |
| SCL, CLD, dataMask | Scene classification, cloud probability in percent, coverage mask |

B10 is not provided in L2A. All channels share the manifest's 20 m UTM grid.
Their source resolutions differ: 10/20/60 m spectral bands and approximately
5 km angle grids. Exporting angles at 20 m does not create 20 m angle information.
This test explicitly uses nearest-neighbour up/downsampling, preserving QA codes
and matching the old downloader's default sampling. Consequently, it is not
claimed to reproduce a full SAFE workflow's interpolation pixel-for-pixel.
CDSE converts reflectance to physical units with baseline harmonization enabled;
the script does not apply an additional factor of 1/10000.

The eight spectral inputs and four angle inputs required by SNAP are exported
to a named-band ENVI intermediate. The generated XML graph executes the installed
`BiophysicalOp`, selecting the S2A or S2B model from acquisition metadata and
requesting FAPAR only. S2C is not silently substituted with an S2A/S2B model.
The original SNAP DIMAP output, flags and execution log remain available.

## Quality and temporal selection

FAPAR is computed separately for each acquisition. A pixel is eligible for the
final composite when:

- Required reflectance and angle values are finite and present; zenith and azimuth
  values pass basic physical-range checks.
- SCL is 4 (vegetation) or 5 (bare soil), with cloud probability 0-40%.
- FAPAR is finite and between 0 and 1, and SNAP's FAPAR flag value is zero.

This is deliberately stricter than the old reflectance composite: water is
excluded and every SNAP warning flag is rejected for this first test. SNAP flags
are bit fields: 1 = input outside definition domain; 2/4 = output clipped within
tolerance to the minimum/maximum; 8/16 = output too low/high. Flags can combine.
The raw FAPAR and flags are saved so this policy can be assessed after the test.

Eligible pixels are chosen by smallest absolute calendar-day distance, then
future date on an equal-distance tie, then lower cloud probability on a same-day
tie. Remaining exact ties follow the deterministic acquisition order recorded in
the summary. Missing pixels remain NoData (-9999), with no neighbour filling.

## Outputs

All outputs go to:

`data/sentinel2_fapar_snap_one_window_test/<aggregation_id>/`

For each acquisition, under `acquisitions/<timestamp>_<sensor>_<scene_id>/`:

- `l2a_inputs_20m.tif`: all 19 input channels, including cloudy observations.
- `fapar_raw_20m.tif`: unfiltered network output where inputs are numerically valid.
- `fapar_quality_20m.tif`: SNAP flags, SCL, CLD, clear-land and usable-FAPAR masks.
- `fapar_clear_20m.tif`: only eligible FAPAR pixels.
- `snap_input.img/.hdr`, `snap_graph.xml`, `snap_fapar.dim/.data`, `snap_gpt.log`.
- Source-scene and returned-download metadata JSON files.

At the window level:

- `fapar_nearest_clear_20m.tif`: final one-band FAPAR composite.
- `fapar_nearest_clear_provenance_20m.tif`: source index, day offset, SCL, cloud
  probability, SNAP flags and validity mask. Source indices are one-based and
  map to `source_index` in the acquisition summary.
- `l2a_inputs_matching_fapar_selection_20m.tif`: all 19 inputs selected from the
  same observations as the final FAPAR, for consistent downstream CNN channels.
- `all_scene_metadata.json` and `scene_selection_audit.json`: selection evidence.
- `acquisition_summary.json`: per-acquisition status, coverage and flag counts.
- `run_configuration.json` and `snap_operator_help.txt`: settings and installed
  operator information.
- `completed.json`: written only after every selected acquisition and the final
  composite succeed. A failure stops the workflow rather than composing a
  silently incomplete set of observations. On reruns, old final TIFFs may still
  exist until replaced; use this completion marker to identify a successful run.

For debugging, share the console error plus `acquisition_summary.json` and the
failing acquisition's `snap_gpt.log`. For scene/date questions, share
`scene_selection_audit.json`. These logs do not intentionally record credentials.

## Method references and validation status

- [CDSE L2A input bands, units and resampling](https://documentation.dataspace.copernicus.eu/APIs/SentinelHub/Data/S2L2A.html)
- [SNAP BiophysicalOp input names and parameters](https://github.com/senbox-org/s2tbx/blob/master/s2tbx-biophysical/src/main/java/org/esa/s2tbx/biophysical/BiophysicalOp.java)
- [SNAP FAPAR flags](https://github.com/senbox-org/s2tbx/blob/master/s2tbx-biophysical/src/main/java/org/esa/s2tbx/biophysical/BiophysicalFlag.java)

Code was reviewed statically. It has not been executed or used to download data
by Codex, in accordance with this project's instructions. The smoke test on
your machine will verify CDSE metadata and the installed SNAP reader/operator.
