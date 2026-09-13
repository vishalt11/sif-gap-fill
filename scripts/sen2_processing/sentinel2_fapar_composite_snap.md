# Local SNAP FAPAR: small-batch composite test

This script reads the existing packed composites. It does **not** download
anything, alter those inputs, or import the downloader scripts. Requirements:
the existing Python environment with NumPy and rasterio, and ESA SNAP with the
Optical Toolbox / BiophysicalOp installed.

## Recommended first test

```powershell
python "D:/UE/4_Semester/code/scripts/sen2_processing/sentinel2_fapar_composite_snap.py" --gpt "C:/Program Files/esa-snap/bin/gpt.exe" --limit 5 --batch-size 5 --workers 1 --compare-individual
```

This selects the first five TIFFs in filename order (not a random sample):

1. Check the packed band layout, metadata scales, grid and sensor IDs.
2. Decode physical inputs once using rasterio's scale metadata. Apply NoData
   masking before scaling: raw 2500 becomes reflectance 0.25, not 0.000025.
3. Prepare named ENVI inputs and run the five windows in one multi-branch graph.
   For each window only the required S2A/S2B models are included. A mixed-sensor
   window has two branches in the same graph, not two extra program launches.
4. Rerun those five windows as individual graphs and compare raw FAPAR and
   sensor-selected flags before output filtering. Flags must be identical and
   FAPAR must match within absolute tolerance 0.000001 (paired NaNs agree).
5. Only after all comparisons in the batch pass, publish the filtered output
   pairs. No new pair from that batch is published if comparison fails.

The test therefore normally launches GPT once for operator help, once for the
batch, and five times for the individual reference graphs. This deliberate test
overhead is **not** the full-production processing rate. The report separates
batch GPT time from reference GPT time. Empty windows need no FAPAR execution.
Equivalence tests check batching, not absolute scientific accuracy or every
possible input/sensor combination. The console reports which models are present;
if none of the five windows is mixed-sensor, also test a known mixed-sensor window
using `--aggregation-id WINDOW_ID`.

Share the `run_summary.json` path printed at the end, or its contents. On error,
also share the failing batch log. Failed ENVI/SNAP intermediates are retained
for diagnosis. `--batch-size 1` is available to isolate a multi-window graph
problem without changing the scientific processing.

## Output

Default folder:
`D:/UE/4_Semester/code/scripts/data/sentinel2_fapar_composite_snap`.

Each window has two files directly in that folder:

- `WINDOW_ID_fapar.tif`: physical Float32 FAPAR, units fraction, NoData=-9999.
- `WINDOW_ID_snap_flag.tif`: UInt8 flags 0,1,2,3 at retained pixels; NoData=255
  elsewhere. Zero is a valid flag and must never be used as NoData.

Both retain the original 200 x 200, 20 m UTM grid and transform. There is no
reprojection, resampling, clipping of invalid FAPAR to [0,1], spatial gap filling,
or additional temporal selection. Validity requires original input availability,
dataMask=1, SCL 4/5, CLD <=40%, physically valid angles, and finite FAPAR in [0,1]
with flag exactly 0,1,2,3. Thus flag 4 is rejected despite being a clipping flag.
Flags 1 and 3 remain warning-bearing estimates, not unflagged retrievals.

For mixed-sensor windows, each model is evaluated on the coherent input image;
only S2A results at sensor_id=1 and S2B results at sensor_id=2 are retained. Values
calculated by the other model at those locations are discarded. Placeholder
inputs for missing pixels are excluded by the original validity mask.

Raw flag counts are saved in the report and TIFF metadata, including rejected
flags. Rejected pixels' spatial flag values are not kept in the final flag TIFF.
`valid_fraction` uses the full 40,000-pixel window as its denominator.

A shared `snap_processing_checklist.csv` tracks completion. `logs/RUN_ID/`
contains graphs, logs, configuration and JSON reports; these are diagnostic
files, not separate permanent acquisition folders. `_work/` contains temporary
ENVI and SNAP outputs. Successful temporary batch folders are deleted unless
`--keep-work` is set; failed ones remain. Retained graphs reference temporary
inputs, so use `--keep-work` if you need to rerun a graph manually afterward.

## Resume and subsequent full run

Without selection options the default is **five windows**, never all 8,347.
`--compare-individual` deliberately recomputes selected windows even if final
outputs exist; it is limited to 20 test windows. Without that switch, validated
matching pairs are skipped. Input content, processing version and GPT identity
are included in output fingerprints. Use `--overwrite` after toolbox/model
updates that leave the executable and operator help unchanged.

After confirming the small test, a full local run is explicitly enabled with:

```powershell
python "D:/UE/4_Semester/code/scripts/sen2_processing/sentinel2_fapar_composite_snap.py" --gpt "C:/Program Files/esa-snap/bin/gpt.exe" --all --batch-size 5 --workers 1
```

Try `--workers 2` only after assessing RAM and CPU. `--snap-threads 2` and
`--cache 512M` apply **per GPT process**; cache is not a cap on total JVM memory.
Multiple processes can use considerably more memory than their tile caches.
Do not start two script instances against the same output directory/checklist.
Changing batch size or workers does not invalidate scientific output caches.

The comparison runs after the batch and can benefit from warm disk/OS caches;
the timing comparison is approximate, not a controlled benchmark. A five-window
test does not establish the runtime for the full set. Logs and failures are
isolated per batch; rerunning retries incomplete pairs. If interrupted, a valid
pair may already exist even before its checklist record was saved; actual TIFF
validation takes precedence over the checklist on resume.

Implementation was reviewed without executing Python, SNAP, or a live test,
as required by this project's instructions. API references used for the graph:
[GPT options](https://step.esa.int/main/wp-content/help/versions/9.0.0/snap/org.esa.snap.snap.gpf.ui/gpf/GraphProcessingTool.html),
[SNAP graph processor](https://github.com/senbox-org/snap-engine/blob/master/snap-gpf/src/main/java/org/esa/snap/core/gpf/graph/GraphProcessor.java).
