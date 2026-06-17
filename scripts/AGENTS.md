# Project Context

This project gap-fills spatially and temporally sparse SIF observations using spatially and temporally complete Sentinel-2 bands at 10 m resolution.

The downstream goal is early-warning crop yield prediction for winter wheat, focusing on the top three winter-wheat-producing German states.

Because complete state-wide multi-band Sentinel-2 processing is very large, the workflow should prioritize high-density winter wheat cluster regions. These clusters are derived from crop-type rasters and used to identify relevant Sentinel-2 tile/subtile regions for download and processing.

Core modeling idea:

SIF polygon mean ~ Sentinel-2 Band 1 + Sentinel-2 Band 2 + Sentinel-2 Band 3 + ...

Sentinel-2 10 m pixels are cropped to SIF polygons and mean-aggregated per polygon. The learned relationship is then used to estimate/fill missing SIF values in places/times where SIF is unavailable but Sentinel-2 bands exist.

Prefer workflow choices that reduce raster size early: crop/mask by wheat-density regions, administrative state boundaries, Sentinel-2 tile chunks, and SIF polygons before expensive operations.


#------------------------------------




#------------------------------------

Chat: geodes_wasp_download summary:

We worked on Sentinel-2 L3A WASP downloads for SIF/Sentinel matching.
Main outcomes:
DLR STAC had historical S2_L3A_WASP metadata, but old V1-2 asset URLs returned 404, both via browser/curl and rasterio. DLR map service can visualize old data, but public download links were unusable for historical files.
We switched to CNES GEODES / pygeodes, collection:
THEIA_REFLECTANCE_SENTINEL2_L3A
GEODES uses:
grid:code = T32UNC style tile codes,
opaque item IDs like URN:FEATURE:DATA:gdh:...,
useful product name in identifier, e.g. SENTINEL2X_20200415-000000-000_L3A_T32UNC_C,
product date in start_datetime, usually YYYY-MM-15T00:00:00Z.

We confirmed GEODES zip archives contain the needed bands:

FRC_B4
FRC_B8
FRC_B11
plus other FRC_* bands, masks, metadata, and quicklook.
Scripts created/modified:
geodes_download_32unc_zip.py: one-off GEODES zip downloader, initially for 32UNC / 2020-04.
geodes_download_32unc_zip_curl.py: batch downloader using pygeodes for search and curl for downloading. It reads a SIF CSV, filters by mgrs_tile, gets unique Delta_Date month-year combos, queries GEODES, downloads full monthly tile zips, and maintains a checklist CSV.
Important checklist behavior:
Completed files are skipped on rerun if file exists and size matches.
Missing products are skipped on rerun if RECHECK_MISSING = False.
Failed and pending rows are retried.
Download observations:
GEODES downloads are real but variable speed, roughly <1 MB/s to 12 MB/s.
Some downloads fail midway with IncompleteRead / ChunkedEncodingError.
GEODES API download endpoint does not appear to support HTTP byte-range resume, so curl -C - cannot truly resume partial files. Retrying can still succeed if a full attempt completes.
Full zip downloads are around 1.1-1.2 GB each.
State/checklists:
32UNC checklist finished: 27 completed, 8 missing, no failed/pending.
32UQV checklist finished: 28 completed, 10 missing, no failed/pending.
User then changed geodes_download_32unc_zip_curl.py to:SIF_CSV = data/sa_sif_mgrs_crop_composition.csv
TILE = 32UPC
GRID_CODE = T32UPC
output/checklist under data/geodes_wasp_zips/32UPC/


#------------------------------------

Chat: niedersachsen_test crop composition summary:

We worked on `niedersachsen_test.R` to calculate crop-type composition for Niedersachsen, Bavaria, Sachsen-Anhalt SIF observations. (Manually changed the filters for each region)

Main task:
Keep only SIF center points that fall inside MGRS tiles:

niedersachsen: 32UND 32UNC 32UPC 
bavaria: 32UQV 32UNA 32UPA 32UPV 32UQU 32UPU
sachsen-anhalt: 32UPB 32UQB 32UQC 32UPC

then build each SIF footprint polygon from its four lat/lon corners and compute raw crop-type raster pixel counts inside that polygon.

Important implementation details:
Use `mgrs_de` from `data/mgrs_de.rds` and crop rasters from `data/crop_type_tif/croptypes_<year>.tif`.
Use the SIF row `Delta_Date` year to select the matching crop raster. Years after 2024 are excluded because crop rasters only exist through 2024.
Use the raster legend from `data/crop_type_tif/LEGEND_CropTypes.txt`.
After loading each yearly raster, set levels with:
`levels(crop_raster) <- data.frame(value = crop_classes$code, crop = crop_classes$label)`
The crop rasters are in EPSG:32632, so SIF polygons are transformed to the raster CRS before extraction.
The extraction uses cell-center/count behavior, not area weighting.

Output:
The final objects are saved as:
`data/ns_sif_mgrs_crop_composition.rds`
`data/ba_sif_mgrs_crop_composition.rds`
`data/sa_sif_mgrs_crop_composition.rds`
It keeps the SIF center geometry and adds raw count columns such as:
`crop_count_winter_wheat`, `crop_count_winter_barley`, `crop_count_maize`, etc.
It also has `crop_pixel_count`, which is the row sum of all crop count columns.

Bug fixed:
`pivot_wider()` originally errored with `Can't convert fill <double> to <list>`.
The fix was to build an explicit `sif_id x crop_code` template, left join observed counts, fill missing counts with `0L`, and then widen with `values_fill = list(n = 0L)` and `values_fn = list(n = sum)`.
The code also handles `terra::extract()` returning either numeric crop codes or crop labels.

Interpretation notes:
Some SIF polygons have `crop_pixel_count == 0`.
After visual inspection, this appears correct because the crop-type raster has no raster values over some non-crop areas; those areas are not represented as class `0 / no_data`.
This is also why `sum(ns_sif_mgrs_crop_composition$crop_count_no_data)` was `0`.
So zero pixel count generally means no crop raster pixels were extracted for that SIF polygon, not that the SIF polygon code failed.

Winter wheat share:
To add winter wheat fraction without multiplying by 100:
`ns_sif_mgrs_crop_composition <- ns_sif_mgrs_crop_composition %>% mutate(ww_pct = if_else(crop_pixel_count > 0, crop_count_winter_wheat / crop_pixel_count, NA_real_))`
Rows with `crop_pixel_count == 0` should get `NA_real_` because `0 / 0` is undefined.
