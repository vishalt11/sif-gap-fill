# Project Context

This project gap-fills spatially and temporally sparse oco-2 SIF observations using spatially and temporally complete Sentinel-2 bands at 10 m resolution.

Main goal is to be able to get a SIF value given the pixel's ndvi, nirv ... spectral indices values and band values.

We can then construct oco-2 SIF level resoluton sif polygons and mean the 10m x 10m gap filled sif values falling in it. The smallest resolution will be 10m x 10m cause that is what the sentinel - 2 bands are at.

The downstream application is early-warning crop yield prediction for crops (starting off with winter wheat), focusing on the top three winter-wheat-producing German states. SIF has been shown to be the single biggest factor in predicting crop yield, due it's close direct relationship with gpp by past studies. 

We will use the gap filled SIF data to find the respective SIF for the growing months of the wheat fields and derive the crop yield / gpp.
Because complete state-wide multi-band Sentinel-2 processing is very large, the workflow should prioritize high-density winter wheat cluster regions. These clusters are derived from crop-type rasters and used to identify relevant Sentinel-2 tile/subtile regions for download and processing.

Core modeling idea:

SIF polygon mean ~ Sentinel-2 Band 1 + Sentinel-2 Band 2 + Sentinel-2 Band 3 + ...

Sentinel-2 10 m pixels are cropped to SIF polygons and mean-aggregated per polygon. The learned relationship is then used to estimate/fill missing SIF values in places/times where SIF is unavailable but Sentinel-2 bands exist.

Prefer workflow choices that reduce raster size early: crop/mask by wheat-density regions, administrative state boundaries, Sentinel-2 tile chunks, and SIF polygons before expensive operations.


# Code Running Clarification

Never run code as the path's for python and R executables you see are not the ones I am running. Just review the code and make sure the logic is correct. I will run the code on my end and report any issues / results back to you.

Assume that I have all packages for python and R installed.

Don't mess around with git for this project. I am doing git on my end.


# Base Data forms

Bavaria only SIF data at data/ba_sif_mgrs_crop_composition.csv (there are other files for niedersachsen, sachsen-anhalt), but you get the general structure from here.

Sentinel-2 bands in data/geodes_wasp_zips
