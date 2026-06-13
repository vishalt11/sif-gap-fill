# Project Context

This project gap-fills spatially and temporally sparse SIF observations using spatially and temporally complete Sentinel-2 bands at 10 m resolution.

The downstream goal is early-warning crop yield prediction for winter wheat, focusing on the top three winter-wheat-producing German states.

Because complete state-wide multi-band Sentinel-2 processing is very large, the workflow should prioritize high-density winter wheat cluster regions. These clusters are derived from crop-type rasters and used to identify relevant Sentinel-2 tile/subtile regions for download and processing.

Core modeling idea:

SIF polygon mean ~ Sentinel-2 Band 1 + Sentinel-2 Band 2 + Sentinel-2 Band 3 + ...

Sentinel-2 10 m pixels are cropped to SIF polygons and mean-aggregated per polygon. The learned relationship is then used to estimate/fill missing SIF values in places/times where SIF is unavailable but Sentinel-2 bands exist.

Prefer workflow choices that reduce raster size early: crop/mask by wheat-density regions, administrative state boundaries, Sentinel-2 tile chunks, and SIF polygons before expensive operations.