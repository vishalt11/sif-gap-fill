from contextlib import ExitStack
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
from rasterio.io import MemoryFile
from rasterio.mask import mask
from rasterio.merge import merge
from rasterio.warp import transform_geom
from shapely.geometry import box, mapping

ROOT = Path(__file__).resolve().parent
WINDOW_ID = "s2_density_4000m_T32ULA_20220629_m1_s00036_db001_w006_lc10"
manifest = ROOT / "data/density_aggregation" / (
    "sentinel2_spatial_aggregation_density_4000m_landcover_combined_"
    "mask60_min4_s2valid95_par_available/density_cluster_4000m_aggregate_manifest.csv"
)
row = pd.read_csv(manifest, low_memory=False).set_index("aggregation_id").loc[WINDOW_ID]
date = pd.Timestamp(row["Delta_Date"])
doy = 1 + 8 * ((date.dayofyear - 1) // 8)  # Same containing-composite rule as prep.

# Find both tiles for exactly the same containing 8-day composite.
paths = []
for tile in ("h18v03", "h18v04"):
    folder = ROOT / "data/glass_geotiff/fapar" / tile / str(date.year)
    matches = sorted(folder.glob(f"*A{date.year}{doy:03d}.{tile}.*.tif"))
    if len(matches) != 1:
        raise ValueError(
            f"Expected exactly one GLASS FAPAR file for {tile}, "
            f"A{date.year}{doy:03d}; found {len(matches)} in {folder}: {matches}"
        )
    paths.append(matches[0])

with ExitStack() as stack:
    sources = [stack.enter_context(rasterio.open(path)) for path in paths]
    first = sources[0]
    if first.crs is None or any(
        src.crs != first.crs
        or not np.allclose(src.res, first.res, rtol=0, atol=1e-6)
        for src in sources
    ):
        raise ValueError("The GLASS tiles must have the same CRS and resolution to merge")

    print(f"Merging A{date.year}{doy:03d}: h18v03 + h18v04", flush=True)
    # Merge first, then crop. Preserve valid zero FAPAR by using NaN as NoData.
    mosaic, mosaic_transform = merge(
        sources, indexes=[1], dtype="float32", nodata=np.nan, method="first"
    )
    outline = mapping(box(row.cell_xmin, row.cell_ymin, row.cell_xmax, row.cell_ymax))
    outline = transform_geom(row.window_crs, first.crs, outline)
    # An in-memory dataset lets mask() crop the merged raster without creating
    # a permanent full-tile mosaic or changing either downloaded source file.
    with MemoryFile() as memory:
        with memory.open(
            driver="GTiff", height=mosaic.shape[1], width=mosaic.shape[2],
            count=1, dtype="float32", crs=first.crs,
            transform=mosaic_transform, nodata=np.nan,
        ) as merged:
            merged.write(mosaic)
            data, _ = mask(merged, [outline], crop=True, filled=False, indexes=1)
    fapar = data.astype("float32")  # Converted GeoTIFF already stores physical FAPAR.
    fapar = np.ma.masked_invalid(np.ma.masked_outside(fapar, 0, 1))

plt.figure(figsize=(6, 5))
plt.imshow(fapar, cmap="YlGn", vmin=0, vmax=1, interpolation="nearest")
plt.colorbar(label="FAPAR")
plt.title(f"GLASS FAPAR 250 m | {date.date()} | composite DOY {doy:03d}")
plt.axis("off")
plt.tight_layout()
plt.show()
