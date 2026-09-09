from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
from rasterio.mask import mask
from rasterio.warp import transform_geom
from shapely.geometry import box, mapping

ROOT = Path(__file__).resolve().parent
WINDOW_ID = "s2_density_4000m_T33UUV_20190212_m1_s00001_db001_w001_lc60"
manifest = ROOT / "data/density_aggregation" / (
    "sentinel2_spatial_aggregation_density_4000m_landcover_combined_"
    "mask60_min4_s2valid95_par_available/density_cluster_4000m_aggregate_manifest.csv"
)
row = pd.read_csv(manifest, low_memory=False).set_index("aggregation_id").loc[WINDOW_ID]
date = pd.Timestamp(row["Delta_Date"])
doy = 1 + 8 * ((date.dayofyear - 1) // 8)  # Same containing-composite rule as prep.

# This window lies in MODIS tile h18v03.
path, = (ROOT / "data/glass_geotiff/fapar/h18v03" / str(date.year)).glob(
    f"*A{date.year}{doy:03d}.h18v03.*.tif"
)

with rasterio.open(path) as src:
    outline = mapping(box(row.cell_xmin, row.cell_ymin, row.cell_xmax, row.cell_ymax))
    outline = transform_geom(row.window_crs, src.crs, outline)
    data, _ = mask(src, [outline], crop=True, filled=False, indexes=1)
    fapar = data.astype("float32")  # Converted GeoTIFF already stores physical FAPAR.
    fapar = np.ma.masked_invalid(np.ma.masked_outside(fapar, 0, 1))

plt.figure(figsize=(6, 5))
plt.imshow(fapar, cmap="YlGn", vmin=0, vmax=1, interpolation="nearest")
plt.colorbar(label="FAPAR")
plt.title(f"GLASS FAPAR 250 m | {date.date()} | composite DOY {doy:03d}")
plt.axis("off")
plt.tight_layout()
plt.show()
