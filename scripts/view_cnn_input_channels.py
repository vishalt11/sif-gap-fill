from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

CHIP_DIR = Path(__file__).resolve().parent / (
    "data/cnn_sentinel2_chips/spatial_aggregate_density_4km_20m_snap_fapar"
)
SHARD_NAME = "chips_00001.npz"
SAMPLE_INDEX = 5  # Sample within this shard.

with np.load(CHIP_DIR / SHARD_NAME, allow_pickle=False) as shard:
    x = shard["X"][SAMPLE_INDEX].astype(np.float32)
    names = shard["channel_names"].astype(str)
    window_id = str(shard["aggregation_id"][SAMPLE_INDEX])
    target = float(shard["y_aggregate"][SAMPLE_INDEX])

fig, axes = plt.subplots(4, 5, figsize=(20, 14), constrained_layout=True)
for ax, values, name in zip(axes.flat, x, names):
    limits = {}
    if name == "fapar" or name.endswith("_fraction"):
        limits = {"vmin": 0, "vmax": 1}
    elif name in ("doy_sin", "doy_cos"):
        limits = {"vmin": -1, "vmax": 1}
    image = ax.imshow(
        np.ma.masked_invalid(values), cmap="viridis", interpolation="nearest", **limits
    )
    ax.set_title(name.replace("_", " "), fontsize=11)
    ax.axis("off")
    fig.colorbar(image, ax=ax, fraction=0.046, pad=0.02)

fig.suptitle(f"{window_id}\nCNN inputs before normalization | aggregate SIF = {target:.4f}")
plt.show()
