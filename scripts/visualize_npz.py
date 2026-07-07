import numpy as np
import matplotlib.pyplot as plt

npz_path = "data/cnn_modis_chips/multisif_24x24_5to10_active_crop/chips_00000.npz"
d = np.load(npz_path, allow_pickle=True)

X = d["X"]
masks = d["footprint_masks"]
y = d["y_targets"]
valid = d["footprint_valid"]
sif_row_ids = d["sif_row_ids"]
channel_names = d["channel_names"]

chip_i = 0
target_i = 2  # 0=757, 1=771, 2=target_modis_sif

valid_slots = np.where(valid[chip_i] == 1)[0]
mask_union = masks[chip_i, valid_slots].sum(axis=0)
mask_union = np.clip(mask_union, 0, 1)

evi_i = list(channel_names).index("evi")
ndvi_i = list(channel_names).index("ndvi")
fapar_i = list(channel_names).index("fapar")

fig, axes = plt.subplots(1, 4, figsize=(14, 4))

axes[0].imshow(X[chip_i, fapar_i], cmap="viridis")
axes[0].set_title("FAPAR")

axes[1].imshow(X[chip_i, evi_i], cmap="viridis")
axes[1].set_title("EVI")

axes[2].imshow(X[chip_i, ndvi_i], cmap="viridis")
axes[2].set_title("NDVI")

axes[3].imshow(X[chip_i, evi_i], cmap="gray")
axes[3].imshow(mask_union, cmap="Reds", alpha=0.5)
axes[3].set_title("Footprint masks on EVI")

for ax in axes:
    ax.axis("off")

plt.tight_layout()
plt.show()

print("valid slots:", valid_slots)
print("sif row ids:", sif_row_ids[chip_i, valid_slots])
print("targets:", y[chip_i, valid_slots, target_i])