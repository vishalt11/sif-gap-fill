"""Visualize predictor channels and SIF footprint masks from one NPZ chip."""

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np


# ---------------------------------------------------------------------------
# Config

DATA_DIR = Path(
    "data/cnn_sentinel2_chips/multisif_6km_20m_indices_fapar_active_crop"
)

# Leave as None to use the first chips_*.npz shard in DATA_DIR.
SHARD_PATH: Path | None = None
CHIP_INDEX = 0

OUTPUT_DIR = DATA_DIR / "visualizations"
SHOW_FIGURES = True


# ---------------------------------------------------------------------------
# Plot helpers


FRACTION_CHANNELS = {
    "fapar",
    "winter_wheat_fraction",
    "winter_barley_fraction",
    "winter_rye_fraction",
    "maize_fraction",
    "grass_forage_fraction",
    "woody_fraction",
    "other_crop_fraction",
    "active_crop_fraction",
    "non_crop_fraction",
}

SIGNED_INDEX_CHANNELS = {"ndmi", "ndvi", "evi", "ndre"}
MONTH_CHANNELS = {"month_sin", "month_cos"}


def choose_shard() -> Path:
    if SHARD_PATH is not None:
        path = Path(SHARD_PATH)
        if not path.exists():
            raise FileNotFoundError(path)
        return path

    matches = sorted(DATA_DIR.glob("chips_*.npz"))
    if not matches:
        raise FileNotFoundError(f"No chips_*.npz files found in {DATA_DIR}")
    return matches[0]


def decode_strings(values: np.ndarray) -> list[str]:
    decoded = []
    for value in values:
        if isinstance(value, bytes):
            decoded.append(value.decode("utf-8"))
        else:
            decoded.append(str(value))
    return decoded


def finite_limits(array: np.ndarray) -> tuple[float, float] | None:
    finite = array[np.isfinite(array)]
    if finite.size == 0:
        return None

    low, high = np.percentile(finite, [2, 98])
    if np.isclose(low, high):
        padding = max(abs(float(low)) * 0.05, 0.01)
        low -= padding
        high += padding
    return float(low), float(high)


def channel_style(name: str, array: np.ndarray) -> tuple[str, float, float] | None:
    if not np.isfinite(array).any():
        return None

    if name in FRACTION_CHANNELS:
        return "viridis", 0.0, 1.0
    if name in SIGNED_INDEX_CHANNELS:
        return "RdYlGn", -1.0, 1.0
    if name in MONTH_CHANNELS:
        return "coolwarm", -1.0, 1.0

    limits = finite_limits(array)
    if limits is None:
        return None
    return "viridis", limits[0], limits[1]


def add_mask_overlay(
    axis,
    background: np.ndarray,
    masks: np.ndarray,
    footprint_valid: np.ndarray,
) -> None:
    background_values = background.astype(np.float32, copy=True)
    limits = finite_limits(background_values)
    if limits is None:
        background_values = np.zeros(background.shape, dtype=np.float32)
        limits = (0.0, 1.0)

    axis.imshow(
        background_values,
        cmap="gray",
        vmin=limits[0],
        vmax=limits[1],
        interpolation="nearest",
    )

    colors = plt.get_cmap("tab10")(np.arange(10))
    for slot in np.flatnonzero(footprint_valid):
        mask = masks[slot]
        rgba = np.zeros((*mask.shape, 4), dtype=np.float32)
        rgba[..., :3] = colors[slot, :3]
        rgba[..., 3] = np.clip(mask, 0.0, 1.0) * 0.80
        axis.imshow(rgba, interpolation="nearest")

    axis.set_title("SIF masks over NDVI", fontsize=10)
    axis.set_xticks([])
    axis.set_yticks([])


def plot_channels(
    x: np.ndarray,
    channel_names: list[str],
    masks: np.ndarray,
    footprint_valid: np.ndarray,
    chip_id: str,
    output_path: Path,
) -> None:
    n_panels = len(channel_names) + 1
    n_cols = 4
    n_rows = int(np.ceil(n_panels / n_cols))

    figure, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(18, 4.2 * n_rows),
        constrained_layout=True,
    )
    axes = np.asarray(axes).reshape(-1)

    for channel_index, (name, channel) in enumerate(zip(channel_names, x)):
        axis = axes[channel_index]
        style = channel_style(name, channel)

        if style is None:
            axis.set_facecolor("#eeeeee")
            axis.text(
                0.5,
                0.5,
                "No valid pixels",
                ha="center",
                va="center",
                transform=axis.transAxes,
            )
        else:
            cmap, vmin, vmax = style
            image = axis.imshow(
                channel,
                cmap=cmap,
                vmin=vmin,
                vmax=vmax,
                interpolation="nearest",
            )
            figure.colorbar(image, ax=axis, fraction=0.046, pad=0.03)

        valid_fraction = float(np.isfinite(channel).mean())
        axis.set_title(f"{channel_index}: {name}\nvalid={valid_fraction:.3f}", fontsize=10)
        axis.set_xticks([])
        axis.set_yticks([])

    ndvi_index = channel_names.index("ndvi")
    add_mask_overlay(
        axes[len(channel_names)],
        x[ndvi_index],
        masks,
        footprint_valid,
    )

    for axis in axes[n_panels:]:
        axis.axis("off")

    figure.suptitle(f"Sentinel-2 CNN predictors: {chip_id}", fontsize=16)
    figure.savefig(output_path, dpi=180, bbox_inches="tight")

    if SHOW_FIGURES:
        plt.show()
    else:
        plt.close(figure)


def plot_footprint_masks(
    masks: np.ndarray,
    footprint_valid: np.ndarray,
    sif_row_ids: np.ndarray,
    y_targets: np.ndarray,
    target_accept: np.ndarray,
    chip_id: str,
    output_path: Path,
) -> None:
    figure, axes = plt.subplots(2, 5, figsize=(18, 7.5), constrained_layout=True)
    axes = axes.reshape(-1)

    colors = plt.get_cmap("tab10")(np.arange(10))
    valid_slots = set(np.flatnonzero(footprint_valid).tolist())

    for slot, axis in enumerate(axes):
        if slot not in valid_slots:
            axis.axis("off")
            continue

        mask = masks[slot]
        cmap = LinearSegmentedColormap.from_list(
            f"mask_{slot}",
            [(1.0, 1.0, 1.0), colors[slot, :3]],
        )
        image = axis.imshow(mask, cmap=cmap, vmin=0.0, vmax=1.0, interpolation="nearest")
        figure.colorbar(image, ax=axis, fraction=0.046, pad=0.03)

        target = float(y_targets[slot])
        accepted = bool(target_accept[slot])
        axis.set_title(
            f"slot={slot}, row={int(sif_row_ids[slot])}\n"
            f"SIF={target:.4f}, accept={accepted}",
            fontsize=10,
        )
        axis.set_xticks([])
        axis.set_yticks([])

    figure.suptitle(f"Fractional OCO-2 footprint masks: {chip_id}", fontsize=16)
    figure.savefig(output_path, dpi=180, bbox_inches="tight")

    if SHOW_FIGURES:
        plt.show()
    else:
        plt.close(figure)


# ---------------------------------------------------------------------------
# Main


def main() -> None:
    shard_path = choose_shard()

    with np.load(shard_path, allow_pickle=False) as shard:
        n_chips = shard["X"].shape[0]
        if CHIP_INDEX < 0 or CHIP_INDEX >= n_chips:
            raise IndexError(
                f"CHIP_INDEX={CHIP_INDEX} is outside this shard's 0-{n_chips - 1} range"
            )

        x = shard["X"][CHIP_INDEX].astype(np.float32)
        masks = shard["footprint_masks"][CHIP_INDEX].astype(np.float32)
        y_targets = shard["y_targets"][CHIP_INDEX].astype(np.float32)
        target_accept = shard["target_accept"][CHIP_INDEX].astype(bool)
        footprint_valid = shard["footprint_valid"][CHIP_INDEX].astype(bool)
        sif_row_ids = shard["sif_row_ids"][CHIP_INDEX]
        channel_names = decode_strings(shard["channel_names"])
        chip_id = decode_strings(shard["chip_id"])[CHIP_INDEX]

    if x.shape[0] != len(channel_names):
        raise ValueError(
            f"X has {x.shape[0]} channels but channel_names has {len(channel_names)}"
        )

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    safe_chip_id = "".join(
        character if character.isalnum() or character in "-_" else "_"
        for character in chip_id
    )

    channel_output = OUTPUT_DIR / f"{safe_chip_id}_channels.png"
    mask_output = OUTPUT_DIR / f"{safe_chip_id}_footprint_masks.png"

    plot_channels(
        x,
        channel_names,
        masks,
        footprint_valid,
        chip_id,
        channel_output,
    )
    plot_footprint_masks(
        masks,
        footprint_valid,
        sif_row_ids,
        y_targets,
        target_accept,
        chip_id,
        mask_output,
    )

    print(f"Shard: {shard_path}")
    print(f"Chip index: {CHIP_INDEX}")
    print(f"Chip ID: {chip_id}")
    print(f"Valid footprints: {int(footprint_valid.sum())}")
    print(f"Wrote {channel_output}")
    print(f"Wrote {mask_output}")


if __name__ == "__main__":
    main()
