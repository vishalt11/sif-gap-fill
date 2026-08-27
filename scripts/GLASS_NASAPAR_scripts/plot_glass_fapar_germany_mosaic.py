from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import rasterio
from shapely.geometry import box
from rasterio.merge import merge


INPUT_DIR = Path("temp")
YEAR = 2024
DOY = 209
TILES = ["h18v03", "h18v04", "h19v03", "h19v04"]
TILE_COLORS = {
    "h18v03": "red",
    "h18v04": "blue",
    "h19v03": "yellow",
    "h19v04": "green",
}

OUTPUT_DIR = Path("eda_images")
OUTPUT_PNG = OUTPUT_DIR / f"glass_fapar_germany_mosaic_{YEAR}_{DOY:03d}.png"

NATURAL_EARTH_COUNTRIES_URL = (
    "https://naturalearth.s3.amazonaws.com/10m_cultural/"
    "ne_10m_admin_0_countries.zip"
)


def find_tif_files() -> list[Path]:
    tif_files: list[Path] = []

    for tile in TILES:
        matches = sorted(INPUT_DIR.glob(f"GLASS09D01.V60.A{YEAR}{DOY:03d}.{tile}.*.tif"))
        if len(matches) != 1:
            raise FileNotFoundError(
                f"Expected exactly one GeoTIFF for {tile}, found {len(matches)} in {INPUT_DIR}"
            )
        tif_files.append(matches[0])

    return tif_files


def germany_border(target_crs) -> gpd.GeoDataFrame:
    countries = gpd.read_file(NATURAL_EARTH_COUNTRIES_URL)
    germany = countries.loc[countries["ADMIN"].eq("Germany")].copy()

    if germany.empty:
        raise RuntimeError("Germany boundary not found in Natural Earth countries file.")

    return germany.to_crs(target_crs)


def tile_footprints(srcs: list[rasterio.io.DatasetReader]) -> gpd.GeoDataFrame:
    rows = []

    for tile, src in zip(TILES, srcs):
        rows.append(
            {
                "tile": tile,
                "color": TILE_COLORS[tile],
                "geometry": box(*src.bounds),
            }
        )

    return gpd.GeoDataFrame(rows, crs=srcs[0].crs)


def plot_mosaic(
    mosaic: np.ndarray,
    transform,
    crs,
    footprints: gpd.GeoDataFrame,
) -> None:
    fapar = mosaic[0].astype("float32")
    fapar = np.ma.masked_invalid(fapar)

    left = transform.c
    top = transform.f
    right = left + transform.a * fapar.shape[1]
    bottom = top + transform.e * fapar.shape[0]

    germany = germany_border(crs)

    fig, ax = plt.subplots(figsize=(11, 10))
    image = ax.imshow(
        fapar,
        extent=(left, right, bottom, top),
        cmap="YlGn",
        vmin=0,
        vmax=1,
    )
    germany.boundary.plot(ax=ax, color="black", linewidth=1.2)
    for _, row in footprints.iterrows():
        gpd.GeoSeries([row.geometry], crs=footprints.crs).boundary.plot(
            ax=ax,
            color=row.color,
            linewidth=2.0,
            label=row.tile,
        )

    ax.set_title(f"GLASS FAPAR MODIS 250 m mosaic, {YEAR} DOY {DOY:03d}")
    ax.set_xlabel("MODIS sinusoidal x")
    ax.set_ylabel("MODIS sinusoidal y")
    ax.set_aspect("equal")
    ax.legend(title="MODIS tile", loc="upper right")

    colorbar = fig.colorbar(image, ax=ax, fraction=0.035, pad=0.02)
    colorbar.set_label("FAPAR")

    fig.tight_layout()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_PNG, dpi=200)
    plt.show()


def main() -> None:
    tif_files = find_tif_files()
    print("Input files:")
    for path in tif_files:
        print(f"  {path}")

    srcs = [rasterio.open(path) for path in tif_files]

    try:
        crs_set = {src.crs.to_string() if src.crs else None for src in srcs}
        if len(crs_set) != 1 or None in crs_set:
            raise RuntimeError(f"Input rasters must all have the same valid CRS. Found: {crs_set}")

        mosaic, transform = merge(srcs)
        footprints = tile_footprints(srcs)
        plot_mosaic(mosaic, transform, srcs[0].crs, footprints)
    finally:
        for src in srcs:
            src.close()

    print(f"Saved plot to: {OUTPUT_PNG}")


if __name__ == "__main__":
    main()
