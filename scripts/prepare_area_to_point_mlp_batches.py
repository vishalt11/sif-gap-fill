from __future__ import annotations

import argparse
import json
import math
import random
import shutil
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import torch
from sklearn.model_selection import train_test_split


GROUP_COL = "polygon_uid"
TARGET_COL = "Daily_SIF_740nm"
WEIGHT_COL = "pixel_weight_equal"

BASE_FEATURE_COLUMNS = [
    "pixel_ndvi",
    "pixel_ndre",
    "pixel_ndre8a",
    "pixel_psri",
    "pixel_osavi",
    "pixel_ndwi",
    "pixel_nirv",
    "pixel_tcari",
    "pixel_ndmi",
    "pixel_msi",
    "pixel_ndmi_swir2",
    "pixel_msi_swir2",
    "pixel_nmdi",
    "pixel_is_winter_wheat",
    "pixel_is_crop",
]

CROP_CODES = [0, 11, 12, 13, 14, 21, 22, 23, 30, 40, 50, 60, 71, 81, 82, 83, 90, 100, 110, 111]

META_COLUMNS = [
    "sif_id",
    "sif_extract_id",
    "Delta_Date",
    "wasp_year_month",
    "sif_month",
    "sif_year",
    "mgrs_tile",
    "Latitude",
    "Longitude",
    "ww_pct",
    "crop_pixel_count",
    "Quality_Flag",
    "Metadata.MeasurementMode",
    "tile_dataset",
]


@dataclass
class PixelFile:
    path: Path
    tile_dataset: str


@dataclass
class PrepareConfig:
    data_dir: str
    output_dir: str
    tiles: list[str] | None
    remove_extreme_sif_outliers: bool
    outlier_quantile: float
    seed: int
    val_fraction: float
    test_fraction: float
    batch_polygons: int
    max_pixel_files: int | None
    max_polygons: int | None
    overwrite: bool


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare compact polygon-batch .pt files for area-to-point MLP training."
    )
    parser.add_argument(
        "--data-dir",
        default="data/area_to_point_nonveg_masked",
        help="Parent folder containing tile/parquet area-to-point outputs.",
    )
    parser.add_argument(
        "--output-dir",
        default="data/area_to_point_batches/all_tiles_20m",
        help="Folder where prepared .pt batch files and metadata will be written.",
    )
    parser.add_argument(
        "--tiles",
        nargs="+",
        default=None,
        help="Tile folders to include. If omitted, all tile folders with parquet manifests are used.",
    )
    parser.add_argument(
        "--keep-extreme-sif-outliers",
        action="store_true",
        help="Keep extreme SIF outliers. By default, symmetric quantile outliers are removed.",
    )
    parser.add_argument(
        "--outlier-quantile",
        type=float,
        default=0.001,
        help="Two-sided quantile threshold for extreme SIF outlier removal.",
    )
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--val-fraction", type=float, default=0.15)
    parser.add_argument("--test-fraction", type=float, default=0.15)
    parser.add_argument("--batch-polygons", type=int, default=64)
    parser.add_argument(
        "--max-pixel-files",
        type=int,
        default=None,
        help="Optional debug limit on number of monthly pixel parquet files.",
    )
    parser.add_argument(
        "--max-polygons",
        type=int,
        default=None,
        help="Optional debug limit on number of polygons after filtering.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Delete and recreate the output folder if it already exists.",
    )
    return parser.parse_args()


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)


def save_json(path: Path, value: object) -> None:
    with path.open("w", encoding="utf-8") as f:
        json.dump(value, f, indent=2)


def resolve_path(path: str | Path, base: Path) -> Path:
    path = Path(path)
    if path.is_absolute():
        return path
    return base / path


def load_manifest(parquet_dir: Path) -> pd.DataFrame:
    parquet_manifest = parquet_dir / "pixel_table_manifest.parquet"
    csv_manifest = parquet_dir / "pixel_table_manifest.csv"

    if parquet_manifest.exists():
        return pd.read_parquet(parquet_manifest)

    if csv_manifest.exists():
        return pd.read_csv(csv_manifest)

    raise FileNotFoundError(f"No pixel table manifest found in {parquet_dir}")


def resolve_pixel_path(path: str | Path, parquet_dir: Path) -> Path:
    path = resolve_path(path, Path.cwd())

    if path.exists():
        return path

    renamed_folder_path = parquet_dir / path.name
    if renamed_folder_path.exists():
        return renamed_folder_path

    return path


def discover_tile_dirs(data_dir: Path, tiles: list[str] | None) -> list[Path]:
    if tiles is not None:
        tile_dirs = [data_dir / tile / "parquet" for tile in tiles]
    else:
        tile_dirs = sorted(
            path for path in data_dir.glob("*/parquet")
            if (path / "pixel_table_manifest.parquet").exists()
            or (path / "pixel_table_manifest.csv").exists()
        )

    missing = [path for path in tile_dirs if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Missing tile parquet folder(s): {missing}")

    if len(tile_dirs) == 0:
        raise FileNotFoundError(f"No tile parquet folders found under {data_dir}")

    return tile_dirs


def tile_dataset_name(parquet_dir: Path) -> str:
    if parquet_dir.name == "parquet":
        return parquet_dir.parent.name
    return parquet_dir.name


def discover_pixel_files(tile_dirs: list[Path], max_pixel_files: int | None) -> list[PixelFile]:
    pixel_files: list[PixelFile] = []

    for parquet_dir in tile_dirs:
        manifest = load_manifest(parquet_dir)

        if "pixel_parquet" not in manifest.columns:
            raise ValueError(f"Manifest must contain a 'pixel_parquet' column: {parquet_dir}")

        tile_name = tile_dataset_name(parquet_dir)
        tile_file_count = 0

        for raw_path in manifest["pixel_parquet"].dropna():
            if max_pixel_files is not None and len(pixel_files) >= max_pixel_files:
                break

            path = resolve_pixel_path(raw_path, parquet_dir)
            if not path.exists():
                raise FileNotFoundError(f"Pixel parquet file not found: {path}")

            pixel_files.append(PixelFile(path=path, tile_dataset=tile_name))
            tile_file_count += 1

        print(f"Discovered {tile_file_count:,} pixel parquet file(s) for {tile_name}.")

        if max_pixel_files is not None and len(pixel_files) >= max_pixel_files:
            break

    if len(pixel_files) == 0:
        raise ValueError("No pixel parquet files found.")

    print(f"Discovered {len(pixel_files):,} pixel parquet file(s) total.")
    return pixel_files


def read_parquet_selected(path: Path, columns: list[str]) -> pd.DataFrame:
    available_columns = set(pq.read_schema(path).names)
    selected_columns = [col for col in columns if col in available_columns]

    if len(selected_columns) == 0:
        raise ValueError(f"None of the requested columns are present in {path}")

    return pd.read_parquet(path, columns=selected_columns)


def add_polygon_uid(data: pd.DataFrame, tile_dataset: str) -> pd.DataFrame:
    if "sif_extract_id" not in data.columns:
        raise ValueError("Data must contain 'sif_extract_id'.")

    if "mgrs_tile" in data.columns:
        tile_part = data["mgrs_tile"].fillna(tile_dataset).astype(str)
    else:
        tile_part = pd.Series(tile_dataset, index=data.index, dtype=str)

    data["tile_dataset"] = tile_dataset
    data[GROUP_COL] = tile_dataset + "__" + tile_part + "__" + data["sif_extract_id"].astype(str)
    return data


def add_time_features(data: pd.DataFrame) -> pd.DataFrame:
    if "sif_month" not in data.columns:
        if "Delta_Date" not in data.columns:
            raise ValueError("Need either 'sif_month' or 'Delta_Date' to create time features.")
        data["sif_month"] = pd.to_datetime(data["Delta_Date"]).dt.month

    month_angle = 2 * math.pi * data["sif_month"].astype(float) / 12.0
    data["sif_month_sin"] = np.sin(month_angle).astype(np.float32)
    data["sif_month_cos"] = np.cos(month_angle).astype(np.float32)
    return data


def add_crop_dummies(data: pd.DataFrame) -> pd.DataFrame:
    if "pixel_crop_code" not in data.columns:
        for code in CROP_CODES:
            data[f"pixel_crop_code_{code}"] = np.float32(0.0)
        data["pixel_crop_code_missing"] = np.float32(1.0)
        return data

    crop_code = pd.to_numeric(data["pixel_crop_code"], errors="coerce").fillna(-1).astype(int)

    for code in CROP_CODES:
        data[f"pixel_crop_code_{code}"] = (crop_code == code).astype(np.float32)

    data["pixel_crop_code_missing"] = (crop_code == -1).astype(np.float32)
    return data


def feature_columns() -> list[str]:
    return BASE_FEATURE_COLUMNS + ["sif_month_sin", "sif_month_cos"]


def print_feature_columns(features: list[str]) -> None:
    print(f"Using {len(features)} feature columns:")
    for i, col in enumerate(features, start=1):
        print(f"  {i:02d}. {col}")


def pixel_columns_for_read() -> list[str]:
    raw_feature_inputs = [
        col for col in BASE_FEATURE_COLUMNS
        if col not in {"pixel_is_winter_wheat", "pixel_is_crop"}
    ]

    return list(
        dict.fromkeys(
            [
                "sif_extract_id",
                TARGET_COL,
                WEIGHT_COL,
                "pixel_index_in_polygon",
                "pixel_crop_code",
                *raw_feature_inputs,
                "pixel_is_winter_wheat",
                "pixel_is_crop",
                *META_COLUMNS,
            ]
        )
    )


def load_polygon_targets(tile_dirs: list[Path]) -> pd.DataFrame:
    frames = []

    for parquet_dir in tile_dirs:
        tile_name = tile_dataset_name(parquet_dir)
        target_path = parquet_dir / "polygon_targets.parquet"

        if not target_path.exists():
            raise FileNotFoundError(f"Polygon target file not found: {target_path}")

        data = pd.read_parquet(target_path)
        data = add_polygon_uid(data, tile_name)
        frames.append(data)

    targets = pd.concat(frames, ignore_index=True)
    print(f"Loaded {len(targets):,} polygon target rows.")
    return targets


def filter_polygon_targets(
    targets: pd.DataFrame,
    remove_extreme_sif_outliers: bool,
    outlier_quantile: float,
) -> pd.DataFrame:
    polygon_index = targets.drop_duplicates(subset=[GROUP_COL]).copy()
    before_polygons = len(polygon_index)

    polygon_index = polygon_index.dropna(subset=[GROUP_COL, TARGET_COL]).copy()

    if "Quality_Flag" in polygon_index.columns:
        polygon_index["Quality_Flag"] = pd.to_numeric(polygon_index["Quality_Flag"], errors="coerce")
        polygon_index = polygon_index[polygon_index["Quality_Flag"].isin([0, 1])].copy()
    else:
        print("Warning: Quality_Flag column not found; no quality filter applied.")

    if "ww_pct" in polygon_index.columns:
        polygon_index["ww_pct"] = polygon_index["ww_pct"].fillna(0)

    if remove_extreme_sif_outliers:
        lower = float(polygon_index[TARGET_COL].quantile(outlier_quantile))
        upper = float(polygon_index[TARGET_COL].quantile(1 - outlier_quantile))
        polygon_index = polygon_index[
            polygon_index[TARGET_COL].between(lower, upper, inclusive="both")
        ].copy()
        print(
            f"Removed extreme SIF outliers using {outlier_quantile:.4f}/"
            f"{1 - outlier_quantile:.4f} polygon quantiles: [{lower:.6f}, {upper:.6f}]"
        )

    polygon_index = polygon_index.drop_duplicates(subset=[GROUP_COL]).reset_index(drop=True)
    print(f"Polygon filters kept {len(polygon_index):,}/{before_polygons:,} polygons.")
    return polygon_index


def split_polygon_ids(
    polygon_ids: np.ndarray,
    val_fraction: float,
    test_fraction: float,
    seed: int,
    max_polygons: int | None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    polygon_ids = np.asarray(polygon_ids)

    if max_polygons is not None and max_polygons < len(polygon_ids):
        rng = np.random.default_rng(seed)
        polygon_ids = rng.choice(polygon_ids, size=max_polygons, replace=False)

    train_val_ids, test_ids = train_test_split(
        polygon_ids,
        test_size=test_fraction,
        random_state=seed,
        shuffle=True,
    )

    val_fraction_of_train_val = val_fraction / (1 - test_fraction)
    train_ids, val_ids = train_test_split(
        train_val_ids,
        test_size=val_fraction_of_train_val,
        random_state=seed,
        shuffle=True,
    )

    return np.asarray(train_ids), np.asarray(val_ids), np.asarray(test_ids)


def prepare_pixel_features(data: pd.DataFrame, features: list[str]) -> pd.DataFrame:
    data = add_time_features(data)
    data = add_crop_dummies(data)

    if "ww_pct" in data.columns:
        data["ww_pct"] = data["ww_pct"].fillna(0)

    for col in ["pixel_is_winter_wheat", "pixel_is_crop"]:
        if col in data.columns:
            data[col] = data[col].fillna(False).astype(np.float32)

    for col in features:
        if col not in data.columns:
            data[col] = np.nan

    return data


def load_pixel_file_for_polygons(
    pixel_file: PixelFile,
    polygon_id_set: set[str],
    features: list[str],
) -> pd.DataFrame:
    data = read_parquet_selected(pixel_file.path, pixel_columns_for_read())
    data = add_polygon_uid(data, pixel_file.tile_dataset)
    data = data[data[GROUP_COL].isin(polygon_id_set)].copy()

    if len(data) == 0:
        return data

    return prepare_pixel_features(data, features)


def compute_feature_stats(
    pixel_files: list[PixelFile],
    train_ids: np.ndarray,
    features: list[str],
) -> tuple[dict[str, float], dict[str, float]]:
    train_id_set = set(map(str, train_ids))
    sums = {col: 0.0 for col in features}
    sums_sq = {col: 0.0 for col in features}
    counts = {col: 0 for col in features}

    print("Computing feature scaling statistics from training pixels.")

    for i, pixel_file in enumerate(pixel_files, start=1):
        data = load_pixel_file_for_polygons(pixel_file, train_id_set, features)

        if len(data) == 0:
            print(f"  stats pass {i:,}/{len(pixel_files):,}: {pixel_file.path.name} (no train polygons)")
            continue

        for col in features:
            values = pd.to_numeric(data[col], errors="coerce").to_numpy(dtype=np.float64)
            mask = np.isfinite(values)
            if not np.any(mask):
                continue
            finite_values = values[mask]
            sums[col] += float(finite_values.sum())
            sums_sq[col] += float((finite_values * finite_values).sum())
            counts[col] += int(mask.sum())

        print(f"  stats pass {i:,}/{len(pixel_files):,}: {pixel_file.path.name}")

    means: dict[str, float] = {}
    stds: dict[str, float] = {}

    for col in features:
        if counts[col] == 0:
            means[col] = 0.0
            stds[col] = 1.0
            continue

        mean = sums[col] / counts[col]
        variance = max((sums_sq[col] / counts[col]) - mean * mean, 0.0)
        std = math.sqrt(variance)

        if not np.isfinite(std) or std == 0:
            std = 1.0

        means[col] = float(mean)
        stds[col] = float(std)

    return means, stds


def standardize_features(
    data: pd.DataFrame,
    features: list[str],
    means: dict[str, float],
    stds: dict[str, float],
) -> pd.DataFrame:
    for col in features:
        values = pd.to_numeric(data[col], errors="coerce").astype(float)
        values = values.fillna(means[col])
        data[col] = ((values - means[col]) / stds[col]).astype(np.float32)
    return data


def make_equal_weights_if_needed(data: pd.DataFrame) -> pd.DataFrame:
    if WEIGHT_COL in data.columns:
        data[WEIGHT_COL] = pd.to_numeric(data[WEIGHT_COL], errors="coerce")
    else:
        data[WEIGHT_COL] = np.nan

    missing_weight = data[WEIGHT_COL].isna()
    if missing_weight.any():
        group_size = data.groupby(GROUP_COL)[GROUP_COL].transform("size")
        data.loc[missing_weight, WEIGHT_COL] = 1.0 / group_size.loc[missing_weight]

    data[WEIGHT_COL] = data[WEIGHT_COL].astype(np.float32)
    return data


def serializable_metadata(meta: pd.DataFrame) -> dict[str, list[object]]:
    result: dict[str, list[object]] = {}

    for col in meta.columns:
        values = meta[col]
        if pd.api.types.is_datetime64_any_dtype(values):
            result[col] = values.astype(str).tolist()
        else:
            result[col] = values.where(pd.notna(values), None).tolist()

    return result


def write_batch(
    split_name: str,
    batch_number: int,
    polygon_frames: list[tuple[str, pd.DataFrame]],
    output_dir: Path,
    features: list[str],
    target_by_group: pd.Series,
    polygon_index: pd.DataFrame,
) -> dict[str, object]:
    group_ids = [group_id for group_id, _ in polygon_frames]
    frame_parts = []
    local_group_parts = []

    for local_index, (_, group_data) in enumerate(polygon_frames):
        frame_parts.append(group_data)
        local_group_parts.append(np.full(len(group_data), local_index, dtype=np.int64))

    batch = pd.concat(frame_parts, ignore_index=True)
    local_group_index = np.concatenate(local_group_parts)

    x = batch[features].to_numpy(dtype=np.float32, copy=True)
    weights = batch[WEIGHT_COL].to_numpy(dtype=np.float32, copy=True)
    y = target_by_group.loc[group_ids].to_numpy(dtype=np.float32)

    available_meta = [col for col in META_COLUMNS if col in polygon_index.columns and col != GROUP_COL]
    meta = polygon_index.set_index(GROUP_COL).loc[group_ids, available_meta].reset_index()

    split_dir = output_dir / split_name
    split_dir.mkdir(parents=True, exist_ok=True)
    batch_path = split_dir / f"batch_{batch_number:06d}.pt"

    torch.save(
        {
            "x": torch.from_numpy(x),
            "weights": torch.from_numpy(weights),
            "group_index": torch.from_numpy(local_group_index),
            "y": torch.from_numpy(y),
            "polygon_uid": group_ids,
            "feature_columns": features,
            "metadata": serializable_metadata(meta),
        },
        batch_path,
    )

    return {
        "split": split_name,
        "batch_number": batch_number,
        "batch_path": str(batch_path),
        "n_polygons": len(group_ids),
        "n_pixels": len(batch),
    }


def flush_ready_batches(
    split_name: str,
    buffers: dict[str, list[tuple[str, pd.DataFrame]]],
    batch_counters: dict[str, int],
    batch_manifest_rows: list[dict[str, object]],
    output_dir: Path,
    features: list[str],
    target_by_group: pd.Series,
    polygon_index: pd.DataFrame,
    batch_polygons: int,
    force: bool = False,
) -> None:
    while len(buffers[split_name]) >= batch_polygons or (force and len(buffers[split_name]) > 0):
        if force:
            batch_frames = buffers[split_name][:]
            buffers[split_name].clear()
        else:
            batch_frames = buffers[split_name][:batch_polygons]
            del buffers[split_name][:batch_polygons]

        batch_counters[split_name] += 1
        row = write_batch(
            split_name=split_name,
            batch_number=batch_counters[split_name],
            polygon_frames=batch_frames,
            output_dir=output_dir,
            features=features,
            target_by_group=target_by_group,
            polygon_index=polygon_index,
        )
        batch_manifest_rows.append(row)

        print(
            f"  wrote {split_name} batch {batch_counters[split_name]:,}: "
            f"{row['n_polygons']:,} polygons, {row['n_pixels']:,} pixels"
        )


def write_batches(
    pixel_files: list[PixelFile],
    polygon_index: pd.DataFrame,
    split_lookup: dict[str, str],
    features: list[str],
    means: dict[str, float],
    stds: dict[str, float],
    output_dir: Path,
    batch_polygons: int,
) -> pd.DataFrame:
    selected_ids = set(split_lookup)
    target_by_group = polygon_index.set_index(GROUP_COL)[TARGET_COL]

    buffers = {
        "train": [],
        "validation": [],
        "test": [],
    }
    batch_counters = {
        "train": 0,
        "validation": 0,
        "test": 0,
    }
    batch_manifest_rows: list[dict[str, object]] = []

    print("Writing prepared .pt polygon batches.")

    for i, pixel_file in enumerate(pixel_files, start=1):
        data = load_pixel_file_for_polygons(pixel_file, selected_ids, features)

        if len(data) == 0:
            print(f"  batch pass {i:,}/{len(pixel_files):,}: {pixel_file.path.name} (no selected polygons)")
            continue

        data = standardize_features(data, features, means, stds)
        data = make_equal_weights_if_needed(data)

        sort_columns = [GROUP_COL]
        if "pixel_index_in_polygon" in data.columns:
            sort_columns.append("pixel_index_in_polygon")
        data = data.sort_values(sort_columns, kind="stable").reset_index(drop=True)

        group_map = {
            str(group_id): group_data.copy()
            for group_id, group_data in data.groupby(GROUP_COL, sort=False)
        }
        group_ids = list(group_map)
        random.shuffle(group_ids)

        for group_id in group_ids:
            split_name = split_lookup.get(str(group_id))
            if split_name is None:
                continue

            group_data = group_map[group_id]
            buffers[split_name].append((str(group_id), group_data))

            flush_ready_batches(
                split_name=split_name,
                buffers=buffers,
                batch_counters=batch_counters,
                batch_manifest_rows=batch_manifest_rows,
                output_dir=output_dir,
                features=features,
                target_by_group=target_by_group,
                polygon_index=polygon_index,
                batch_polygons=batch_polygons,
                force=False,
            )

        print(f"  batch pass {i:,}/{len(pixel_files):,}: {pixel_file.path.name}")

    for split_name in ["train", "validation", "test"]:
        flush_ready_batches(
            split_name=split_name,
            buffers=buffers,
            batch_counters=batch_counters,
            batch_manifest_rows=batch_manifest_rows,
            output_dir=output_dir,
            features=features,
            target_by_group=target_by_group,
            polygon_index=polygon_index,
            batch_polygons=batch_polygons,
            force=True,
        )

    return pd.DataFrame(batch_manifest_rows)


def prepare_output_dir(output_dir: Path, overwrite: bool) -> None:
    if output_dir.exists() and any(output_dir.glob("*")):
        if not overwrite:
            raise FileExistsError(
                f"Output folder already exists and is not empty: {output_dir}. "
                "Use --overwrite to recreate it."
            )
        shutil.rmtree(output_dir)

    output_dir.mkdir(parents=True, exist_ok=True)


def main() -> None:
    args = parse_args()
    config = PrepareConfig(
        data_dir=args.data_dir,
        output_dir=args.output_dir,
        tiles=args.tiles,
        remove_extreme_sif_outliers=not args.keep_extreme_sif_outliers,
        outlier_quantile=args.outlier_quantile,
        seed=args.seed,
        val_fraction=args.val_fraction,
        test_fraction=args.test_fraction,
        batch_polygons=args.batch_polygons,
        max_pixel_files=args.max_pixel_files,
        max_polygons=args.max_polygons,
        overwrite=args.overwrite,
    )

    set_seed(config.seed)

    data_dir = Path(config.data_dir)
    output_dir = Path(config.output_dir)
    prepare_output_dir(output_dir, config.overwrite)

    tile_dirs = discover_tile_dirs(data_dir, config.tiles)
    pixel_files = discover_pixel_files(tile_dirs, config.max_pixel_files)
    features = feature_columns()
    print_feature_columns(features)

    polygon_targets = load_polygon_targets(tile_dirs)
    polygon_index = filter_polygon_targets(
        polygon_targets,
        remove_extreme_sif_outliers=config.remove_extreme_sif_outliers,
        outlier_quantile=config.outlier_quantile,
    )

    polygon_ids = polygon_index[GROUP_COL].drop_duplicates().to_numpy()
    train_ids, val_ids, test_ids = split_polygon_ids(
        polygon_ids,
        val_fraction=config.val_fraction,
        test_fraction=config.test_fraction,
        seed=config.seed,
        max_polygons=config.max_polygons,
    )

    selected_ids = set(map(str, np.concatenate([train_ids, val_ids, test_ids])))
    polygon_index = polygon_index[polygon_index[GROUP_COL].isin(selected_ids)].copy()

    print(
        f"Split polygons: train={len(train_ids):,}, "
        f"validation={len(val_ids):,}, test={len(test_ids):,}"
    )

    means, stds = compute_feature_stats(
        pixel_files=pixel_files,
        train_ids=train_ids,
        features=features,
    )

    split_lookup = {str(group_id): "train" for group_id in train_ids}
    split_lookup.update({str(group_id): "validation" for group_id in val_ids})
    split_lookup.update({str(group_id): "test" for group_id in test_ids})

    batch_manifest = write_batches(
        pixel_files=pixel_files,
        polygon_index=polygon_index,
        split_lookup=split_lookup,
        features=features,
        means=means,
        stds=stds,
        output_dir=output_dir,
        batch_polygons=config.batch_polygons,
    )

    save_json(output_dir / "config.json", asdict(config))
    save_json(output_dir / "feature_columns.json", features)
    save_json(output_dir / "feature_means.json", means)
    save_json(output_dir / "feature_stds.json", stds)
    save_json(output_dir / "feature_impute_values.json", means)
    save_json(
        output_dir / "split_ids.json",
        {
            "train": [str(x) for x in train_ids],
            "validation": [str(x) for x in val_ids],
            "test": [str(x) for x in test_ids],
        },
    )

    polygon_index.to_csv(output_dir / "polygon_index_filtered.csv", index=False)
    batch_manifest.to_csv(output_dir / "batch_manifest.csv", index=False)

    print(f"Prepared {len(batch_manifest):,} batch file(s).")
    print(f"Saved prepared batches to {output_dir}")


if __name__ == "__main__":
    main()
