from __future__ import annotations

import argparse
import json
import math
import random
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import torch
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import train_test_split
from torch import nn


GROUP_COL = "polygon_uid"
TARGET_COL = "Daily_SIF_740nm"
WEIGHT_COL = "pixel_weight_equal"

#DEFAULT_TILES = ["32UNA", "32UNC", "32UPC"]
DEFAULT_TILES = ["32UNC", "32UPC"]
BASE_FEATURE_COLUMNS = [
    "pixel_b2",
    "pixel_b3",
    "pixel_b4",
    "pixel_b8",
    "pixel_b5",
    "pixel_b6",
    "pixel_b7",
    "pixel_b8a",
    "pixel_b11",
    "pixel_b12",
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
class TrainConfig:
    data_dir: str
    output_dir: str
    tiles: list[str]
    remove_extreme_sif_outliers: bool
    outlier_quantile: float
    seed: int
    val_fraction: float
    test_fraction: float
    batch_polygons: int
    epochs: int
    patience: int
    learning_rate: float
    weight_decay: float
    dropout: float
    hidden_layers: list[int]
    loss: str
    max_pixel_files: int | None
    max_polygons: int | None


@dataclass
class InMemoryDataset:
    features: np.ndarray
    weights: np.ndarray
    group_bounds: dict[str, tuple[int, int]]
    polygon_index: pd.DataFrame


class PixelToSifMLP(nn.Module):
    def __init__(self, n_features: int, hidden_layers: Iterable[int], dropout: float):
        super().__init__()

        layers: list[nn.Module] = []
        in_features = n_features

        for hidden in hidden_layers:
            layers.append(nn.Linear(in_features, hidden))
            layers.append(nn.ReLU())
            if dropout > 0:
                layers.append(nn.Dropout(dropout))
            in_features = hidden

        layers.append(nn.Linear(in_features, 1))
        self.network = nn.Sequential(*layers)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.network(x).squeeze(-1)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Train a 20 m area-to-point MLP from in-memory pixel tables."
    )
    parser.add_argument(
        "--data-dir",
        default="data/area_to_point_nonveg_masked",
        help="Parent folder containing tile area-to-point parquet folders.",
    )
    parser.add_argument(
        "--output-dir",
        default="data/area_to_point_models/three_tiles_mlp_area_to_point_20m",
        help="Folder where model outputs will be written.",
    )
    parser.add_argument(
        "--tiles",
        nargs="+",
        default=DEFAULT_TILES,
        help="Tile folders to load from data-dir. Default: 32UNA 32UNC 32UPC.",
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
    parser.add_argument("--batch-polygons", type=int, default=32)
    parser.add_argument("--epochs", type=int, default=150)
    parser.add_argument("--patience", type=int, default=20)
    parser.add_argument("--learning-rate", type=float, default=1e-3)
    parser.add_argument("--weight-decay", type=float, default=1e-4)
    parser.add_argument("--dropout", type=float, default=0.1)
    parser.add_argument("--hidden-layers", type=int, nargs="+", default=[64, 32])
    parser.add_argument("--loss", choices=["huber", "mse"], default="huber")
    parser.add_argument(
        "--max-pixel-files",
        type=int,
        default=None,
        help="Optional debug limit on number of monthly pixel parquet files to use.",
    )
    parser.add_argument(
        "--max-polygons",
        type=int,
        default=None,
        help="Optional debug limit on number of SIF polygons to keep after filtering.",
    )
    return parser.parse_args()


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)

    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


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


def discover_pixel_files(data_dir: Path, tiles: list[str], max_pixel_files: int | None) -> list[PixelFile]:
    pixel_files: list[PixelFile] = []

    for tile in tiles:
        parquet_dir = data_dir / tile / "parquet"

        if not parquet_dir.exists():
            raise FileNotFoundError(f"Tile parquet folder not found: {parquet_dir}")

        manifest = load_manifest(parquet_dir)

        if "pixel_parquet" not in manifest.columns:
            raise ValueError(f"Manifest must contain a 'pixel_parquet' column: {parquet_dir}")

        tile_file_count = 0
        for raw_path in manifest["pixel_parquet"].dropna():
            if max_pixel_files is not None and len(pixel_files) >= max_pixel_files:
                break

            path = resolve_pixel_path(raw_path, parquet_dir)
            if not path.exists():
                raise FileNotFoundError(f"Pixel parquet file not found: {path}")

            pixel_files.append(PixelFile(path=path, tile_dataset=tile))
            tile_file_count += 1

        print(f"Discovered {tile_file_count:,} pixel parquet file(s) for {tile}.")

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
        raise ValueError("Pixel data must contain 'sif_extract_id'.")

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
    crop_dummy_columns = [f"pixel_crop_code_{code}" for code in CROP_CODES]
    crop_dummy_columns.append("pixel_crop_code_missing")
    return BASE_FEATURE_COLUMNS + crop_dummy_columns + ["sif_month_sin", "sif_month_cos"]


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


def load_pixel_data(pixel_files: list[PixelFile]) -> pd.DataFrame:
    frames = []

    for i, pixel_file in enumerate(pixel_files, start=1):
        print(f"Loading {i:,}/{len(pixel_files):,}: {pixel_file.path.name}")
        data = read_parquet_selected(pixel_file.path, pixel_columns_for_read())
        data = add_polygon_uid(data, pixel_file.tile_dataset)
        frames.append(data)

    combined = pd.concat(frames, ignore_index=True)
    print(f"Loaded {len(combined):,} pixel rows into memory.")
    return combined


def prepare_model_features(data: pd.DataFrame, features: list[str]) -> pd.DataFrame:
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


def build_polygon_index(
    data: pd.DataFrame,
    remove_extreme_sif_outliers: bool,
    outlier_quantile: float,
) -> pd.DataFrame:
    polygon_index = data.drop_duplicates(subset=[GROUP_COL]).copy()
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


def compute_feature_stats(data: pd.DataFrame, train_ids: np.ndarray, features: list[str]) -> tuple[dict[str, float], dict[str, float]]:
    train_id_set = set(map(str, train_ids))
    train_data = data[data[GROUP_COL].isin(train_id_set)]

    means: dict[str, float] = {}
    stds: dict[str, float] = {}

    print("Computing feature scaling statistics from training pixels.")

    for col in features:
        values = pd.to_numeric(train_data[col], errors="coerce").to_numpy(dtype=np.float64)
        finite_values = values[np.isfinite(values)]

        if len(finite_values) == 0:
            means[col] = 0.0
            stds[col] = 1.0
            continue

        mean = float(finite_values.mean())
        std = float(finite_values.std())

        if not np.isfinite(std) or std == 0:
            std = 1.0

        means[col] = mean
        stds[col] = std

    return means, stds


def standardize_features(data: pd.DataFrame, features: list[str], means: dict[str, float], stds: dict[str, float]) -> pd.DataFrame:
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


def build_in_memory_dataset(data: pd.DataFrame, features: list[str], polygon_index: pd.DataFrame) -> InMemoryDataset:
    sort_columns = [GROUP_COL]
    if "pixel_index_in_polygon" in data.columns:
        sort_columns.append("pixel_index_in_polygon")

    data = data.sort_values(sort_columns, kind="stable").reset_index(drop=True)

    group_values = data[GROUP_COL].astype(str).to_numpy()
    change_positions = np.flatnonzero(group_values[1:] != group_values[:-1]) + 1
    starts = np.r_[0, change_positions]
    ends = np.r_[change_positions, len(group_values)]
    group_bounds = {
        str(group_values[start]): (int(start), int(end))
        for start, end in zip(starts, ends)
    }

    feature_array = data[features].to_numpy(dtype=np.float32, copy=True)
    weight_array = data[WEIGHT_COL].to_numpy(dtype=np.float32, copy=True)

    return InMemoryDataset(
        features=feature_array,
        weights=weight_array,
        group_bounds=group_bounds,
        polygon_index=polygon_index.copy(),
    )


def iter_polygon_batches(
    dataset: InMemoryDataset,
    polygon_ids: np.ndarray,
    batch_polygons: int,
    shuffle: bool,
    device: torch.device,
):
    ids = [str(x) for x in polygon_ids if str(x) in dataset.group_bounds]

    if shuffle:
        random.shuffle(ids)

    target_by_group = dataset.polygon_index.set_index(GROUP_COL)[TARGET_COL]

    for start in range(0, len(ids), batch_polygons):
        batch_ids = ids[start : start + batch_polygons]
        if len(batch_ids) == 0:
            continue

        index_parts = []
        local_group_parts = []
        kept_group_ids = []

        for local_index, group_id in enumerate(batch_ids):
            group_start, group_end = dataset.group_bounds[group_id]
            index_parts.append(np.arange(group_start, group_end, dtype=np.int64))
            local_group_parts.append(np.full(group_end - group_start, local_index, dtype=np.int64))
            kept_group_ids.append(group_id)

        pixel_indices = np.concatenate(index_parts)
        local_group_index_np = np.concatenate(local_group_parts)

        x = torch.as_tensor(dataset.features[pixel_indices], dtype=torch.float32, device=device)
        weights = torch.as_tensor(dataset.weights[pixel_indices], dtype=torch.float32, device=device)
        local_group_index = torch.as_tensor(local_group_index_np, dtype=torch.long, device=device)
        y = torch.as_tensor(
            target_by_group.loc[kept_group_ids].to_numpy(dtype=np.float32),
            dtype=torch.float32,
            device=device,
        )

        yield x, weights, local_group_index, y, np.asarray(kept_group_ids, dtype=object)


def aggregate_polygon_predictions(
    pixel_predictions: torch.Tensor,
    weights: torch.Tensor,
    local_group_index: torch.Tensor,
    n_polygons: int,
) -> torch.Tensor:
    weighted_sum = torch.zeros(n_polygons, dtype=torch.float32, device=pixel_predictions.device)
    weight_sum = torch.zeros(n_polygons, dtype=torch.float32, device=pixel_predictions.device)

    weighted_sum.scatter_add_(0, local_group_index, pixel_predictions * weights)
    weight_sum.scatter_add_(0, local_group_index, weights)

    return weighted_sum / weight_sum.clamp_min(1e-8)


def train_one_epoch(
    model: nn.Module,
    dataset: InMemoryDataset,
    train_ids: np.ndarray,
    optimizer: torch.optim.Optimizer,
    loss_fn: nn.Module,
    batch_polygons: int,
    device: torch.device,
) -> float:
    model.train()
    total_loss = 0.0
    total_polygons = 0

    for x, weights, local_group_index, y, _ in iter_polygon_batches(
        dataset=dataset,
        polygon_ids=train_ids,
        batch_polygons=batch_polygons,
        shuffle=True,
        device=device,
    ):
        optimizer.zero_grad(set_to_none=True)
        pixel_pred = model(x)
        polygon_pred = aggregate_polygon_predictions(
            pixel_pred,
            weights,
            local_group_index,
            n_polygons=len(y),
        )
        loss = loss_fn(polygon_pred, y)
        loss.backward()
        optimizer.step()

        total_loss += float(loss.detach().cpu()) * len(y)
        total_polygons += len(y)

    if total_polygons == 0:
        raise RuntimeError("No training batches were produced.")

    return total_loss / total_polygons


@torch.no_grad()
def predict_polygons(
    model: nn.Module,
    dataset: InMemoryDataset,
    polygon_ids: np.ndarray,
    batch_polygons: int,
    device: torch.device,
) -> pd.DataFrame:
    model.eval()
    rows = []

    for x, weights, local_group_index, y, group_ids in iter_polygon_batches(
        dataset=dataset,
        polygon_ids=polygon_ids,
        batch_polygons=batch_polygons,
        shuffle=False,
        device=device,
    ):
        pixel_pred = model(x)
        polygon_pred = aggregate_polygon_predictions(
            pixel_pred,
            weights,
            local_group_index,
            n_polygons=len(y),
        )

        rows.append(
            pd.DataFrame(
                {
                    GROUP_COL: group_ids,
                    "observed_sif": y.detach().cpu().numpy(),
                    "predicted_sif": polygon_pred.detach().cpu().numpy(),
                }
            )
        )

    if len(rows) == 0:
        raise RuntimeError("No prediction batches were produced.")

    predictions = pd.concat(rows, ignore_index=True)

    available_meta = [col for col in META_COLUMNS if col in dataset.polygon_index.columns and col != GROUP_COL]
    meta = dataset.polygon_index[[GROUP_COL, *available_meta]].drop_duplicates(subset=[GROUP_COL])
    predictions = meta.merge(predictions, on=GROUP_COL, how="right")
    return predictions


def compute_metrics(predictions: pd.DataFrame, split: str) -> dict[str, float | str | int]:
    y_true = predictions["observed_sif"].to_numpy(dtype=float)
    y_pred = predictions["predicted_sif"].to_numpy(dtype=float)
    residual = y_pred - y_true
    mse = float(mean_squared_error(y_true, y_pred))

    if len(y_true) > 1 and np.std(y_true) > 0 and np.std(y_pred) > 0:
        pearson_r = float(np.corrcoef(y_true, y_pred)[0, 1])
    else:
        pearson_r = np.nan

    return {
        "split": split,
        "n_polygons": int(len(predictions)),
        "rmse": float(math.sqrt(mse)),
        "mae": float(mean_absolute_error(y_true, y_pred)),
        "bias": float(np.mean(residual)),
        "r2": float(r2_score(y_true, y_pred)) if len(y_true) > 1 else np.nan,
        "pearson_r": pearson_r,
    }


def save_json(path: Path, value: object) -> None:
    with path.open("w", encoding="utf-8") as f:
        json.dump(value, f, indent=2)


def main() -> None:
    args = parse_args()
    config = TrainConfig(
        data_dir=args.data_dir,
        output_dir=args.output_dir,
        tiles=args.tiles,
        remove_extreme_sif_outliers=not args.keep_extreme_sif_outliers,
        outlier_quantile=args.outlier_quantile,
        seed=args.seed,
        val_fraction=args.val_fraction,
        test_fraction=args.test_fraction,
        batch_polygons=args.batch_polygons,
        epochs=args.epochs,
        patience=args.patience,
        learning_rate=args.learning_rate,
        weight_decay=args.weight_decay,
        dropout=args.dropout,
        hidden_layers=args.hidden_layers,
        loss=args.loss,
        max_pixel_files=args.max_pixel_files,
        max_polygons=args.max_polygons,
    )

    set_seed(config.seed)

    data_dir = Path(config.data_dir)
    output_dir = Path(config.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    save_json(output_dir / "config.json", asdict(config))

    pixel_files = discover_pixel_files(
        data_dir=data_dir,
        tiles=config.tiles,
        max_pixel_files=config.max_pixel_files,
    )
    features = feature_columns()

    data = load_pixel_data(pixel_files)
    data = prepare_model_features(data, features)

    polygon_index = build_polygon_index(
        data,
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
    data = data[data[GROUP_COL].isin(selected_ids)].copy()

    print(
        f"Split polygons: train={len(train_ids):,}, "
        f"validation={len(val_ids):,}, test={len(test_ids):,}"
    )

    means, stds = compute_feature_stats(data, train_ids, features)
    data = standardize_features(data, features, means, stds)
    data = make_equal_weights_if_needed(data)

    dataset = build_in_memory_dataset(
        data=data,
        features=features,
        polygon_index=polygon_index,
    )

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

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")
    print(f"Using {len(features)} feature columns.")

    model = PixelToSifMLP(
        n_features=len(features),
        hidden_layers=config.hidden_layers,
        dropout=config.dropout,
    ).to(device)

    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=config.learning_rate,
        weight_decay=config.weight_decay,
    )

    if config.loss == "huber":
        loss_fn: nn.Module = nn.HuberLoss(delta=0.1)
    else:
        loss_fn = nn.MSELoss()

    history = []
    best_val_rmse = math.inf
    best_epoch = 0
    best_model_path = output_dir / "best_model.pt"

    for epoch in range(1, config.epochs + 1):
        train_loss = train_one_epoch(
            model=model,
            dataset=dataset,
            train_ids=train_ids,
            optimizer=optimizer,
            loss_fn=loss_fn,
            batch_polygons=config.batch_polygons,
            device=device,
        )

        val_predictions = predict_polygons(
            model=model,
            dataset=dataset,
            polygon_ids=val_ids,
            batch_polygons=config.batch_polygons,
            device=device,
        )
        val_metrics = compute_metrics(val_predictions, "validation")
        val_rmse = float(val_metrics["rmse"])

        history_row = {
            "epoch": epoch,
            "train_loss": train_loss,
            **{f"val_{k}": v for k, v in val_metrics.items() if k != "split"},
        }
        history.append(history_row)

        print(
            f"epoch={epoch:03d} train_loss={train_loss:.6f} "
            f"val_rmse={val_rmse:.6f} val_r={val_metrics['pearson_r']:.3f}"
        )

        if val_rmse < best_val_rmse:
            best_val_rmse = val_rmse
            best_epoch = epoch
            torch.save(
                {
                    "model_state_dict": model.state_dict(),
                    "n_features": len(features),
                    "feature_columns": features,
                    "hidden_layers": config.hidden_layers,
                    "dropout": config.dropout,
                    "feature_means": means,
                    "feature_stds": stds,
                    "feature_impute_values": means,
                    "config": asdict(config),
                    "best_epoch": best_epoch,
                    "best_val_rmse": best_val_rmse,
                },
                best_model_path,
            )

        if epoch - best_epoch >= config.patience:
            print(f"Stopping early at epoch {epoch}; best epoch was {best_epoch}.")
            break

    pd.DataFrame(history).to_csv(output_dir / "training_history.csv", index=False)

    checkpoint = torch.load(best_model_path, map_location=device)
    model.load_state_dict(checkpoint["model_state_dict"])

    split_predictions = []
    split_metrics = []

    for split_name, ids in [
        ("train", train_ids),
        ("validation", val_ids),
        ("test", test_ids),
    ]:
        predictions = predict_polygons(
            model=model,
            dataset=dataset,
            polygon_ids=ids,
            batch_polygons=config.batch_polygons,
            device=device,
        )
        predictions["split"] = split_name
        predictions.to_csv(output_dir / f"polygon_predictions_{split_name}.csv", index=False)

        split_predictions.append(predictions)
        split_metrics.append(compute_metrics(predictions, split_name))

    pd.concat(split_predictions, ignore_index=True).to_csv(
        output_dir / "polygon_predictions_all_splits.csv",
        index=False,
    )
    pd.DataFrame(split_metrics).to_csv(output_dir / "metrics.csv", index=False)

    print(f"Best epoch: {best_epoch}")
    print(f"Saved model and outputs to {output_dir}")


if __name__ == "__main__":
    main()
