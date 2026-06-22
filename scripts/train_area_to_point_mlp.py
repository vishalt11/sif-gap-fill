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
import torch
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import train_test_split
from torch import nn


GROUP_COL = "sif_extract_id"
TARGET_COL = "Daily_SIF_740nm"
WEIGHT_COL = "pixel_weight_equal"

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
]


@dataclass
class TrainConfig:
    data_dir: str
    output_dir: str
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


class PolygonBatcher:
    def __init__(
        self,
        data: pd.DataFrame,
        polygon_ids: np.ndarray,
        feature_columns: list[str],
        device: torch.device,
    ):
        self.device = device
        self.feature_columns = feature_columns

        subset = data[data[GROUP_COL].isin(polygon_ids)].copy()
        subset = subset.sort_values([GROUP_COL, "pixel_index_in_polygon"], kind="stable")
        subset = subset.reset_index(drop=True)

        self.group_ids = subset[GROUP_COL].drop_duplicates().to_numpy()
        group_to_index = {group_id: i for i, group_id in enumerate(self.group_ids)}

        self.x = subset[feature_columns].to_numpy(dtype=np.float32, copy=True)
        self.weights = subset[WEIGHT_COL].to_numpy(dtype=np.float32, copy=True)
        self.row_group_index = subset[GROUP_COL].map(group_to_index).to_numpy(dtype=np.int64)

        grouped_indices = subset.groupby(GROUP_COL, sort=False).indices
        self.indices_by_group = [
            np.asarray(grouped_indices[group_id], dtype=np.int64)
            for group_id in self.group_ids
        ]

        target_values = subset.groupby(GROUP_COL, sort=False)[TARGET_COL].first()
        self.y = target_values.loc[self.group_ids].to_numpy(dtype=np.float32)

        available_meta = [col for col in META_COLUMNS if col in subset.columns]
        self.meta = (
            subset.groupby(GROUP_COL, sort=False)[available_meta]
            .first()
            .reindex(self.group_ids)
            .reset_index()
        )

    def __len__(self) -> int:
        return len(self.group_ids)

    def iter_batches(self, batch_polygons: int, shuffle: bool) -> Iterable[tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]]:
        group_order = np.arange(len(self.group_ids))

        if shuffle:
            np.random.shuffle(group_order)

        for start in range(0, len(group_order), batch_polygons):
            batch_group_indices = group_order[start : start + batch_polygons]
            row_indices = np.concatenate([self.indices_by_group[i] for i in batch_group_indices])

            global_to_local = {global_i: local_i for local_i, global_i in enumerate(batch_group_indices)}
            local_group_index = np.fromiter(
                (global_to_local[i] for i in self.row_group_index[row_indices]),
                dtype=np.int64,
                count=len(row_indices),
            )

            x = torch.as_tensor(self.x[row_indices], dtype=torch.float32, device=self.device)
            weights = torch.as_tensor(self.weights[row_indices], dtype=torch.float32, device=self.device)
            local_group_index_t = torch.as_tensor(local_group_index, dtype=torch.long, device=self.device)
            y = torch.as_tensor(self.y[batch_group_indices], dtype=torch.float32, device=self.device)

            yield x, weights, local_group_index_t, y


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Train a 20 m area-to-point MLP with OCO polygon-aggregated loss."
    )
    parser.add_argument(
        "--data-dir",
        default="data/area_to_point_nonveg_masked/ba_sif_32UQV_wasp_area_to_point_20m/parquet",
        help="Folder containing pixel_table_manifest and pixel parquet files.",
    )
    parser.add_argument(
        "--output-dir",
        default="data/area_to_point_models/ba_sif_32UQV_mlp_area_to_point_20m",
        help="Folder where model outputs will be written.",
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
        help="Optional debug limit on number of monthly pixel parquet files to load.",
    )
    parser.add_argument(
        "--max-polygons",
        type=int,
        default=None,
        help="Optional debug limit on number of SIF polygons to keep after loading.",
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


def load_manifest(data_dir: Path) -> pd.DataFrame:
    parquet_manifest = data_dir / "pixel_table_manifest.parquet"
    csv_manifest = data_dir / "pixel_table_manifest.csv"

    if parquet_manifest.exists():
        return pd.read_parquet(parquet_manifest)

    if csv_manifest.exists():
        return pd.read_csv(csv_manifest)

    raise FileNotFoundError(f"No pixel table manifest found in {data_dir}")


def read_parquet_selected(path: Path, columns: list[str]) -> pd.DataFrame:
    try:
        return pd.read_parquet(path, columns=columns)
    except Exception:
        data = pd.read_parquet(path)
        available = [col for col in columns if col in data.columns]
        return data[available]


def load_pixel_tables(data_dir: Path, max_pixel_files: int | None) -> pd.DataFrame:
    manifest = load_manifest(data_dir)

    if "pixel_parquet" not in manifest.columns:
        raise ValueError("Manifest must contain a 'pixel_parquet' column.")

    pixel_paths = [resolve_path(path, Path.cwd()) for path in manifest["pixel_parquet"].dropna()]

    if max_pixel_files is not None:
        pixel_paths = pixel_paths[:max_pixel_files]

    if len(pixel_paths) == 0:
        raise ValueError("No pixel parquet files listed in manifest.")

    read_columns = list(
        dict.fromkeys(
            [
                GROUP_COL,
                TARGET_COL,
                WEIGHT_COL,
                "pixel_index_in_polygon",
                "pixel_crop_code",
                *BASE_FEATURE_COLUMNS,
                *META_COLUMNS,
            ]
        )
    )

    frames = []
    for path in pixel_paths:
        if not path.exists():
            raise FileNotFoundError(f"Pixel parquet file not found: {path}")
        print(f"Loading {path}")
        frames.append(read_parquet_selected(path, read_columns))

    data = pd.concat(frames, ignore_index=True)
    print(f"Loaded {len(data):,} pixel rows from {len(pixel_paths):,} parquet files.")
    return data


def add_time_features(data: pd.DataFrame) -> pd.DataFrame:
    if "sif_month" not in data.columns:
        if "Delta_Date" not in data.columns:
            raise ValueError("Need either 'sif_month' or 'Delta_Date' to create time features.")
        data["sif_month"] = pd.to_datetime(data["Delta_Date"]).dt.month

    month_angle = 2 * math.pi * data["sif_month"].astype(float) / 12.0
    data["sif_month_sin"] = np.sin(month_angle)
    data["sif_month_cos"] = np.cos(month_angle)
    return data


def add_crop_dummies(data: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    if "pixel_crop_code" not in data.columns:
        return data, []

    crop_code = data["pixel_crop_code"].fillna(-1).astype(int)
    dummy_columns = []

    for code in CROP_CODES:
        col = f"pixel_crop_code_{code}"
        data[col] = (crop_code == code).astype(np.float32)
        dummy_columns.append(col)

    data["pixel_crop_code_missing"] = (crop_code == -1).astype(np.float32)
    dummy_columns.append("pixel_crop_code_missing")
    return data, dummy_columns


def prepare_features(data: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    data = data.dropna(subset=[GROUP_COL, TARGET_COL]).copy()
    data = add_time_features(data)
    data, crop_dummy_columns = add_crop_dummies(data)

    for col in ["pixel_is_winter_wheat", "pixel_is_crop"]:
        if col in data.columns:
            data[col] = data[col].fillna(False).astype(np.float32)

    feature_columns = [
        col for col in BASE_FEATURE_COLUMNS
        if col in data.columns and col not in {"pixel_is_winter_wheat", "pixel_is_crop"}
    ]

    for col in ["pixel_is_winter_wheat", "pixel_is_crop"]:
        if col in data.columns:
            feature_columns.append(col)

    feature_columns.extend(crop_dummy_columns)
    feature_columns.extend(["sif_month_sin", "sif_month_cos"])

    if len(feature_columns) == 0:
        raise ValueError("No feature columns found.")

    data = data.dropna(subset=[GROUP_COL, TARGET_COL])
    data = data[data[feature_columns].notna().any(axis=1)].copy()

    print(f"Using {len(feature_columns)} feature columns.")
    print(f"Kept {len(data):,} pixel rows after target/feature filtering.")
    return data, feature_columns


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
    )

    val_fraction_of_train_val = val_fraction / (1 - test_fraction)
    train_ids, val_ids = train_test_split(
        train_val_ids,
        test_size=val_fraction_of_train_val,
        random_state=seed,
    )

    return np.asarray(train_ids), np.asarray(val_ids), np.asarray(test_ids)


def fit_preprocessing(
    data: pd.DataFrame,
    feature_columns: list[str],
    train_ids: np.ndarray,
) -> tuple[pd.DataFrame, dict[str, float], dict[str, float], dict[str, float], list[str]]:
    train_mask = data[GROUP_COL].isin(train_ids)

    medians: dict[str, float] = {}
    means: dict[str, float] = {}
    stds: dict[str, float] = {}
    kept_features: list[str] = []

    for col in feature_columns:
        train_values = pd.to_numeric(data.loc[train_mask, col], errors="coerce")
        median = float(train_values.median(skipna=True))

        if not np.isfinite(median):
            median = 0.0

        data[col] = pd.to_numeric(data[col], errors="coerce").fillna(median)
        train_values = data.loc[train_mask, col].astype(float)

        mean = float(train_values.mean())
        std = float(train_values.std(ddof=0))

        if not np.isfinite(std) or std == 0:
            std = 1.0

        medians[col] = median
        means[col] = mean
        stds[col] = std
        kept_features.append(col)

    for col in kept_features:
        data[col] = ((data[col].astype(float) - means[col]) / stds[col]).astype(np.float32)

    return data, medians, means, stds, kept_features


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
    batcher: PolygonBatcher,
    optimizer: torch.optim.Optimizer,
    loss_fn: nn.Module,
    batch_polygons: int,
) -> float:
    model.train()
    losses = []

    for x, weights, local_group_index, y in batcher.iter_batches(batch_polygons, shuffle=True):
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
        losses.append(float(loss.detach().cpu()))

    return float(np.mean(losses))


@torch.no_grad()
def predict_polygons(
    model: nn.Module,
    batcher: PolygonBatcher,
    batch_polygons: int,
) -> pd.DataFrame:
    model.eval()
    predictions = []
    observed = []
    group_ids = []

    for x, weights, local_group_index, y in batcher.iter_batches(batch_polygons, shuffle=False):
        pixel_pred = model(x)
        polygon_pred = aggregate_polygon_predictions(
            pixel_pred,
            weights,
            local_group_index,
            n_polygons=len(y),
        )

        start = len(group_ids)
        end = start + len(y)
        group_ids.extend(batcher.group_ids[start:end])
        observed.extend(y.detach().cpu().numpy())
        predictions.extend(polygon_pred.detach().cpu().numpy())

    pred_df = pd.DataFrame(
        {
            GROUP_COL: group_ids,
            "observed_sif": np.asarray(observed, dtype=np.float32),
            "predicted_sif": np.asarray(predictions, dtype=np.float32),
        }
    )

    return batcher.meta.merge(pred_df, on=GROUP_COL, how="right")


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

    data = load_pixel_tables(data_dir, config.max_pixel_files)
    data, feature_columns = prepare_features(data)

    polygon_ids = data[GROUP_COL].drop_duplicates().to_numpy()
    train_ids, val_ids, test_ids = split_polygon_ids(
        polygon_ids,
        val_fraction=config.val_fraction,
        test_fraction=config.test_fraction,
        seed=config.seed,
        max_polygons=config.max_polygons,
    )

    selected_ids = np.concatenate([train_ids, val_ids, test_ids])
    data = data[data[GROUP_COL].isin(selected_ids)].copy()

    data, medians, means, stds, feature_columns = fit_preprocessing(
        data=data,
        feature_columns=feature_columns,
        train_ids=train_ids,
    )

    save_json(output_dir / "feature_columns.json", feature_columns)
    save_json(output_dir / "feature_medians.json", medians)
    save_json(output_dir / "feature_means.json", means)
    save_json(output_dir / "feature_stds.json", stds)
    save_json(
        output_dir / "split_ids.json",
        {
            "train": [int(x) if isinstance(x, (np.integer, int)) else str(x) for x in train_ids],
            "validation": [int(x) if isinstance(x, (np.integer, int)) else str(x) for x in val_ids],
            "test": [int(x) if isinstance(x, (np.integer, int)) else str(x) for x in test_ids],
        },
    )

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")

    train_batcher = PolygonBatcher(data, train_ids, feature_columns, device)
    val_batcher = PolygonBatcher(data, val_ids, feature_columns, device)
    test_batcher = PolygonBatcher(data, test_ids, feature_columns, device)

    model = PixelToSifMLP(
        n_features=len(feature_columns),
        hidden_layers=config.hidden_layers,
        dropout=config.dropout,
    ).to(device)

    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=config.learning_rate,
        weight_decay=config.weight_decay,
    )

    loss_fn: nn.Module
    if config.loss == "huber":
        loss_fn = nn.HuberLoss(delta=0.1)
    else:
        loss_fn = nn.MSELoss()

    history = []
    best_val_rmse = math.inf
    best_epoch = 0
    best_model_path = output_dir / "best_model.pt"

    for epoch in range(1, config.epochs + 1):
        train_loss = train_one_epoch(
            model,
            train_batcher,
            optimizer,
            loss_fn,
            batch_polygons=config.batch_polygons,
        )

        val_predictions = predict_polygons(model, val_batcher, config.batch_polygons)
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
                    "n_features": len(feature_columns),
                    "feature_columns": feature_columns,
                    "hidden_layers": config.hidden_layers,
                    "dropout": config.dropout,
                    "feature_medians": medians,
                    "feature_means": means,
                    "feature_stds": stds,
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

    for split_name, batcher in [
        ("train", train_batcher),
        ("validation", val_batcher),
        ("test", test_batcher),
    ]:
        predictions = predict_polygons(model, batcher, config.batch_polygons)
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
