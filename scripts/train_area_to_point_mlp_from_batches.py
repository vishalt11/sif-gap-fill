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
from torch import nn


GROUP_COL = "polygon_uid"


@dataclass
class TrainConfig:
    batch_dir: str
    output_dir: str
    seed: int
    epochs: int
    patience: int
    learning_rate: float
    weight_decay: float
    dropout: float
    hidden_layers: list[int]
    loss: str
    val_every: int
    max_train_batches: int | None


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
        description="Train a 20 m area-to-point MLP from prepared .pt polygon batches."
    )
    parser.add_argument(
        "--batch-dir",
        default="data/area_to_point_batches/all_tiles_20m",
        help="Folder created by prepare_area_to_point_mlp_batches.py.",
    )
    parser.add_argument(
        "--output-dir",
        default="data/area_to_point_models/all_tiles_mlp_area_to_point_20m_from_batches_17_features",
        help="Folder where model outputs will be written.",
    )
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--epochs", type=int, default=150)
    parser.add_argument("--patience", type=int, default=20)
    parser.add_argument("--learning-rate", type=float, default=1e-3)
    parser.add_argument("--weight-decay", type=float, default=1e-4)
    parser.add_argument("--dropout", type=float, default=0.1)
    parser.add_argument("--hidden-layers", type=int, nargs="+", default=[64, 32])
    parser.add_argument("--loss", choices=["huber", "mse"], default="huber")
    parser.add_argument(
        "--val-every",
        type=int,
        default=1,
        help="Run validation every N epochs, plus epoch 1.",
    )
    parser.add_argument(
        "--max-train-batches",
        type=int,
        default=100,
        help="Random number of train batches to use per epoch. Use 0 to use all train batches.",
    )
    return parser.parse_args()


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)

    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


def save_json(path: Path, value: object) -> None:
    with path.open("w", encoding="utf-8") as f:
        json.dump(value, f, indent=2)


def load_json(path: Path) -> object:
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def resolve_batch_path(path_value: str, batch_dir: Path) -> Path:
    path = Path(path_value)
    if path.is_absolute():
        return path

    cwd_path = Path.cwd() / path
    if cwd_path.exists():
        return cwd_path

    batch_dir_path = batch_dir / path.name
    if batch_dir_path.exists():
        return batch_dir_path

    split_path = batch_dir / path.parent.name / path.name
    if split_path.exists():
        return split_path

    return cwd_path


def load_batch_manifest(batch_dir: Path) -> pd.DataFrame:
    manifest_path = batch_dir / "batch_manifest.csv"

    if not manifest_path.exists():
        raise FileNotFoundError(f"Batch manifest not found: {manifest_path}")

    manifest = pd.read_csv(manifest_path)

    required = {"split", "batch_path", "n_polygons", "n_pixels"}
    missing = required.difference(manifest.columns)
    if missing:
        raise ValueError(f"Batch manifest is missing required columns: {sorted(missing)}")

    manifest["resolved_batch_path"] = manifest["batch_path"].apply(
        lambda value: str(resolve_batch_path(str(value), batch_dir))
    )

    missing_files = [
        path for path in manifest["resolved_batch_path"]
        if not Path(path).exists()
    ]
    if missing_files:
        raise FileNotFoundError(f"Some batch files listed in manifest are missing: {missing_files[:5]}")

    return manifest


def load_torch_file(path: Path, map_location: str | torch.device = "cpu") -> dict:
    try:
        return torch.load(path, map_location=map_location, weights_only=False)
    except TypeError:
        return torch.load(path, map_location=map_location)


def load_batch(path: Path, device: torch.device) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor, list[str], dict]:
    batch = load_torch_file(path, map_location="cpu")

    x = batch["x"].to(device=device, dtype=torch.float32, non_blocking=True)
    weights = batch["weights"].to(device=device, dtype=torch.float32, non_blocking=True)
    group_index = batch["group_index"].to(device=device, dtype=torch.long, non_blocking=True)
    y = batch["y"].to(device=device, dtype=torch.float32, non_blocking=True)
    polygon_uid = [str(x) for x in batch["polygon_uid"]]
    metadata = batch.get("metadata", {})

    return x, weights, group_index, y, polygon_uid, metadata


def print_feature_columns(features: list[str]) -> None:
    print(f"Using {len(features)} feature columns:")
    for i, col in enumerate(features, start=1):
        print(f"  {i:02d}. {col}")


def aggregate_polygon_predictions(
    pixel_predictions: torch.Tensor,
    weights: torch.Tensor,
    group_index: torch.Tensor,
    n_polygons: int,
) -> torch.Tensor:
    weighted_sum = torch.zeros(n_polygons, dtype=torch.float32, device=pixel_predictions.device)
    weight_sum = torch.zeros(n_polygons, dtype=torch.float32, device=pixel_predictions.device)

    weighted_sum.scatter_add_(0, group_index, pixel_predictions * weights)
    weight_sum.scatter_add_(0, group_index, weights)

    return weighted_sum / weight_sum.clamp_min(1e-8)


def train_one_epoch(
    model: nn.Module,
    batch_paths: list[Path],
    optimizer: torch.optim.Optimizer,
    loss_fn: nn.Module,
    device: torch.device,
    max_train_batches: int | None,
) -> float:
    model.train()
    total_loss = 0.0
    total_polygons = 0

    epoch_paths = list(batch_paths)
    random.shuffle(epoch_paths)

    if max_train_batches is not None and max_train_batches > 0:
        epoch_paths = epoch_paths[:max_train_batches]

    for batch_path in epoch_paths:
        x, weights, group_index, y, _, _ = load_batch(batch_path, device)

        optimizer.zero_grad(set_to_none=True)
        pixel_pred = model(x)
        polygon_pred = aggregate_polygon_predictions(
            pixel_pred,
            weights,
            group_index,
            n_polygons=len(y),
        )
        loss = loss_fn(polygon_pred, y)
        loss.backward()
        optimizer.step()

        total_loss += float(loss.detach().cpu()) * len(y)
        total_polygons += len(y)

        del x, weights, group_index, y, pixel_pred, polygon_pred, loss

    if total_polygons == 0:
        raise RuntimeError("No training batches were produced.")

    return total_loss / total_polygons


@torch.no_grad()
def predict_batches(
    model: nn.Module,
    batch_paths: list[Path],
    split: str,
    device: torch.device,
) -> pd.DataFrame:
    model.eval()
    rows = []

    for batch_path in batch_paths:
        x, weights, group_index, y, polygon_uid, metadata = load_batch(batch_path, device)
        pixel_pred = model(x)
        polygon_pred = aggregate_polygon_predictions(
            pixel_pred,
            weights,
            group_index,
            n_polygons=len(y),
        )

        pred_df = pd.DataFrame(
            {
                GROUP_COL: polygon_uid,
                "observed_sif": y.detach().cpu().numpy(),
                "predicted_sif": polygon_pred.detach().cpu().numpy(),
                "split": split,
            }
        )

        if isinstance(metadata, dict) and len(metadata) > 0:
            meta_df = pd.DataFrame(metadata)
            if GROUP_COL in meta_df.columns:
                pred_df = meta_df.merge(pred_df, on=GROUP_COL, how="right")

        rows.append(pred_df)
        del x, weights, group_index, y, pixel_pred, polygon_pred

    if len(rows) == 0:
        raise RuntimeError(f"No prediction batches were produced for split: {split}")

    return pd.concat(rows, ignore_index=True)


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


def paths_for_split(manifest: pd.DataFrame, split: str) -> list[Path]:
    rows = manifest[manifest["split"] == split].copy()
    return [Path(path) for path in rows["resolved_batch_path"].tolist()]


def main() -> None:
    args = parse_args()
    config = TrainConfig(
        batch_dir=args.batch_dir,
        output_dir=args.output_dir,
        seed=args.seed,
        epochs=args.epochs,
        patience=args.patience,
        learning_rate=args.learning_rate,
        weight_decay=args.weight_decay,
        dropout=args.dropout,
        hidden_layers=args.hidden_layers,
        loss=args.loss,
        val_every=args.val_every,
        max_train_batches=args.max_train_batches,
    )

    set_seed(config.seed)

    batch_dir = Path(config.batch_dir)
    output_dir = Path(config.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    feature_columns_path = batch_dir / "feature_columns.json"
    if not feature_columns_path.exists():
        raise FileNotFoundError(f"Feature column file not found: {feature_columns_path}")

    features = load_json(feature_columns_path)
    if not isinstance(features, list):
        raise ValueError(f"Expected feature_columns.json to contain a list: {feature_columns_path}")

    manifest = load_batch_manifest(batch_dir)

    train_paths = paths_for_split(manifest, "train")
    val_paths = paths_for_split(manifest, "validation")
    test_paths = paths_for_split(manifest, "test")

    if len(train_paths) == 0 or len(val_paths) == 0 or len(test_paths) == 0:
        raise ValueError(
            "Need at least one batch for each split. "
            f"Found train={len(train_paths)}, validation={len(val_paths)}, test={len(test_paths)}."
        )

    save_json(output_dir / "config.json", asdict(config))
    save_json(output_dir / "feature_columns.json", features)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")
    print_feature_columns(features)
    print(
        f"Batch files: train={len(train_paths):,}, "
        f"validation={len(val_paths):,}, test={len(test_paths):,}"
    )
    if config.max_train_batches is not None and config.max_train_batches > 0:
        print(
            f"Each epoch will use {min(config.max_train_batches, len(train_paths)):,} "
            f"random train batch(es) from the {len(train_paths):,} available."
        )
    else:
        print(f"Each epoch will use all {len(train_paths):,} train batch(es).")

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
            batch_paths=train_paths,
            optimizer=optimizer,
            loss_fn=loss_fn,
            device=device,
            max_train_batches=config.max_train_batches,
        )

        should_validate = epoch == 1 or epoch % config.val_every == 0

        if should_validate:
            val_predictions = predict_batches(
                model=model,
                batch_paths=val_paths,
                split="validation",
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
                        "config": asdict(config),
                        "best_epoch": best_epoch,
                        "best_val_rmse": best_val_rmse,
                    },
                    best_model_path,
                )
        else:
            history.append(
                {
                    "epoch": epoch,
                    "train_loss": train_loss,
                    "val_n_polygons": np.nan,
                    "val_rmse": np.nan,
                    "val_mae": np.nan,
                    "val_bias": np.nan,
                    "val_r2": np.nan,
                    "val_pearson_r": np.nan,
                }
            )
            print(f"epoch={epoch:03d} train_loss={train_loss:.6f} val_skipped=1")

        if best_epoch > 0 and epoch - best_epoch >= config.patience:
            print(f"Stopping early at epoch {epoch}; best epoch was {best_epoch}.")
            break

    if best_epoch == 0:
        raise RuntimeError("No validation pass completed; cannot save a best model.")

    pd.DataFrame(history).to_csv(output_dir / "training_history.csv", index=False)

    checkpoint = load_torch_file(best_model_path, map_location=device)
    model.load_state_dict(checkpoint["model_state_dict"])

    split_predictions = []
    split_metrics = []

    for split_name, paths in [
        ("train", train_paths),
        ("validation", val_paths),
        ("test", test_paths),
    ]:
        predictions = predict_batches(
            model=model,
            batch_paths=paths,
            split=split_name,
            device=device,
        )
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
