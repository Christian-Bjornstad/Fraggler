from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, average_precision_score, roc_auc_score


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Offline evaluation scaffold for clonality candidate artifacts.")
    parser.add_argument("--ladder-csv", type=Path, required=True, help="Path to clonality_ladder_candidates.csv")
    parser.add_argument("--pk-csv", type=Path, required=True, help="Path to clonality_pk_candidates.csv")
    parser.add_argument("--labels-csv", type=Path, default=None, help="Optional gold-label CSV")
    parser.add_argument("--output-dir", type=Path, required=True, help="Directory for metrics and prediction outputs")
    parser.add_argument("--task", choices=("ladder", "pk", "both"), default="both")
    parser.add_argument("--group-by", type=str, default="source_run_dir")
    parser.add_argument("--holdout-pattern", type=str, default=None)
    parser.add_argument("--min-labels", type=int, default=10)
    parser.add_argument("--include-sl", action="store_true")
    parser.add_argument("--write-predictions", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    output_dir = Path(args.output_dir).expanduser()
    output_dir.mkdir(parents=True, exist_ok=True)

    ladder = _read_csv(args.ladder_csv)
    pk = _read_csv(args.pk_csv)
    labels = _read_csv(args.labels_csv) if args.labels_csv else pd.DataFrame()
    if not args.include_sl:
        if "assay" in ladder.columns:
            ladder = ladder[ladder["assay"].astype(str).str.upper() != "SL"].copy()
        if "assay" in pk.columns:
            pk = pk[pk["assay"].astype(str).str.upper() != "SL"].copy()

    tasks = ["ladder", "pk"] if args.task == "both" else [args.task]
    metrics: dict[str, dict[str, object]] = {}
    outputs: dict[str, str] = {}
    for task in tasks:
        table = ladder if task == "ladder" else pk
        task_metrics, predictions, weights = _evaluate_task(
            task,
            table,
            labels,
            group_by=args.group_by,
            holdout_pattern=args.holdout_pattern,
            min_labels=args.min_labels,
        )
        metrics[task] = task_metrics
        if args.write_predictions and predictions is not None:
            prediction_path = output_dir / f"{task}_predictions.csv"
            predictions.to_csv(prediction_path, index=False)
            outputs[f"{task}_predictions"] = str(prediction_path)
            
        if weights is not None:
            weight_path = output_dir / f"{task}_weights.json"
            weight_path.write_text(json.dumps(weights, indent=2, sort_keys=True) + "\n", encoding="utf-8")
            outputs[f"{task}_weights"] = str(weight_path)

    metrics_path = output_dir / "metrics.json"
    manifest_path = output_dir / "model_manifest.json"
    metrics_path.write_text(json.dumps(metrics, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    manifest_path.write_text(
        json.dumps(
            {
                "task": args.task,
                "group_by": args.group_by,
                "holdout_pattern": args.holdout_pattern,
                "include_sl": bool(args.include_sl),
                "inputs": {
                    "ladder_csv": str(Path(args.ladder_csv).expanduser()),
                    "pk_csv": str(Path(args.pk_csv).expanduser()),
                    "labels_csv": str(Path(args.labels_csv).expanduser()) if args.labels_csv else "",
                },
                "outputs": outputs,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"metrics: {metrics_path}")
    print(f"manifest: {manifest_path}")
    for name, path in outputs.items():
        print(f"{name}: {path}")
    return 0


def _read_csv(path: Path | None) -> pd.DataFrame:
    if path is None:
        return pd.DataFrame()
    return pd.read_csv(Path(path).expanduser()).fillna("")


def _evaluate_task(
    task: str,
    table: pd.DataFrame,
    labels: pd.DataFrame,
    *,
    group_by: str,
    holdout_pattern: str | None,
    min_labels: int,
) -> tuple[dict[str, object], pd.DataFrame | None, dict[str, Any] | None]:
    metrics: dict[str, object] = {
        "row_count": int(len(table)),
        "group_count": int(table[group_by].nunique()) if group_by in table.columns and not table.empty else 0,
    }
    if table.empty:
        metrics["status"] = "empty_table"
        return metrics, None, None
    if labels.empty:
        metrics["status"] = "missing_labels"
        return metrics, None, None

    label_view = labels.copy()
    label_view["artifact_table"] = label_view.get("artifact_table", pd.Series("", index=label_view.index)).astype(str)
    label_view["artifact_row_key"] = label_view.get("artifact_row_key", pd.Series("", index=label_view.index)).astype(str)
    if "label_source" in label_view.columns:
        label_source = label_view["label_source"].astype(str).str.strip().str.lower()
        auto_mask = label_source.isin({"auto_generated_from_pipeline", "pipeline_pseudo_label"})
        if auto_mask.any():
            metrics["excluded_auto_label_rows"] = int(auto_mask.sum())
            label_view = label_view.loc[~auto_mask].copy()
    task_table_name = "ladder_candidates" if task == "ladder" else "pk_candidates"
    label_view = label_view[label_view["artifact_table"] == task_table_name]
    if label_view.empty:
        metrics["status"] = "missing_labels"
        return metrics, None, None

    merged = table.merge(label_view, how="inner", on=["artifact_table", "artifact_row_key"])
    metrics["labeled_row_count"] = int(len(merged))
    if len(merged) < min_labels:
        metrics["status"] = "insufficient_labels"
        return metrics, merged, None

    y = _label_target(task, merged["label"].astype(str))
    if y.nunique() < 2:
        metrics["status"] = "single_class_labels"
        return metrics, merged, None

    holdout_mask = _build_holdout_mask(merged, group_by=group_by, holdout_pattern=holdout_pattern)
    if not holdout_mask.any():
        metrics["status"] = "train_only"
        train_mask = np.ones(len(merged), dtype=bool)
        eval_mask = train_mask
    else:
        metrics["status"] = "ok"
        train_mask = ~holdout_mask
        eval_mask = holdout_mask
        if y[train_mask].nunique() < 2:
            metrics["status"] = "holdout_single_class_train"
            train_mask = np.ones(len(merged), dtype=bool)
            eval_mask = train_mask

    X = _prepare_features(task, merged)
    model = LogisticRegression(max_iter=500)
    model.fit(X.loc[train_mask], y.loc[train_mask])
    eval_probs = model.predict_proba(X.loc[eval_mask])[:, 1]
    eval_true = y.loc[eval_mask]
    eval_pred = (eval_probs >= 0.5).astype(int)

    metrics["train_rows"] = int(train_mask.sum())
    metrics["eval_rows"] = int(eval_mask.sum())
    metrics["positive_rate"] = float(eval_true.mean())
    metrics["accuracy"] = float(accuracy_score(eval_true, eval_pred))
    if eval_true.nunique() > 1:
        metrics["roc_auc"] = float(roc_auc_score(eval_true, eval_probs))
        metrics["average_precision"] = float(average_precision_score(eval_true, eval_probs))

    # Save model weights for production use
    weights = {
        "columns": X.columns.tolist(),
        "coef": model.coef_[0].tolist(),
        "intercept": float(model.intercept_[0]),
    }
    
    predictions = merged.loc[eval_mask, ["artifact_table", "artifact_row_key", "label"]].copy()
    predictions["predicted_probability"] = eval_probs
    predictions["predicted_label"] = eval_pred
    return metrics, predictions, weights


def _label_target(task: str, labels: pd.Series) -> pd.Series:
    normalized = labels.astype(str).str.strip().str.lower()
    if task == "ladder":
        positives = {"accept", "selected_correct"}
    else:
        positives = {"selected_correct"}
    return normalized.isin(positives).astype(int)


def _build_holdout_mask(frame: pd.DataFrame, *, group_by: str, holdout_pattern: str | None) -> np.ndarray:
    if group_by not in frame.columns:
        return np.zeros(len(frame), dtype=bool)
    groups = frame[group_by].astype(str)
    if holdout_pattern:
        return groups.str.contains(holdout_pattern, regex=True, na=False).to_numpy()
    unique_groups = sorted(group for group in groups.unique().tolist() if group)
    if len(unique_groups) < 2:
        return np.zeros(len(frame), dtype=bool)
    return (groups == unique_groups[-1]).to_numpy()


def _prepare_features(task: str, frame: pd.DataFrame) -> pd.DataFrame:
    if task == "ladder":
        numeric_cols = [
            "candidate_index",
            "candidate_time",
            "candidate_intensity",
            "selected_step_bp",
            "ladder_r2",
            "ladder_review_required",
        ]
        categorical_cols = ["candidate_source", "assay", "control", "sample_kind", "ladder", "ladder_fit_strategy"]
    else:
        numeric_cols = [
            "expected_bp",
            "window_bp",
            "search_window_bp",
            "ok",
            "found_bp",
            "delta_bp",
            "height",
            "area",
            "selection_score",
            "fallback_from_window_bp",
        ]
        categorical_cols = ["marker_name", "kind", "channel", "search_mode", "assay", "control", "sample_kind"]

    numeric = frame.reindex(columns=numeric_cols, fill_value=0).copy()
    for col in numeric.columns:
        numeric[col] = pd.to_numeric(numeric[col], errors="coerce").fillna(0.0)
    categorical = pd.get_dummies(frame.reindex(columns=categorical_cols, fill_value="").astype(str), dtype=float)
    return pd.concat([numeric, categorical], axis=1)


if __name__ == "__main__":
    raise SystemExit(main())
