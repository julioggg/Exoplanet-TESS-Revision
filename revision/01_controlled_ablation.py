#!/usr/bin/env python3
"""
Controlled feature-ablation experiment for the submitted TESS Random Forest.

Purpose
-------
Quantify the contribution of n_points and feature families while holding constant:
- the exact train/test membership,
- Random Forest hyperparameters,
- class weighting,
- random seed,
- decision threshold.

This is a baseline sensitivity experiment on the SUBMITTED feature table. It does
not solve the known class-dependent cadence/download/BLS-grid confounding; those
will be addressed by the symmetric reconstruction pipeline.

Run from repository root:
    python revision/01_controlled_ablation.py

Outputs:
    results/ablations/
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    confusion_matrix,
    f1_score,
    precision_score,
    recall_score,
    roc_auc_score,
)
from sklearn.model_selection import train_test_split


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATASET = REPO_ROOT / "data/dataset/final_dataset.csv"
DEFAULT_OUTDIR = REPO_ROOT / "results/ablations"

RANDOM_STATE = 42
TEST_SIZE = 0.30
N_ESTIMATORS = 500
CLASS_WEIGHT = "balanced"
THRESHOLD = 0.50
N_BOOTSTRAP = 2000
BOOTSTRAP_SEED = 20261027

FEATURE_SETS = {
    # Submitted feature vector.
    "FULL": [
        "mean_flux",
        "median_flux",
        "std_flux",
        "mad_flux",
        "amplitude",
        "n_points",
        "best_period",
        "transit_duration",
        "transit_depth",
        "bls_power",
        "snr_estimate",
    ],

    # Reviewer-requested controlled ablation.
    "NO_N_POINTS": [
        "mean_flux",
        "median_flux",
        "std_flux",
        "mad_flux",
        "amplitude",
        "best_period",
        "transit_duration",
        "transit_depth",
        "bls_power",
        "snr_estimate",
    ],

    # Reviewer-requested coverage proxy baseline.
    "N_POINTS_ONLY": [
        "n_points",
    ],

    # Reviewer-requested BLS-only model.
    "BLS_ONLY": [
        "best_period",
        "transit_duration",
        "transit_depth",
        "bls_power",
        "snr_estimate",
    ],

    # Transit parameters without detection-strength proxies.
    "TRANSIT_PARAMETERS_ONLY": [
        "best_period",
        "transit_duration",
        "transit_depth",
    ],

    # Light-curve descriptive statistics, excluding coverage.
    "STATISTICAL_ONLY": [
        "mean_flux",
        "median_flux",
        "std_flux",
        "mad_flux",
        "amplitude",
    ],

    # Noise/variability descriptors only.
    "VARIABILITY_ONLY": [
        "std_flux",
        "mad_flux",
        "amplitude",
    ],

    # Detection-strength proxies only.
    "DETECTION_STRENGTH_ONLY": [
        "bls_power",
        "snr_estimate",
    ],
}


def dataframe_to_markdown(df: pd.DataFrame) -> str:
    """Minimal Markdown table writer with no optional dependencies."""
    cols = [str(c) for c in df.columns]
    lines = [
        "| " + " | ".join(cols) + " |",
        "| " + " | ".join(["---"] * len(cols)) + " |",
    ]
    for _, row in df.iterrows():
        vals = []
        for col in df.columns:
            value = row[col]
            if isinstance(value, float):
                vals.append(f"{value:.6g}")
            else:
                vals.append(str(value))
        lines.append("| " + " | ".join(vals) + " |")
    return "\\n".join(lines)


def specificity_score(y_true, y_pred):
    cm = confusion_matrix(y_true, y_pred, labels=[0, 1])
    tn, fp, fn, tp = cm.ravel()
    return tn / (tn + fp) if (tn + fp) else np.nan


def metrics_from_predictions(y_true, y_prob, threshold=THRESHOLD):
    y_pred = (np.asarray(y_prob) >= threshold).astype(int)
    y_true = np.asarray(y_true).astype(int)

    cm = confusion_matrix(y_true, y_pred, labels=[0, 1])
    tn, fp, fn, tp = cm.ravel()

    return {
        "accuracy": float(accuracy_score(y_true, y_pred)),
        "precision": float(precision_score(y_true, y_pred, zero_division=0)),
        "recall": float(recall_score(y_true, y_pred, zero_division=0)),
        "specificity": float(specificity_score(y_true, y_pred)),
        "f1": float(f1_score(y_true, y_pred, zero_division=0)),
        "roc_auc": float(roc_auc_score(y_true, y_prob)),
        "average_precision": float(average_precision_score(y_true, y_prob)),
        "tn": int(tn),
        "fp": int(fp),
        "fn": int(fn),
        "tp": int(tp),
    }


def bootstrap_metrics(y_true, y_prob, n_bootstrap=N_BOOTSTRAP, seed=BOOTSTRAP_SEED):
    rng = np.random.default_rng(seed)
    y_true = np.asarray(y_true).astype(int)
    y_prob = np.asarray(y_prob, dtype=float)
    n = len(y_true)
    records = []

    for b in range(n_bootstrap):
        idx = rng.integers(0, n, size=n)
        yt = y_true[idx]
        yp = y_prob[idx]

        # ROC-AUC requires both classes.
        if len(np.unique(yt)) < 2:
            continue

        row = metrics_from_predictions(yt, yp)
        row["bootstrap_id"] = b
        records.append(row)

    return pd.DataFrame(records)


def ci_summary(bootstrap_df, model_name):
    metrics = [
        "accuracy",
        "precision",
        "recall",
        "specificity",
        "f1",
        "roc_auc",
        "average_precision",
    ]
    rows = []
    for metric in metrics:
        x = pd.to_numeric(bootstrap_df[metric], errors="coerce").dropna()
        rows.append({
            "model": model_name,
            "metric": metric,
            "n_bootstrap_valid": int(len(x)),
            "bootstrap_mean": float(x.mean()),
            "ci95_low": float(x.quantile(0.025)),
            "ci95_high": float(x.quantile(0.975)),
        })
    return rows


def paired_bootstrap_delta(
    y_true,
    prob_model,
    prob_reference,
    model_name,
    reference_name="FULL",
    n_bootstrap=N_BOOTSTRAP,
    seed=BOOTSTRAP_SEED + 100,
):
    rng = np.random.default_rng(seed)
    y_true = np.asarray(y_true).astype(int)
    p_model = np.asarray(prob_model, dtype=float)
    p_ref = np.asarray(prob_reference, dtype=float)
    n = len(y_true)

    metrics = ["roc_auc", "average_precision", "f1", "accuracy"]
    deltas = {m: [] for m in metrics}

    for _ in range(n_bootstrap):
        idx = rng.integers(0, n, size=n)
        yt = y_true[idx]
        if len(np.unique(yt)) < 2:
            continue

        mm = metrics_from_predictions(yt, p_model[idx])
        rr = metrics_from_predictions(yt, p_ref[idx])

        for metric in metrics:
            deltas[metric].append(mm[metric] - rr[metric])

    rows = []
    for metric, vals in deltas.items():
        x = pd.Series(vals, dtype=float)
        rows.append({
            "model": model_name,
            "reference": reference_name,
            "metric": metric,
            "delta_mean": float(x.mean()),
            "delta_median": float(x.median()),
            "ci95_low": float(x.quantile(0.025)),
            "ci95_high": float(x.quantile(0.975)),
            "p_delta_gt_0": float((x > 0).mean()),
            "p_delta_lt_0": float((x < 0).mean()),
            "n_bootstrap_valid": int(len(x)),
        })
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", type=Path, default=DEFAULT_DATASET)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--n-bootstrap", type=int, default=N_BOOTSTRAP)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.dataset).copy()
    required = {"tic_id", "label"}
    missing_required = required - set(df.columns)
    if missing_required:
        raise ValueError(f"Missing required columns: {sorted(missing_required)}")

    # Validate all requested features.
    all_requested = sorted({f for fs in FEATURE_SETS.values() for f in fs})
    missing_features = [f for f in all_requested if f not in df.columns]
    if missing_features:
        raise ValueError(f"Missing requested features: {missing_features}")

    # Keep exact row indices so every model sees the exact same stars.
    all_idx = np.arange(len(df))
    y = pd.to_numeric(df["label"], errors="raise").astype(int).to_numpy()

    train_idx, test_idx = train_test_split(
        all_idx,
        test_size=TEST_SIZE,
        random_state=RANDOM_STATE,
        stratify=y,
    )

    split = pd.DataFrame({
        "row_index": np.concatenate([train_idx, test_idx]),
        "tic_id": pd.concat(
            [
                df.iloc[train_idx]["tic_id"],
                df.iloc[test_idx]["tic_id"],
            ],
            ignore_index=True,
        ),
        "label": np.concatenate([y[train_idx], y[test_idx]]),
        "split": ["train"] * len(train_idx) + ["test"] * len(test_idx),
    })
    split.to_csv(args.outdir / "split_membership.csv", index=False)

    metric_rows = []
    prediction_frames = []
    importance_frames = []
    ci_rows = []
    probabilities = {}

    for model_name, features in FEATURE_SETS.items():
        X = df[features].apply(pd.to_numeric, errors="raise")

        X_train = X.iloc[train_idx]
        X_test = X.iloc[test_idx]
        y_train = y[train_idx]
        y_test = y[test_idx]

        model = RandomForestClassifier(
            n_estimators=N_ESTIMATORS,
            class_weight=CLASS_WEIGHT,
            random_state=RANDOM_STATE,
            n_jobs=-1,
        )
        model.fit(X_train, y_train)

        prob = model.predict_proba(X_test)[:, 1]
        pred = (prob >= THRESHOLD).astype(int)
        probabilities[model_name] = prob

        m = metrics_from_predictions(y_test, prob)
        m.update({
            "model": model_name,
            "n_features": len(features),
            "features": "|".join(features),
            "n_train": len(train_idx),
            "n_test": len(test_idx),
            "n_estimators": N_ESTIMATORS,
            "class_weight": CLASS_WEIGHT,
            "random_state": RANDOM_STATE,
            "test_size": TEST_SIZE,
            "threshold": THRESHOLD,
        })
        metric_rows.append(m)

        prediction_frames.append(pd.DataFrame({
            "model": model_name,
            "row_index": test_idx,
            "tic_id": df.iloc[test_idx]["tic_id"].to_numpy(),
            "y_true": y_test,
            "y_prob": prob,
            "y_pred": pred,
        }))

        importance_frames.append(pd.DataFrame({
            "model": model_name,
            "feature": features,
            "importance": model.feature_importances_,
        }))

        bdf = bootstrap_metrics(
            y_test,
            prob,
            n_bootstrap=args.n_bootstrap,
            seed=BOOTSTRAP_SEED,
        )
        ci_rows.extend(ci_summary(bdf, model_name))

    metrics_df = pd.DataFrame(metric_rows)
    metrics_df = metrics_df[
        [
            "model",
            "n_features",
            "features",
            "accuracy",
            "precision",
            "recall",
            "specificity",
            "f1",
            "roc_auc",
            "average_precision",
            "tn",
            "fp",
            "fn",
            "tp",
            "n_train",
            "n_test",
            "n_estimators",
            "class_weight",
            "random_state",
            "test_size",
            "threshold",
        ]
    ]
    metrics_df.to_csv(args.outdir / "controlled_ablation_metrics.csv", index=False)

    predictions_df = pd.concat(prediction_frames, ignore_index=True)
    predictions_df.to_csv(args.outdir / "controlled_ablation_predictions.csv", index=False)

    importances_df = pd.concat(importance_frames, ignore_index=True)
    importances_df.to_csv(args.outdir / "controlled_ablation_importances.csv", index=False)

    ci_df = pd.DataFrame(ci_rows)
    ci_df.to_csv(args.outdir / "controlled_ablation_bootstrap_ci.csv", index=False)

    paired_rows = []
    y_test = y[test_idx]
    ref_prob = probabilities["FULL"]
    for model_name, prob in probabilities.items():
        if model_name == "FULL":
            continue
        paired_rows.extend(
            paired_bootstrap_delta(
                y_test,
                prob,
                ref_prob,
                model_name=model_name,
                reference_name="FULL",
                n_bootstrap=args.n_bootstrap,
            )
        )
    paired_df = pd.DataFrame(paired_rows)
    paired_df.to_csv(args.outdir / "paired_bootstrap_deltas_vs_full.csv", index=False)

    config = {
        "dataset": str(args.dataset),
        "n_rows": int(len(df)),
        "n_unique_tics": int(df["tic_id"].nunique()),
        "test_size": TEST_SIZE,
        "random_state": RANDOM_STATE,
        "n_estimators": N_ESTIMATORS,
        "class_weight": CLASS_WEIGHT,
        "threshold": THRESHOLD,
        "n_bootstrap": int(args.n_bootstrap),
        "bootstrap_seed": BOOTSTRAP_SEED,
        "feature_sets": FEATURE_SETS,
        "interpretation_guardrail": (
            "These are sensitivity analyses on the submitted feature table. "
            "Because positive/negative light curves were generated with different "
            "cadences/download strategies and different BLS grids, these results "
            "must not be interpreted as the final deconfounded performance."
        ),
    }
    (args.outdir / "controlled_ablation_config.json").write_text(
        json.dumps(config, indent=2),
        encoding="utf-8",
    )

    # Compact Markdown report.
    display_cols = [
        "model",
        "accuracy",
        "precision",
        "recall",
        "specificity",
        "f1",
        "roc_auc",
        "average_precision",
        "tn",
        "fp",
        "fn",
        "tp",
    ]
    report = []
    report.append("# Controlled ablation on the submitted feature table")
    report.append("")
    report.append(
        "All models use the exact same train/test stars and identical Random Forest "
        "hyperparameters (500 trees, balanced class weights, random_state=42). "
        "Only the feature set changes."
    )
    report.append("")
    report.append(
        "**Important:** this experiment quantifies baseline feature dependence but "
        "does not remove the known class-dependent cadence/download/BLS-grid confounding."
    )
    report.append("")
    report.append("## Hold-out metrics")
    report.append("")
    report.append(dataframe_to_markdown(metrics_df[display_cols].round(6)))
    report.append("")
    report.append("## Key reviewer comparisons")
    report.append("")
    full = metrics_df.set_index("model").loc["FULL"]
    no_n = metrics_df.set_index("model").loc["NO_N_POINTS"]
    np_only = metrics_df.set_index("model").loc["N_POINTS_ONLY"]
    bls = metrics_df.set_index("model").loc["BLS_ONLY"]
    report.append(
        f"- FULL ROC-AUC={full.roc_auc:.4f}, F1={full.f1:.4f}."
    )
    report.append(
        f"- NO_N_POINTS ROC-AUC={no_n.roc_auc:.4f}, F1={no_n.f1:.4f}."
    )
    report.append(
        f"- N_POINTS_ONLY ROC-AUC={np_only.roc_auc:.4f}, F1={np_only.f1:.4f}."
    )
    report.append(
        f"- BLS_ONLY ROC-AUC={bls.roc_auc:.4f}, F1={bls.f1:.4f}."
    )
    report.append("")
    report.append(
        "Interpretation must be deferred until the symmetric reconstruction is run, "
        "because n_points is strongly coupled to class-dependent cadence and the BLS "
        "features were computed using different positive/negative grids."
    )
    report.append("")
    (args.outdir / "controlled_ablation_report.md").write_text(
        "\n".join(report) + "\n",
        encoding="utf-8",
    )

    print("=" * 78)
    print("CONTROLLED ABLATION — SUBMITTED FEATURE TABLE")
    print("=" * 78)
    print(metrics_df[display_cols].round(4).to_string(index=False))
    print()
    print("Saved:", args.outdir)
    print("NOTE: baseline sensitivity only; not yet deconfounded.")
    print("=" * 78)


if __name__ == "__main__":
    main()
