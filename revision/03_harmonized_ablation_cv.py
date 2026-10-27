#!/usr/bin/env python3
"""Controlled RF ablations on 30-minute harmonized TESS data.

Identical 70/30 stratified holdout and hyperparameters for every feature set.
Repeated stratified 5-fold CV uses ONLY the 70% training partition.
Unique TIC per row is asserted; injection experiments later require group-by-base-TIC.

Run: python revision/03_harmonized_ablation_cv.py --n-jobs 8
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import ks_2samp, mannwhitneyu
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (accuracy_score, average_precision_score,
    balanced_accuracy_score, confusion_matrix, f1_score, precision_score,
    recall_score, roc_auc_score)
from sklearn.model_selection import RepeatedStratifiedKFold, train_test_split

ROOT = Path(__file__).resolve().parents[1]
FEATURES = ROOT / "results/harmonized/harmonized_features.csv"
MANIFEST = ROOT / "results/harmonized/harmonization_manifest.csv"
OUTPUT = ROOT / "results/harmonized_eval"
SEED, TEST_FRAC, N_TREES, N_FOLDS = 42, 0.30, 500, 5

STATS = ["mean_flux", "median_flux", "std_flux", "mad_flux", "amplitude"]
BLS = ["best_period", "transit_duration", "transit_depth", "bls_power", "snr_estimate"]
BLS_SCALED = ["best_period", "transit_duration", "transit_depth",
              "bls_power_noise_scaled", "bls_depth_snr_noise_scaled"]
MODEL_SETS = {
    "FULL_RAW": STATS + ["n_points"] + BLS,
    "NO_N_POINTS_RAW": STATS + BLS,
    "N_POINTS_ONLY": ["n_points"],
    "BLS_ONLY_RAW": BLS,
    "FULL_NOISE_SCALED": ["mean_flux", "std_flux", "mad_flux",
                          "amplitude", "n_points"] + BLS_SCALED,
    "NO_N_POINTS_NOISE_SCALED": ["mean_flux", "std_flux",
                                 "mad_flux", "amplitude"] + BLS_SCALED,
    "BLS_ONLY_NOISE_SCALED": BLS_SCALED,
    "OBSERVATIONAL_ONLY": ["n_points", "baseline_days_after",
                            "median_cadence_minutes_after"],
    "STATISTICAL_ONLY": ["mean_flux", "std_flux", "mad_flux", "amplitude"],
    "TRANSIT_PARAMETERS_ONLY": ["best_period", "transit_duration", "transit_depth"],
}


def score(y, proba):
    y = np.asarray(y, dtype=int)
    p = np.asarray(proba, dtype=float)
    pred = (p >= 0.5).astype(int)
    tn, fp, fn, tp = confusion_matrix(y, pred, labels=[0, 1]).ravel()
    return {
        "roc_auc": float(roc_auc_score(y, p)),
        "average_precision": float(average_precision_score(y, p)),
        "accuracy": float(accuracy_score(y, pred)),
        "balanced_accuracy": float(balanced_accuracy_score(y, pred)),
        "precision": float(precision_score(y, pred, zero_division=0)),
        "recall": float(recall_score(y, pred, zero_division=0)),
        "specificity": float(tn / (tn + fp)) if tn + fp else np.nan,
        "f1": float(f1_score(y, pred, zero_division=0)),
        "tn": int(tn), "fp": int(fp), "fn": int(fn), "tp": int(tp)
    }


def classifier(n_jobs):
    return RandomForestClassifier(
        n_estimators=N_TREES, class_weight="balanced",
        random_state=SEED, n_jobs=n_jobs
    )


def observational_balance(df):
    out = []
    for col in ["n_points", "baseline_days_after",
                "median_cadence_minutes_after"]:
        a = df.loc[df.label == 1, col].to_numpy(dtype=float)
        b = df.loc[df.label == 0, col].to_numpy(dtype=float)
        gt = sum(int(np.sum(v > b)) for v in a)
        lt = sum(int(np.sum(v < b)) for v in a)
        ks = ks_2samp(a, b)
        mw = mannwhitneyu(a, b, alternative="two-sided")
        out.append({
            "feature": col, "negative_median": float(np.median(b)),
            "positive_median": float(np.median(a)),
            "cliffs_delta_pos_minus_neg": float((gt - lt) / (len(a) * len(b))),
            "ks_statistic": float(ks.statistic),
            "ks_pvalue": float(ks.pvalue),
            "mannwhitney_pvalue": float(mw.pvalue)
        })
    return pd.DataFrame(out)


def paired_bootstrap(y, p1, p0, candidate, reference, count):
    """Paired class-stratified holdout bootstrap: test-sample uncertainty only."""
    y = np.asarray(y)
    pos = np.flatnonzero(y == 1)
    neg = np.flatnonzero(y == 0)
    rng = np.random.default_rng(20261027)
    metrics = ["roc_auc", "average_precision", "f1", "balanced_accuracy"]
    values = {m: [] for m in metrics}
    point1, point0 = score(y, p1), score(y, p0)
    for _ in range(count):
        ii = np.concatenate([
            rng.choice(pos, size=len(pos), replace=True),
            rng.choice(neg, size=len(neg), replace=True)
        ])
        r1, r0 = score(y[ii], p1[ii]), score(y[ii], p0[ii])
        for m in metrics:
            values[m].append(r1[m] - r0[m])
    rows = []
    for m in metrics:
        v = np.asarray(values[m])
        rows.append({
            "model": candidate, "reference": reference, "metric": m,
            "observed_delta": point1[m] - point0[m],
            "ci95_low": float(np.quantile(v, 0.025)),
            "ci95_high": float(np.quantile(v, 0.975)),
            "bootstrap_delta_mean": float(np.mean(v)),
            "n_bootstrap": count
        })
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-jobs", type=int, default=8)
    parser.add_argument("--n-repeats", type=int, default=3)
    parser.add_argument("--n-bootstrap", type=int, default=2000)
    parser.add_argument("--outdir", type=Path, default=OUTPUT)
    args = parser.parse_args()
    if min(args.n_jobs, args.n_repeats, args.n_bootstrap) <= 0:
        raise ValueError("Counts must be positive")
    args.outdir.mkdir(parents=True, exist_ok=True)

    raw = pd.read_csv(FEATURES)
    manifest = pd.read_csv(MANIFEST)
    if raw.duplicated(subset=["tic_id"]).any():
        raise ValueError("Duplicate TICs: star grouping required")
    if manifest.duplicated(subset=["tic_id", "label"]).any():
        raise ValueError("Duplicate harmonization manifest TIC/label pairs")
    if not manifest["status"].eq("ok").all():
        raise ValueError("Harmonization manifest contains failed targets")
    obs_cols = ["tic_id", "label", "baseline_days_after",
                "median_cadence_minutes_after"]
    df = raw.merge(manifest[obs_cols], on=["tic_id", "label"], how="inner",
                   validate="one_to_one")
    if len(df) != len(raw):
        raise ValueError("Incomplete feature/manifest join")
    if set(df["label"]) != {0, 1}:
        raise ValueError("Both classes required")

    # BLS received no dy, so raw 'depth_snr' is not physically calibrated.
    # A constant per-curve noise approximation rescales the original signal:
    # sigma = 1.4826 * MAD(normalized flux). These are sensitivity proxies.
    sigma = 1.4826 * df["mad_flux"].to_numpy(dtype=float)
    if not np.isfinite(sigma).all() or (sigma <= 0).any():
        raise ValueError("Nonpositive/invalid robust noise scale")
    df["noise_sigma_mad"] = sigma
    df["bls_depth_snr_noise_scaled"] = df["bls_depth_snr"] / sigma
    df["bls_power_noise_scaled"] = df["bls_power"] / sigma**2

    feature_cols = sorted({x for feats in MODEL_SETS.values() for x in feats})
    all_values = df[feature_cols].to_numpy(dtype=float)
    if not np.isfinite(all_values).all():
        raise ValueError("Nonfinite model inputs")
    y = df["label"].to_numpy(dtype=int)
    idx = np.arange(len(df))
    train, test = train_test_split(idx, test_size=TEST_FRAC, random_state=SEED,
                                   stratify=y)
    y_train, y_test = y[train], y[test]
    pd.DataFrame({
        "tic_id": df.tic_id, "label": y,
        "split": np.where(np.isin(idx, test), "test", "train")
    }).to_csv(args.outdir / "holdout_split_membership.csv", index=False)

    balance = observational_balance(df)
    balance.to_csv(args.outdir / "observational_balance_tests.csv", index=False)
    df.groupby("label")["noise_sigma_mad"].agg(
        ["count", "median", "mean", "min", "max"]
    ).to_csv(args.outdir / "noise_proxy_by_class.csv")

    cv = RepeatedStratifiedKFold(n_splits=N_FOLDS,
        n_repeats=args.n_repeats, random_state=SEED)
    folds = list(cv.split(np.zeros(len(train)), y_train))
    holdout, cv_scores, pred_rows, importances, proba = [], [], [], [], {}

    for name, features in MODEL_SETS.items():
        print(f"[RUN] {name}: {len(features)} variables", flush=True)
        X = df[features].to_numpy(dtype=float)
        model = classifier(args.n_jobs)
        model.fit(X[train], y_train)
        p = model.predict_proba(X[test])[:, 1]
        proba[name] = p
        holdout.append({"model": name, "n_features": len(features),
                        **score(y_test, p)})
        importances.extend({
            "model": name, "feature": f, "importance": float(imp)
        } for f, imp in zip(features, model.feature_importances_))
        pred_rows.extend({
            "model": name, "tic_id": str(tic), "label": int(label),
            "probability": float(pp), "prediction": int(pp >= 0.5)
        } for tic, label, pp in zip(df.iloc[test].tic_id, y_test, p))
        for k, (tr, va) in enumerate(folds):
            clf = classifier(args.n_jobs)
            clf.fit(X[train[tr]], y_train[tr])
            v = clf.predict_proba(X[train[va]])[:, 1]
            cv_scores.append({
                "model": name, "repeat": int(k // N_FOLDS + 1),
                "fold": int(k % N_FOLDS + 1),
                "n_validation": int(len(va)), **score(y_train[va], v)
            })
        print(f"[DONE] {name}", flush=True)

    holdout = pd.DataFrame(holdout)
    cv_scores = pd.DataFrame(cv_scores)
    holdout.to_csv(args.outdir / "holdout_metrics.csv", index=False)
    cv_scores.to_csv(args.outdir / "train_repeated_cv_folds.csv", index=False)
    pd.DataFrame(pred_rows).to_csv(args.outdir / "holdout_predictions.csv",
                                    index=False)
    pd.DataFrame(importances).to_csv(args.outdir / "feature_importances.csv",
                                      index=False)
    val = cv_scores.melt(
        id_vars=["model", "repeat", "fold"],
        value_vars=["roc_auc", "average_precision", "f1", "balanced_accuracy"],
        var_name="metric", value_name="value"
    )
    cv_summary = val.groupby(["model", "metric"])["value"].agg(
        ["count", "mean", "std", "min", "max"]
    ).reset_index()
    cv_summary.to_csv(args.outdir / "train_repeated_cv_summary.csv", index=False)

    comparisons = [
        ("NO_N_POINTS_RAW", "FULL_RAW"),
        ("N_POINTS_ONLY", "FULL_RAW"),
        ("BLS_ONLY_RAW", "FULL_RAW"),
        ("NO_N_POINTS_NOISE_SCALED", "FULL_NOISE_SCALED"),
        ("BLS_ONLY_NOISE_SCALED", "FULL_NOISE_SCALED"),
        ("OBSERVATIONAL_ONLY", "FULL_NOISE_SCALED"),
    ]
    paired = []
    for candidate, reference in comparisons:
        paired.extend(paired_bootstrap(y_test, proba[candidate], proba[reference],
                                        candidate, reference, args.n_bootstrap))
    pd.DataFrame(paired).to_csv(
        args.outdir / "paired_holdout_bootstrap_deltas.csv", index=False
    )
    cfg = {
        "n_stars": int(len(df)), "n_unique_tics": int(df.tic_id.nunique()),
        "train_size": len(train), "test_size": len(test),
        "seed": SEED, "test_fraction": TEST_FRAC,
        "trees": N_TREES, "class_weight": "balanced",
        "n_repeats": args.n_repeats, "n_folds": N_FOLDS,
        "cv_training_only": True, "n_jobs": args.n_jobs,
        "n_bootstrap": args.n_bootstrap, "feature_sets": MODEL_SETS,
        "noise_scaled_columns_are_proxies_not_calibrated_snr": True,
        "limitations": [
            "Processed light curves were not regenerated from raw FITS.",
            "No confirmed planet is not a verified true-negative label.",
            "Injection/recovery and eclipsing-binary challenges are not yet done.",
            "BLS period grid extends beyond 27-day time baseline.",
            "Holdout bootstrap excludes retraining, selection and label uncertainty."
        ],
        "python_version": sys.version.split()[0],
        "numpy_version": np.__version__, "pandas_version": pd.__version__
    }
    (args.outdir / "evaluation_config.json").write_text(
        json.dumps(cfg, indent=2), encoding="utf-8"
    )
    report = [
        "# Harmonized TESS: controlled RF ablation and repeated CV",
        "",
        "Each TIC appears once. Same 70/30 holdout and RF hyperparameters.",
        "Repeated 5-fold CV is restricted to the 70% training set.",
        "",
        "HOLDOUT METRICS",
        holdout.round(5).to_string(index=False),
        "",
        "POST-HARMONIZATION OBSERVATIONAL BALANCE",
        balance.round(6).to_string(index=False),
        "",
        "REPEATED CV FOLD SUMMARY (NOT CONFIDENCE INTERVALS)",
        cv_summary.round(5).to_string(index=False),
        "",
        "CAUTIONS",
        "Class-dependent preprocessing upstream of the harmonized curves remains.",
        "Negative classes denote catalogue status, not proven absence of transits.",
        "BLS noise-scaled variables are uncalibrated constant-noise proxies.",
        "BLS searches periods up to 30 days even though coverage is <=27 days."
    ]
    (args.outdir / "evaluation_report.md").write_text(
        "\n".join(report) + "\n", encoding="utf-8"
    )
    print("\nHARMONIZED HOLDOUT METRICS", flush=True)
    print(holdout[
        ["model", "roc_auc", "average_precision", "f1", "balanced_accuracy",
         "tn", "fp", "fn", "tp"]
    ].round(5).to_string(index=False), flush=True)
    print("\nOBSERVATIONAL BALANCE", flush=True)
    print(balance[
        ["feature", "negative_median", "positive_median",
         "cliffs_delta_pos_minus_neg", "ks_pvalue"]
    ].round(6).to_string(index=False), flush=True)
    print(f"\nSaved: {args.outdir}", flush=True)


if __name__ == "__main__":
    main()
