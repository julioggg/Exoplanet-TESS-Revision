#!/usr/bin/env python3
"""Forensic negative controls for acquisition/coverage shortcuts in TESS classification.

Reuses the EXACT star-level train/test membership of experiment 03.
Tests whether almost-identical 30-minute cadence *jitter* and observation-window
geometry still encode the original class assignment after harmonization.

Also audits raw FITS headers (on CEDIA) for differences in product provenance.

Run:
  python revision/04_observational_shortcuts.py --n-jobs 8
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
from scipy.stats import ks_2samp
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.model_selection import StratifiedKFold

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "results/harmonized/harmonized_features.csv"
MANIFEST = ROOT / "results/harmonized/harmonization_manifest.csv"
SPLITS = ROOT / "results/harmonized_eval/holdout_split_membership.csv"
ORIGINAL_ROOT = (
    Path.home() / "exoplanet_revision/scratch/original/Aticulo CVV G3/"
    "exo_moo_tess_data_module"
)
OUTPUT = ROOT / "results/provenance_diagnostics"
SEED = 42
TREE_COUNT = 500

MODELS = {
    "OBS_ALL_PRECISE": [
        "n_points", "baseline_days_after", "median_cadence_minutes_after"
    ],
    "COUNT_PLUS_BASELINE": ["n_points", "baseline_days_after"],
    "COUNT_PLUS_CADENCE": ["n_points", "median_cadence_minutes_after"],
    "BASELINE_PLUS_CADENCE": [
        "baseline_days_after", "median_cadence_minutes_after"
    ],
    "COUNT_ONLY": ["n_points"],
    "BASELINE_ONLY": ["baseline_days_after"],
    "CADENCE_ONLY_PRECISE": ["median_cadence_minutes_after"],
    "OBS_CADENCE_ROUNDED_1MIN": [
        "n_points", "baseline_days_after", "cadence_rounded_minute"
    ],
    "OBS_CADENCE_FIXED_30MIN": [
        "n_points", "baseline_days_after", "cadence_exact_30min"
    ],
    "OCCUPANCY_ONLY": ["occupancy_fraction"],
    "COUNT_BASELINE_OCCUPANCY": [
        "n_points", "baseline_days_after", "occupancy_fraction"
    ],
}


def build_rf(n_jobs):
    return RandomForestClassifier(
        n_estimators=TREE_COUNT,
        class_weight="balanced",
        random_state=SEED,
        n_jobs=n_jobs,
    )


def calculate_metric(y, p):
    from sklearn.metrics import (
        confusion_matrix, f1_score, balanced_accuracy_score
    )
    y = np.asarray(y, dtype=int)
    p = np.asarray(p, dtype=float)
    pred = (p >= .5).astype(int)
    tn, fp, fn, tp = confusion_matrix(y, pred, labels=[0, 1]).ravel()
    return {
        "roc_auc": float(roc_auc_score(y, p)),
        "average_precision": float(average_precision_score(y, p)),
        "f1": float(f1_score(y, pred, zero_division=0)),
        "balanced_accuracy": float(balanced_accuracy_score(y, pred)),
        "tn": int(tn), "fp": int(fp), "fn": int(fn), "tp": int(tp),
    }


def audit_original_fits(original_root):
    rows = []
    for label, category in [(0, "negatives"), (1, "positives")]:
        folder = original_root / "data" / "raw" / category
        if not folder.exists():
            continue
        for path in sorted(folder.glob("*.fits")):
            row = {
                "tic_id": path.stem, "label": label, "fits_file": str(path),
                "fits_size_kb": round(path.stat().st_size / 1024.0, 2),
                "error": "",
            }
            try:
                with fits.open(path, memmap=False, lazy_load_hdus=True) as h:
                    primary = h[0].header
                    extension = h[1].header if len(h) > 1 else {}
                    for key in [
                        "TELESCOP", "ORIGIN", "CREATOR", "PROCVER",
                        "SECTOR", "CAMERA", "CCD", "AUTHOR", "OBJECT",
                        "TICID", "TIMESYS", "TIMEDEL", "BJDREFI",
                        "CADENCE", "EXTNAME", "NAXIS2", "DATE-OBS"
                    ]:
                        value = extension.get(key, primary.get(key, None))
                        row[key.lower()] = (
                            str(value)[:100] if value is not None else ""
                        )
            except Exception as exc:
                row["error"] = repr(exc)
            rows.append(row)
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-jobs", type=int, default=8)
    parser.add_argument("--outdir", type=Path, default=OUTPUT)
    parser.add_argument("--original-root", type=Path, default=ORIGINAL_ROOT)
    parser.add_argument("--cv-folds", type=int, default=5)
    args = parser.parse_args()
    if args.n_jobs <= 0 or args.cv_folds < 2:
        raise ValueError("Invalid n-jobs or cv-folds")
    args.outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(INPUT)
    man = pd.read_csv(MANIFEST)
    saved = pd.read_csv(SPLITS)
    if df["tic_id"].duplicated().any() or saved["tic_id"].duplicated().any():
        raise ValueError("TICs must be unique")
    if not saved["split"].isin(["train", "test"]).all():
        raise ValueError("Unexpected saved split value")
    if len(saved) != len(df):
        raise ValueError("Split membership does not cover all rows")
    if man.duplicated(subset=["tic_id", "label"]).any():
        raise ValueError("Duplicate manifest labels/TICs")
    df = df.merge(
        man[["tic_id", "label", "baseline_days_after",
             "median_cadence_minutes_after"]],
        on=["tic_id", "label"], how="left", validate="one_to_one"
    )
    df = df.merge(saved[["tic_id", "label", "split"]],
                  on=["tic_id", "label"], how="left", validate="one_to_one")
    if df["split"].isna().any():
        raise ValueError("Missing split membership in joined data")
    if not np.isfinite(df["baseline_days_after"]).all():
        raise ValueError("Missing coverage")
    if not np.isfinite(df["median_cadence_minutes_after"]).all():
        raise ValueError("Missing cadence")

    df["cadence_rounded_minute"] = np.round(
        df["median_cadence_minutes_after"]
    )
    df["cadence_exact_30min"] = np.full(len(df), 30.0)
    df["cadence_offset_seconds"] = (
        df["median_cadence_minutes_after"] - 30.0
    ) * 60.0
    df["occupancy_fraction"] = (
        df["n_points"] /
        (df["baseline_days_after"] * 24.0 * 60.0 / 30.0)
    )

    y = df["label"].to_numpy(dtype=int)
    train = np.flatnonzero(df["split"].to_numpy() == "train")
    test = np.flatnonzero(df["split"].to_numpy() == "test")
    if set(y[train]) != {0, 1} or set(y[test]) != {0, 1}:
        raise ValueError("Both classes required in train and test")

    print("=" * 75)
    print("OBSERVATIONAL SHORTCUT FORENSICS")
    print("=" * 75)
    print(f"Targets: {len(df)}, train: {len(train)}, test: {len(test)}", flush=True)

    scores = []
    importances = []
    predictions = []
    cvrows = []

    # Same folds are used for all diagnostic feature combinations.
    splitter = StratifiedKFold(
        n_splits=args.cv_folds, shuffle=True, random_state=SEED
    )
    cv_splits = list(splitter.split(np.zeros(len(train)), y[train]))
    for name, features in MODELS.items():
        x = df[features].to_numpy(dtype=float)
        if not np.isfinite(x).all():
            raise ValueError(f"Nonfinite features in {name}")
        clf = build_rf(args.n_jobs)
        clf.fit(x[train], y[train])
        p = clf.predict_proba(x[test])[:, 1]
        scores.append({
            "model": name, "n_features": len(features),
            **calculate_metric(y[test], p)
        })
        importances.extend({
            "model": name, "feature": col, "importance": float(imp)
        } for col, imp in zip(features, clf.feature_importances_))
        predictions.extend({
            "model": name, "tic_id": str(tic), "label": int(label),
            "probability": float(pp)
        } for tic, label, pp in zip(df.iloc[test]["tic_id"], y[test], p))
        for fold, (a, b) in enumerate(cv_splits, 1):
            model = build_rf(args.n_jobs)
            model.fit(x[train[a]], y[train[a]])
            yp = model.predict_proba(x[train[b]])[:, 1]
            cvrows.append({
                "model": name, "fold": fold,
                **calculate_metric(y[train[b]], yp)
            })
        print(
            f"{name:<30} holdout AUC={scores[-1]['roc_auc']:.4f} "
            f"AP={scores[-1]['average_precision']:.4f}", flush=True
        )

    pd.DataFrame(scores).to_csv(
        args.outdir / "observational_shortcut_models.csv", index=False
    )
    pd.DataFrame(importances).to_csv(
        args.outdir / "observational_shortcut_importances.csv", index=False
    )
    pd.DataFrame(predictions).to_csv(
        args.outdir / "observational_shortcut_holdout_predictions.csv",
        index=False
    )
    pd.DataFrame(cvrows).to_csv(
        args.outdir / "observational_shortcut_train_cv.csv", index=False
    )

    class_balance = []
    for col in [
        "n_points", "baseline_days_after", "median_cadence_minutes_after",
        "cadence_offset_seconds", "occupancy_fraction",
        "cadence_rounded_minute"
    ]:
        pos = df.loc[df.label == 1, col].to_numpy(dtype=float)
        neg = df.loc[df.label == 0, col].to_numpy(dtype=float)
        test_result = ks_2samp(pos, neg)
        class_balance.append({
            "feature": col,
            "negative_median": float(np.median(neg)),
            "positive_median": float(np.median(pos)),
            "negative_q1": float(np.quantile(neg, 0.25)),
            "negative_q3": float(np.quantile(neg, 0.75)),
            "positive_q1": float(np.quantile(pos, 0.25)),
            "positive_q3": float(np.quantile(pos, 0.75)),
            "ks_statistic": float(test_result.statistic),
            "ks_pvalue": float(test_result.pvalue)
        })
    pd.DataFrame(class_balance).to_csv(
        args.outdir / "observational_proxy_distributions.csv", index=False
    )

    fits_inventory = audit_original_fits(args.original_root)
    fits_inventory.to_csv(args.outdir / "original_fits_provenance.csv",
                          index=False)
    if not fits_inventory.empty:
        categories = []
        for key in ["creator", "origin", "procver", "sector", "timedel"]:
            if key not in fits_inventory:
                continue
            for (label, value), grp in fits_inventory.groupby(
                ["label", key], dropna=False
            ):
                categories.append({
                    "field": key, "label": label, "value": value,
                    "n_files": len(grp)
                })
        pd.DataFrame(categories).to_csv(
            args.outdir / "fits_provenance_counts.csv", index=False
        )

    config = {
        "source_dataset": str(INPUT),
        "heldout_membership_reused": str(SPLITS),
        "original_fits_root": str(args.original_root),
        "n_unique_tics": int(df.tic_id.nunique()),
        "n_train": int(len(train)), "n_test": int(len(test)),
        "tree_count": TREE_COUNT, "seed": SEED, "cv_folds": args.cv_folds,
        "cv_on_training_only": True,
        "diagnostic_sets": MODELS,
        "warning": (
            "Predictive observational metadata cannot establish transit "
            "detectability. This test diagnoses acquisition/provenance shortcuts. "
            "Raw FITS provenance fields can be missing or altered on export."
        )
    }
    (args.outdir / "diagnostic_config.json").write_text(
        json.dumps(config, indent=2), encoding="utf-8"
    )
    summary = pd.DataFrame(scores)
    lines = [
        "# TESS observational shortcut diagnostics", "",
        "Repeated the experiment 03 holdout; only observational variables",
        "are included. A true 30-min cadence, rounded to one minute,",
        "eliminates sub-second timestamp fluctuations that could serve",
        "as product or preprocessing fingerprints.", "",
        summary.round(5).to_string(index=False),
        "", "FITS headers audited: " + str(len(fits_inventory)),
        "",
        "Do not use these models for exoplanet detection or planet validation.",
    ]
    (args.outdir / "diagnostic_report.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )
    print("\nFITS headers inventoried:", len(fits_inventory))
    print("Output:", args.outdir)


if __name__ == "__main__":
    main()
