#!/usr/bin/env python3
"""Sector + coverage common-support study for TESS catalogue classification.

Exploratory, NOT definitive detection validation. Matched pairs are constructed
using only observational covariates: recorded FITS sector, n_points, baseline,
and occupancy. Each positive and negative TIC is used at most once.

Paired-star grouped 5-fold CV holds BOTH members of each matched pair together
in training or validation, preventing matched-pair leakage.

Run from repository root:
    python revision/05_sector_matched_evaluation.py --n-jobs 8
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from scipy.stats import ks_2samp
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (
    roc_auc_score, average_precision_score, f1_score,
    balanced_accuracy_score, confusion_matrix,
)
from sklearn.model_selection import StratifiedGroupKFold

ROOT = Path(__file__).resolve().parents[1]
FEATURES_FILE = ROOT / "results/harmonized/harmonized_features.csv"
MANIFEST_FILE = ROOT / "results/harmonized/harmonization_manifest.csv"
FITS_FILE = ROOT / "results/provenance_diagnostics/original_fits_provenance.csv"
OUTDIR = ROOT / "results/sector_matched"

SEED = 42
N_TREES = 500
N_FOLDS = 5
REPEATS = 3

STAT = ["mean_flux", "std_flux", "mad_flux", "amplitude"]
BLS = [
    "best_period", "transit_duration", "transit_depth",
    "bls_power", "snr_estimate",
]
FEATURE_SETS = {
    "OBSERVATIONAL_ONLY": ["n_points", "baseline_days_after", "occupancy_fraction"],
    "BLS_ONLY": BLS,
    "NO_N_POINTS": STAT + BLS,
    "FULL": STAT + BLS + ["n_points"],
}


def clf(n_jobs):
    return RandomForestClassifier(
        n_estimators=N_TREES, class_weight="balanced",
        random_state=SEED, n_jobs=n_jobs,
    )


def perf(y, p):
    y = np.asarray(y, dtype=int)
    pred = (np.asarray(p) >= 0.5).astype(int)
    tn, fp, fn, tp = confusion_matrix(y, pred, labels=[0, 1]).ravel()
    return {
        "roc_auc": float(roc_auc_score(y, p)),
        "average_precision": float(average_precision_score(y, p)),
        "f1": float(f1_score(y, pred, zero_division=0)),
        "balanced_accuracy": float(balanced_accuracy_score(y, pred)),
        "tn": int(tn), "fp": int(fp), "fn": int(fn), "tp": int(tp),
    }


def load():
    feat = pd.read_csv(FEATURES_FILE)
    man = pd.read_csv(MANIFEST_FILE)
    fits_data = pd.read_csv(FITS_FILE)
    if feat["tic_id"].duplicated().any():
        raise ValueError("Duplicate feature TIC")
    if man.duplicated(["tic_id", "label"]).any():
        raise ValueError("Duplicate manifest TIC/label")
    if fits_data.duplicated(["tic_id", "label"]).any():
        raise ValueError("Duplicate FITS provenance TIC/label")
    for needed in ("sector", "error"):
        if needed not in fits_data:
            raise ValueError("Missing FITS provenance column " + needed)

    df = feat.merge(
        man[["tic_id", "label", "baseline_days_after"]],
        on=["tic_id", "label"], how="left", validate="one_to_one"
    ).merge(
        fits_data[["tic_id", "label", "sector", "error"]],
        on=["tic_id", "label"], how="left", validate="one_to_one"
    )
    df["sector"] = pd.to_numeric(df["sector"], errors="coerce")
    if df["sector"].isna().any():
        raise ValueError("Missing/invalid FITS sector values")
    if df["baseline_days_after"].isna().any():
        raise ValueError("Missing baseline values")
    if df["error"].fillna("").str.len().gt(0).any():
        raise ValueError("Some FITS inventories contain errors")
    if len(df) != len(feat):
        raise ValueError("Incomplete metadata join")
    df["sector"] = df["sector"].astype(int)
    df["occupancy_fraction"] = (
        df["n_points"] / (48.0 * df["baseline_days_after"])
    )
    if not np.isfinite(df["occupancy_fraction"]).all():
        raise ValueError("Invalid occupancy")
    return df


def match(df, max_dpoints, max_baseline, max_occupancy):
    """
    Hungarian within-sector pairing. Hard caliper violations receive a large
    penalty in the assignment; rejected pairs remain excluded from tight set.
    This is a design-restriction diagnostic, not random sampling from TESS.
    """
    rows = []
    summary = []
    pair_counter = 0
    for sector in sorted(df.sector.unique()):
        sub = df[df.sector == sector]
        neg = sub[sub.label == 0].sort_values("tic_id")
        pos = sub[sub.label == 1].sort_values("tic_id")
        summary.append({
            "sector": int(sector),
            "n_negative": int(len(neg)),
            "n_positive": int(len(pos)),
            "maximum_possible_pairs": int(min(len(neg), len(pos))),
        })
        if len(neg) == 0 or len(pos) == 0:
            continue
        n0 = neg["n_points"].to_numpy(dtype=float)[:, None]
        n1 = pos["n_points"].to_numpy(dtype=float)[None, :]
        b0 = neg["baseline_days_after"].to_numpy(dtype=float)[:, None]
        b1 = pos["baseline_days_after"].to_numpy(dtype=float)[None, :]
        o0 = neg["occupancy_fraction"].to_numpy(dtype=float)[:, None]
        o1 = pos["occupancy_fraction"].to_numpy(dtype=float)[None, :]
        dn = np.abs(n0 - n1)
        db = np.abs(b0 - b1)
        do = np.abs(o0 - o1)
        good = (dn <= max_dpoints) & (db <= max_baseline) & (
            do <= max_occupancy
        )
        # Positive penalties prioritize in-caliper assignment where possible,
        # while still reporting all available sector-matched pairings.
        cost = dn / 150.0 + db / 1.5 + do / 0.10
        cost = cost + (~good) * 10000.0
        ri, ci = linear_sum_assignment(cost)
        for a, b in zip(ri, ci):
            pair_counter += 1
            rows.append({
                "pair_id": f"P{pair_counter:04d}",
                "sector": int(sector),
                "negative_tic": str(neg.iloc[a]["tic_id"]),
                "positive_tic": str(pos.iloc[b]["tic_id"]),
                "n_points_delta": float(dn[a, b]),
                "baseline_days_delta": float(db[a, b]),
                "occupancy_delta": float(do[a, b]),
                "match_cost_unpenalized": float(
                    dn[a, b] / 150.0 + db[a, b] / 1.5 + do[a, b] / 0.10
                ),
                "passes_calipers": bool(good[a, b]),
            })
    return pd.DataFrame(rows), pd.DataFrame(summary)


def matched_dataset(df, pairs):
    if len(pairs) == 0:
        return pd.DataFrame()
    neg = pairs[["pair_id", "negative_tic"]].rename(
        columns={"negative_tic": "tic_id"}
    )
    pos = pairs[["pair_id", "positive_tic"]].rename(
        columns={"positive_tic": "tic_id"}
    )
    membership = pd.concat([neg, pos], ignore_index=True)
    if membership["tic_id"].duplicated().any():
        raise ValueError("A TIC was reused in matched sample")
    out = df.merge(membership, on="tic_id", how="inner", validate="one_to_one")
    if len(out) != 2 * len(pairs):
        raise ValueError("Pair membership incomplete")
    checks = out.groupby("pair_id")["label"].agg(["count", "sum"])
    if not ((checks["count"] == 2) & (checks["sum"] == 1)).all():
        raise ValueError("Matched pairs must contain 1 positive and 1 negative")
    return out.sort_values(["pair_id", "label"]).reset_index(drop=True)


def balance(df, design):
    if len(df) == 0:
        return pd.DataFrame()
    rows = []
    for col in ["n_points", "baseline_days_after", "occupancy_fraction"]:
        a = df.loc[df.label == 1, col].to_numpy(dtype=float)
        b = df.loc[df.label == 0, col].to_numpy(dtype=float)
        k = ks_2samp(a, b)
        rows.append({
            "design": design, "feature": col,
            "negative_median": float(np.median(b)),
            "positive_median": float(np.median(a)),
            "ks_statistic": float(k.statistic),
            "ks_pvalue": float(k.pvalue),
            "mean_difference_positive_minus_negative": float(
                np.mean(a) - np.mean(b)
            ),
        })
    return pd.DataFrame(rows)


def evaluate_grouped(df, design, n_jobs, repeats):
    if df.empty:
        return [], []
    n_pairs = df.pair_id.nunique()
    if n_pairs < N_FOLDS:
        print(f"[SKIP] {design}: only {n_pairs} pairs", flush=True)
        return [], []

    y = df["label"].to_numpy(dtype=int)
    groups = df["pair_id"].to_numpy()
    output = []
    predictions = []

    for repeat in range(repeats):
        splitter = StratifiedGroupKFold(
            n_splits=N_FOLDS, shuffle=True,
            random_state=SEED + repeat * 101,
        )
        splits = list(splitter.split(
            np.zeros((len(df), 1)), y, groups=groups
        ))
        for name, features in FEATURE_SETS.items():
            X = df[features].to_numpy(dtype=float)
            if not np.isfinite(X).all():
                raise ValueError(f"Nonfinite feature in {name}")
            for fold, (tr, te) in enumerate(splits, 1):
                if set(groups[tr]) & set(groups[te]):
                    raise AssertionError("Pair leakage across folds")
                if len(np.unique(y[te])) != 2:
                    raise ValueError("Fold lacks both classes")
                model = clf(n_jobs)
                model.fit(X[tr], y[tr])
                p = model.predict_proba(X[te])[:, 1]
                output.append({
                    "design": design, "model": name,
                    "repeat": repeat + 1, "fold": fold,
                    "n_pairs_total": n_pairs,
                    "n_train": len(tr), "n_valid": len(te),
                    **perf(y[te], p),
                })
                predictions.extend({
                    "design": design, "model": name,
                    "repeat": repeat + 1, "fold": fold,
                    "pair_id": str(pair_id), "tic_id": str(tic),
                    "label": int(yi), "probability": float(pi),
                } for pair_id, tic, yi, pi in zip(
                    df.iloc[te]["pair_id"], df.iloc[te]["tic_id"],
                    y[te], p
                ))
            print(f"[DONE] {design} / {name}", flush=True)
    return output, predictions


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--n-jobs", type=int, default=8)
    p.add_argument("--repeats", type=int, default=REPEATS)
    p.add_argument("--max-points-diff", type=float, default=150.0)
    p.add_argument("--max-baseline-diff", type=float, default=1.5)
    p.add_argument("--max-occupancy-diff", type=float, default=0.10)
    p.add_argument("--outdir", type=Path, default=OUTDIR)
    args = p.parse_args()
    if args.n_jobs < 1 or args.repeats < 1:
        raise ValueError("Positive counts required")
    args.outdir.mkdir(parents=True, exist_ok=True)

    df = load()
    pairs, sectors = match(
        df, args.max_points_diff, args.max_baseline_diff,
        args.max_occupancy_diff
    )
    sectors.to_csv(args.outdir / "sector_class_counts.csv", index=False)
    pairs.to_csv(args.outdir / "sector_matched_pairs.csv", index=False)

    all_matched = matched_dataset(df, pairs)
    tight = matched_dataset(df, pairs[pairs.passes_calipers].copy())
    all_matched[["pair_id", "tic_id", "label", "sector",
                 "n_points", "baseline_days_after",
                 "occupancy_fraction"]].to_csv(
        args.outdir / "all_sector_matched_membership.csv", index=False
    )
    if not tight.empty:
        tight[["pair_id", "tic_id", "label", "sector",
               "n_points", "baseline_days_after",
               "occupancy_fraction"]].to_csv(
            args.outdir / "caliper_matched_membership.csv", index=False
        )
    else:
        pd.DataFrame(columns=[
            "pair_id", "tic_id", "label", "sector", "n_points",
            "baseline_days_after", "occupancy_fraction"
        ]).to_csv(
            args.outdir / "caliper_matched_membership.csv", index=False
        )

    balance_df = pd.concat([
        balance(df, "ALL_TARGETS"),
        balance(all_matched, "SECTOR_ONLY_MATCHED"),
        balance(tight, "SECTOR_AND_CALIPER_MATCHED"),
    ], ignore_index=True)
    balance_df.to_csv(args.outdir / "matched_covariate_balance.csv",
                      index=False)

    records, predictions = [], []
    for name, data in [
        ("SECTOR_ONLY_MATCHED", all_matched),
        ("SECTOR_AND_CALIPER_MATCHED", tight),
    ]:
        r, pred = evaluate_grouped(data, name, args.n_jobs, args.repeats)
        records.extend(r)
        predictions.extend(pred)

    if records:
        cv = pd.DataFrame(records)
        cv.to_csv(args.outdir / "pair_grouped_cv_folds.csv", index=False)
        metrics = ["roc_auc", "average_precision", "f1", "balanced_accuracy"]
        summary = cv.melt(
            id_vars=["design", "model", "repeat", "fold"],
            value_vars=metrics, var_name="metric", value_name="value"
        ).groupby(["design", "model", "metric"])["value"].agg(
            ["count", "mean", "std", "min", "max"]
        ).reset_index()
    else:
        cv = pd.DataFrame()
        summary = pd.DataFrame()
        cv.to_csv(args.outdir / "pair_grouped_cv_folds.csv", index=False)
    summary.to_csv(args.outdir / "pair_grouped_cv_summary.csv", index=False)
    pd.DataFrame(predictions).to_csv(
        args.outdir / "pair_grouped_cv_predictions.csv", index=False
    )

    config = {
        "full_dataset_size": len(df),
        "recorded_sector_overlap": int(
            ((sectors.n_negative > 0) & (sectors.n_positive > 0)).sum()
        ),
        "maximum_sector_pairs": int(sectors.maximum_possible_pairs.sum()),
        "matched_sector_pairs": int(len(pairs)),
        "tight_pairs": int(len(tight) // 2),
        "calipers": {
            "max_points_delta": args.max_points_diff,
            "max_baseline_days_delta": args.max_baseline_diff,
            "max_occupancy_fraction_delta": args.max_occupancy_diff,
        },
        "pair_grouped_cross_validation": True,
        "n_folds": N_FOLDS, "n_repeats": args.repeats,
        "n_trees": N_TREES, "random_seed": SEED,
        "models": FEATURE_SETS,
        "limitations": [
            "Sector from lightkurve-exported FITS may not represent every stitched sector.",
            "Matching uses labels to select common-support pairs.",
            "Group CV fold metrics are not independent confidence intervals.",
            "Small retained sample; not independently representative of TESS.",
            "The original positive and negative upstream pipelines differ.",
            "Catalogue-negative class is not a verified astrophysical true negative."
        ]
    }
    (args.outdir / "match_config.json").write_text(
        json.dumps(config, indent=2), encoding="utf-8"
    )

    report = [
        "# Within-sector and observational-covariate matched diagnosis", "",
        f"Original targets: {len(df)}",
        f"Sector-overlap count: {config['recorded_sector_overlap']}",
        f"Maximum pairs by recorded sector: {config['maximum_sector_pairs']}",
        f"Pairs matched within sector: {len(pairs)}",
        f"Pairs passing all strict calipers: {len(tight) // 2}", "",
        "## Covariate balance", "",
        balance_df.round(6).to_string(index=False), "",
        "## Pair-grouped cross-validation summaries", "",
        summary.round(5).to_string(index=False) if not summary.empty else (
            "No evaluable matched pairs"
        ), "",
        "## Limitations", "",
        "Matching is only a sensitivity check on a restricted support region.",
        "FITS-recorded sector may not describe all observations in stitched curves.",
        "No confirmed planet is not the same as a verified non-transiting target.",
        "Rebuild from homogeneous source products and run injection/recovery next.",
    ]
    (args.outdir / "match_report.md").write_text(
        "\n".join(report) + "\n", encoding="utf-8"
    )

    print("=" * 78)
    print("SECTOR AND COVERAGE MATCHED ANALYSIS")
    print("=" * 78)
    print(f"Sector common-support maximum pairs: {len(pairs)}")
    print(f"Caliper-matched pairs: {len(tight) // 2}")
    print("Matched CV result summary:")
    print(summary.round(4).to_string(index=False) if not summary.empty
          else "No evaluable matched pairs")
    print(f"Saved to {args.outdir}")


if __name__ == "__main__":
    main()
