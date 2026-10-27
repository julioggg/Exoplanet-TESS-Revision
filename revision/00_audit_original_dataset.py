#!/usr/bin/env python3
"""
Audit the submitted TESS exoplanet dataset and its original preprocessing pipeline.

This script is intentionally diagnostic: it does NOT alter the submitted baseline
datasets. It writes compact audit artifacts under results/audit/ so they can be
versioned in GitHub while the large raw/processed light curves remain on CEDIA.

Run from the repository root:
    python revision/00_audit_original_dataset.py

Optional:
    python revision/00_audit_original_dataset.py \
      --original-root "$HOME/exoplanet_revision/scratch/original/Aticulo CVV G3/exo_moo_tess_data_module"
"""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATASET = REPO_ROOT / "data/dataset/final_dataset.csv"
DEFAULT_ORIGINAL_ROOT = (
    Path.home()
    / "exoplanet_revision/scratch/original/Aticulo CVV G3/exo_moo_tess_data_module"
)
OUTDIR = REPO_ROOT / "results/audit"


def tic_normalize(value) -> str:
    s = str(value).strip()
    s = s.replace(" ", "_")
    if not s.startswith("TIC_"):
        if s.startswith("TIC"):
            s = s.replace("TIC", "TIC_", 1).replace("__", "_")
        else:
            s = f"TIC_{s}"
    return s


def safe_float(value):
    try:
        x = float(value)
        return x if np.isfinite(x) else None
    except Exception:
        return None


def summary_row(values: pd.Series) -> dict:
    x = pd.to_numeric(values, errors="coerce").dropna().astype(float)
    if len(x) == 0:
        return {
            "n": 0,
            "mean": np.nan,
            "std": np.nan,
            "median": np.nan,
            "q1": np.nan,
            "q3": np.nan,
            "iqr": np.nan,
            "min": np.nan,
            "max": np.nan,
        }
    q1 = x.quantile(0.25)
    q3 = x.quantile(0.75)
    return {
        "n": int(len(x)),
        "mean": float(x.mean()),
        "std": float(x.std(ddof=1)) if len(x) > 1 else 0.0,
        "median": float(x.median()),
        "q1": float(q1),
        "q3": float(q3),
        "iqr": float(q3 - q1),
        "min": float(x.min()),
        "max": float(x.max()),
    }


def cliffs_delta(a, b) -> float:
    """Cliff's delta: positive means values in a tend to exceed values in b."""
    a = np.asarray(pd.to_numeric(pd.Series(a), errors="coerce").dropna(), dtype=float)
    b = np.asarray(pd.to_numeric(pd.Series(b), errors="coerce").dropna(), dtype=float)
    if len(a) == 0 or len(b) == 0:
        return np.nan
    gt = 0
    lt = 0
    for x in a:
        gt += int(np.sum(x > b))
        lt += int(np.sum(x < b))
    return (gt - lt) / (len(a) * len(b))


def standardized_mean_difference(a, b) -> float:
    a = np.asarray(pd.to_numeric(pd.Series(a), errors="coerce").dropna(), dtype=float)
    b = np.asarray(pd.to_numeric(pd.Series(b), errors="coerce").dropna(), dtype=float)
    if len(a) < 2 or len(b) < 2:
        return np.nan
    va = np.var(a, ddof=1)
    vb = np.var(b, ddof=1)
    pooled = math.sqrt(((len(a) - 1) * va + (len(b) - 1) * vb) / (len(a) + len(b) - 2))
    if pooled == 0:
        return 0.0 if np.mean(a) == np.mean(b) else np.inf
    return float((np.mean(a) - np.mean(b)) / pooled)


def distribution_tests(df: pd.DataFrame, features: list[str]) -> pd.DataFrame:
    rows = []
    have_scipy = False
    try:
        from scipy.stats import ks_2samp, mannwhitneyu
        have_scipy = True
    except Exception:
        ks_2samp = mannwhitneyu = None

    for feature in features:
        if feature not in df.columns:
            continue
        neg = pd.to_numeric(df.loc[df["label"] == 0, feature], errors="coerce").dropna()
        pos = pd.to_numeric(df.loc[df["label"] == 1, feature], errors="coerce").dropna()
        row = {
            "feature": feature,
            "n_negative": int(len(neg)),
            "n_positive": int(len(pos)),
            "negative_median": float(neg.median()) if len(neg) else np.nan,
            "positive_median": float(pos.median()) if len(pos) else np.nan,
            "median_ratio_pos_over_neg": (
                float(pos.median() / neg.median())
                if len(pos) and len(neg) and neg.median() != 0
                else np.nan
            ),
            "standardized_mean_difference_pos_minus_neg": standardized_mean_difference(pos, neg),
            "cliffs_delta_pos_minus_neg": cliffs_delta(pos, neg),
            "ks_statistic": np.nan,
            "ks_pvalue": np.nan,
            "mannwhitney_u": np.nan,
            "mannwhitney_pvalue": np.nan,
        }
        if have_scipy and len(neg) and len(pos):
            ks = ks_2samp(pos, neg, alternative="two-sided", mode="auto")
            mw = mannwhitneyu(pos, neg, alternative="two-sided")
            row.update({
                "ks_statistic": float(ks.statistic),
                "ks_pvalue": float(ks.pvalue),
                "mannwhitney_u": float(mw.statistic),
                "mannwhitney_pvalue": float(mw.pvalue),
            })
        rows.append(row)

    return pd.DataFrame(rows)


def audit_processed_curves(original_root: Path, dataset: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Measure per-target coverage/cadence from processed curves kept on CEDIA."""
    records = []
    class_paths = {
        0: original_root / "data/processed/negatives",
        1: original_root / "data/processed/positives",
    }

    for label, directory in class_paths.items():
        if not directory.exists():
            continue

        for path in sorted(directory.glob("*_processed.csv")):
            tic_id = tic_normalize(path.stem.replace("_processed", ""))
            try:
                usecols = None
                probe = pd.read_csv(path, nrows=1)
                if "TIME" not in probe.columns:
                    raise ValueError("TIME column missing")
                usecols = ["TIME"]
                d = pd.read_csv(path, usecols=usecols)
                t = pd.to_numeric(d["TIME"], errors="coerce").dropna().to_numpy(dtype=float)
                t = np.sort(np.unique(t))
                n = len(t)
                if n >= 2:
                    dt = np.diff(t)
                    dt = dt[np.isfinite(dt) & (dt > 0)]
                else:
                    dt = np.array([], dtype=float)

                cadence_days = float(np.median(dt)) if len(dt) else np.nan
                cadence_min = cadence_days * 24.0 * 60.0 if np.isfinite(cadence_days) else np.nan
                baseline = float(t[-1] - t[0]) if n >= 2 else 0.0
                large_gap_threshold = 5.0 * cadence_days if np.isfinite(cadence_days) else np.nan
                large_gaps = int(np.sum(dt > large_gap_threshold)) if len(dt) and np.isfinite(large_gap_threshold) else 0

                records.append({
                    "tic_id": tic_id,
                    "label_from_directory": label,
                    "n_points_processed": int(len(d)),
                    "n_unique_times": int(n),
                    "baseline_days": baseline,
                    "median_cadence_minutes": cadence_min,
                    "cadence_p10_minutes": float(np.quantile(dt, 0.10) * 1440.0) if len(dt) else np.nan,
                    "cadence_p90_minutes": float(np.quantile(dt, 0.90) * 1440.0) if len(dt) else np.nan,
                    "large_gap_count_gt_5xcadence": large_gaps,
                    "processed_file": str(path),
                    "read_error": "",
                })
            except Exception as exc:
                records.append({
                    "tic_id": tic_id,
                    "label_from_directory": label,
                    "n_points_processed": np.nan,
                    "n_unique_times": np.nan,
                    "baseline_days": np.nan,
                    "median_cadence_minutes": np.nan,
                    "cadence_p10_minutes": np.nan,
                    "cadence_p90_minutes": np.nan,
                    "large_gap_count_gt_5xcadence": np.nan,
                    "processed_file": str(path),
                    "read_error": repr(exc),
                })

    coverage = pd.DataFrame(records)
    if coverage.empty:
        return coverage, pd.DataFrame()

    dsmall = dataset[["tic_id", "label", "n_points"]].copy()
    dsmall["tic_id"] = dsmall["tic_id"].map(tic_normalize)
    coverage = coverage.merge(dsmall, on="tic_id", how="left")
    coverage["n_points_difference_dataset_minus_processed"] = (
        pd.to_numeric(coverage["n_points"], errors="coerce")
        - pd.to_numeric(coverage["n_points_processed"], errors="coerce")
    )
    coverage["label_matches_directory"] = coverage["label"] == coverage["label_from_directory"]

    metrics = [
        "n_points_processed",
        "baseline_days",
        "median_cadence_minutes",
        "cadence_p10_minutes",
        "cadence_p90_minutes",
        "large_gap_count_gt_5xcadence",
    ]
    rows = []
    for label in [0, 1]:
        sub = coverage[coverage["label_from_directory"] == label]
        for metric in metrics:
            s = summary_row(sub[metric])
            s.update({"label": label, "metric": metric})
            rows.append(s)
    class_summary = pd.DataFrame(rows)

    return coverage, class_summary


def parse_bls_grid(text: str) -> dict:
    lin = re.search(
        r"np\.linspace\s*\(\s*([0-9.]+)\s*,\s*([0-9.]+)\s*,\s*([0-9]+)\s*\)",
        text,
        flags=re.S,
    )
    dur = re.search(
        r"durations\s*=\s*np\.array\s*\(\s*\[([^\]]+)\]",
        text,
        flags=re.S,
    )
    out = {
        "period_min": None,
        "period_max": None,
        "period_grid_points": None,
        "durations": None,
    }
    if lin:
        out.update({
            "period_min": float(lin.group(1)),
            "period_max": float(lin.group(2)),
            "period_grid_points": int(lin.group(3)),
        })
    if dur:
        vals = re.findall(r"[-+]?(?:\d*\.\d+|\d+)", dur.group(1))
        out["durations"] = [float(x) for x in vals]
    return out


def source_code_audit(repo_root: Path) -> dict:
    src = repo_root / "src"
    files = {
        "positive_download": src / "03c_download_200_positive_targets.py",
        "negative_download": src / "12_download_negative_fits.py",
        "positive_preprocess": src / "08_preprocess_positive_dataset.py",
        "negative_preprocess": src / "13_preprocess_negatives.py",
        "positive_bls": src / "09_extract_bls_positive.py",
        "negative_bls": src / "14_extract_negative_bls_features.py",
        "positive_stats": src / "10_extract_statistical_features.py",
        "negative_stats": src / "15_extract_negative_statistical_features.py",
        "rf": src / "18_train_random_forest.py",
        "rf_without_npoints": src / "19_train_without_npoints.py",
    }
    text = {k: p.read_text(encoding="utf-8", errors="replace") for k, p in files.items()}

    pos_bls = parse_bls_grid(text["positive_bls"])
    neg_bls = parse_bls_grid(text["negative_bls"])

    return {
        "positive_download_uses_first_search_result_only": "search[0].download()" in text["positive_download"],
        "negative_download_uses_download_all": "download_all()" in text["negative_download"],
        "negative_download_stitches_collection": ".stitch()" in text["negative_download"],
        "positive_preprocess_minimum_500_points": bool(
            re.search(r"len\(df\)\s*<\s*500", text["positive_preprocess"])
        ),
        "negative_preprocess_has_same_minimum_500_points": bool(
            re.search(r"len\(df\)\s*<\s*500", text["negative_preprocess"])
        ),
        "positive_bls_grid": pos_bls,
        "negative_bls_grid": neg_bls,
        "bls_grids_identical": pos_bls == neg_bls,
        "positive_std_uses_numpy_default_ddof0": "np.std(flux)" in text["positive_stats"],
        "negative_std_uses_pandas_series_default_ddof1": "flux.std()" in text["negative_stats"],
        "submitted_rf_n_estimators": 500 if "n_estimators=500" in text["rf"].replace(" ", "") else None,
        "submitted_rf_balanced_class_weight": 'class_weight="balanced"' in text["rf"].replace(" ", ""),
        "without_npoints_rf_n_estimators": 200 if "n_estimators=200" in text["rf_without_npoints"].replace(" ", "") else None,
        "without_npoints_rf_balanced_class_weight": 'class_weight="balanced"' in text["rf_without_npoints"].replace(" ", ""),
        "clean_ablation_comparison": (
            "n_estimators=500" in text["rf"].replace(" ", "")
            and "n_estimators=500" in text["rf_without_npoints"].replace(" ", "")
            and ('class_weight="balanced"' in text["rf"].replace(" ", ""))
            == ('class_weight="balanced"' in text["rf_without_npoints"].replace(" ", ""))
        ),
    }


def write_markdown_report(
    dataset: pd.DataFrame,
    feature_summary: pd.DataFrame,
    tests: pd.DataFrame,
    coverage: pd.DataFrame,
    coverage_summary: pd.DataFrame,
    source_audit: dict,
    outpath: Path,
) -> None:
    lines = []
    lines.append("# Original dataset and pipeline audit")
    lines.append("")
    lines.append("This report describes the submitted baseline. It does not modify the original dataset.")
    lines.append("")
    lines.append("## Dataset")
    lines.append("")
    lines.append(f"- Rows: {len(dataset)}")
    lines.append(f"- Unique TIC IDs: {dataset['tic_id'].nunique()}")
    lines.append(f"- Duplicate TIC IDs: {int(dataset['tic_id'].duplicated(keep=False).sum())}")
    lines.append(f"- Positive rows: {int((dataset['label'] == 1).sum())}")
    lines.append(f"- Negative rows: {int((dataset['label'] == 0).sum())}")
    lines.append("")

    if "n_points" in dataset.columns:
        for label, name in [(0, "negative"), (1, "positive")]:
            x = dataset.loc[dataset["label"] == label, "n_points"]
            lines.append(
                f"- n_points {name}: median={x.median():.1f}, "
                f"mean={x.mean():.1f}, IQR=[{x.quantile(.25):.1f}, {x.quantile(.75):.1f}], "
                f"range=[{x.min():.0f}, {x.max():.0f}]"
            )
        lines.append("")

    lines.append("## Source-code asymmetries detected")
    lines.append("")
    lines.append(f"- Positive download uses only first search result: {source_audit['positive_download_uses_first_search_result_only']}")
    lines.append(f"- Negative download uses download_all(): {source_audit['negative_download_uses_download_all']}")
    lines.append(f"- Negative download stitches products: {source_audit['negative_download_stitches_collection']}")
    lines.append(f"- Positive preprocessing enforces >=500 points: {source_audit['positive_preprocess_minimum_500_points']}")
    lines.append(f"- Negative preprocessing enforces same threshold: {source_audit['negative_preprocess_has_same_minimum_500_points']}")
    lines.append(f"- Positive BLS grid: {source_audit['positive_bls_grid']}")
    lines.append(f"- Negative BLS grid: {source_audit['negative_bls_grid']}")
    lines.append(f"- BLS grids identical: {source_audit['bls_grids_identical']}")
    lines.append(f"- Positive std uses NumPy ddof=0: {source_audit['positive_std_uses_numpy_default_ddof0']}")
    lines.append(f"- Negative std uses pandas ddof=1: {source_audit['negative_std_uses_pandas_series_default_ddof1']}")
    lines.append(f"- Existing no-n_points comparison is controlled: {source_audit['clean_ablation_comparison']}")
    lines.append("")

    if "median_flux" in dataset.columns:
        mf = pd.to_numeric(dataset["median_flux"], errors="coerce")
        lines.append("## Median-flux diagnostic")
        lines.append("")
        lines.append(f"- Overall median(median_flux): {mf.median():.12g}")
        lines.append(f"- SD(median_flux): {mf.std(ddof=1):.12g}")
        lines.append(f"- Fraction within 1e-3 of unity: {float((mf.sub(1).abs() <= 1e-3).mean()):.4f}")
        lines.append("")

    if not coverage.empty:
        lines.append("## Processed-curve observational coverage")
        lines.append("")
        lines.append(f"- Processed curves audited: {len(coverage)}")
        lines.append(f"- Curves linked to final dataset: {int(coverage['label'].notna().sum())}")
        for metric in ["n_points_processed", "baseline_days", "median_cadence_minutes"]:
            lines.append(f"- {metric}:")
            for label, name in [(0, "negative"), (1, "positive")]:
                row = coverage_summary[
                    (coverage_summary["label"] == label)
                    & (coverage_summary["metric"] == metric)
                ]
                if len(row):
                    r = row.iloc[0]
                    lines.append(
                        f"  - {name}: median={r['median']:.6g}, "
                        f"IQR=[{r['q1']:.6g}, {r['q3']:.6g}], n={int(r['n'])}"
                    )
        lines.append("")

    lines.append("## Interpretation")
    lines.append("")
    lines.append(
        "The submitted dataset must be treated as potentially confounded by observational coverage "
        "until a symmetric positive/negative light-curve construction and feature extraction pipeline "
        "is rerun. RF impurity importance and SHAP agreement alone cannot rule out this confounding."
    )
    lines.append(
        "The existing script that removes n_points is not a controlled ablation because its Random "
        "Forest hyperparameters differ from the submitted model."
    )
    lines.append(
        "All revised ablations should therefore use identical splits and identical classifier "
        "hyperparameters, changing only the feature set."
    )
    lines.append("")

    outpath.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", type=Path, default=DEFAULT_DATASET)
    parser.add_argument("--original-root", type=Path, default=DEFAULT_ORIGINAL_ROOT)
    parser.add_argument("--outdir", type=Path, default=OUTDIR)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    if not args.dataset.exists():
        raise FileNotFoundError(f"Dataset not found: {args.dataset}")

    df = pd.read_csv(args.dataset)
    required = {"tic_id", "label"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    df["tic_id"] = df["tic_id"].map(tic_normalize)
    df["label"] = pd.to_numeric(df["label"], errors="raise").astype(int)

    numeric_features = [
        c for c in df.columns
        if c not in {"tic_id", "label"} and pd.api.types.is_numeric_dtype(df[c])
    ]

    duplicate_tics = df[df["tic_id"].duplicated(keep=False)].sort_values("tic_id")
    duplicate_tics.to_csv(args.outdir / "duplicate_tics.csv", index=False)

    conflicts = (
        df.groupby("tic_id")["label"]
        .nunique()
        .rename("n_unique_labels")
        .reset_index()
    )
    conflicts = conflicts[conflicts["n_unique_labels"] > 1]
    conflicts.to_csv(args.outdir / "cross_label_tic_conflicts.csv", index=False)

    feature_rows = []
    for feature in numeric_features:
        for label in [0, 1]:
            s = summary_row(df.loc[df["label"] == label, feature])
            s.update({"feature": feature, "label": label})
            feature_rows.append(s)
    feature_summary = pd.DataFrame(feature_rows)
    feature_summary.to_csv(args.outdir / "feature_summary_by_class.csv", index=False)

    tests = distribution_tests(df, numeric_features)
    tests.to_csv(args.outdir / "feature_distribution_tests.csv", index=False)

    coverage, coverage_summary = audit_processed_curves(args.original_root, df)
    coverage.to_csv(args.outdir / "observational_per_target.csv", index=False)
    coverage_summary.to_csv(args.outdir / "observational_summary_by_class.csv", index=False)

    if not coverage.empty:
        obs_metrics = [
            c for c in [
                "n_points_processed",
                "baseline_days",
                "median_cadence_minutes",
                "cadence_p10_minutes",
                "cadence_p90_minutes",
                "large_gap_count_gt_5xcadence",
            ]
            if c in coverage.columns
        ]
        ctest = coverage.rename(columns={"label_from_directory": "label"})
        obs_tests = distribution_tests(ctest, obs_metrics)
    else:
        obs_tests = pd.DataFrame()
    obs_tests.to_csv(args.outdir / "observational_distribution_tests.csv", index=False)

    source_audit = source_code_audit(REPO_ROOT)
    (args.outdir / "source_code_audit.json").write_text(
        json.dumps(source_audit, indent=2, sort_keys=True),
        encoding="utf-8",
    )

    missingness = df.isna().sum().rename("missing_count").to_frame()
    missingness["missing_fraction"] = missingness["missing_count"] / len(df)
    missingness.to_csv(args.outdir / "missingness.csv")

    finite_rows = []
    for feature in numeric_features:
        x = pd.to_numeric(df[feature], errors="coerce").to_numpy(dtype=float)
        finite_rows.append({
            "feature": feature,
            "n_total": len(x),
            "n_finite": int(np.isfinite(x).sum()),
            "n_nonfinite": int((~np.isfinite(x)).sum()),
        })
    pd.DataFrame(finite_rows).to_csv(args.outdir / "finite_value_audit.csv", index=False)

    overview = {
        "dataset_path": str(args.dataset),
        "original_root": str(args.original_root),
        "n_rows": int(len(df)),
        "n_unique_tics": int(df["tic_id"].nunique()),
        "n_duplicate_rows_by_tic_membership": int(df["tic_id"].duplicated(keep=False).sum()),
        "n_cross_label_tic_conflicts": int(len(conflicts)),
        "class_counts": {str(k): int(v) for k, v in df["label"].value_counts().sort_index().items()},
        "features": numeric_features,
        "processed_curve_audit_available": bool(not coverage.empty),
        "n_processed_curves_audited": int(len(coverage)),
        "source_code_audit": source_audit,
    }
    if "median_flux" in df.columns:
        mf = pd.to_numeric(df["median_flux"], errors="coerce")
        overview["median_flux"] = {
            "median": safe_float(mf.median()),
            "std": safe_float(mf.std(ddof=1)),
            "fraction_within_1e-3_of_one": safe_float((mf.sub(1).abs() <= 1e-3).mean()),
        }
    if "n_points" in df.columns:
        overview["n_points_by_class"] = {}
        for label in [0, 1]:
            x = pd.to_numeric(df.loc[df["label"] == label, "n_points"], errors="coerce").dropna()
            overview["n_points_by_class"][str(label)] = {
                "n": int(len(x)),
                "mean": safe_float(x.mean()),
                "median": safe_float(x.median()),
                "q1": safe_float(x.quantile(0.25)),
                "q3": safe_float(x.quantile(0.75)),
                "min": safe_float(x.min()),
                "max": safe_float(x.max()),
            }

    (args.outdir / "audit_overview.json").write_text(
        json.dumps(overview, indent=2, sort_keys=True),
        encoding="utf-8",
    )

    write_markdown_report(
        df,
        feature_summary,
        tests,
        coverage,
        coverage_summary,
        source_audit,
        args.outdir / "audit_report.md",
    )

    print("=" * 72)
    print("ORIGINAL TESS DATASET / PIPELINE AUDIT")
    print("=" * 72)
    print(f"Rows: {len(df)}")
    print(f"Unique TICs: {df['tic_id'].nunique()}")
    print(f"Class counts: {df['label'].value_counts().sort_index().to_dict()}")
    print(f"Duplicate TIC membership rows: {df['tic_id'].duplicated(keep=False).sum()}")
    print(f"Cross-label TIC conflicts: {len(conflicts)}")
    if "n_points" in df.columns:
        for label, name in [(0, "negative"), (1, "positive")]:
            x = df.loc[df["label"] == label, "n_points"]
            print(
                f"n_points {name}: median={x.median():.1f}, "
                f"mean={x.mean():.1f}, IQR=({x.quantile(.25):.1f}, {x.quantile(.75):.1f})"
            )
    print(f"Processed curves audited: {len(coverage)}")
    print(f"Source asymmetries written to: {args.outdir / 'source_code_audit.json'}")
    print(f"Full report: {args.outdir / 'audit_report.md'}")
    print("=" * 72)


if __name__ == "__main__":
    main()
