#!/usr/bin/env python3
"""Pre-specified validation of the expanded, paired TESS injection study.

Runs AFTER revision/06_injection_recovery.py --tag full. The 12 pilot source
stars are excluded from the primary validation cohort to avoid reporting
objects already inspected during iterative method development.

All splitting and bootstrap resampling use SOURCE TIC as the independent unit.
The threshold is calibrated ONLY on training control stars. Reported 95% CIs
condition on that learned threshold; they do not incorporate calibration-stage
or target-selection uncertainty.

Run:
    python -u revision/08_validate_expanded_injections.py --bootstrap 3000
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import ks_2samp
from sklearn.model_selection import train_test_split

ROOT = Path(__file__).resolve().parents[1]
FULL = ROOT / "results/injection_full"
PILOT = ROOT / "results/injection_pilot"
OUT = ROOT / "results/injection_validation"

SEED = 42
TEST_SIZE = 0.30
BOOT_SEED = 24011027
TARGET_FPRS = (0.01, 0.05, 0.10)
PRIMARY_TARGET_FPR = 0.05


def validate_input(truth, features, pilot):
    if truth.empty or features.empty or pilot.empty:
        raise ValueError("Empty full/pilot injection results")
    for file_name, df in [
        ("full truth", truth),
        ("full paired features", features),
        ("pilot truth", pilot),
    ]:
        if "tic_id" not in df.columns:
            raise ValueError(f"{file_name} lacks tic_id")
    if truth.pair_id.duplicated().any():
        raise ValueError("Duplicated injection pair IDs")
    if features.duplicated(["pair_id", "label"]).any():
        raise ValueError("Repeated label within one pair")
    observed = features.groupby("pair_id").label.agg(["size", "sum"])
    if not ((observed["size"] == 2) & (observed["sum"] == 1)).all():
        raise ValueError("Each injection must have one control and one injected curve")
    if set(features.pair_id) != set(truth.pair_id):
        raise ValueError("Feature-pair set differs from injection-truth pair set")
    labels = features.groupby("pair_id").tic_id.nunique()
    if not labels.eq(1).all():
        raise ValueError("One pair contains more than one source star")
    control = features[features.label == 0].set_index("pair_id")
    injected = features[features.label == 1].set_index("pair_id")
    for col in ("n_points", "baseline_days", "noise_sigma_baseline"):
        if not np.allclose(
            control.loc[truth.pair_id, col].to_numpy(dtype=float),
            injected.loc[truth.pair_id, col].to_numpy(dtype=float),
            rtol=0, atol=1e-9,
        ):
            raise ValueError(f"Control/injection not matched on {col}")
    if not set(features.label.unique()) == {0, 1}:
        raise ValueError("Both binary classes required")


def build_stars(truth, features):
    control = features[features.label == 0][
        ["pair_id", "tic_id", "bls_depth_snr", "n_points", "baseline_days"]
    ].rename(columns={"bls_depth_snr": "control_snr"})
    injected = features[features.label == 1][
        ["pair_id", "tic_id", "bls_depth_snr"]
    ].rename(columns={"bls_depth_snr": "injected_snr"})
    records = truth.merge(control, on=["pair_id", "tic_id"],
                          how="inner", validate="one_to_one").merge(
        injected, on=["pair_id", "tic_id"], how="inner", validate="one_to_one"
    )
    if len(records) != len(truth):
        raise ValueError("Incomplete pair-feature join")
    if not np.isfinite(records[[
        "control_snr", "injected_snr", "true_period_days", "radius_ratio"
    ]].to_numpy(float)).all():
        raise ValueError("Non-finite validation variables")
    if (records["control_snr"] < 0).any() or (records["injected_snr"] < 0).any():
        raise ValueError("Negative BLS depth_snr")
    return records


def get_unique_controls(records):
    # If one TIC has multiple injections, the unmodified control must yield
    # exactly the same SNR in every pair (not independent controls).
    grouped = records.groupby("tic_id")["control_snr"]
    span = grouped.max() - grouped.min()
    if (span > 1e-7).any():
        first = span[span > 1e-7].index[:5].tolist()
        raise ValueError("Control features drift across injections: " + str(first))
    return grouped.first().sort_index()


def interval_by_cluster(test_records, test_controls, threshold, n_boot, seed):
    rng = np.random.default_rng(seed)
    by_star = {
        str(tic): sub for tic, sub in test_records.groupby("tic_id")
    }
    tics = np.array(sorted(by_star))
    ctrl = test_controls.to_dict()

    # Metric definitions averaged equally over independent source stars.
    def measure(stars):
        false_alarms = []
        detections = []
        direct = []
        alias_inclusive = []
        joint = []
        delta_positive = []
        for tic in stars:
            tic = str(tic)
            block = by_star[tic]
            false_alarms.append(float(ctrl[tic] >= threshold))
            pred = block.injected_snr.to_numpy(float) >= threshold
            recovered = block.period_recovered_direct.to_numpy(bool)
            harmonic = block.period_recovered_alias.to_numpy(bool)
            detections.append(float(np.mean(pred)))
            direct.append(float(np.mean(recovered)))
            alias_inclusive.append(float(np.mean(recovered | harmonic)))
            joint.append(float(np.mean(pred & recovered)))
            delta_positive.append(float(np.mean(
                block.injected_snr.to_numpy(float) > 
                block.control_snr.to_numpy(float)
            )))
        return {
            "control_fpr": float(np.mean(false_alarms)),
            "injected_detection": float(np.mean(detections)),
            "direct_period_recovery": float(np.mean(direct)),
            "period_or_alias_recovery": float(np.mean(alias_inclusive)),
            "joint_detection_and_direct_period": float(np.mean(joint)),
            "paired_snr_improvement_rate": float(np.mean(delta_positive)),
        }

    point = measure(tics)
    boot = {key: [] for key in point}
    for _ in range(n_boot):
        draw = rng.choice(tics, size=len(tics), replace=True)
        value = measure(draw)
        for key in boot:
            boot[key].append(value[key])

    rows = []
    for name, v in boot.items():
        arr = np.asarray(v, dtype=float)
        rows.append({
            "metric": name,
            "point_estimate": point[name],
            "ci95_low": float(np.quantile(arr, 0.025)),
            "ci95_high": float(np.quantile(arr, 0.975)),
            "n_independent_test_stars": len(tics),
            "bootstrap_unit": "TIC",
            "n_cluster_bootstrap": n_boot,
        })
    return pd.DataFrame(rows)


def summarize_completeness(records, threshold):
    data = records.copy()
    data["detected"] = data.injected_snr >= threshold
    data["direct_recovered"] = data.period_recovered_direct.astype(bool)
    data["joint"] = data["detected"] & data["direct_recovered"]
    data["period_group"] = pd.cut(
        data.true_period_days,
        [1.5, 3.5, 6.5, 9.5],
        include_lowest=True, right=True,
    ).astype(str)
    data["radius_group"] = pd.cut(
        data.radius_ratio,
        [0.035, 0.07, 0.11, 0.15],
        include_lowest=True, right=True,
    ).astype(str)
    rows = []
    for col in ("period_group", "radius_group"):
        for category, sub in data.groupby(col, dropna=False):
            rows.append({
                "dimension": col,
                "bin": str(category),
                "n_injections": len(sub),
                "n_unique_source_stars": sub.tic_id.nunique(),
                "detection_rate": float(sub.detected.mean()),
                "direct_period_recovery": float(sub.direct_recovered.mean()),
                "joint_detection_period_rate": float(sub.joint.mean()),
            })
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bootstrap", type=int, default=3000)
    ap.add_argument("--outdir", type=Path, default=OUT)
    args = ap.parse_args()
    if args.bootstrap < 100:
        raise ValueError("At least 100 clustered bootstraps recommended")
    args.outdir.mkdir(parents=True, exist_ok=True)

    truth = pd.read_csv(FULL / "injection_truth_and_recovery.csv")
    features = pd.read_csv(FULL / "paired_features.csv")
    pilot = pd.read_csv(PILOT / "injection_truth_and_recovery.csv")
    validate_input(truth, features, pilot)

    pilot_tics = set(pilot.tic_id.astype(str))
    raw_stars = set(truth.tic_id.astype(str))
    nonpilot = raw_stars - pilot_tics
    if len(nonpilot) < 40:
        raise ValueError(
            f"Only {len(nonpilot)} uninspected full-cohort TICs; "
            "cannot report meaningful expanded validation"
        )
    selected_truth = truth[
        ~truth.tic_id.astype(str).isin(pilot_tics)
    ].copy()
    records = build_stars(selected_truth, features)
    controls = get_unique_controls(records)
    star_ids = sorted(controls.index.astype(str))
    train_ids, test_ids = train_test_split(
        star_ids, test_size=TEST_SIZE, random_state=SEED
    )
    train_ids, test_ids = set(train_ids), set(test_ids)
    split = pd.DataFrame({
        "tic_id": star_ids,
        "split": [
            "calibration" if t in train_ids else "evaluation"
            for t in star_ids
        ],
        "was_in_pilot": [False] * len(star_ids),
    })
    split.to_csv(args.outdir / "source_star_split.csv", index=False)

    calib_ctrl = controls.loc[sorted(train_ids)]
    eval_ctrl = controls.loc[sorted(test_ids)]
    test_records = records[
        records.tic_id.astype(str).isin(test_ids)
    ].copy()
    # Leakage checks.
    if set(test_records.tic_id.astype(str)) & train_ids:
        raise AssertionError("Star-level split leakage")
    if set(test_records.tic_id.astype(str)) & pilot_tics:
        raise AssertionError("Pilot star included in expanded evaluation")

    per_threshold = []
    primary_ci = None
    primary_bins = None
    for nominal_fpr in TARGET_FPRS:
        threshold = float(np.quantile(calib_ctrl.values, 1 - nominal_fpr))
        ci = interval_by_cluster(
            test_records, eval_ctrl, threshold,
            args.bootstrap, BOOT_SEED + int(nominal_fpr*1000)
        )
        ci.insert(0, "nominal_control_fpr", nominal_fpr)
        ci.insert(1, "calibrated_snr_threshold", threshold)
        per_threshold.append(ci)
        if abs(nominal_fpr - PRIMARY_TARGET_FPR) < 1e-10:
            primary_ci = ci
            primary_bins = summarize_completeness(test_records, threshold)
    metrics = pd.concat(per_threshold, ignore_index=True)
    metrics.to_csv(args.outdir / "cluster_bootstrap_results.csv", index=False)
    primary_bins.to_csv(args.outdir / "period_radius_completeness.csv",
                        index=False)

    source_diag = test_records.copy()
    source_diag["control_fpr_flag_5pct"] = (
        source_diag.control_snr >= float(primary_ci.calibrated_snr_threshold.iloc[0])
    )
    source_diag["injection_detected_5pct"] = (
        source_diag.injected_snr >= float(primary_ci.calibrated_snr_threshold.iloc[0])
    )
    source_diag.to_csv(args.outdir / "evaluation_pairs.csv", index=False)

    config = {
        "full_input": str(FULL),
        "pilot_input": str(PILOT),
        "original_tics_in_full": len(raw_stars),
        "pilot_tics_excluded": len(raw_stars & pilot_tics),
        "n_remaining_tics": len(nonpilot),
        "n_calibration_tics": len(train_ids),
        "n_evaluation_tics": len(test_ids),
        "n_evaluation_injections": len(test_records),
        "star_split_seed": SEED,
        "nominal_false_positive_targets": list(TARGET_FPRS),
        "primary_nominal_fpr": PRIMARY_TARGET_FPR,
        "bootstrap_clusters": args.bootstrap,
        "bootstrap_unit": "source TIC",
        "mask_and_photometry_match_in_each_pair": True,
        "precautions": [
            "Pilot TICs were inspected during development, excluded from primary analysis.",
            "Catalogue negatives are not proven non-transiting stars.",
            "Threshold learned using CALIBRATION control TICs only.",
            "CIs reflect sampling of EVALUATION stars with fixed threshold.",
            "Per-source replicated injections are not independent.",
            "Assumed solar stellar parameters; injection completeness is conditional.",
            "Small-bin summaries are descriptive, not confidence intervals.",
            "This evaluates BLS detectability and period recovery, not an RF exoplanet detector.",
        ],
    }
    (args.outdir / "validation_config.json").write_text(
        json.dumps(config, indent=2), encoding="utf-8"
    )
    report = [
        "# Expanded paired transit injection/recovery validation", "",
        f"Full cohort with successful injections: {len(raw_stars)} distinct source stars",
        f"Pilot stars excluded: {len(raw_stars & pilot_tics)}",
        f"Independent calibration stars: {len(train_ids)}",
        f"Independent evaluation stars: {len(test_ids)}",
        f"Evaluated injected trials: {len(test_records)}",
        "", "## Primary, calibrated 5%-nominal-FPR results (star-cluster CIs)",
        "",
        primary_ci.round(5).to_string(index=False),
        "", "## Sensitivity by fixed target-FPR levels (not test-tuned)",
        "",
        metrics.round(5).to_string(index=False),
        "", "## Completeness by injection parameters (descriptive)", "",
        primary_bins.round(5).to_string(index=False),
        "", "## Interpretation safeguards", "",
        *["- " + v for v in config["precautions"]], "",
    ]
    (args.outdir / "validation_report.md").write_text(
        "\n".join(report), encoding="utf-8"
    )
    print("\n" + "\n".join(report[:18]), flush=True)
    print(f"\nSaved to: {args.outdir}", flush=True)


if __name__ == "__main__":
    main()
