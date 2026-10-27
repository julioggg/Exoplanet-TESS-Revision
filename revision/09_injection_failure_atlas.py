#!/usr/bin/env python3
"""Post-hoc failure atlas and TWO-STAGE uncertainty for TESS injection study.

Input results are held fixed. Never uses evaluation targets to tune pipeline
hyperparameters or choose a new threshold. The 5%-nominal control threshold
was PRE-SPECIFIED. This audit quantifies the extra sampling uncertainty
introduced by threshold estimation using only 95 calibration stars.

Calibrator bootstrap samples source TICs with replacement, recalibrates the
nominal 95th-percentile threshold, and independently resamples EVALUATION TICs
(including their complete triples of injections) with replacement.

These are exploratory, not independent confirmatory results.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
TRUTH = ROOT / "results/injection_full/injection_truth_and_recovery.csv"
FEATURES = ROOT / "results/injection_full/paired_features.csv"
SPLIT = ROOT / "results/injection_validation/source_star_split.csv"
OUT = ROOT / "results/injection_failure_atlas"
SEED = 20261027
Q = 0.95


def load_data():
    truth = pd.read_csv(TRUTH)
    # The initial full-run split is not the independent, pilot-excluded split.
    # Use ONLY the separately frozen source-star split from experiment 08.
    truth = truth.drop(columns=["split"], errors="ignore")
    feats = pd.read_csv(FEATURES)
    split = pd.read_csv(SPLIT)
    if feats.duplicated(["pair_id", "label"]).any():
        raise ValueError("Repeated feature member")
    c = feats.loc[feats.label == 0, [
        "tic_id", "pair_id", "bls_depth_snr", "noise_sigma_baseline",
        "baseline_days", "n_points",
    ]].rename(columns={"bls_depth_snr": "control_snr"})
    i = feats.loc[feats.label == 1, [
        "tic_id", "pair_id", "bls_depth_snr",
    ]].rename(columns={"bls_depth_snr": "injected_snr"})
    df = truth.merge(c, on=["tic_id", "pair_id"], validate="one_to_one")
    df = df.merge(i, on=["tic_id", "pair_id"], validate="one_to_one")
    df = df.merge(split, on="tic_id", validate="many_to_one")
    if df.split.isna().any():
        raise ValueError("Null split: check pilot exclusion")
    if set(df.split) != {"calibration", "evaluation"}:
        raise ValueError("Split column unexpected")
    df["delta_snr"] = df.injected_snr - df.control_snr
    df["nominal_depth"] = 1.0 - df.raw_model_min_flux
    df["nominal_depth_to_mad"] = (
        df.nominal_depth / df.noise_sigma_baseline
    )
    if not np.isfinite(df[[
        "control_snr", "injected_snr", "nominal_depth_to_mad"
    ]].to_numpy(float)).all():
        raise ValueError("Unexpected nonfinite signal")
    for tic, group in df.groupby("tic_id"):
        if group.control_snr.max() - group.control_snr.min() > 1e-7:
            raise ValueError(f"Control SNR differs across injections for {tic}")
        if len(group) != 3:
            raise ValueError(f"Expected exactly 3 injections for TIC {tic}")
    return df


def to_stars(df):
    by = df.groupby("tic_id", sort=True)
    tics = sorted(by.groups.keys())
    snr_c = np.array([by.get_group(t).control_snr.iloc[0] for t in tics])
    snr_i = np.array([
        by.get_group(t).sort_values("pair_id").injected_snr.to_numpy(float)
        for t in tics
    ])
    recovered = np.array([
        by.get_group(t).sort_values("pair_id").period_recovered_direct.to_numpy(bool)
        for t in tics
    ])
    ratios = np.array([
        by.get_group(t).sort_values("pair_id").radius_ratio.to_numpy(float)
        for t in tics
    ])
    return np.asarray(tics), snr_c, snr_i, recovered, ratios


def two_stage(calib_ctrl, evaluation, nboot, rng):
    _, c, _, _, _ = to_stars(calib_ctrl)
    tics, ec, ei, ep, er = to_stars(evaluation)
    threshold = float(np.quantile(c, Q))
    detected = ei >= threshold
    point = {
        "threshold": threshold,
        "false_positive_rate": float(np.mean(ec >= threshold)),
        "injected_detection_rate": float(np.mean(detected)),
        "joint_detection_and_direct_period_rate":
            float(np.mean(detected & ep)),
    }
    tr = np.empty((nboot, 4), dtype=float)
    for b in range(nboot):
        cc = c[rng.integers(0, len(c), len(c))]
        sampled_tics = rng.integers(0, len(tics), len(tics))
        q = float(np.quantile(cc, Q))
        ecb = ec[sampled_tics]
        eib = ei[sampled_tics]
        epb = ep[sampled_tics]
        hits = eib >= q
        tr[b, :] = [
            q,
            np.mean(ecb >= q),
            np.mean(hits),
            np.mean(hits & epb),
        ]
    metrics = list(point)
    results = []
    for index, name in enumerate(metrics):
        vals = tr[:, index]
        results.append({
            "metric": name,
            "point_estimate": point[name],
            "two_stage_ci95_low": float(np.quantile(vals, 0.025)),
            "two_stage_ci95_high": float(np.quantile(vals, 0.975)),
            "bootstrap_median": float(np.median(vals)),
            "n_calibration_stars": len(c),
            "n_evaluation_stars": len(tics),
            "bootstrap_unit": "TIC both stages",
            "n_bootstrap": nboot,
        })
    return pd.DataFrame(results), tr


def per_star_diagnostics(test, threshold):
    rows = []
    for tic, group in test.groupby("tic_id", sort=True):
        ctrl = float(group.control_snr.iloc[0])
        p = group.injected_snr.to_numpy(float)
        recovered = group.period_recovered_direct.to_numpy(bool)
        rows.append({
            "tic_id": tic,
            "n_injections": len(group),
            "control_snr": ctrl,
            "control_false_alarm": bool(ctrl >= threshold),
            "injected_snr_max": float(p.max()),
            "injected_snr_mean": float(p.mean()),
            "mean_delta_snr": float(group.delta_snr.mean()),
            "fraction_snr_improved": float((group.delta_snr > 0).mean()),
            "fraction_period_recovered": float(recovered.mean()),
            "detection_rate": float(np.mean(p >= threshold)),
            "joint_detection_period": float(np.mean(
                (p >= threshold) & recovered
            )),
            "median_nominal_depth": float(group.nominal_depth.median()),
            "median_depth_to_robust_sigma": float(
                group.nominal_depth_to_mad.median()
            ),
            "median_noise_sigma": float(
                group.noise_sigma_baseline.median()
            ),
        })
    return pd.DataFrame(rows)


def plot_figures(calib, test, results, draws, output):
    figures = output / "figures"
    figures.mkdir(exist_ok=True, parents=True)
    control = calib.groupby("tic_id").control_snr.first().to_numpy(float)
    evalctrl = test.groupby("tic_id").control_snr.first().to_numpy(float)
    threshold = float(results.loc[
        results.metric == "threshold", "point_estimate"
    ].iloc[0])
    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.hist(control, bins=30, alpha=.65, label="Calibration controls")
    ax.hist(evalctrl, bins=30, alpha=.5, label="Held-out controls")
    ax.axvline(threshold, linestyle="--", label="95th percentile (calibration)")
    ax.set(
        xlabel="Maximum BLS depth SNR (dimensionless white-noise proxy)",
        ylabel="Number of independent source stars",
        title="Control-score heavy tail and threshold calibration"
    )
    ax.legend()
    fig.tight_layout()
    fig.savefig(figures / "control_snr_calibration.png", dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(test.control_snr, test.injected_snr, s=16, alpha=.50)
    upper = float(np.nanpercentile(np.r_[
        test.control_snr.values, test.injected_snr.values], 99.5))
    upper = max(upper, threshold * 1.1)
    ax.plot([0, upper], [0, upper], "--", linewidth=1,
            label="No change in BLS SNR")
    ax.axhline(threshold, linestyle=":", label="Detection threshold")
    ax.set(xlabel="Paired original control BLS SNR",
           ylabel="Paired injected BLS SNR",
           title="BLS peak response to physically modeled transit")
    ax.set_xlim(0, upper)
    ax.set_ylim(0, upper)
    ax.legend()
    fig.tight_layout()
    fig.savefig(figures / "control_vs_injected_snr.png", dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.5, 5))
    ax.hist(draws[:, 0], bins=35, alpha=.8)
    ax.axvline(threshold, linestyle="--", label="Original calibration threshold")
    ax.set(xlabel="Bootstrap-reestimated 95th-percentile threshold",
           ylabel="Number of two-stage bootstrap draws",
           title="Threshold uncertainty across calibration star resamples")
    ax.legend()
    fig.tight_layout()
    fig.savefig(figures / "threshold_stability.png", dpi=160)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bootstrap", type=int, default=3000)
    ap.add_argument("--outdir", type=Path, default=OUT)
    args = ap.parse_args()
    if args.bootstrap < 100:
        raise ValueError("Choose at least 100 bootstrap replications")
    args.outdir.mkdir(parents=True, exist_ok=True)
    df = load_data()
    cal = df.loc[df.split == "calibration"].copy()
    eval_data = df.loc[df.split == "evaluation"].copy()
    rng = np.random.default_rng(SEED)
    intervals, draws = two_stage(cal, eval_data, args.bootstrap, rng)
    intervals.to_csv(args.outdir / "two_stage_cluster_intervals.csv",
                     index=False)
    threshold = float(intervals.loc[
        intervals.metric == "threshold", "point_estimate"
    ].iloc[0])
    table = per_star_diagnostics(eval_data, threshold)
    table.to_csv(args.outdir / "heldout_star_failures.csv", index=False)
    # High-control-SNR tail is an audit finding, not astrophysical truth.
    high = table.sort_values("control_snr", ascending=False).head(15)
    high.to_csv(args.outdir / "highest_control_snr_stars.csv", index=False)
    plot_figures(cal, eval_data, intervals, draws, args.outdir)
    comp = {
        "n_source_stars_calibration": int(cal.tic_id.nunique()),
        "n_source_stars_evaluation": int(eval_data.tic_id.nunique()),
        "n_injections_evaluation": int(len(eval_data)),
        "fraction_injected_max_BLS_SNR_higher_than_control": float(
            np.mean(eval_data.delta_snr > 0)
        ),
        "fraction_injected_period_recovered": float(
            np.mean(eval_data.period_recovered_direct.astype(bool))
        ),
        "report_is_post_hoc": True,
        "do_not_reoptimize_on_current_heldout_stars": True,
        "two_stage_bootstrap": (
            "Each replicate reestimates the 95th percentile on a "
            "source-star bootstrap of calibration controls, and "
            "independently resamples the heldout source stars (keeping "
            "all three injections together)."
        ),
    }
    (args.outdir / "audit_config.json").write_text(
        json.dumps(comp, indent=2), encoding="utf-8"
    )
    report = [
        "# TESS injection failure atlas and uncertainty audit", "",
        "This is a post-hoc descriptive analysis of the FIXED expanded sample.",
        "No thresholds or model parameters have been changed.", "",
        "## Two-stage bootstrap: threshold AND heldout TIC uncertainty", "",
        intervals.round(5).to_string(index=False), "",
        "## Baseline SNR tail: highest-scoring heldout control stars", "",
        high.round(5).to_string(index=False), "",
        "## Interpretive limits", "",
        "Nominal 5% calibration does not guarantee observed 5% false positives.",
        "Correlated stellar and instrumental noise invalidates simple white-noise significance interpretations.",
        "The baseline source stars are catalogue negatives, not confirmed transit-free objects.",
        "All exploratory findings here must be confirmed using new, independent targets.",
        "Never tune using these heldout stars and then reuse the same data as test.",
    ]
    (args.outdir / "audit_report.md").write_text(
        "\n".join(report) + "\n", encoding="utf-8"
    )
    print("\n" + "\n".join(report[:11]))
    print(f"\nSaved: {args.outdir}")


if __name__ == "__main__":
    main()
