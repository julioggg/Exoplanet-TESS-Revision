#!/usr/bin/env python3
"""Forensic review of paired transit injection/recovery (diagnostic only).

Reconstructs each pilot pair from its original FITS and saved injection truth.
Compares BLIND BLS results with ORACLE searches at the KNOWN injection period,
without treating the latter as a deployable detector. Generates phased
control/injection plots to inspect intrinsic variability and competing dips.

The paired tests must group by TIC, not injection. None of these diagnostics
are intended as final estimates of exoplanet detection completeness.

Run from repository root:
  python -u revision/07_injection_forensics.py --tag pilot
"""
from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.timeseries import BoxLeastSquares

ROOT = Path(__file__).resolve().parents[1]
PREVIOUS = ROOT / "revision/06_injection_recovery.py"
ORIGINAL = (
    Path.home() / "exoplanet_revision/scratch/original/"
    "Aticulo CVV G3/exo_moo_tess_data_module/data/raw/negatives"
)


def load_pipeline():
    """Use the exact same preprocessing/transit model as the pilot."""
    spec = importlib.util.spec_from_file_location("injection_pipeline", PREVIOUS)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def get_true_period_bls(time, flux, period, sigma, durations):
    """ORACLE search: period is supplied from injection truth."""
    valid_durations = np.asarray(durations, dtype=float)
    valid_durations = valid_durations[valid_durations < period]
    if len(valid_durations) == 0:
        raise ValueError("No duration shorter than injected period")
    bls = BoxLeastSquares(
        time, flux, dy=np.full_like(time, sigma, dtype=float)
    )
    model = bls.power(np.array([period]), valid_durations)
    return {
        "oracle_period_power": float(np.asarray(model.power).reshape(-1)[0]),
        "oracle_period_depth_snr": float(np.asarray(model.depth_snr).reshape(-1)[0]),
        "oracle_period_duration_days": float(np.asarray(model.duration).reshape(-1)[0]),
        "oracle_period_transit_time": float(
            np.asarray(model.transit_time).reshape(-1)[0]
        )
    }


def phase_curve(time, flux, t0, period):
    return ((time - t0 + 0.5 * period) % period) / period - 0.5, flux


def folded_binned(phases, flux, nbins=65):
    boundaries = np.linspace(-0.5, 0.5, nbins + 1)
    ind = np.digitize(phases, boundaries) - 1
    xx, yy = [], []
    for b in range(nbins):
        values = flux[ind == b]
        if len(values) >= 2:
            xx.append(float((boundaries[b] + boundaries[b+1]) / 2))
            yy.append(float(np.median(values)))
    return np.asarray(xx), np.asarray(yy)


def make_plot(path, time, base, injected, truth, featctrl, featinj):
    t0 = float(truth["t0"])
    p = float(truth["true_period_days"])
    phase_c, _ = phase_curve(time, base, t0, p)
    phase_i, _ = phase_curve(time, injected, t0, p)
    xctrl, yctrl = folded_binned(phase_c, base)
    xinj, yinj = folded_binned(phase_i, injected)

    fig, axs = plt.subplots(2, 1, figsize=(10, 7), layout="constrained")
    axs[0].plot(time - time[0], base, ".", ms=2.0, alpha=0.60,
                label="Original control")
    axs[0].plot(time - time[0], injected, ".", ms=2.0, alpha=0.45,
                label="Same source + injection")
    axs[0].set(xlabel="Days since first valid observation",
               ylabel="Normalized flux",
               title=f"{truth['tic_id']} | P true={p:.4f} d")
    axs[0].legend(loc="best", fontsize=8)
    axs[1].plot(phase_c, base, ".", ms=2.0, alpha=.18, label="Control")
    axs[1].plot(phase_i, injected, ".", ms=2.0, alpha=.18, label="Injected")
    if len(xctrl):
        axs[1].plot(xctrl, yctrl, "-", lw=1.7)
    if len(xinj):
        axs[1].plot(xinj, yinj, "-", lw=1.7)
    axs[1].axvline(0, linewidth=.8, linestyle="--", color="black")
    axs[1].set(xlabel="Phase at KNOWN injected ephemeris",
               ylabel="Normalized flux",
               title=(
                   f"Blind best P control={featctrl['best_period']:.3f} d, "
                   f"injected={featinj['best_period']:.3f} d"
               ))
    axrng = np.percentile(
        np.concatenate([base, injected]), [0.3, 99.7]
    )
    if np.all(np.isfinite(axrng)) and axrng[1] > axrng[0]:
        pad = 0.1 * (axrng[1] - axrng[0])
        axs[1].set_ylim(axrng[0] - pad, axrng[1] + pad)
    fig.savefig(path, dpi=125)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--tag", choices=["pilot", "full"], default="pilot")
    p.add_argument("--limit", type=int, default=0,
                   help="0 examines all saved successful pairs")
    p.add_argument("--no-plots", action="store_true")
    args = p.parse_args()
    indir = ROOT / "results" / f"injection_{args.tag}"
    outdir = ROOT / "results" / f"injection_{args.tag}_forensics"
    outdir.mkdir(parents=True, exist_ok=True)
    plotdir = outdir / "figures"
    if not args.no_plots:
        plotdir.mkdir(parents=True, exist_ok=True)

    truth = pd.read_csv(indir / "injection_truth_and_recovery.csv")
    features = pd.read_csv(indir / "paired_features.csv")
    if truth.empty or features.empty:
        raise ValueError("No saved injection pairs to diagnose")
    if features.duplicated(["pair_id", "label"]).any():
        raise ValueError("Duplicate pair member detected")
    pipeline = load_pipeline()
    rows, fails = [], []
    selected = truth.head(args.limit) if args.limit > 0 else truth

    for k, trow in enumerate(selected.to_dict("records"), 1):
        tic = str(trow["tic_id"])
        pair = str(trow["pair_id"])
        try:
            matches = features[features["pair_id"] == pair]
            if len(matches) != 2 or set(matches["label"]) != {0, 1}:
                raise ValueError("Control/injection feature pair missing")
            fc = matches[matches.label == 0].iloc[0]
            fi = matches[matches.label == 1].iloc[0]
            time, flux = pipeline.load_original_curve(ORIGINAL / (tic + ".fits"))
            period = float(trow["true_period_days"])
            phase = (float(trow["t0"]) - float(time[0])) / period
            model, params = pipeline.transit_model(
                time, period, float(trow["radius_ratio"]),
                float(trow["impact_parameter"]), phase
            )
            tt, control, injected, n_after = pipeline.shared_preprocessing(
                time, flux, flux * model
            )
            if len(tt) != int(fc["n_points"]) or len(tt) != int(fi["n_points"]):
                raise AssertionError("Reconstructed number of points differs from pilot")
            raw_diff = float(np.max(np.abs(injected - control)))
            if abs(raw_diff - float(trow["max_injected_flux_difference"])) > 1e-7:
                raise AssertionError("Reconstructed injected signal differs from saved pilot")
            sigma = float(fc["noise_sigma_baseline"])
            oracle_ctrl = get_true_period_bls(
                tt, control, period, sigma, pipeline.DURATIONS
            )
            oracle_inj = get_true_period_bls(
                tt, injected, period, sigma, pipeline.DURATIONS
            )
            row = {
                "tic_id": tic, "pair_id": pair, "split": str(trow["split"]),
                "period_true": period,
                "period_blind_injected": float(fi["best_period"]),
                "period_blind_control": float(fc["best_period"]),
                "radius_ratio": float(trow["radius_ratio"]),
                "expected_max_transit_depth": (
                    1.0 - float(trow["raw_model_min_flux"])
                ),
                "intrinsic_control_flux_std": float(np.std(control, ddof=0)),
                "intrinsic_control_robust_sigma": sigma,
                "transit_depth_over_baseline_sigma": (
                    (1.0 - float(trow["raw_model_min_flux"])) / sigma
                ),
                "blind_control_depth_snr": float(fc["bls_depth_snr"]),
                "blind_injected_depth_snr": float(fi["bls_depth_snr"]),
                "blind_delta_depth_snr": (
                    float(fi["bls_depth_snr"]) - float(fc["bls_depth_snr"])
                ),
                "blind_direct_period_recovery": bool(
                    trow["period_recovered_direct"]
                ),
                "oracle_control_power": oracle_ctrl["oracle_period_power"],
                "oracle_injected_power": oracle_inj["oracle_period_power"],
                "oracle_delta_power": (
                    oracle_inj["oracle_period_power"]
                    - oracle_ctrl["oracle_period_power"]
                ),
                "oracle_control_depth_snr": oracle_ctrl["oracle_period_depth_snr"],
                "oracle_injected_depth_snr": oracle_inj["oracle_period_depth_snr"],
                "oracle_delta_depth_snr": (
                    oracle_inj["oracle_period_depth_snr"]
                    - oracle_ctrl["oracle_period_depth_snr"]
                ),
                "n_points_same_as_control": int(len(tt)),
                "control_flux_amplitude": float(np.ptp(control)),
            }
            rows.append(row)
            if not args.no_plots:
                make_plot(
                    plotdir / f"{tic}__{pair.split('__')[-1]}.png",
                    tt, control, injected, trow, fc, fi
                )
            print(f"[{k}/{len(selected)}] {tic}: blind={row['blind_direct_period_recovery']}, "
                  f"oracle Δpower={row['oracle_delta_power']:.2f}", flush=True)
        except Exception as exc:
            fails.append({
                "tic_id": tic, "pair_id": pair, "error": repr(exc)
            })
            print(f"[{k}/{len(selected)}] {tic}: FAILED {exc}", flush=True)

    pd.DataFrame(rows).to_csv(outdir / "per_injection_diagnostics.csv",
                              index=False)
    pd.DataFrame(
        fails, columns=["tic_id", "pair_id", "error"]
    ).to_csv(outdir / "diagnostic_failures.csv", index=False)
    if not rows:
        raise RuntimeError("All diagnostics failed")
    frame = pd.DataFrame(rows)
    report = [
        f"# Injection diagnostics: {args.tag}", "",
        f"Pairs audited: {len(rows)}",
        f"Failures: {len(fails)}",
        f"Direct blind period recovery: {int(frame.blind_direct_period_recovery.sum())}/{len(frame)}",
        f"Oracle period power increased after injection: "
        f"{int((frame.oracle_delta_power > 0).sum())}/{len(frame)}",
        "", "Important: 'oracle' means true period was supplied. This is only "
        "a diagnostic of potential detectability at a known ephemeris.",
        "High BLS depth_snr in control stars may indicate stellar variability, "
        "eclipses, systematics or correlated photometric noise; it is not proof "
        "that the catalogue-negative star contains a planet.",
        "", "## Per-injection results", "",
        frame[[
            "tic_id", "split", "period_true", "period_blind_injected",
            "blind_direct_period_recovery", "radius_ratio",
            "blind_control_depth_snr", "blind_injected_depth_snr",
            "oracle_delta_power", "oracle_delta_depth_snr",
        ]].round(4).to_string(index=False), ""
    ]
    (outdir / "diagnostic_report.md").write_text(
        "\n".join(report) + "\n", encoding="utf-8"
    )
    (outdir / "config.json").write_text(json.dumps({
        "tag": args.tag, "n_pairs": len(rows), "n_failures": len(fails),
        "input_truth": str(indir / "injection_truth_and_recovery.csv"),
        "raw_source": str(ORIGINAL),
        "source_pipeline": str(PREVIOUS),
        "blind_search_changed": False,
        "oracle_knows_true_period": True,
        "notes": [
            "This is a post hoc failure analysis, not independent validation.",
            "No tuning or sample selection should be based on held-out outcomes.",
            "The noise uses a constant MAD estimate, not full correlated-noise modeling."
        ]
    }, indent=2), encoding="utf-8")
    print("\n" + "\n".join(report[:10]), flush=True)


if __name__ == "__main__":
    main()
