#!/usr/bin/env python3
"""Paired transit injection/recovery in REAL TESS light curves.

Controls and injections derive from the SAME original negative-catalogue FITS
file. Transit signals are added to flux BEFORE shared quality/outlier masking,
normalization, 30-minute binning and BLS. The outlier mask is computed only
from the unmodified base flux and applied identically to both members of a
pair, preventing label-dependent cadence removal.

PyTransit RoadRunner is used with quadratic limb darkening and finite-exposure
supersampling. Solar-mass/radius stellar parameters are explicitly ASSUMED,
not inferred; this is a controlled detectability test, not a physical planet
population study.

Run smoke/pilot:
  python -u revision/06_injection_recovery.py --max-stars 12 \
    --injections-per-star 1 --workers 1 --tag pilot

Run full:
  python -u revision/06_injection_recovery.py --max-stars 0 \
    --injections-per-star 3 --workers 4 --tag full
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

# Prevent oversubscribed Numba workers in the CEDIA process pool.
os.environ.setdefault("NUMBA_NUM_THREADS", "1")

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.stats import sigma_clip
from astropy.timeseries import BoxLeastSquares
from sklearn.model_selection import train_test_split

ROOT = Path(__file__).resolve().parents[1]
ORIGINAL = (
    Path.home()
    / "exoplanet_revision/scratch/original/Aticulo CVV G3/"
    "exo_moo_tess_data_module"
)
RAW = ORIGINAL / "data/raw/negatives"
DATASET = ROOT / "data/dataset/final_dataset.csv"
CADENCE_DAYS = 30.0 / 1440.0
WINDOW_DAYS = 27.0
PERIODS = np.linspace(0.5, 15.0, 4000)
DURATIONS = np.array([0.05, 0.10, 0.15, 0.20, 0.25])
MIN_POINTS = 300
SEED = 20261027


def seed_for_tic(tic: str) -> int:
    digest = hashlib.sha256(f"{SEED}|{tic}".encode()).digest()
    return int.from_bytes(digest[:8], "little")


def load_original_curve(path: Path):
    with fits.open(path, memmap=False) as h:
        if len(h) < 2:
            raise ValueError("FITS has no light-curve extension")
        dat = h[1].data
        for col in ("TIME", "FLUX"):
            if col not in dat.names:
                raise ValueError(f"FITS is missing {col}")
        t = np.asarray(dat["TIME"], dtype=np.float64)
        f = np.asarray(dat["FLUX"], dtype=np.float64)
        q = (np.asarray(dat["QUALITY"]) == 0
             if "QUALITY" in dat.names else np.ones(len(t), dtype=bool))
        valid = np.isfinite(t) & np.isfinite(f) & (f > 0) & q
        if "FLUX_ERR" in dat.names:
            valid &= np.isfinite(np.asarray(dat["FLUX_ERR"], dtype=float))
    t, f = t[valid], f[valid]
    if len(t) < MIN_POINTS:
        raise ValueError(f"Only {len(t)} valid input cadences")
    order = np.argsort(t)
    t, f = t[order], f[order]
    within = t <= t[0] + WINDOW_DAYS
    t, f = t[within], f[within]
    if len(t) < MIN_POINTS:
        raise ValueError(f"Only {len(t)} input cadences in first 27 days")
    return t, f


def shared_preprocessing(t, raw_control, raw_injected):
    """Identical time, common baseline outlier mask and normalization."""
    baseline = float(np.median(raw_control))
    if baseline <= 0 or not np.isfinite(baseline):
        raise ValueError("Invalid baseline median")
    norm_base = raw_control / baseline
    clipped = sigma_clip(norm_base, sigma=5.0, maxiters=3)
    good = ~np.ma.getmaskarray(clipped)
    tt = t[good]
    c = raw_control[good] / baseline
    s = raw_injected[good] / baseline
    if len(tt) < MIN_POINTS:
        raise ValueError("Too few observations after shared clipping")

    # The binning origin comes from the first common observation, not labels.
    bins = np.floor((tt - tt[0]) / CADENCE_DAYS + 1e-12).astype("int64")
    frame = pd.DataFrame({"bin": bins, "time": tt, "base": c, "inj": s})
    b = frame.groupby("bin", sort=True).agg(
        time=("time", "median"),
        base=("base", "median"),
        inj=("inj", "median"),
        nraw=("base", "size"),
    )
    if len(b) < MIN_POINTS:
        raise ValueError(f"Only {len(b)} common 30-min bins")
    bt = b["time"].to_numpy(float)
    ctrl = b["base"].to_numpy(float)
    injected = b["inj"].to_numpy(float)
    if not (np.isfinite(ctrl).all() and np.isfinite(injected).all()):
        raise ValueError("Nonfinite binned flux")
    return bt, ctrl, injected, int(len(tt))


def transit_model(t, p, k, impact, phase):
    """Supersampled circular-orbit transit with assumed solar-like host."""
    from pytransit import RoadRunnerModel
    # Kepler's third law for M*=1 Msun and R*=1 Rsun:
    # a/Rstar ~ (1 AU/Rsun)*(P/365.25 days)**(2/3)
    a_rs = 215.032 * (p / 365.25)**(2.0 / 3.0)
    if a_rs <= impact:
        raise ValueError("Impossible impact parameter")
    inclination = float(np.arccos(impact / a_rs))
    dtime = np.diff(np.unique(t))
    dtime = dtime[dtime > 0]
    exptime = (
        float(min(np.median(dtime), CADENCE_DAYS))
        if len(dtime) else CADENCE_DAYS
    )
    tm = RoadRunnerModel("quadratic")
    tm.set_data(t, exptimes=exptime, nsamples=7)
    t0 = float(t[0] + phase * p)
    flux = np.asarray(
        tm.evaluate(
            k=float(k), ldc=[0.3, 0.2], t0=t0, p=float(p),
            a=float(a_rs), i=inclination
        ), dtype=float
    ).reshape(-1)
    if len(flux) != len(t) or not np.isfinite(flux).all():
        raise ValueError("Transit model yielded invalid flux array")
    if np.min(flux) >= 1.0:
        raise ValueError("No in-transit samples in source window")
    return flux, {
        "true_period_days": float(p),
        "radius_ratio": float(k),
        "impact_parameter": float(impact),
        "t0": t0,
        "a_over_rstar": float(a_rs),
        "inclination_rad": inclination,
        "exptime_days": exptime,
        "raw_model_min_flux": float(np.min(flux)),
    }


def bls_features(time, flux, base_noise):
    bls = BoxLeastSquares(
        time, flux, dy=np.full_like(time, base_noise, dtype=float)
    )
    scan = bls.power(PERIODS, DURATIONS)
    pw = np.asarray(scan.power, dtype=float)
    idx = int(np.nanargmax(pw))
    return {
        "mean_flux": float(np.mean(flux)),
        "median_flux": float(np.median(flux)),
        "std_flux": float(np.std(flux, ddof=0)),
        "mad_flux": float(np.median(np.abs(flux - np.median(flux)))),
        "amplitude": float(np.ptp(flux)),
        "n_points": int(len(flux)),
        "best_period": float(np.asarray(scan.period)[idx]),
        "transit_duration": float(np.asarray(scan.duration)[idx]),
        "transit_depth": float(np.asarray(scan.depth)[idx]),
        "bls_power": float(pw[idx]),
        "bls_depth_snr": float(np.asarray(scan.depth_snr)[idx]),
        "noise_sigma_baseline": float(base_noise),
        "baseline_days": float(time[-1] - time[0]),
    }


def period_recovery(measured, true):
    """Strict period match or integer 2:1 alias (both reported separately)."""
    ratio = measured / true
    direct = bool(abs(ratio - 1.0) <= 0.02)
    harmonic = bool(
        abs(ratio - 0.5) <= 0.01 or abs(ratio - 2.0) <= 0.04
    )
    return direct, harmonic


def one_source(task):
    tic, n_injections, rawdir = task
    out_features, out_truth, errors = [], [], []
    try:
        t, base = load_original_curve(Path(rawdir) / f"{tic}.fits")
        rng = np.random.default_rng(seed_for_tic(tic))
        # Base light curve and its outlier rejection depend only on original data.
        # All injected and control counterparts use the SAME common masks.
        for attempt in range(n_injections):
            pair_id = f"{tic}__inj{attempt:02d}"
            try:
                period = float(rng.uniform(1.5, 9.5))
                k = float(rng.uniform(0.035, 0.15))
                impact = float(rng.uniform(0.0, 0.80))
                phase = float(rng.uniform(0.0, 1.0))
                model, truth = transit_model(t, period, k, impact, phase)
                bt, ctrl, inj, n_after = shared_preprocessing(t, base, base * model)
                if not np.array_equal(np.isfinite(ctrl), np.isfinite(inj)):
                    raise AssertionError("Unequal finite observation masks")
                diff = float(np.max(np.abs(inj - ctrl)))
                if diff < 1e-7:
                    raise ValueError("No measurable injected signal remains")
                sigma = float(
                    1.4826 * np.median(np.abs(ctrl - np.median(ctrl)))
                )
                if not np.isfinite(sigma) or sigma <= 0:
                    sigma = float(np.std(ctrl, ddof=0))
                if not np.isfinite(sigma) or sigma <= 0:
                    raise ValueError("Invalid common photometric noise scale")
                cfeat = bls_features(bt, ctrl, sigma)
                sfeat = bls_features(bt, inj, sigma)
                matched, alias = period_recovery(
                    sfeat["best_period"], period
                )
                # Twin records: a repeated control per injection (independent
                # rows are NOT claimed; all evaluation must group by TIC).
                for label, features in [(0, cfeat), (1, sfeat)]:
                    out_features.append({
                        "tic_id": tic, "pair_id": pair_id,
                        "label": label, **features,
                    })
                out_truth.append({
                    "tic_id": tic, "pair_id": pair_id, **truth,
                    "n_native_cadences": int(len(t)),
                    "n_after_shared_clipping": n_after,
                    "n_binned_cadences": int(len(bt)),
                    "max_injected_flux_difference": diff,
                    "control_bls_snr": cfeat["bls_depth_snr"],
                    "injected_bls_snr": sfeat["bls_depth_snr"],
                    "delta_bls_snr": (
                        sfeat["bls_depth_snr"] - cfeat["bls_depth_snr"]
                    ),
                    "control_best_period": cfeat["best_period"],
                    "injected_best_period": sfeat["best_period"],
                    "period_recovered_direct": matched,
                    "period_recovered_alias": alias,
                })
            except Exception as ex:
                errors.append({
                    "tic_id": tic, "pair_id": pair_id,
                    "stage": "injection_or_extraction",
                    "error": repr(ex),
                })
    except Exception as ex:
        errors.append({
            "tic_id": tic, "pair_id": "",
            "stage": "load_raw_fits", "error": repr(ex),
        })
    return out_features, out_truth, errors


def summarize(frame, truth, n_stars, outdir):
    if frame.empty or truth.empty:
        return {"status": "no_successful_pairs"}
    grouped = frame.groupby("pair_id")["label"].agg(["count", "sum"])
    if not ((grouped["count"] == 2) & (grouped["sum"] == 1)).all():
        raise AssertionError("Broken paired dataset")
    star_ids = sorted(frame.tic_id.unique())
    train, test = train_test_split(
        star_ids, test_size=0.30, random_state=42
    )
    assignment = {tic: ("train" if tic in set(train) else "test")
                  for tic in star_ids}
    frame["split"] = frame["tic_id"].map(assignment)
    truth["split"] = truth["tic_id"].map(assignment)
    # Threshold calibrated using DISTINCT control stars; repeated control
    # rows for the same star must not inflate calibration sample sizes.
    controls = frame[frame.label == 0].drop_duplicates("tic_id")
    calibration = controls[controls.split == "train"]
    held_ctrl = controls[controls.split == "test"]
    held_inj = frame[(frame.label == 1) & (frame.split == "test")]
    threshold = float(np.quantile(calibration.bls_depth_snr, 0.95))
    fpr = float((held_ctrl.bls_depth_snr >= threshold).mean())
    tpr = float((held_inj.bls_depth_snr >= threshold).mean())
    t_truth = truth[truth.split == "test"]
    recovery = float(t_truth["period_recovered_direct"].mean())
    rec_alias = float(
        (t_truth["period_recovered_direct"] |
         t_truth["period_recovered_alias"]).mean()
    )
    detections = t_truth["injected_bls_snr"] >= threshold
    direct_and_detect = (
        detections & t_truth["period_recovered_direct"]
    )
    stats = {
        "status": "completed",
        "n_original_stars_requested": n_stars,
        "n_unique_stars_usable": len(star_ids),
        "n_pairs": len(truth),
        "n_star_train": len(train),
        "n_star_test": len(test),
        "n_test_injections": int(len(t_truth)),
        "n_test_control_stars": int(len(held_ctrl)),
        "threshold_calibrated_on_training_controls_95pct": threshold,
        "heldout_control_false_positive_rate": fpr,
        "heldout_injected_detection_rate": tpr,
        "heldout_direct_period_recovery_fraction": recovery,
        "heldout_period_or_2to1_alias_recovery_fraction": rec_alias,
        "heldout_detection_and_direct_period_recovery_fraction": float(
            direct_and_detect.mean()
        ),
        "note": (
            "Detection uses training-calibrated 95th percentile of control "
            "BLS depth_snr. Period recovery requires relative error <=2%, "
            "separately allowing 2:1 aliases. These are descriptive point "
            "estimates; stars (not injections) are the independent units."
        ),
    }
    frame.to_csv(outdir / "paired_features.csv", index=False)
    truth.to_csv(outdir / "injection_truth_and_recovery.csv", index=False)
    return stats


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-stars", type=int, default=12)
    parser.add_argument("--injections-per-star", type=int, default=1)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--tag", choices=["pilot", "full"], default="pilot")
    args = parser.parse_args()
    if args.max_stars < 0 or args.injections_per_star < 1 or args.workers < 1:
        raise ValueError("Invalid number of stars, injections or workers")
    outdir = ROOT / "results" / f"injection_{args.tag}"
    outdir.mkdir(parents=True, exist_ok=True)
    catalogue = pd.read_csv(DATASET)
    target_tics = sorted(catalogue.loc[
        catalogue.label == 0, "tic_id"
    ].astype(str).unique())
    target_tics = [
        tic for tic in target_tics if (RAW / f"{tic}.fits").exists()
    ]
    if args.max_stars:
        target_tics = target_tics[:args.max_stars]
    if len(target_tics) < 5:
        raise ValueError("Too few baseline curves for grouped split")
    print(
        f"INJECTION/RECOVERY: {len(target_tics)} source stars, "
        f"{args.injections_per_star} injections/star, workers={args.workers}",
        flush=True
    )
    tasks = [
        (tic, args.injections_per_star, str(RAW)) for tic in target_tics
    ]
    rows, truths, failures = [], [], []
    if args.workers == 1:
        iterator = map(one_source, tasks)
        for k, (f, t, e) in enumerate(iterator, 1):
            rows.extend(f); truths.extend(t); failures.extend(e)
            print(f"[{k:03d}/{len(tasks):03d}] pairs={len(t)} "
                  f"errors={len(e)}", flush=True)
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            futures = {pool.submit(one_source, task): task[0]
                       for task in tasks}
            for k, completed in enumerate(as_completed(futures), 1):
                f, t, e = completed.result()
                rows.extend(f); truths.extend(t); failures.extend(e)
                print(f"[{k:03d}/{len(tasks):03d}] "
                      f"{futures[completed]} pairs={len(t)} "
                      f"errors={len(e)}", flush=True)
    frame = pd.DataFrame(rows)
    truth = pd.DataFrame(truths)
    failure_df = pd.DataFrame(
        failures, columns=["tic_id", "pair_id", "stage", "error"]
    )
    failure_df.to_csv(outdir / "failures.csv", index=False)
    if not frame.empty:
        frame.sort_values(["tic_id", "pair_id", "label"],
                          inplace=True, ignore_index=True)
        truth.sort_values(["tic_id", "pair_id"], inplace=True,
                          ignore_index=True)
    stats = summarize(frame, truth, len(target_tics), outdir)
    config = {
        "experiment": f"injection_{args.tag}",
        "raw_sources": str(RAW),
        "target_catalogue": str(DATASET),
        "n_stars_requested": len(target_tics),
        "n_injections_per_star_requested": args.injections_per_star,
        "n_workers": args.workers,
        "fixed_seed": SEED,
        "pytransit_model": "RoadRunnerModel('quadratic')",
        "limb_darkening": [0.3, 0.2],
        "assumed_host_mass_msun": 1.0,
        "assumed_host_radius_rsun": 1.0,
        "period_range_days": [1.5, 9.5],
        "radius_ratio_range": [0.035, 0.15],
        "impact_parameter_range": [0.0, 0.8],
        "exposure_supersamples": 7,
        "baseline_max_days": WINDOW_DAYS,
        "cadence_minutes": 30.0,
        "bls_period_grid": [0.5, 15.0, len(PERIODS)],
        "bls_durations_days": DURATIONS.tolist(),
        "shared_control_derived_clip_mask": True,
        "shared_control_derived_flux_normalization": True,
        "split_group": "tic_id",
        "limitations": [
            "The host-star mass/radius/limb darkening are assumed Solar-like.",
            "Catalogued negative sources may contain unknown real transits.",
            "The noise dy is a per-star MAD-based white-noise surrogate.",
            "This is a controlled injection diagnostic, not planet validation.",
            "Repeated controls/injections from one TIC are not independent."
        ],
        "summary": stats,
        "number_of_failed_pairs_or_files": len(failures),
    }
    (outdir / "config_and_summary.json").write_text(
        json.dumps(config, indent=2), encoding="utf-8"
    )
    lines = [
        f"# TESS paired transit injection/recovery ({args.tag})", "",
        f"Stars requested: {len(target_tics)}",
        f"Successful injection pairs: {len(truth)}",
        f"Failures: {len(failures)}",
        "", "## Held-out recovery metrics", "",
    ]
    lines += [f"- {key}: {value}" for key, value in stats.items()]
    lines += ["", "## Protocol notes", "",
        "A synthetic transit was multiplied into each original Lightkurve-",
        "exported raw negative-catalogue FITS before common preprocessing.",
        "The SAME base-derived clipping mask, normalization and time bins",
        "were applied to injected/control members of each pair.",
        "The split is by TIC; no source star crosses train/test.",
        "The detection threshold is calibrated on training controls only.",
        "Recovered period and thresholded detection are both reported.",
        "Only pilot feasibility should be inferred from small pilot results.",
    ]
    (outdir / "report.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )
    print("\n" + "\n".join(lines), flush=True)
    if len(truth) < 1:
        raise RuntimeError("No successful injections; inspect failures.csv")


if __name__ == "__main__":
    main()
