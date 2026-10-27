#!/usr/bin/env python3
"""
Build a cadence-matched, symmetrically processed TESS feature table.

Scientific goal
---------------
Remove the strongest class-dependent observational confounding identified in the
submitted baseline:
- positives were predominantly ~2-min cadence;
- negatives were predominantly ~30-min cadence;
- some negatives contained stitched multi-sector baselines;
- positive/negative BLS grids and statistical definitions differed.

This script starts from the original processed light curves kept on CEDIA and
applies the SAME transformation to both classes:
1. finite TIME/FLUX_NORM selection;
2. optional truncation to the first 27 days;
3. 30-minute time binning by median flux;
4. re-normalization by the binned median;
5. identical 5-sigma clipping;
6. identical statistical feature definitions;
7. identical BLS period/duration grid.

Large cadence-matched curves are written OUTSIDE the Git repository by default.
Only compact feature/manifest outputs are written under results/harmonized/.

Run from repository root:
    python revision/02_build_harmonized_dataset.py

Recommended HPC run:
    python revision/02_build_harmonized_dataset.py --n-workers 16
"""

from __future__ import annotations

import argparse
import json
import math
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.stats import sigma_clip
from astropy.timeseries import BoxLeastSquares


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ORIGINAL_ROOT = (
    Path.home()
    / "exoplanet_revision/scratch/original/Aticulo CVV G3/exo_moo_tess_data_module"
)
DEFAULT_WORKDIR = Path.home() / "exoplanet_revision/data/harmonized_30min"
DEFAULT_OUTDIR = REPO_ROOT / "results/harmonized"

CADENCE_MINUTES = 30.0
CADENCE_DAYS = CADENCE_MINUTES / (24.0 * 60.0)
MAX_BASELINE_DAYS = 27.0
MIN_POINTS_AFTER_BINNING = 300

PERIOD_MIN_DAYS = 0.5
PERIOD_MAX_DAYS = 30.0
PERIOD_GRID_POINTS = 5000
DURATIONS_DAYS = np.array([0.05, 0.10, 0.15, 0.20, 0.25], dtype=float)

SIGMA_CLIP_SIGMA = 5.0
SIGMA_CLIP_MAXITERS = 3


def tic_normalize(value: str) -> str:
    s = str(value).strip().replace(" ", "_")
    if not s.startswith("TIC_"):
        if s.startswith("TIC"):
            s = s.replace("TIC", "TIC_", 1).replace("__", "_")
        else:
            s = f"TIC_{s}"
    return s


def safe_float(x):
    try:
        x = float(x)
        return x if np.isfinite(x) else np.nan
    except Exception:
        return np.nan


def basic_time_characteristics(time: np.ndarray) -> dict:
    t = np.asarray(time, dtype=float)
    t = t[np.isfinite(t)]
    t = np.sort(np.unique(t))
    if len(t) < 2:
        return {
            "n_points": int(len(t)),
            "baseline_days": 0.0,
            "median_cadence_minutes": np.nan,
        }
    dt = np.diff(t)
    dt = dt[np.isfinite(dt) & (dt > 0)]
    return {
        "n_points": int(len(t)),
        "baseline_days": float(t[-1] - t[0]),
        "median_cadence_minutes": float(np.median(dt) * 1440.0) if len(dt) else np.nan,
    }


def harmonize_curve(path: Path, max_baseline_days: float) -> tuple[pd.DataFrame, dict]:
    df = pd.read_csv(path)

    required = {"TIME", "FLUX_NORM"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path.name}: missing columns {sorted(missing)}")

    time = pd.to_numeric(df["TIME"], errors="coerce").to_numpy(dtype=float)
    flux = pd.to_numeric(df["FLUX_NORM"], errors="coerce").to_numpy(dtype=float)

    mask = np.isfinite(time) & np.isfinite(flux)
    time = time[mask]
    flux = flux[mask]

    if len(time) < 2:
        raise ValueError("fewer than 2 finite observations")

    order = np.argsort(time)
    time = time[order]
    flux = flux[order]

    before = basic_time_characteristics(time)

    # Remove multi-sector/extreme coverage asymmetry while preserving a standard
    # TESS-sector-length window.
    t0 = float(time[0])
    keep = time <= (t0 + max_baseline_days)
    time = time[keep]
    flux = flux[keep]

    if len(time) < 2:
        raise ValueError("fewer than 2 observations after baseline truncation")

    # Common 30-minute bins, anchored independently at the first valid cadence.
    bin_index = np.floor((time - time[0]) / CADENCE_DAYS + 1e-12).astype(np.int64)

    tmp = pd.DataFrame({
        "bin": bin_index,
        "time": time,
        "flux": flux,
    })

    binned = (
        tmp.groupby("bin", sort=True)
        .agg(TIME=("time", "median"), FLUX_NORM=("flux", "median"), N_RAW=("flux", "size"))
        .reset_index(drop=True)
    )

    binned = binned.replace([np.inf, -np.inf], np.nan).dropna()

    if len(binned) < MIN_POINTS_AFTER_BINNING:
        raise ValueError(
            f"only {len(binned)} points after common-cadence binning "
            f"(<{MIN_POINTS_AFTER_BINNING})"
        )

    median_flux = float(np.median(binned["FLUX_NORM"]))
    if not np.isfinite(median_flux) or median_flux == 0:
        raise ValueError("invalid binned median flux")

    binned["FLUX_NORM"] = binned["FLUX_NORM"] / median_flux

    clipped = sigma_clip(
        binned["FLUX_NORM"].to_numpy(dtype=float),
        sigma=SIGMA_CLIP_SIGMA,
        maxiters=SIGMA_CLIP_MAXITERS,
    )
    binned = binned.loc[~np.asarray(clipped.mask)].copy()
    binned = binned.reset_index(drop=True)

    if len(binned) < MIN_POINTS_AFTER_BINNING:
        raise ValueError(
            f"only {len(binned)} points after common clipping "
            f"(<{MIN_POINTS_AFTER_BINNING})"
        )

    after = basic_time_characteristics(binned["TIME"].to_numpy(dtype=float))

    manifest = {
        "source_file": str(path),
        "n_points_before": before["n_points"],
        "baseline_days_before": before["baseline_days"],
        "median_cadence_minutes_before": before["median_cadence_minutes"],
        "n_points_after": after["n_points"],
        "baseline_days_after": after["baseline_days"],
        "median_cadence_minutes_after": after["median_cadence_minutes"],
        "raw_points_represented_after_binning": int(binned["N_RAW"].sum()),
    }
    return binned, manifest


def extract_features(curve: pd.DataFrame) -> dict:
    time = curve["TIME"].to_numpy(dtype=float)
    flux = curve["FLUX_NORM"].to_numpy(dtype=float)

    mean_flux = float(np.mean(flux))
    median_flux = float(np.median(flux))
    std_flux = float(np.std(flux, ddof=0))
    mad_flux = float(np.median(np.abs(flux - median_flux)))
    amplitude = float(np.max(flux) - np.min(flux))

    if std_flux <= 0 or not np.isfinite(std_flux):
        raise ValueError("non-positive/invalid flux scatter")

    bls = BoxLeastSquares(time, flux)
    periods = np.linspace(
        PERIOD_MIN_DAYS,
        PERIOD_MAX_DAYS,
        PERIOD_GRID_POINTS,
        dtype=float,
    )
    result = bls.power(periods, DURATIONS_DAYS)

    power_arr = np.asarray(result.power, dtype=float)
    if not np.isfinite(power_arr).any():
        raise ValueError("BLS returned no finite power")

    idx = int(np.nanargmax(power_arr))

    best_period = float(np.asarray(result.period)[idx])
    best_duration = float(np.asarray(result.duration)[idx])
    best_depth = float(np.asarray(result.depth)[idx])
    best_power = float(np.asarray(result.power)[idx])

    legacy_snr = float(abs(best_depth) / std_flux)

    depth_snr = np.nan
    if hasattr(result, "depth_snr"):
        try:
            depth_snr = float(np.asarray(result.depth_snr)[idx])
        except Exception:
            depth_snr = np.nan

    # Prefer the BLS-reported depth S/N. Fall back explicitly to the legacy proxy
    # only if the installed Astropy version does not provide depth_snr.
    snr_estimate = depth_snr if np.isfinite(depth_snr) else legacy_snr
    snr_definition = (
        "astropy_bls_depth_snr"
        if np.isfinite(depth_snr)
        else "legacy_abs_depth_over_flux_std"
    )

    return {
        "mean_flux": mean_flux,
        "median_flux": median_flux,
        "std_flux": std_flux,
        "mad_flux": mad_flux,
        "amplitude": amplitude,
        "n_points": int(len(flux)),
        "best_period": best_period,
        "transit_duration": best_duration,
        "transit_depth": best_depth,
        "bls_power": best_power,
        "snr_estimate": snr_estimate,
        "bls_depth_snr": depth_snr,
        "legacy_snr_estimate": legacy_snr,
        "snr_definition": snr_definition,
    }


def process_one(task):
    label, path_str, workdir_str, max_baseline_days, save_curves = task
    path = Path(path_str)
    workdir = Path(workdir_str)
    tic_id = tic_normalize(path.stem.replace("_processed", ""))

    try:
        curve, manifest = harmonize_curve(path, max_baseline_days=max_baseline_days)
        features = extract_features(curve)

        curve_path = ""
        if save_curves:
            class_dir = workdir / ("positives" if label == 1 else "negatives")
            class_dir.mkdir(parents=True, exist_ok=True)
            out = class_dir / f"{tic_id}_harmonized.csv"
            curve.to_csv(out, index=False)
            curve_path = str(out)

        row = {
            "tic_id": tic_id,
            **features,
            "label": int(label),
        }
        manifest_row = {
            "tic_id": tic_id,
            "label": int(label),
            **manifest,
            "harmonized_curve_file": curve_path,
            "status": "ok",
            "error": "",
        }
        return row, manifest_row

    except Exception as exc:
        manifest_row = {
            "tic_id": tic_id,
            "label": int(label),
            "source_file": str(path),
            "status": "failed",
            "error": repr(exc),
        }
        return None, manifest_row


def class_summary(manifest: pd.DataFrame) -> pd.DataFrame:
    ok = manifest[manifest["status"] == "ok"].copy()
    metrics = [
        "n_points_before",
        "baseline_days_before",
        "median_cadence_minutes_before",
        "n_points_after",
        "baseline_days_after",
        "median_cadence_minutes_after",
    ]
    rows = []
    for label in [0, 1]:
        sub = ok[ok["label"] == label]
        for metric in metrics:
            if metric not in sub.columns:
                continue
            x = pd.to_numeric(sub[metric], errors="coerce").dropna()
            rows.append({
                "label": label,
                "metric": metric,
                "n": int(len(x)),
                "mean": float(x.mean()) if len(x) else np.nan,
                "std": float(x.std(ddof=1)) if len(x) > 1 else np.nan,
                "median": float(x.median()) if len(x) else np.nan,
                "q1": float(x.quantile(0.25)) if len(x) else np.nan,
                "q3": float(x.quantile(0.75)) if len(x) else np.nan,
                "min": float(x.min()) if len(x) else np.nan,
                "max": float(x.max()) if len(x) else np.nan,
            })
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--original-root", type=Path, default=DEFAULT_ORIGINAL_ROOT)
    parser.add_argument("--workdir", type=Path, default=DEFAULT_WORKDIR)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument(
        "--n-workers",
        type=int,
        default=min(16, max(1, os.cpu_count() or 1)),
    )
    parser.add_argument("--max-baseline-days", type=float, default=MAX_BASELINE_DAYS)
    parser.add_argument(
        "--no-save-curves",
        action="store_true",
        help="Do not save cadence-matched per-target curves outside the repository.",
    )
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    args.workdir.mkdir(parents=True, exist_ok=True)

    class_dirs = {
        0: args.original_root / "data/processed/negatives",
        1: args.original_root / "data/processed/positives",
    }

    tasks = []
    for label, directory in class_dirs.items():
        if not directory.exists():
            raise FileNotFoundError(f"Missing processed-curve directory: {directory}")
        for path in sorted(directory.glob("*_processed.csv")):
            tasks.append(
                (
                    label,
                    str(path),
                    str(args.workdir),
                    float(args.max_baseline_days),
                    not args.no_save_curves,
                )
            )

    print("=" * 80)
    print("SYMMETRIC 30-MIN HARMONIZATION + UNIFIED BLS")
    print("=" * 80)
    print(f"Targets discovered: {len(tasks)}")
    print(f"Workers: {args.n_workers}")
    print(f"Cadence: {CADENCE_MINUTES:.1f} min")
    print(f"Maximum baseline: {args.max_baseline_days:.1f} d")
    print(
        f"BLS grid: {PERIOD_MIN_DAYS}-{PERIOD_MAX_DAYS} d, "
        f"{PERIOD_GRID_POINTS} periods, durations={DURATIONS_DAYS.tolist()}"
    )
    print("=" * 80)

    feature_rows = []
    manifest_rows = []

    if args.n_workers == 1:
        iterator = map(process_one, tasks)
        for i, (row, manifest) in enumerate(iterator, start=1):
            manifest_rows.append(manifest)
            if row is not None:
                feature_rows.append(row)
            print(
                f"[{i:03d}/{len(tasks):03d}] "
                f"{manifest['tic_id']} label={manifest['label']} {manifest['status']}"
            )
    else:
        with ProcessPoolExecutor(max_workers=args.n_workers) as ex:
            futures = {ex.submit(process_one, task): task for task in tasks}
            completed = 0
            for fut in as_completed(futures):
                completed += 1
                row, manifest = fut.result()
                manifest_rows.append(manifest)
                if row is not None:
                    feature_rows.append(row)
                print(
                    f"[{completed:03d}/{len(tasks):03d}] "
                    f"{manifest['tic_id']} label={manifest['label']} {manifest['status']}"
                )

    features = pd.DataFrame(feature_rows)
    manifest = pd.DataFrame(manifest_rows)

    if not features.empty:
        features = features.sort_values(["label", "tic_id"]).reset_index(drop=True)
    if not manifest.empty:
        manifest = manifest.sort_values(["label", "tic_id"]).reset_index(drop=True)

    features.to_csv(args.outdir / "harmonized_features.csv", index=False)
    manifest.to_csv(args.outdir / "harmonization_manifest.csv", index=False)

    failures = manifest[manifest["status"] != "ok"].copy()
    failures.to_csv(args.outdir / "harmonization_failures.csv", index=False)

    summary = class_summary(manifest)
    summary.to_csv(args.outdir / "harmonization_summary_by_class.csv", index=False)

    config = {
        "scientific_status": "cadence-matched sensitivity dataset",
        "input_root": str(args.original_root),
        "curve_workdir": str(args.workdir),
        "cadence_minutes": CADENCE_MINUTES,
        "max_baseline_days": float(args.max_baseline_days),
        "min_points_after_binning": MIN_POINTS_AFTER_BINNING,
        "sigma_clip_sigma": SIGMA_CLIP_SIGMA,
        "sigma_clip_maxiters": SIGMA_CLIP_MAXITERS,
        "period_min_days": PERIOD_MIN_DAYS,
        "period_max_days": PERIOD_MAX_DAYS,
        "period_grid_points": PERIOD_GRID_POINTS,
        "durations_days": DURATIONS_DAYS.tolist(),
        "n_workers": int(args.n_workers),
        "n_targets_discovered": int(len(tasks)),
        "n_success": int((manifest["status"] == "ok").sum()) if len(manifest) else 0,
        "n_failed": int((manifest["status"] != "ok").sum()) if len(manifest) else 0,
        "n_positive_success": int(((manifest["status"] == "ok") & (manifest["label"] == 1)).sum()) if len(manifest) else 0,
        "n_negative_success": int(((manifest["status"] == "ok") & (manifest["label"] == 0)).sum()) if len(manifest) else 0,
        "snr_primary": "Astropy BoxLeastSquares result.depth_snr where available",
        "legacy_snr_retained": "abs(transit_depth)/std_flux",
        "guardrail": (
            "This harmonizes cadence, coverage truncation, clipping, statistical "
            "definitions and BLS grid. It does not by itself convert catalogue "
            "negatives into astrophysically verified true negatives."
        ),
    }
    (args.outdir / "harmonization_config.json").write_text(
        json.dumps(config, indent=2),
        encoding="utf-8",
    )

    report = []
    report.append("# Cadence-matched symmetric reconstruction")
    report.append("")
    report.append(
        f"Both classes were rebinned to {CADENCE_MINUTES:.0f}-minute cadence, "
        f"limited to the first {args.max_baseline_days:.0f} days, and processed "
        "with identical clipping, statistical definitions and BLS grids."
    )
    report.append("")
    report.append(f"- Targets discovered: {len(tasks)}")
    report.append(f"- Successful: {config['n_success']}")
    report.append(f"- Failed: {config['n_failed']}")
    report.append(f"- Positive successful: {config['n_positive_success']}")
    report.append(f"- Negative successful: {config['n_negative_success']}")
    report.append("")
    if len(summary):
        for metric in ["n_points_after", "baseline_days_after", "median_cadence_minutes_after"]:
            report.append(f"## {metric}")
            report.append("")
            for label, name in [(0, "negative"), (1, "positive")]:
                x = summary[(summary["label"] == label) & (summary["metric"] == metric)]
                if len(x):
                    r = x.iloc[0]
                    report.append(
                        f"- {name}: median={r['median']:.6g}, "
                        f"IQR=[{r['q1']:.6g}, {r['q3']:.6g}], "
                        f"range=[{r['min']:.6g}, {r['max']:.6g}], n={int(r['n'])}"
                    )
            report.append("")
    report.append(
        "The next step is to statistically test the post-harmonization observational "
        "distributions and then train controlled ablations on harmonized features."
    )
    report.append("")

    (args.outdir / "harmonization_report.md").write_text(
        "\n".join(report),
        encoding="utf-8",
    )

    print()
    print("=" * 80)
    print("HARMONIZATION COMPLETE")
    print("=" * 80)
    print(f"Successful: {config['n_success']}/{len(tasks)}")
    print(f"Failures: {config['n_failed']}")
    if len(summary):
        print()
        print(summary[
            summary["metric"].isin(
                ["n_points_after", "baseline_days_after", "median_cadence_minutes_after"]
            )
        ].to_string(index=False))
    print()
    print("Outputs:", args.outdir)
    print("Large harmonized curves:", args.workdir)
    print("=" * 80)


if __name__ == "__main__":
    main()
