#!/usr/bin/env python3
"""Frozen, label-blind preparation of original TGLC CAL_APER_FLUX curves.

This engineering pilot follows the failed PSF channel in one original FITS.
The policy is fixed now (before any NEW heldout TICs are selected):
  - ONLY original, SHA256-verified MAST TGLC products;
  - ONLY calibrated APERTURE flux (cal_aper_flux) for EVERY target;
  - common TESS default-quality mask, finite BTJD and finite positive aperture;
  - TGLC_flags inventoried but NOT filtered in the main stream;
  - no sigma clipping, no detrending, no per-class/sector exceptions;
  - one normalization by a target's accepted-flux median;
  - retain native BTJD; align sources by cadence_num only for audit.

The original, previously studied stars are ENGINEERING SAMPLES; no BLS,
train/test split, planet classification or detection performance is measured.

Use:
  python -u revision/14_build_uniform_aperture_curves.py

Later, supply an externally verified prospective manifest via --manifest.
The policy applies unchanged to that manifest.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits
from lightkurve.utils import TessQualityFlags

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / "results/mast_original_pilot/verified_original_fits.csv"
DEFAULT_OUT = ROOT / "results/uniform_aperture_pilot"
DEFAULT_PRIVATE = Path.home() / "exoplanet_revision/data/uniform_aperture_curves"
FLUX_COLUMN = "cal_aper_flux"
QUALITY_BITMASK = "default"
MIN_CADENCES = 300
POLICY_NAME = "TGLC_APERTURE_1800S_V1"


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for part in iter(lambda: f.read(1024 * 1024), b""):
            h.update(part)
    return h.hexdigest()


def stable_robust_sigma(a):
    a = np.asarray(a, dtype=float)
    return float(1.4826 * np.median(np.abs(a - np.median(a))))


def load_one(item):
    path = Path(str(item["local_path"])).expanduser()
    if not path.is_file():
        raise FileNotFoundError(path)
    original_hash = digest(path)
    if original_hash != str(item["sha256"]):
        raise ValueError("Original FITS SHA256 does not match MAST manifest")
    if str(item["author"]).upper() != "TGLC":
        raise ValueError("Non-TGLC product not supported by policy")
    if abs(float(item["exptime_seconds"]) - 1800.0) > 3:
        raise ValueError("Policy requires original nominal 1800-s exposure")

    with fits.open(path, memmap=False) as hdul:
        candidates = []
        for h in hdul:
            columns = getattr(h, "columns", None)
            names = list(columns.names or []) if columns is not None else []
            norm = {str(n).strip().lower(): str(n) for n in names}
            required = {
                "time", "cadence_num", "tess_flags", FLUX_COLUMN
            }
            if required.issubset(norm):
                candidates.append((h, norm))
        if len(candidates) != 1:
            raise ValueError(f"Expected 1 TGLC time series, found {len(candidates)}")
        ext, col = candidates[0]

        for h in hdul:
            if "SECTOR" in h.header and int(h.header["SECTOR"]) != int(item["sector"]):
                raise ValueError("Sector mismatch in original FITS header")
            if "TICID" in h.header and str(int(h.header["TICID"])) != "".join(
                c for c in str(item["tic_id"]) if c.isdigit()
            ):
                raise ValueError("TIC mismatch in original FITS header")

        data = ext.data
        t = np.asarray(data[col["time"]], dtype=np.float64)
        aperture = np.asarray(data[col[FLUX_COLUMN]], dtype=np.float64)
        cadence = np.asarray(data[col["cadence_num"]], dtype=np.int64)
        quality = np.asarray(data[col["tess_flags"]], dtype=np.int64)
        tglc_flag = (
            np.asarray(data[col["tglc_flags"]], dtype=np.int64)
            if "tglc_flags" in col else np.zeros(len(t), dtype=np.int64)
        )

    if len({len(t), len(aperture), len(cadence), len(quality), len(tglc_flag)}) != 1:
        raise ValueError("Photometry columns have unequal row counts")
    t_good = np.isfinite(t)
    f_good = np.isfinite(aperture) & (aperture > 0)
    q_good = TessQualityFlags.create_quality_mask(
        quality, bitmask=QUALITY_BITMASK
    )
    accepted = t_good & f_good & q_good
    if int(accepted.sum()) < MIN_CADENCES:
        raise ValueError(
            f"Only {accepted.sum()} aperture-valid good cadences out of {len(t)}"
        )
    idx = np.flatnonzero(accepted)
    idx = idx[np.argsort(t[idx], kind="mergesort")]
    time = t[idx]
    if not (np.diff(time) > 0).all():
        raise ValueError("Accepted BTJD is not strictly increasing")
    if len(np.unique(cadence[idx])) != len(idx):
        raise ValueError("Duplicate cadence_num after quality filtering")

    value = aperture[idx]
    median = float(np.median(value))
    if median <= 0 or not np.isfinite(median):
        raise ValueError("Invalid normalized flux reference")
    normalized = value / median
    if not np.isfinite(normalized).all():
        raise ValueError("Normalized flux is nonfinite")
    curve = pd.DataFrame({
        "time_btjd": time,
        "cadence_num": cadence[idx],
        "flux_cal_aper_normalized": normalized,
        "tess_flags": quality[idx],
        "tglc_flags": tglc_flag[idx],
    })
    delta = np.diff(time)
    p2p = stable_robust_sigma(np.diff(normalized))/np.sqrt(2.0)
    info = {
        "tic_id": str(item["tic_id"]),
        "sector": int(item["sector"]),
        "pair_id": str(item.get("pair_id", "")),
        "n_raw": int(len(t)),
        "n_good_TESS_default": int(q_good.sum()),
        "n_positive_finite_aperture": int(f_good.sum()),
        "n_accepted": len(curve),
        "fraction_accepted": float(accepted.mean()),
        "n_TGLC_flag_nonzero_raw": int((tglc_flag != 0).sum()),
        "n_TGLC_flag_nonzero_accepted": int((tglc_flag[idx] != 0).sum()),
        "n_TGLC_strict_zero_if_sensitivity": int((tglc_flag[idx] == 0).sum()),
        "median_cadence_min": float(np.median(delta) * 1440),
        "baseline_days": float(time[-1]-time[0]),
        "n_gaps_over_3h": int(np.sum(delta > 3/24)),
        "flux_median_raw": median,
        "flux_mad_sigma_relative": stable_robust_sigma(normalized),
        "flux_p2p_sigma_ppm": float(1e6 * p2p),
        "normalized_flux_p01": float(np.quantile(normalized, 0.01)),
        "normalized_flux_p99": float(np.quantile(normalized, 0.99)),
        "fits_sha256": original_hash,
        "policy_id": POLICY_NAME,
    }
    return curve, info


def assess_pair(a, b, tic_a, tic_b, sector, pair):
    shared = a[["cadence_num", "time_btjd"]].merge(
        b[["cadence_num", "time_btjd"]],
        on="cadence_num", how="inner", validate="one_to_one",
        suffixes=("_a", "_b"),
    )
    offset_seconds = (
        (shared.time_btjd_a - shared.time_btjd_b).to_numpy(float) * 86400
    )
    n_a = len(a)
    n_b = len(b)
    return {
        "pair_id": pair,
        "sector": int(sector),
        "tic_a": tic_a, "tic_b": tic_b,
        "n_valid_a": n_a, "n_valid_b": n_b,
        "n_joint_cadence_nums": len(shared),
        "jaccard_cadence_nums": (
            len(shared) / (n_a + n_b - len(shared))
            if n_a+n_b-len(shared)>0 else np.nan
        ),
        "median_abs_BTJD_difference_seconds": (
            float(np.median(np.abs(offset_seconds)))
            if len(shared) else np.nan
        ),
        "same_source_timestamps_forced": False,
    }


def figure(tic, sector, curve, outdir):
    t = curve.time_btjd.to_numpy(dtype=float)
    flux = curve.flux_cal_aper_normalized.to_numpy(dtype=float)
    fig, axs = plt.subplots(2, 1, figsize=(10, 6.5),
                             layout="constrained")
    axs[0].plot(t-t[0], flux, ".", markersize=2.8)
    axs[0].set(
        title=f"{tic} | TESS sector {sector} | {POLICY_NAME}",
        xlabel="Days since first accepted original cadence",
        ylabel="Normalized calibrated aperture flux",
    )
    axs[1].plot(t-t[0], np.r_[np.nan, np.diff(flux)], ".", markersize=2.8)
    axs[1].set(
        xlabel="Days since first accepted original cadence",
        ylabel="Consecutive flux difference",
    )
    fig.savefig(
        outdir / f"{tic}_sector_{sector:02d}_aperture.png", dpi=140
    )
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    ap.add_argument("--outdir", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--private-dir", type=Path, default=DEFAULT_PRIVATE)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    args.private_dir.mkdir(parents=True, exist_ok=True)
    # Prevent accidental git commits of all time-series data.
    if args.private_dir.resolve().is_relative_to(ROOT.resolve()):
        raise ValueError("Private cleaned time-series must NOT live inside repo")
    figures = args.outdir / "figures"
    figures.mkdir(parents=True, exist_ok=True)

    man = pd.read_csv(args.manifest, dtype={"tic_id": str, "pair_id": str})
    required = {
        "tic_id", "sector", "sha256", "local_path",
        "author", "exptime_seconds"
    }
    if not required.issubset(man):
        raise ValueError("Input manifest lacks: " + str(required - set(man)))
    if man.duplicated(["tic_id", "sector"]).any():
        raise ValueError("Duplicate source TIC+sector")
    if man.empty:
        raise ValueError("Empty original-FITS manifest")

    curves, summary, errors = {}, [], []
    for item in man.to_dict("records"):
        tic = str(item["tic_id"])
        sector = int(item["sector"])
        try:
            curve, row = load_one(item)
            path = args.private_dir / (
                f"{tic}_sector_{sector:02d}_{POLICY_NAME}.csv"
            )
            curve.to_csv(path, index=False)
            curves[(tic,sector)] = curve
            row["private_cleaned_file"] = str(path)
            summary.append(row)
            figure(tic, sector, curve, figures)
            print(
                f"[OK] {tic} s{sector}: {row['n_accepted']}/"
                f"{row['n_raw']} valid; aperture p2p "
                f"{row['flux_p2p_sigma_ppm']:.0f} ppm",
                flush=True,
            )
        except Exception as exc:
            errors.append({
                "tic_id": tic, "sector": sector,
                "error": repr(exc)
            })
            print(f"[FAIL] {tic} s{sector}: {exc}", flush=True)

    joined = []
    if "pair_id" in man.columns:
        for pair, group in man.groupby("pair_id"):
            if len(group) != 2 or group.sector.nunique() != 1:
                continue
            a, b = group.iloc[0], group.iloc[1]
            ka = (str(a.tic_id), int(a.sector))
            kb = (str(b.tic_id), int(b.sector))
            if ka in curves and kb in curves:
                joined.append(assess_pair(
                    curves[ka], curves[kb],
                    ka[0], kb[0], ka[1], str(pair)
                ))

    stat = pd.DataFrame(summary)
    stat.to_csv(args.outdir/"aperture_curve_metrics.csv", index=False)
    pd.DataFrame(joined).to_csv(
        args.outdir/"same_sector_cadence_alignment.csv", index=False
    )
    pd.DataFrame(errors, columns=["tic_id", "sector", "error"]).to_csv(
        args.outdir/"processing_errors.csv", index=False
    )
    policy = {
        "id": POLICY_NAME,
        "status": "FROZEN FOR NEXT PROSPECTIVE COHORT",
        "reason_for_choice": (
            "A developmental original MAST TGLC FITS had all PSF "
            "measurements missing but usable calibrated aperture photometry. "
            "Selection does NOT depend on original class label."
        ),
        "product_author": "TGLC",
        "nominal_exposure_seconds": 1800,
        "flux_for_ALL_targets": FLUX_COLUMN,
        "flux_selection_per_source_or_class": False,
        "time_field": "time (BTJD)",
        "quality_flags": "TessQualityFlags default",
        "TGLC_flags_policy": "report only; do not exclude in main pipeline",
        "raw_cadence_filter": "finite time AND finite positive cal_aper_flux",
        "normalization": "single median of per-source accepted flux",
        "sigma_clip": False,
        "detrending": False,
        "rebinning": False,
        "minimum_accepted_cadences": MIN_CADENCES,
        "pair_alignment_for_audit": "cadence_num, NOT rounded BTJD",
        "pair_alignment_used_for_processing": False,
        "raw_source_files_sha256_verified": True,
        "n_sources_requested": len(man),
        "n_sources_succeeded": len(summary),
        "n_sources_failed": len(errors),
        "new_stars_required_for_independent_validation": True,
        "warning": (
            "This policy was chosen on previously examined stars. "
            "No performance from these sources is independent, and "
            "TGLC calibrated aperture depths may be diluted by contamination."
        ),
    }
    (args.outdir/"frozen_aperture_policy.json").write_text(
        json.dumps(policy, indent=2), encoding="utf-8"
    )
    report = [
        "# TGLC original uniform calibrated-aperture photometry pilot", "",
        f"Policy: {POLICY_NAME}",
        f"Processed: {len(summary)}/{len(man)}",
        f"Errors: {len(errors)}",
        "Each target uses the identical flux selection and quality rules.",
        "No transit features were extracted; no labels were used.", "",
        "## Per-star flux QA", "",
        stat.drop(columns=["fits_sha256","private_cleaned_file"],
                  errors="ignore").round(6).to_string(index=False)
        if not stat.empty else "(none)",
        "", "## Shared sector cadence numbers (audit ONLY)", "",
        pd.DataFrame(joined).round(5).to_string(index=False)
        if joined else "(none)",
        "", "## Interpretation", "",
        "Previous nonfinite PSF values do not reject a valid aperture target.",
        "No per-target fallback is applied; aperture is mandatory for ALL.",
        "Barycentric timestamps are NEVER replaced with a matched source's times.",
        "TGLC_flags remain recorded to support future sensitivity checks.",
        "The four development stars must not serve as a new heldout test set.",
    ]
    (args.outdir/"processing_report.md").write_text(
        "\n".join(report)+"\n", encoding="utf-8"
    )
    print("\n"+"\n".join(report[:8]), flush=True)
    if errors:
        raise RuntimeError(
            f"Uniform aperture pilot had {len(errors)} failures; "
            "inspect processing_errors.csv"
        )


if __name__ == "__main__":
    main()
