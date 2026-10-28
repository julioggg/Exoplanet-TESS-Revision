#!/usr/bin/env python3
"""Provenance-checked, label-blind TGLC photometry QA for the 4 original FITS.

This script does NOT fit or select a classifier, BLS hyperparameters, or transit
signals. Flux choice is predeclared before seeing class labels:
  primary = calibrated PSF (cal_psf_flux);
  sensitivity = calibrated aperture (cal_aper_flux).

Both use the SAME observation timestamps and the SAME joint validity mask:
finite positive corrected fluxes in BOTH channels, finite time, and the TESS
quality default bitmask. TGLC_flags are separately inventoried and compared
with a strict-zero TGLC-flag sensitivity mask, but not automatically rejected.
No sigma-clipping/detrending is performed here, to avoid erasing transits.

Cleaned, aligned photometry goes to the HPC private data directory (not Git).
Git stores only summary tables, provenance checks, and diagnostic figures.

Run from repository root:
    python -u revision/12_original_tglc_photometry_qa.py
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits
from lightkurve.utils import TessQualityFlags

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "results/mast_original_pilot/verified_original_fits.csv"
RESULT_DIR = ROOT / "results/mast_photometry_qa"
CLEANED_DIR = Path.home() / "exoplanet_revision/data/tglc_qa_cleaned"
DEFAULT_EXP_DAYS = 1800.0 / 86400.0

COLUMNS = {
    "primary_psf": "cal_psf_flux",
    "sensitivity_aperture": "cal_aper_flux",
}

def checksum(path):
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def robust_sigma(v):
    x = np.asarray(v, dtype=float)
    x = x[np.isfinite(x)]
    if not len(x):
        return float("nan")
    return float(1.4826 * np.median(np.abs(x - np.median(x))))


def metric(t, f):
    x = np.asarray(f, dtype=float)
    tt = np.asarray(t, dtype=float)
    assert len(x) == len(tt) and len(x) > 0
    baseline = float(np.median(x))
    norm = x / baseline
    dt = np.diff(tt)
    increasing = dt > 0
    in_order = bool(np.all(increasing))
    if not in_order:
        raise ValueError("Non-increasing or duplicate observation times")
    p2p = (
        float(robust_sigma(np.diff(norm)) / np.sqrt(2.0))
        if len(norm) > 1 else float("nan")
    )
    return {
        "n_valid_cadences": int(len(x)),
        "normalized_mean": float(np.mean(norm)),
        "normalized_median": float(np.median(norm)),
        "normalized_std": float(np.std(norm, ddof=0)),
        "normalized_robust_sigma": robust_sigma(norm),
        "normalized_p2p_sigma": p2p,
        "normalized_p2p_ppm": 1e6 * p2p,
        "normalized_amplitude_p99_p01": float(
            np.quantile(norm, .99) - np.quantile(norm, .01)
        ),
        "raw_median_flux": baseline,
        "min_normalized_flux": float(np.min(norm)),
        "max_normalized_flux": float(np.max(norm)),
        "baseline_days": float(tt[-1] - tt[0]),
        "median_cadence_minutes": float(np.median(dt)*1440) if len(dt) else np.nan,
        "n_gaps_over_1h": int(np.sum(dt > 1.0 / 24.0)),
        "n_gaps_over_3h": int(np.sum(dt > 3.0 / 24.0)),
    }


def hdu_header_value(hdul, key, default=""):
    for h in hdul:
        value = h.header.get(key)
        if value not in (None, ""):
            return value
    return default


def extract(path, info):
    with fits.open(path, memmap=False) as hdul:
        candidates = []
        for idx, h in enumerate(hdul):
            cols = getattr(h, "columns", None)
            names = [str(n) for n in cols.names] if cols is not None else []
            norm = {n.lower(): n for n in names}
            if ("time" in norm and all(c in norm for c in COLUMNS.values())
                    and "tess_flags" in norm):
                candidates.append((idx, h, norm))
        if len(candidates) != 1:
            raise ValueError(
                f"Expected one unambiguous TGLC table, found {len(candidates)}"
            )
        idx, h, norm = candidates[0]
        data = h.data
        t = np.asarray(data[norm["time"]], dtype=np.float64)
        psf = np.asarray(data[norm["cal_psf_flux"]], dtype=np.float64)
        aper = np.asarray(data[norm["cal_aper_flux"]], dtype=np.float64)
        q0 = np.asarray(data[norm["tess_flags"]], dtype=np.int64)
        if "tglc_flags" in norm:
            q1 = np.asarray(data[norm["tglc_flags"]], dtype=np.int64)
        else:
            q1 = np.zeros(len(t), dtype=np.int64)
        cadence_num = (
            np.asarray(data[norm["cadence_num"]])
            if "cadence_num" in norm else np.arange(len(t))
        )
        if len({len(t), len(psf), len(aper), len(q0), len(q1),
                len(cadence_num)}) != 1:
            raise ValueError("TGLC column row count mismatch")
        sector = hdu_header_value(hdul, "SECTOR")
        if sector != "" and int(sector) != int(info["sector"]):
            raise ValueError("Sector mismatch against verified FITS manifest")
        tic = hdu_header_value(hdul, "TICID")
        if tic != "" and re.sub(r"\D", "", str(tic)) != re.sub(
            r"\D", "", str(info["tic_id"])
        ):
            raise ValueError("TIC header mismatch")
    return t, psf, aper, q0, q1, cadence_num


def make_plot(tic, sector, time, psf, aper, outdir, q0, q1):
    fig, axes = plt.subplots(3, 1, figsize=(10.5, 8.0),
                             sharex=True, layout="constrained")
    d = time-time[0]
    axes[0].plot(d, psf, ".", ms=2.5, alpha=.8)
    axes[0].set(ylabel="Normalized PSF flux",
                title=f"{tic} | Sector {sector} | TGLC corrected flux")
    axes[1].plot(d, aper, ".", ms=2.5, alpha=.8)
    axes[1].set(ylabel="Normalized aperture flux")
    axes[2].plot(d, psf-aper, ".", ms=2.5, alpha=.8)
    axes[2].axhline(0, color="black", ls=":", lw=.7)
    axes[2].set(xlabel="Days from first accepted exposure",
                ylabel="PSF - aperture")
    fig.savefig(outdir / f"{tic}_s{sector}_psf_aperture_qa.png",
                dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=Path, default=RESULT_DIR)
    parser.add_argument("--cleaned-dir", type=Path, default=CLEANED_DIR)
    args = parser.parse_args()
    args.outdir.mkdir(exist_ok=True, parents=True)
    args.cleaned_dir.mkdir(exist_ok=True, parents=True)
    figures = args.outdir / "figures"
    figures.mkdir(exist_ok=True, parents=True)

    manifest = pd.read_csv(MANIFEST, dtype={"tic_id": str})
    if manifest.empty or manifest.tic_id.duplicated().any():
        raise ValueError("Empty or duplicated verified FITS manifest")
    if len(manifest) != 4:
        print(f"WARNING: pilot manifest contains {len(manifest)} files, not 4",
              flush=True)

    metrics = []
    audit = []
    failure = []
    for item in manifest.to_dict("records"):
        tic = str(item["tic_id"])
        pair = str(item["pair_id"])
        sector = int(item["sector"])
        try:
            path = Path(str(item["local_path"]))
            if not path.is_file():
                raise FileNotFoundError(path)
            digest = checksum(path)
            if digest != str(item["sha256"]):
                raise ValueError("SHA256 mismatch from original verified file")
            t, psf, aper, q0, q1, cadno = extract(path, item)
            if len(t) < 300:
                raise ValueError("Too few original observations")

            finite_time = np.isfinite(t)
            finite_psf = np.isfinite(psf)
            finite_aper = np.isfinite(aper)
            psf_pos = finite_psf & (psf > 0)
            aper_pos = finite_aper & (aper > 0)
            tess_mask = TessQualityFlags.create_quality_mask(
                q0, bitmask="default"
            )
            # Common eligibility independent of positive/negative catalogue
            # membership, chosen to make corrected PSF/aperture comparable.
            mask = finite_time & psf_pos & aper_pos & tess_mask
            strict_mask = mask & (q1 == 0)
            kept = np.flatnonzero(mask)
            if len(kept) < 300:
                raise ValueError(f"Only {len(kept)} jointly valid cadences")
            order = np.argsort(t[kept], kind="mergesort")
            indices = kept[order]
            tt = t[indices]
            f1 = psf[indices]
            f2 = aper[indices]
            if len(np.unique(tt)) != len(tt):
                raise ValueError("Duplicate time stamps among accepted cadences")
            f1n = f1 / np.median(f1)
            f2n = f2 / np.median(f2)
            assert np.isfinite(f1n).all() and np.isfinite(f2n).all()
            flags_detail = {
                "n_total": len(t),
                "n_finite_time": int(finite_time.sum()),
                "n_nonfinite_cal_psf_flux": int((~finite_psf).sum()),
                "n_nonfinite_cal_aper_flux": int((~finite_aper).sum()),
                "n_nonpositive_cal_psf_flux": int(
                    (finite_psf & ~psf_pos).sum()
                ),
                "n_nonpositive_cal_aper_flux": int(
                    (finite_aper & ~aper_pos).sum()
                ),
                "n_nonzero_TESS_flags": int((q0 != 0).sum()),
                "n_failed_TESS_default": int((~tess_mask).sum()),
                "n_nonzero_TGLC_flags": int((q1 != 0).sum()),
                "n_common_quality_kept": len(indices),
                "n_strict_TGLC_zero_kept": int(strict_mask.sum()),
                "fraction_kept": len(indices) / len(t),
                "tess_flags_unique": ";".join(
                    f"{int(x)}:{int(c)}" for x, c in
                    zip(*np.unique(q0, return_counts=True))
                )[:400],
                "tglc_flags_unique": ";".join(
                    f"{int(x)}:{int(c)}" for x, c in
                    zip(*np.unique(q1, return_counts=True))
                )[:400],
            }
            audit.append({"tic_id": tic, "pair_id": pair,
                          "sector": sector, **flags_detail})
            for role, flux in [
                ("primary_psf", f1), ("sensitivity_aperture", f2)
            ]:
                vals = metric(tt, flux)
                vals.update({
                    "tic_id": tic, "pair_id": pair, "sector": sector,
                    "method": role,
                    "n_valid_common": len(tt),
                    "fraction_kept": len(indices) / len(t),
                    "tess_default_mask": True,
                    "tglc_flag_zero_required": False,
                    "no_sigma_clipping": True,
                    "no_detrending": True,
                })
                # Quantify strict-TGLC-flag sensitivity on IDENTICAL masks.
                strict_local = q1[indices] == 0
                vals["n_tglc_strict_valid"] = int(strict_local.sum())
                vals["strict_tglc_robust_sigma"] = (
                    metric(tt[strict_local], flux[strict_local])[
                        "normalized_robust_sigma"]
                    if strict_local.sum() >= 300 else np.nan
                )
                metrics.append(vals)
            make_plot(tic, sector, tt, f1n, f2n, figures, q0, q1)
            out = args.cleaned_dir / f"{tic}_sector_{sector:02d}.csv"
            pd.DataFrame({
                "time_btjd": tt, "cal_psf_flux_normalized": f1n,
                "cal_aper_flux_normalized": f2n,
                "cadence_num": cadno[indices],
                "TESS_flags": q0[indices],
                "TGLC_flags": q1[indices],
            }).to_csv(out, index=False)
            print(f"[OK] {tic} s{sector}: "
                  f"{len(indices)}/{len(t)} common cadences; "
                  f"PSF p2p={metric(tt, f1)['normalized_p2p_ppm']:.0f} ppm",
                  flush=True)
        except Exception as exc:
            failure.append({
                "tic_id": tic, "sector": sector, "pair_id": pair,
                "error": repr(exc),
            })
            print(f"[FAIL] {tic} s{sector}: {exc}", flush=True)

    stat_table = pd.DataFrame(metrics)
    audit_table = pd.DataFrame(audit)
    stat_table.to_csv(args.outdir / "flux_method_metrics.csv", index=False)
    audit_table.to_csv(args.outdir / "cadence_quality_audit.csv", index=False)
    pd.DataFrame(failure, columns=[
        "tic_id", "sector", "pair_id", "error"
    ]).to_csv(args.outdir / "qa_failures.csv", index=False)

    # Paired timestamp comparisons are independent of catalogue labels.
    pairs = []
    for (pair, sector), rows in audit_table.groupby(
        ["pair_id", "sector"]
    ) if not audit_table.empty else []:
        tics = sorted(rows.tic_id.unique())
        if len(tics) != 2:
            continue
        a = pd.read_csv(
            args.cleaned_dir / f"{tics[0]}_sector_{sector:02d}.csv"
        ).time_btjd.to_numpy(float)
        b = pd.read_csv(
            args.cleaned_dir / f"{tics[1]}_sector_{sector:02d}.csv"
        ).time_btjd.to_numpy(float)
        # Compare timestamps at a strict 30-second quantization; no claim of
        # perfect pair alignment if time samples differ by >30 seconds.
        ka = set(np.rint(a * 2880).astype("int64"))
        kb = set(np.rint(b * 2880).astype("int64"))
        overlap = len(ka & kb)
        names = {
            str(entry["tic_id"]): str(entry.get("product_filename", ""))
            for entry in manifest.to_dict("records")
            if str(entry["pair_id"]) == str(pair)
        }
        def camera_ccd(tic):
            m = re.search(r"cam([1-4])-ccd([1-4])", names.get(tic, ""))
            return (int(m.group(1)), int(m.group(2))) if m else None
        loc_a = camera_ccd(tics[0])
        loc_b = camera_ccd(tics[1])
        pairs.append({
            "pair_id": pair, "sector": sector,
            "tic_a": tics[0], "tic_b": tics[1],
            "camera_a": loc_a[0] if loc_a else np.nan,
            "camera_b": loc_b[0] if loc_b else np.nan,
            "ccd_a": loc_a[1] if loc_a else np.nan,
            "ccd_b": loc_b[1] if loc_b else np.nan,
            "same_camera_and_ccd": (
                bool(loc_a == loc_b) if loc_a and loc_b else None
            ),
            "n_a": len(ka), "n_b": len(kb),
            "n_matched_half_minute_bins": overlap,
            "jaccard_half_minute_bins": (
                overlap / len(ka | kb) if len(ka | kb) else 0.0
            ),
        })
    pd.DataFrame(pairs).to_csv(
        args.outdir / "paired_time_grid_overlap.csv", index=False
    )

    config = {
        "data_source": "Original MAST TGLC FITS sector 2 and 7 pilot",
        "flux_primary": "cal_psf_flux",
        "flux_secondary_sensitivity": "cal_aper_flux",
        "flux_choice_made_without_labels": True,
        "common_mask": (
            "finite time and finite POSITIVE calibrated PSF/aperture "
            "flux AND Lightkurve default TESS quality mask"
        ),
        "tglc_flags_policy": (
            "Inventory only in primary; strict 0 separately analyzed "
            "as sensitivity"
        ),
        "time_frame": "TGLC TIME (BTJD; no conversion needed for relative days)",
        "no_clipping_and_no_detrending": True,
        "raw_FITS_never_uploaded_to_GitHub": True,
        "private_cleaned_data": str(args.cleaned_dir),
        "n_fits_requested": len(manifest),
        "n_fits_success": len(audit),
        "n_fits_failed": len(failure),
        "limitations": [
            "Four previously inspected targets: engineering feasibility, not independent validation.",
            "Different camera/CCD combinations may cause remaining differences even within sector.",
            "Primary-mission calibrated TGLC flux can be corrupted for faint/contaminated sources.",
            "No transit-injection/recovery or false-positive calibration is attempted here.",
            "Point-to-point scatter is a descriptive proxy, not a rigorous CDPP calculation.",
            "Rejecting nonpositive calibrated flux may bias rare, strongly variable targets.",
            "Manual exclusion of TGLC-flagged cadences would require separate validation.",
        ],
    }
    (args.outdir / "qa_config.json").write_text(
        json.dumps(config, indent=2), encoding="utf-8"
    )
    lines = [
        "# Original MAST TGLC photometric QA (label-blind)", "",
        f"FITS audited: {len(audit)}/{len(manifest)}",
        f"Failed: {len(failure)}",
        f"Predeclared primary: {config['flux_primary']}",
        f"Secondary sensitivity: {config['flux_secondary_sensitivity']}",
        "", "## Shared sampling and flag breakdown", "",
        (audit_table.round(5).to_string(index=False)
         if not audit_table.empty else "(none)"),
        "", "## Per-curve photometric scatter (no clipping/detrending)", "",
        (stat_table.round(6).to_string(index=False)
         if not stat_table.empty else "(none)"),
        "", "## Pair timestamp support overlap", "",
        (pd.DataFrame(pairs).round(5).to_string(index=False)
         if pairs else "(none)"),
        "", "## Cautions", "",
        *["- "+s for s in config["limitations"]],
    ]
    (args.outdir / "qa_report.md").write_text(
        "\n".join(lines)+"\n", encoding="utf-8"
    )
    print("\nORIGINAL TGLC QA")
    print(f"Audited: {len(audit)}/{len(manifest)}")
    print(f"Failures: {len(failure)}")
    print(stat_table[
        ["tic_id", "method", "n_valid_cadences",
         "normalized_robust_sigma", "normalized_p2p_ppm"]
    ].round(6).to_string(index=False) if len(stat_table) else "No results")
    print(f"\nSaved report: {args.outdir / 'qa_report.md'}")


if __name__ == "__main__":
    main()
