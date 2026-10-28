#!/usr/bin/env python3
"""TGLC original-product forensic audit WITHOUT changing the QA criteria.

After experiment 12 found a target with zero jointly valid cadences and
zero overlapping rounded BTJD timestamps in one pair, this audit determines:
  1) exactly which photometry channel / value signs / flags rejected the star;
  2) which alternative channels contain usable samples;
  3) whether matched sector pairs share TESS 'cadence_num' indexes despite
     per-target barycentric timestamp offsets.

No flux fallback, manual exclusion, new model, transit search or class-based
decision is made here. Only original MAST FITS bytes are read, verified against
the SHA256 acquisition manifest. All outputs are aggregate tables/figures.

Run:
  python -u revision/13_tglc_quality_and_cadence_forensics.py
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
MANIFEST = ROOT / "results/mast_original_pilot/verified_original_fits.csv"
OUT = ROOT / "results/tglc_integrity_forensics"
CHANNELS = [
    "cal_psf_flux", "cal_aper_flux", "psf_flux", "aperture_flux"
]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fp:
        for part in iter(lambda: fp.read(1024 * 1024), b""):
            h.update(part)
    return h.hexdigest()


def robust_noise(x):
    x = np.asarray(x, dtype=float)
    med = np.median(x)
    return float(1.4826 * np.median(np.abs(x - med)))


def stats(values, common):
    """Channel integrity independently of *other* photometry channels."""
    raw = np.asarray(values, dtype=float)
    fin = np.isfinite(raw)
    pos = fin & (raw > 0)
    eligible = common & pos
    after = raw[eligible]
    base = float(np.median(after)) if len(after) else np.nan
    normalized = after/base if np.isfinite(base) and base > 0 else np.array([])
    return {
        "n_finite": int(fin.sum()),
        "n_nonfinite": int((~fin).sum()),
        "n_positive": int(pos.sum()),
        "n_zero": int((fin & (raw == 0)).sum()),
        "n_negative": int((fin & (raw < 0)).sum()),
        "n_quality_good_positive": int(eligible.sum()),
        "fraction_quality_good_positive": float(
            eligible.sum() / max(1, common.sum())
        ),
        "n_quality_good_negative": int(
            (common & fin & (raw < 0)).sum()
        ),
        "n_quality_good_zero": int(
            (common & fin & (raw == 0)).sum()
        ),
        "median_finite_flux": (
            float(np.median(raw[fin])) if fin.any() else np.nan
        ),
        "median_positive_good_flux": base,
        "finite_flux_q01": (
            float(np.quantile(raw[fin], .01)) if fin.any() else np.nan
        ),
        "finite_flux_q99": (
            float(np.quantile(raw[fin], .99)) if fin.any() else np.nan
        ),
        "normalized_robust_sigma_if_eligible": (
            robust_noise(normalized) if len(normalized)>10 else np.nan
        ),
        "eligible_by_channel_300_points": bool(len(after) >= 300),
    }


def read_one(row):
    filepath = Path(row["local_path"])
    if not filepath.is_file():
        raise FileNotFoundError(filepath)
    if sha256(filepath) != str(row["sha256"]):
        raise ValueError("Original FITS SHA256 mismatch")
    with fits.open(filepath, memmap=False) as h:
        table_matches = []
        for ix, ext in enumerate(h):
            cols = getattr(ext, "columns", None)
            names = list(cols.names or []) if cols is not None else []
            dictionary = {str(n).lower(): str(n) for n in names}
            if set(["time", "tess_flags", "cadence_num"]) <= set(dictionary):
                table_matches.append((ix, ext, dictionary))
        if len(table_matches) != 1:
            raise ValueError(
                f"Expected 1 TGLC time table; got {len(table_matches)}"
            )
        _, ext, cols = table_matches[0]
        data = ext.data
        def array(name, dtype=float):
            if name not in cols:
                return None
            return np.asarray(data[cols[name]], dtype=dtype)
        time = array("time")
        cad = array("cadence_num", np.int64)
        tq = array("tess_flags", np.int64)
        gq = array("tglc_flags", np.int64)
        if gq is None:
            gq = np.zeros(len(time), dtype=np.int64)
        channels = {key: array(key) for key in CHANNELS}
        if any(v is None for v in channels.values()):
            raise ValueError("One or more expected TGLC flux columns missing")
        if any(len(v) != len(time) for v in [
            cad, tq, gq, *channels.values()
        ]):
            raise ValueError("Column length mismatch")
        sector = None
        ticid = None
        for ext2 in h:
            if "SECTOR" in ext2.header:
                sector = int(ext2.header["SECTOR"])
            if "TICID" in ext2.header:
                ticid = str(int(ext2.header["TICID"]))
        if sector != int(row["sector"]):
            raise ValueError("Sector differs from verified metadata")
        if ticid and ticid != "".join(filter(str.isdigit, row["tic_id"])):
            raise ValueError("TIC differs from verified metadata")

    time_good = np.isfinite(time)
    tess_good = TessQualityFlags.create_quality_mask(tq, bitmask="default")
    both_base = time_good & tess_good
    channel_ok = {
        name: both_base & np.isfinite(flux) & (flux > 0)
        for name, flux in channels.items()
    }
    corrected_joint = (
        channel_ok["cal_psf_flux"] & channel_ok["cal_aper_flux"]
    )
    all_four_joint = np.logical_and.reduce(
        [channel_ok[x] for x in CHANNELS]
    )
    info = {
        "tic_id": str(row["tic_id"]),
        "sector": int(row["sector"]),
        "pair_id": str(row["pair_id"]),
        "n_raw": len(time),
        "n_time_finite": int(time_good.sum()),
        "n_tess_default_good": int(tess_good.sum()),
        "n_time_and_tess_good": int(both_base.sum()),
        "n_tglc_flags_nonzero": int(np.count_nonzero(gq)),
        "n_corrected_joint_positive": int(corrected_joint.sum()),
        "n_all_four_joint_positive": int(all_four_joint.sum()),
        "n_corrected_psf_only_good": int(
            (channel_ok["cal_psf_flux"] &
             ~channel_ok["cal_aper_flux"]).sum()
        ),
        "n_corrected_aperture_only_good": int(
            (~channel_ok["cal_psf_flux"] &
             channel_ok["cal_aper_flux"]).sum()
        ),
        "n_unprocessed_psf_good": int(channel_ok["psf_flux"].sum()),
        "n_unprocessed_aperture_good": int(
            channel_ok["aperture_flux"].sum()
        ),
    }
    detail = []
    for name, flux in channels.items():
        detail.append({
            "tic_id": str(row["tic_id"]),
            "sector": int(row["sector"]),
            "pair_id": str(row["pair_id"]),
            "channel": name,
            **stats(flux, both_base),
        })
    keep = both_base & np.isfinite(cad)
    if len(np.unique(cad[keep])) != keep.sum():
        raise ValueError("Repeated cadence_num among quality-good cadences")
    cad_table = pd.DataFrame({
        "cadence_num": cad[keep].astype("int64"),
        "time_btjd": time[keep],
        "quality": tq[keep],
        "cal_psf_valid": channel_ok["cal_psf_flux"][keep],
        "cal_aper_valid": channel_ok["cal_aper_flux"][keep],
    }).sort_values("cadence_num").reset_index(drop=True)
    return info, detail, cad_table, channels, time, both_base


def figure(tic, time, channels, quality, output):
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 6.5),
                             layout="constrained")
    for ax, channel in zip(axes.flat, CHANNELS):
        flux = channels[channel]
        fin = np.isfinite(flux)
        if np.any(fin):
            ti = time[fin]
            val = flux[fin]
            ax.plot(ti-ti[0], val, ".", alpha=.5, markersize=2.6)
            ax.plot(time[fin & quality] - ti[0],
                    flux[fin & quality], ".", alpha=.4,
                    markersize=1.8)
            ax.axhline(0, color="black", linestyle=":", lw=.8)
            ax.set_title(channel)
            ax.set_xlabel("Days from first finite exposure")
            ax.set_ylabel("Original FITS flux")
        else:
            ax.text(.5,.5,"No finite flux",ha="center",va="center",
                    transform=ax.transAxes)
    fig.suptitle(f"TGLC channel integrity: {tic}")
    fig.savefig(output / f"{tic}_all_channels.png", dpi=140)
    plt.close(fig)


def alignment(table_a, table_b, a, b, pair, sector):
    merged = table_a.merge(
        table_b, how="inner", on="cadence_num",
        suffixes=("_a","_b"), validate="one_to_one",
    )
    ta = table_a["time_btjd"].to_numpy(float)
    tb = table_b["time_btjd"].to_numpy(float)
    cad_a = set(table_a.cadence_num.to_list())
    cad_b = set(table_b.cadence_num.to_list())
    dt_s = (
        (merged.time_btjd_a - merged.time_btjd_b).to_numpy(float) * 86400.
        if len(merged) else np.array([])
    )
    qa_good_both = (
        merged.cal_psf_valid_a & merged.cal_psf_valid_b &
        merged.cal_aper_valid_a & merged.cal_aper_valid_b
    )
    return {
        "pair_id": pair, "sector": sector,
        "tic_a": a, "tic_b": b,
        "n_quality_good_a": len(table_a),
        "n_quality_good_b": len(table_b),
        "n_common_cadence_nums": len(merged),
        "cadence_num_jaccard": (
            len(merged) / len(cad_a | cad_b)
            if (cad_a or cad_b) else np.nan
        ),
        "n_common_all_four_corrected_valid": int(qa_good_both.sum()),
        "median_absolute_time_offset_seconds": (
            float(np.median(np.abs(dt_s))) if len(dt_s) else np.nan
        ),
        "median_signed_time_offset_seconds": (
            float(np.median(dt_s)) if len(dt_s) else np.nan
        ),
        "range_time_offset_seconds": (
            float(np.ptp(dt_s)) if len(dt_s) else np.nan
        ),
        "timing_comparison_is": "same TESS cadence identifier; BTJD differs per target",
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", type=Path, default=OUT)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    figures = args.outdir / "figures"
    figures.mkdir(exist_ok=True)

    manifest = pd.read_csv(MANIFEST, dtype={
        "tic_id": str, "pair_id": str
    })
    info, channels, timing, fails = [], [], {}, []
    for row in manifest.to_dict("records"):
        tic = row["tic_id"]
        try:
            one, detail, cadence, fluxes, t, mask = read_one(row)
            info.append(one)
            channels.extend(detail)
            timing[tic] = cadence
            figure(tic, t, fluxes, mask, figures)
            print(
                f"[OK] {tic} s{row['sector']} "
                f"PSF={one['n_corrected_psf_only_good'] + one['n_corrected_joint_positive']}, "
                f"APER={one['n_corrected_aperture_only_good'] + one['n_corrected_joint_positive']}, "
                f"joint={one['n_corrected_joint_positive']}",
                flush=True,
            )
        except Exception as exc:
            fails.append({"tic_id": tic, "error": repr(exc)})
            print(f"[FAIL] {tic}: {exc}", flush=True)

    detail_table = pd.DataFrame(channels)
    totals = pd.DataFrame(info)
    summary = []
    for pair, group in manifest.groupby("pair_id"):
        if len(group) != 2:
            continue
        rowa, rowb = group.iloc[0], group.iloc[1]
        if rowa.tic_id in timing and rowb.tic_id in timing:
            summary.append(alignment(
                timing[rowa.tic_id], timing[rowb.tic_id],
                rowa.tic_id, rowb.tic_id, pair, int(rowa.sector),
            ))
    aligned = pd.DataFrame(summary)
    totals.to_csv(args.outdir / "per_fits_integrity.csv", index=False)
    detail_table.to_csv(args.outdir / "per_flux_column_integrity.csv",
                        index=False)
    aligned.to_csv(args.outdir / "pair_cadence_num_alignment.csv",
                   index=False)
    pd.DataFrame(fails, columns=[
        "tic_id", "error"
    ]).to_csv(args.outdir / "forensics_errors.csv", index=False)
    report = [
        "# Original TGLC channel and cadence forensics", "",
        "Previously obtained FITS rechecked by SHA256.", "",
        "## Which original channels reject each star?", "",
        totals.to_string(index=False) if not totals.empty else "(none)",
        "", "## Quality-good and signed-flux counts", "",
        detail_table.round(5).to_string(index=False) if not detail_table.empty else "(none)",
        "", "## Matching actual observing cadence numbers", "",
        aligned.round(4).to_string(index=False) if not aligned.empty else "(none)",
        "", "## Design constraints", "",
        "- No previous flux-selection or clipping rule was changed.",
        "- Independent per-column validity is inspected BEFORE any new policy.",
        "- A TGLC calibrated flux may be corrupt even if a raw column is finite.",
        "- Matched cadence_num is not automatically identical local barycentric time.",
        "- No astrophysical labels were used; these old TICs are engineering samples.",
        "- Next sample must be newly drawn with an a priori quality policy.",
    ]
    (args.outdir / "forensics_report.md").write_text(
        "\n".join(report) + "\n", encoding="utf-8"
    )
    (args.outdir / "forensics_config.json").write_text(
        json.dumps({
            "source_manifest": str(MANIFEST),
            "n_original_fits": len(manifest),
            "n_success": len(info),
            "n_errors": len(fails),
            "original_fits_bytes_changed": False,
            "flux_fallback_applied": False,
            "quality_masks_changed": False,
            "classification_or_BLS_fitted": False,
            "selection_used_catalogue_label": False,
            "timestamp_matching_method": "cadence_num",
            "recommended_followup": (
                "Predeclare the same photometric channel and source "
                "quality policy for ALL new TICs; use originals only."
            ),
        }, indent=2), encoding="utf-8"
    )
    print("\n" + "\n".join(report[:15]), flush=True)
    if fails:
        raise RuntimeError(
            f"{len(fails)} source FITS forensics failed; "
            "inspect forensics_errors.csv"
        )


if __name__ == "__main__":
    main()
