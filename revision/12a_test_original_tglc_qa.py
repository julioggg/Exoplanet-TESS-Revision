#!/usr/bin/env python3
"""Offline synthetic smoke-test for the label-blind original TGLC QA pipeline.

Tests case handling, common cadence masking, independent PSF/aperture channel
metrics, and the requirement that cleaned raw data stay outside the repo.
No network, original FITS or astrophysical labels required.

Run: python revision/12a_test_original_tglc_qa.py
"""
from __future__ import annotations

import importlib.util
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[1]
MAIN = ROOT / "revision/12_original_tglc_photometry_qa.py"
spec = importlib.util.spec_from_file_location("tglc_labelblind_audit", MAIN)
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


def create_synthetic(path, tic, sector):
    n = 510
    phase = np.arange(n, dtype=float)
    time = 1500.0 + phase * (1800.0 / 86400.0)
    cal_psf = 1.0 + 0.008 * np.sin(phase / 40.0)
    cal_aper = 1.0 + 0.015 * np.sin(phase / 45.0)
    cal_psf[3] = np.nan
    cal_aper[4] = np.nan
    cal_psf[10] = -0.03
    cal_aper[12] = -0.02
    tess = np.zeros(n, dtype=np.int32)
    tess[14] = 16  # default TESS quality bitmask excludes safe mode
    tglc = np.zeros(n, dtype=np.int32)
    tglc[15:25] = 1
    cols = [
        fits.Column(name="time", format="D", array=time),
        fits.Column(name="cal_psf_flux", format="D", array=cal_psf),
        fits.Column(name="cal_aper_flux", format="D", array=cal_aper),
        fits.Column(name="TESS_flags", format="J", array=tess),
        fits.Column(name="TGLC_flags", format="J", array=tglc),
        fits.Column(name="cadence_num", format="J",
                    array=np.arange(n, dtype=np.int32)),
    ]
    pri = fits.PrimaryHDU()
    pri.header["TELESCOP"] = "TESS"
    pri.header["SECTOR"] = sector
    pri.header["TICID"] = tic.replace("TIC_", "")
    tab = fits.BinTableHDU.from_columns(cols, name="LIGHTCURVE")
    fits.HDUList([pri, tab]).writeto(path)


def main():
    with tempfile.TemporaryDirectory(prefix="tglc_qa_test_") as d:
        root = Path(d)
        outdir = root / "qa_report"
        cleaned = root / "private_cleaned"
        records = []
        for tic in ("TIC_90000001", "TIC_90000002"):
            filename = root / f"{tic}.fits"
            create_synthetic(filename, tic, 7)
            records.append({
                "tic_id": tic, "sector": 7, "pair_id": "PAIR_TEST",
                "local_path": str(filename),
                "sha256": audit.checksum(filename),
            })
        manifest = root / "synthetic_manifest.csv"
        pd.DataFrame(records).to_csv(manifest, index=False)
        audit.MANIFEST = manifest
        sys.argv = [
            str(MAIN), "--outdir", str(outdir),
            "--cleaned-dir", str(cleaned),
        ]
        audit.main()
        errors = pd.read_csv(outdir / "qa_failures.csv")
        if len(errors):
            raise AssertionError("Synthetic QA failed: " + errors.to_string())
        summary = pd.read_csv(outdir / "flux_method_metrics.csv")
        cadence = pd.read_csv(outdir / "cadence_quality_audit.csv")
        assert len(summary) == 4
        assert len(cadence) == 2
        assert set(summary.method) == {
            "primary_psf", "sensitivity_aperture"
        }
        assert (cadence.n_common_quality_kept < 510).all()
        assert (cadence.n_common_quality_kept > 300).all()
        assert (cadence.n_nonzero_TGLC_flags == 10).all()
        assert (summary.normalized_p2p_ppm > 0).all()
        assert len(list(cleaned.glob("*.csv"))) == 2
        assert len(list((outdir / "figures").glob("*.png"))) == 2
        overlap = pd.read_csv(outdir / "paired_time_grid_overlap.csv")
        assert len(overlap) == 1
        assert overlap.jaccard_half_minute_bins.iloc[0] == 1.0
        print("[PASS] 2 synthetic original TGLC FITS")
        print("[PASS] Shared quality and channel masks")
        print("[PASS] Label-blind PSF/aperture QA")
        print("[PASS] Pairwise timestamp grid alignment")
        print("[PASS] Private cleaned files, report and diagnostic figures")


if __name__ == "__main__":
    main()
