#!/usr/bin/env python3
"""Offline regression tests for original TGLC independent-channel quality audit.

Two artificial original FITS have the same cadence_num exposures, but a 120-s
different barycentric-like time offset. One target has entirely negative
cal_psf_flux while calibrated aperture remains valid, reproducing the
observed failure of requiring BOTH corrected channels to be positive.

This test ensures the forensic reader reports rather than hides that case.
"""
from __future__ import annotations

import importlib.util
import tempfile
from pathlib import Path

import numpy as np
from astropy.io import fits

PATH = Path(__file__).with_name("13_tglc_quality_and_cadence_forensics.py")
spec = importlib.util.spec_from_file_location("tglc_forensics", PATH)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def synthetic_fits(path, tic, all_psf_negative=False, time_offset_seconds=0):
    n = 510
    idx = np.arange(n)
    time = 1500.0 + idx / 48.0 + time_offset_seconds / 86400.0
    psf = 1.0 + 0.01*np.sin(idx/28)
    cal_psf = -psf if all_psf_negative else psf
    aper = 1.0 + .003*np.sin(idx/17)
    tq = np.zeros(n, dtype=np.int32)
    tq[7] = 16
    q1 = np.zeros(n, dtype=np.int32)
    q1[25] = 1
    table = fits.BinTableHDU.from_columns([
        fits.Column(name="time", format="D", array=time),
        fits.Column(name="cadence_num", format="J", array=idx.astype(np.int32)),
        fits.Column(name="tess_flags", format="J", array=tq),
        fits.Column(name="TGLC_flags", format="J", array=q1),
        fits.Column(name="cal_psf_flux", format="D", array=cal_psf),
        fits.Column(name="cal_aper_flux", format="D", array=aper),
        fits.Column(name="psf_flux", format="D", array=psf*20000.0),
        fits.Column(name="aperture_flux", format="D", array=aper*20000.0),
    ], name="LIGHTCURVE")
    primary = fits.PrimaryHDU()
    primary.header["SECTOR"] = 7
    primary.header["TICID"] = tic.replace("TIC_", "")
    fits.HDUList([primary, table]).writeto(path)


def main():
    with tempfile.TemporaryDirectory(prefix="original_tglc13_") as folder:
        root = Path(folder)
        p1 = root/"a.fits"
        p2 = root/"b.fits"
        synthetic_fits(p1, "TIC_90000101", True, 0)
        synthetic_fits(p2, "TIC_90000102", False, 120)
        row1 = {
            "local_path": str(p1),
            "sha256": mod.sha256(p1),
            "tic_id": "TIC_90000101",
            "sector": 7,
            "pair_id": "PTEST"
        }
        row2 = {
            **row1,
            "local_path": str(p2),
            "sha256": mod.sha256(p2),
            "tic_id": "TIC_90000102"
        }
        inf1, det1, cad1, _, _, _ = mod.read_one(row1)
        inf2, det2, cad2, _, _, _ = mod.read_one(row2)
        assert inf1["n_corrected_joint_positive"] == 0
        assert inf1["n_corrected_aperture_only_good"] == 509
        assert inf2["n_corrected_joint_positive"] == 509
        assert any(
            x["channel"]=="cal_psf_flux"
            and x["n_quality_good_negative"] == 509
            for x in det1
        )
        pair = mod.alignment(
            cad1, cad2, row1["tic_id"], row2["tic_id"], "PTEST", 7
        )
        assert pair["n_common_cadence_nums"] == 509
        assert pair["cadence_num_jaccard"] == 1.0
        assert pair["n_common_all_four_corrected_valid"] == 0
        assert abs(abs(pair["median_signed_time_offset_seconds"])-120.0) < .001
        print("[PASS] Entirely nonpositive calibrated PSF is identified")
        print("[PASS] Corrected aperture stays independently usable")
        print("[PASS] Same 509 cadence numbers match despite 120-s offset")
        print("[PASS] Strict joint mask correctly remains empty")


if __name__ == "__main__":
    main()
