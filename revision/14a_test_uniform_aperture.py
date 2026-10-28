#!/usr/bin/env python3
"""Offline tests for frozen uniform original-TGLC aperture extraction.

The test deliberately supplies an all-NaN PSF channel to prove that the
universal aperture channel still works without ANY class-conditioned fallback.
"""
from __future__ import annotations

import importlib.util
import tempfile
from pathlib import Path

import numpy as np
from astropy.io import fits


SRC = Path(__file__).with_name("14_build_uniform_aperture_curves.py")
spec = importlib.util.spec_from_file_location("aperture_policy_14", SRC)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def write_original(path, tic, sector=7, dt=0.0, nan_psf=False):
    n = 510
    cadence = np.arange(n, dtype=np.int32)
    aperture = 1.0 + 0.002 * np.sin(np.arange(n)/13)
    aperture[5] = np.nan
    aperture[6] = -0.5
    psf = np.full(n, np.nan) if nan_psf else 1.0+np.cos(np.arange(n)/12)*.001
    q = np.zeros(n, dtype=np.int32)
    q[15] = 16
    tg = np.zeros(n, dtype=np.int32)
    tg[12:20] = 1
    curve = fits.BinTableHDU.from_columns([
        fits.Column(name="time", format="D",
            array=1500.0+cadence/48+dt/86400),
        fits.Column(name="cal_aper_flux", format="D", array=aperture),
        fits.Column(name="cal_psf_flux", format="D", array=psf),
        fits.Column(name="cadence_num", format="J", array=cadence),
        fits.Column(name="tess_flags", format="J", array=q),
        fits.Column(name="tglc_flags", format="J", array=tg),
    ], name="LIGHTCURVE")
    primary = fits.PrimaryHDU()
    primary.header["TELESCOP"] = "TESS"
    primary.header["SECTOR"] = sector
    primary.header["TICID"] = tic.replace("TIC_", "")
    fits.HDUList([primary, curve]).writeto(path)


def main():
    with tempfile.TemporaryDirectory(prefix="exoplanet_tglc14_") as folder:
        work = Path(folder)
        a, b = work/"original_a.fits", work/"original_b.fits"
        write_original(a, "TIC_90000301", nan_psf=True)
        write_original(b, "TIC_90000302", dt=80.5)
        items = [
            {
                "tic_id": f"TIC_9000030{i}",
                "sector": 7,
                "pair_id": "P_SYSTEST",
                "author": "TGLC",
                "exptime_seconds": 1800.0,
                "local_path": str(p),
                "sha256": mod.digest(p),
            }
            for i, p in enumerate((a,b), start=1)
        ]
        curve_a, report_a = mod.load_one(items[0])
        curve_b, report_b = mod.load_one(items[1])
        assert len(curve_a) == 507
        assert len(curve_b) == 507
        assert report_a["n_accepted"] == report_b["n_accepted"] == 507
        assert np.isclose(np.median(curve_a.flux_cal_aper_normalized), 1.0)
        assert report_a["n_TGLC_flag_nonzero_accepted"] == 8
        assert report_a["n_good_TESS_default"] == 509
        pair = mod.assess_pair(
            curve_a, curve_b, items[0]["tic_id"],
            items[1]["tic_id"], 7, "P_SYSTEST"
        )
        assert pair["jaccard_cadence_nums"] == 1.0
        assert pair["n_joint_cadence_nums"] == 507
        assert abs(pair["median_abs_BTJD_difference_seconds"]-80.5) < 0.01
        bad_item = dict(items[0], sha256="0"*64)
        try:
            mod.load_one(bad_item)
        except ValueError as exc:
            assert "SHA256" in str(exc)
        else:
            raise AssertionError("Corrupted source SHA not detected")
        other_product = dict(items[0], author="OTHER")
        try:
            mod.load_one(other_product)
        except ValueError as exc:
            assert "Non-TGLC" in str(exc)
        else:
            raise AssertionError("Unexpected product author was accepted")
        print("[PASS] Aperture works with 100% missing PSF")
        print("[PASS] Universal quality mask, 507 accepted cadences")
        print("[PASS] Same cadence numbers despite 80.5-s BTJD offset")
        print("[PASS] TGLC flags are preserved without changing main mask")
        print("[PASS] SHA256 and pipeline provenance enforced")
        print("All uniform original-TGLC aperture tests passed.")


if __name__ == "__main__":
    main()
