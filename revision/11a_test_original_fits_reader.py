#!/usr/bin/env python3
"""Regression tests for original TGLC/SPOC FITS recognition.

No internet connection and no real MAST downloads required.
Run:
    python revision/11a_test_original_fits_reader.py
"""
from __future__ import annotations

import importlib.util
import tempfile
from pathlib import Path

import numpy as np
from astropy.io import fits


SCRIPT = Path(__file__).with_name("11_download_original_mast_pilot.py")
spec = importlib.util.spec_from_file_location("mast_pilot_reader", SCRIPT)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def write_fits(path, time_col, flux_col, sector=7, tic="144415503"):
    data = fits.BinTableHDU.from_columns([
        fits.Column(
            name=time_col, format="D", unit="d",
            array=np.arange(400, dtype=float) / 48.0
        ),
        fits.Column(
            name=flux_col, format="E",
            array=np.ones(400, dtype=np.float32)
        ),
        fits.Column(
            name="tess_flags" if time_col.islower() else "QUALITY",
            format="J", array=np.zeros(400, dtype=np.int32)
        )
    ], name="LIGHTCURVE")
    primary = fits.PrimaryHDU()
    primary.header["TELESCOP"] = "TESS"
    primary.header["SECTOR"] = sector
    primary.header["TICID"] = tic.replace("TIC_", "")
    fits.HDUList([primary, data]).writeto(path, overwrite=True)


def main():
    tic = "TIC_144415503"
    with tempfile.TemporaryDirectory() as d:
        root = Path(d)
        for column, flux in [
            ("time", "cal_psf_flux"),
            ("TIME", "PDCSAP_FLUX"),
        ]:
            path = root / f"{column}.fits"
            write_fits(path, column, flux)
            header, inventory = mod.inspect_fits(
                path, {"tic_id": tic, "sector": 7, "author": "TGLC"}
            )
            assert header["n_science_rows"] == 400
            assert header["time_column_name"] == column
            assert flux in header["flux_columns_available"]
            assert header["science_hdu_index"] == 1
            assert any(x["column"] == flux for x in inventory)
            print(f"[PASS] {column} + {flux}: 400 rows")
        bad = root / "wrong-sector.fits"
        write_fits(bad, "time", "cal_aper_flux", sector=8)
        try:
            mod.inspect_fits(
                bad, {"tic_id": tic, "sector": 7, "author": "TGLC"}
            )
        except ValueError as exc:
            if "sector differs" not in str(exc):
                raise
            print("[PASS] Wrong sector rejected")
        else:
            raise AssertionError("Wrong-sector science FITS not rejected")
    print("All original-FITS reader regression tests passed.")


if __name__ == "__main__":
    main()
