#!/usr/bin/env python3
"""Offline tests for frozen prospective MAST FITS acquisition.

Runs without external MAST access. Confirms that the exact selected
original URI determines the local file identity, validates FITS provenance,
accepts a well-formed TGLC table and rejects wrong TIC/sector/missing flux.
"""
from __future__ import annotations

import importlib.util
import tempfile
from pathlib import Path

import numpy as np
from astropy.io import fits

SRC=Path(__file__).with_name("16_acquire_frozen_new_hosts.py")
spec=importlib.util.spec_from_file_location("frozen_hosts",SRC)
code=importlib.util.module_from_spec(spec)
spec.loader.exec_module(code)


def write_tglc(path, tic="TIC_90000401", sector=4,
               include_aperture=True):
    n=510
    columns=[
        fits.Column(
            name="time",format="D",
            array=1350.0+np.arange(n,dtype=float)/48
        ),
        fits.Column(
            name="TESS_flags",format="J",
            array=np.zeros(n,dtype=np.int32)
        ),
        fits.Column(
            name="cadence_num",format="J",
            array=np.arange(n,dtype=np.int32)
        )
    ]
    if include_aperture:
        columns.append(fits.Column(
            name="cal_aper_flux",format="E",
            array=np.ones(n,dtype=np.float32)
        ))
    tab=fits.BinTableHDU.from_columns(columns,name="LIGHTCURVE")
    primary=fits.PrimaryHDU()
    primary.header["TELESCOP"]="TESS"
    primary.header["ORIGIN"]="UCSB/TGLC"
    primary.header["SECTOR"]=sector
    primary.header["TICID"]=tic.split("_")[1]
    fits.HDUList([primary,tab]).writeto(path)


def main():
    item={
        "tic_id":"TIC_90000401",
        "sector":4,
        "author":"TGLC",
        "exptime_seconds":1800,
        "mast_uri":"mast:HLSP/tglc/s0004/cam2-ccd1/a_original.fits",
    }
    name=code.file_name_for(item)
    assert name.startswith("TIC_90000401_s0004_TGLC1800_")
    assert name.endswith(".fits")
    assert code.file_name_for(item)==name
    changed=dict(item)
    changed["mast_uri"]=changed["mast_uri"].replace("s0004","s0005")
    try:
        code.file_name_for(changed)
    except ValueError as exc:
        assert "sector" in str(exc)
    else:
        raise AssertionError("Expected URI-sector mismatch")

    with tempfile.TemporaryDirectory(prefix="original_mast16_") as work:
        root=Path(work)
        original=root/name
        write_tglc(original)
        info,cols=code.inspect_original(original,item)
        assert info["n_science_rows"]==510
        assert info["flux_column_name"]=="cal_aper_flux"
        assert any(row["column_name"]=="cal_aper_flux" for row in cols)
        assert len(code.sha256(original))==64
        assert code.download_once(item["mast_uri"],original)=="previously_downloaded"

        wrong_sector=root/"wrong_sector.fits"
        write_tglc(wrong_sector,sector=6)
        try:
            code.inspect_original(wrong_sector,item)
        except ValueError as exc:
            assert "sector" in str(exc)
        else:
            raise AssertionError("Wrong-sector FITS was accepted")

        wrong_tic=root/"wrong_tic.fits"
        write_tglc(wrong_tic,tic="TIC_90000409")
        try:
            code.inspect_original(wrong_tic,item)
        except ValueError as exc:
            assert "TIC" in str(exc)
        else:
            raise AssertionError("Wrong-TIC FITS was accepted")

        no_aperture=root/"no_aperture.fits"
        write_tglc(no_aperture,include_aperture=False)
        try:
            code.inspect_original(no_aperture,item)
        except ValueError as exc:
            assert "science extension" in str(exc)
        else:
            raise AssertionError("FITS without cal_aper_flux was accepted")

    print("[PASS] Frozen MAST URI and sector determine safe destination name")
    print("[PASS] Original FITS verified, 510 raw science rows")
    print("[PASS] TIC/sector/flux provenance rejects incorrect products")
    print("[PASS] Existing science FITS are reused, not redownloaded")
    print("All prospective MAST acquisition tests passed.")


if __name__=="__main__":
    main()
