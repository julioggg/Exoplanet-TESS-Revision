#!/usr/bin/env python3
"""Offline regression test for new, leakage-free NASA TESS host inventory.

No NASA/MAST network connection is necessary. Tests deterministic queue,
full TIC exclusion, balanced orbital-period interleaving and exact MAST ID /
author / exposure / sector checks.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd


SCRIPT=Path(__file__).with_name("15_new_confirmed_tess_cohort.py")
spec=importlib.util.spec_from_file_location("new_tess_hosts",SCRIPT)
module=importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class FakeQuantity:
    def __init__(self,value):
        self.value=value


class FakeSearchResult:
    def __init__(self):
        self.target_name=["90001201","90001201","90001201","90001202"]
        self.mission=[
            "TESS Sector 02", "TESS Sector 25",
            "TESS Sector 27", "TESS Sector 02"
        ]
        self.author=["TGLC"]*4
        self.exptime=[FakeQuantity(1800)]*4
        self.table=pd.DataFrame({
            "dataURI":[
                "mast:HLSP/tglc/s0002/one.fits",
                "mast:HLSP/tglc/s0025/two.fits",
                "mast:HLSP/tglc/s0027/three.fits",
                "mast:HLSP/tglc/s0002/wrong-tic.fits"
            ],
            "productFilename":[
                "one.fits","two.fits","three.fits","wrong-tic.fits"
            ],
            "obsid":["1","2","3","4"],
        })
        # Astropy tables expose .colnames but pandas tables do not.
        self.table.colnames=list(self.table.columns)
    def __len__(self):
        return len(self.target_name)


class FakeLightkurve:
    def search_lightcurve(self,*args,**kwargs):
        assert kwargs["mission"]=="TESS"
        assert kwargs["author"]=="TGLC"
        assert kwargs["exptime"]==1800
        return FakeSearchResult()


def main():
    assert module.normalized_tic("TIC 90001201")=="TIC_90001201"
    assert module.normalized_tic("90001201")=="TIC_90001201"
    assert module.normalized_tic("unmapped-source")==""
    records=[]
    for idx in range(24):
        records.append({
            "pl_name":f"New Test Planet {idx}",
            "hostname":f"Test Star {idx}",
            "tic_id":f"TIC {90001000+idx}",
            "pl_orbper": [2.0,4.0,7.0][idx % 3],
            "pl_tranmid":2458000.0,
            "pl_controv_flag":0,
            "default_flag":1,
            "tran_flag":1,
        })
    data=pd.DataFrame(records)
    excluded={"TIC_90001000","TIC_90001002"}
    queue=module.prepare_queue(data,excluded)
    again=module.prepare_queue(data,excluded)
    assert queue.equals(again)
    assert len(queue)==22
    assert not set(queue.tic_norm).intersection(excluded)
    assert queue.queue_rank.is_unique
    assert set(queue.bin)==set(module.PERIOD_LABELS)
    # Exact matching drops wrong TIC and extended mission sector 27.
    rows=module.metadata_rows({"tic_norm":"TIC_90001201"},FakeLightkurve())
    assert len(rows)==2
    assert [row["sector"] for row in rows]==[2,25]
    assert rows[0]["mast_uri"].startswith("mast:")
    assert all(row["tic_id"]=="TIC_90001201" for row in rows)
    print("[PASS] Deterministic period-stratified new-TIC host queue")
    print("[PASS] All excluded prior TICs removed")
    print("[PASS] Exact target, author, 1800-s and primary-sector filtering")
    print("[PASS] Archive product URI inventory without downloading FITS")


if __name__=="__main__":
    main()
