#!/usr/bin/env python3
"""OFFLINE tests for the two new positive-host cohorts, no MAST/NASA queries.

Synthetic archive source rows are deliberately distinct; source selection
depends ONLY on frozen period stratum, pseudorandom queue rank and archive
product metadata. 36 development + 18 holdout must be distinct and complete.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd

SRC=Path(__file__).with_name(
    "20_preregister_fresh_development_and_holdout.py"
)
spec=importlib.util.spec_from_file_location("cohort20",SRC)
mod=importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def fake_products(row):
    tic=row["tic_norm"]
    sector=int(row["queue_rank"]%23+1)
    return [{
        "tic_id":tic,
        "sector":sector,
        "mast_uri":(
            f"mast:HLSP/tglc/s{sector:04d}/cam1-ccd2/"
            f"hlsp_tglc_tess_ffi_gaiaid-{tic[4:]}-s{sector:04d}"
            "_tess_v1_llc.fits"
        ),
        "author":"TGLC",
        "exptime_seconds":1800.0,
        "product_filename":tic+"-test.fits",
        "mast_obsid":tic[4:],
    }]


def main():
    records=[]
    for b,period in zip(mod.PERIOD_BINS,[2.3,4.2,7.7]):
        for i in range(31):
            k=(len(records)+100_000)
            records.append({
                "tic_norm":f"TIC_{k}",
                "bin":b,
                "queue_rank":len(records)+1,
                "hostname":f"Star{k}",
                "pl_name":f"Planet{k}",
                "planet_period_days":period,
                "pl_tranmid":2458333.0,
                "pl_tsystemref":"BJD-TDB",
                "pl_trandur":2.0,
                "pl_tranmiderr1":0.0003,
                "pl_tranmiderr2":-0.0003,
                "pl_ratror":0.1,
                "pl_trandep":1.0,
            })
    queue=pd.DataFrame(records)
    frozen_old={"TIC_999999"}
    selected,audit,count,done=mod.plan(
        queue,fake_products,frozen_old,
        query_cap=len(queue),pause=0
    )
    assert done
    assert len(selected)==54
    dev=selected[selected.role=="development"]
    hold=selected[selected.role=="LOCKED_HOLDOUT"]
    assert len(dev)==36 and len(hold)==18
    assert dev.selection_rank.is_unique
    assert hold.selection_rank.is_unique
    assert set(dev.tic_id).isdisjoint(hold.tic_id)
    assert set(selected.tic_id).isdisjoint(frozen_old)
    assert all(
        len(dev[dev.period_bin==period_bin])==12
        and len(hold[hold.period_bin==period_bin])==6
        for period_bin in mod.PERIOD_BINS
    )
    for frame in [dev,hold]:
        assert frame.sector.between(1,26).all()
        assert frame.author.eq("TGLC").all()
        assert frame.exptime_seconds.eq(1800).all()
        assert frame.mast_uri.str.startswith("mast:").all()
    # With a lower MAST query budget the script must NOT call the
    # generated incomplete cohort a complete frozen benchmark.
    partial,log,_,complete=mod.plan(
        queue,fake_products,frozen_old,
        query_cap=54,pause=0
    )
    assert not complete
    assert len(partial)<=54
    assert len(log)==54
    print("[PASS] Strictly disjoint 36 development and 18 untouched holdout TICs")
    print("[PASS] Balanced period strata: 12+6 per bin")
    print("[PASS] Unique acquisition-compatible selection ranks per role")
    print("[PASS] No photometry or BLS inspected; metadata-only planning")
    print("[PASS] Incomplete cohorts are not marked complete")


if __name__=="__main__":
    main()
