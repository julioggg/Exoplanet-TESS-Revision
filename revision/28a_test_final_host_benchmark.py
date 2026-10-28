#!/usr/bin/env python3
"""Synthetic test for final 18-host TESS evaluation BEFORE protocol freeze.

NO real holdout names, MAST, raw FITS or selected photometry is opened.
Uses fake source TICs and synthetic transit profiles in a temporary folder,
then verifies the original BLS and event score helper interfaces.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd

SRC=Path(__file__).with_name("28_final_host_benchmark.py")
spec=importlib.util.spec_from_file_location("exp28_test",SRC)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def main():
    assert m.match_period(4.0,4.0)==(True,False)
    assert m.match_period(2.0,4.0)==(False,True)
    assert m.match_period(8.0,4.0)==(False,True)
    assert m.match_period(7.,4.)==(False,False)
    assert m.exact_interval(0,18)[0]==0
    assert m.exact_interval(18,18)[1]==1.0

    hosts=[]
    model=[]
    for i in range(18):
        tic=f"TIC_{1000000+i}"
        period=2.5+(i//6)*2
        pred=period if i<6 else period*(1.3 if i%2 else .6)
        hosts.append({
            "tic_id":tic,"planet":f"Synthetic{i} b",
            "planet_period_days":period,
            "period_bin":["1.5-3.5","3.5-6.5","6.5-9.5"][i//6],
        })
        model.append({
            "tic_id":tic,"best_BLS_period_days":pred,
            "raw_log_power":12.0 if i<9 else 5.0,
            "multievent_evidence":6.0 if i<6 else 2.0
        })
    scores={"raw_log_power":10.93563724046706,
            "multievent_evidence":5.329805494852721}
    a,b,r=m.summarise(pd.DataFrame(model),pd.DataFrame(hosts),scores)
    assert len(a)==18 and len(b)==3
    assert r["n_direct_blind_period_agreements"]==6
    assert r["n_correct_and_raw_threshold_accepted"]==6
    assert r["n_correct_and_evidence_threshold_accepted"]==6
    assert r["n_evidence_threshold_accepted"]==6
    assert r["n_raw_threshold_accepted"]==9
    assert r["n_wrong_period_but_raw_threshold_accepted"]==3
    assert r["n_wrong_period_but_evidence_threshold_accepted"]==0
    assert r["NO_estimation_of_detection_sensitivity"]
    assert a.direct_period_match_2pct.sum()==6

    exp23=m.import_module(
        "test_frozen23","revision/23_development_event_evidence.py"
    )
    exp25=m.import_module(
        "test_frozen25","revision/25_development_surrogate_null_calibration.py"
    )
    t=1500.+np.arange(1200)/48.
    p=3.17
    epoch=1501.
    cyc=np.rint((t-epoch)/p)
    in_tr=np.abs(t-(epoch+cyc*p))<.075
    flux=1.+np.random.default_rng(2026).normal(0,.0003,len(t))
    flux[in_tr]-=.015
    with TemporaryDirectory(prefix="tess_fake_final_") as temp:
        root=Path(temp)
        tic="TIC_987654321"
        fname=root/f"{tic}_sector_09_{m.POLICY}.csv"
        pd.DataFrame({
            "time_btjd":t,
            "flux_cal_aper_normalized":flux,
            "tglc_flags":np.zeros(len(t),dtype=int),
            "cadence_num":np.arange(len(t)),
        }).to_csv(fname,index=False)
        row={"tic_id":tic,"sector":9}
        qa={"private_cleaned_file":str(fname),
            "n_accepted":len(t),"policy_id":m.POLICY}
        fits={"sha256":"f"*64}
        # No true period is provided to blind_one.
        result=m.blind_one(row,qa,fits,root,exp23,exp25)
        assert result["tic_id"]==tic
        assert len(result["cleaned_curve_sha256"])==64
        assert result["original_fits_sha256"]=="f"*64
        assert result["raw_log_power"]>0
        assert result["multievent_evidence"]>=0
        assert abs(result["best_BLS_period_days"]-p)/p<.02
    print("[PASS] Blind 18-source period/alias scoring and exact CI")
    print("[PASS] Frozen raw and event-evidence threshold reporting")
    print("[PASS] Original 4000x5 TESS BLS and event score on synthetic data")
    print("[PASS] No NASA period supplied to blind source scanning")
    print("[PASS] No final locked source FITS/flux or real MAST queries")


if __name__=="__main__":
    main()
