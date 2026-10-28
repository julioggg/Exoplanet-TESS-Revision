#!/usr/bin/env python3
"""Offline synthetic regression for experiment 26 paired TIC-cluster audit.

NO MAST queries, NO private light curves, NO experiment or holdout stars.
Ensures 144 null draws are NOT treated as 144 independent host labels.
"""
from __future__ import annotations

import importlib.util
import math
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPT=Path(__file__).with_name("26_paired_development_score_audit.py")
spec=importlib.util.spec_from_file_location("exp26",SCRIPT)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def synthetic():
    ids=[f"TIC_{100000+i}" for i in range(36)]
    cal=set(ids[:24])
    internal=set(ids[24:])
    rows=[]
    inj=[]
    for i,tic in enumerate(ids):
        role="calibration" if tic in cal else "internal_development_audit"
        for block in m.BLOCKS:
            for draw in range(4):
                # Diverse deterministic nondegenerate scores with unique keys.
                key=((i*13+block//24+draw*7)%37)/37
                rows.append({
                    "tic_id":tic,"development_split":role,
                    "block_cadences":block,"draw":draw,
                    "raw_log_power":float(key),
                    "multievent_evidence":float(1.1*key)
                })
        for level in range(3):
            inj.append({
                "tic_id":tic,"development_split":role,
                "pair_id":f"{tic}_physical_{level}",
                "true_injected_period_recovered":(level==0),
                "raw_log_power":(.95 if level==0 and i%2 else .1),
                "multievent_evidence":(1.03 if level==0 else .1),
            })
    nul=pd.DataFrame(rows)
    injected=pd.DataFrame(inj)
    calibration=nul[nul.development_split=="calibration"]
    thresholds={}
    summary=[]
    for score in ["raw_log_power","multievent_evidence"]:
        t=float(np.quantile(calibration[score],.95,method="higher"))
        thresholds[score]=t
        for split in ["calibration","internal_development_audit"]:
            summary.append({
                "score_name":score,"source_split":split,
                "threshold_95pct_calibration_surrogate":t
            })
    rules={
        "calibration_null_tic_set":list(cal),
        "internal_audit_tic_set":list(internal)
    }
    return nul,injected,pd.DataFrame(summary),rules


def main():
    nul,inj,settings,rules=synthetic()
    frozen,tics=m.validate_inputs(nul,inj,settings,rules)
    assert len(tics)==12
    assert len(frozen)==2
    n,j,a,b,block,summary=m.evaluate(nul,inj,frozen,tics,n_boot=500)
    assert len(n)==144 and len(j)==36
    assert len(a)==12 and len(b)==12 and len(block)==3
    assert b.n.eq(3).all() and a.n.eq(12).all()
    assert (block.n_trials==48).all()
    assert (block.n_sources==12).all()
    assert summary["n_host_TIC_clusters"]==12
    assert summary["n_surrogate_null_trials"]==144
    assert summary["n_physical_injection_trials"]==36
    assert summary["physical_injections_with_correct_blind_period"]==12
    for key in (
        "surrogate_difference_paired_cluster",
        "correct_period_and_accepted_difference_paired_cluster"
    ):
        assert summary[key]["n_distinct_host_clusters"]==12
        assert 0<=summary[key]["source_level_signflip_exact_p_exploratory"]<=1
        assert np.isfinite(summary[key]["bootstrap95_difference_low"])
    rng=np.random.default_rng(4321)
    x=rng.uniform(0,1,size=12)
    z=np.zeros(12)
    assert m.source_cluster_bootstrap_ci(z,z,44)["mean_source_cluster_difference"]==0
    assert m.paired_signflip_exact(z,z)==1.
    assert math.isclose(m.paired_signflip_exact(z,x),2/2**12)
    # A single duplicated surrogate must invalidate source-structured data.
    doubled=pd.concat([nul,nul.iloc[[0]]],ignore_index=True)
    try:
        m.validate_inputs(doubled,inj,settings,rules)
    except ValueError as exc:
        assert "432" in str(exc) or "keys" in str(exc)
    else:
        raise AssertionError("Duplicated surrogate not detected")
    # Neither internal-audit data nor injection outcomes can change the
    # threshold, which derives ONLY from the frozen calibration surrogates.
    altered=settings.copy()
    altered.loc[
        altered.source_split=="internal_development_audit",
        "threshold_95pct_calibration_surrogate"
    ]*=1.01
    try:
        m.validate_inputs(nul,inj,altered,rules)
    except ValueError as exc:
        assert "threshold" in str(exc)
    else:
        raise AssertionError("Inconsistent frozen thresholds were accepted")

    print("[PASS] Frozen 24/12 source split and calibrated 95% threshold verified")
    print("[PASS] Paired 144 null / 36 injection data over 12 independent TICs")
    print("[PASS] Source-cluster bootstrap and exact 12-TIC sign-flip test")
    print("[PASS] Block-size sensitivity, duplicate and threshold-drift checks")
    print("[PASS] No MAST, private FITS, or locked holdout data opened")


if __name__=="__main__":
    main()
