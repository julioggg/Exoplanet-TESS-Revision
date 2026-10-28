#!/usr/bin/env python3
"""Offline tests for development-only surrogate-null BLS score calibration.

All source TIC identifiers and fluxes are SYNTHETIC. No real MAST/FITS
files, NASA planet labels, or final holdout files are ever opened.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

SRC=Path(__file__).with_name(
    "25_development_surrogate_null_calibration.py"
)
spec=importlib.util.spec_from_file_location("exp25",SRC)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def main():
    tics=[f"TIC_{i:09d}" for i in range(100_001,100_037)]
    split=m.assign_sources(tics)
    assert len(split)==36
    assert sum(x=="calibration" for x in split.values())==24
    assert sum(x=="internal_development_audit" for x in split.values())==12
    assert m.assign_sources(list(reversed(tics)))==split
    try:
        m.assign_sources(tics[:-1]+[tics[-2]])
    except ValueError as exc:
        assert "distinct" in str(exc)
    else:
        raise AssertionError("Duplicate synthetic source not rejected")

    t=1500.+np.arange(1200)/48.
    f=np.ones(len(t))+np.random.default_rng(123).normal(
        0,.0002,len(t)
    )
    q=np.zeros(len(t),dtype=np.int64)
    q[100:108]=1
    # Artificial >3-hour gap in timestamps: original time/gap grid unchanged
    t[700:]+=1.0
    surrogate,sflags,indices=m.block_surrogate(
        t,f,q,48,np.random.default_rng(456)
    )
    assert len(surrogate)==len(t)
    assert len(sflags)==len(t)
    assert len(indices)==len(t)
    assert np.array_equal(surrogate,f[indices])
    assert np.array_equal(sflags,q[indices])
    assert np.all(indices[:700]<700)
    assert np.all(indices[700:]>=700)
    assert not np.array_equal(indices,np.arange(len(t)))
    again=m.block_surrogate(
        t,f,q,48,np.random.default_rng(456)
    )[2]
    assert np.array_equal(indices,again)
    assert m.draw_seed(tics[0],24,0)==m.draw_seed(tics[0],24,0)
    assert m.draw_seed(tics[0],24,0)!=m.draw_seed(tics[0],24,1)

    sc=m.score_from_features(100.0,4,.2,.75)
    assert np.isclose(sc["raw_log_power"],np.log1p(100.0))
    assert 0<sc["multievent_evidence"]<sc["raw_log_power"]
    assert m.score_from_features(100.0,0,0,1)["multievent_evidence"]==0
    assert m.score_from_features(100.0,4,1,1)["multievent_evidence"]==0
    assert m.score_from_features(100.0,4,0,0)["multievent_evidence"]==0
    assert m.score_from_features(100.0,4,0,10.)["multievent_evidence"]==(
        m.score_from_features(100.0,4,0,1.)["multievent_evidence"]
    )
    try:
        m.score_from_features(1.,3,1.1,.5)
    except ValueError:
        pass
    else:
        raise AssertionError("Invalid flag ratio was accepted")

    rows=[]
    for tic in tics:
        for block in m.BLOCK_SIZES:
            for draw in range(m.DRAW_PER_SIZE):
                seed=m.draw_seed(tic,block,draw)
                rows.append({
                    "tic_id":tic,
                    "block_cadences":block,
                    "draw":draw,
                    "raw_log_power":float(
                        np.random.default_rng(seed).uniform(1,10)),
                    "multievent_evidence":float(
                        np.random.default_rng(seed+1).uniform(0,4))
                })
    null=pd.DataFrame(rows)
    assert len(null)==432
    assert null.groupby("tic_id").size().eq(12).all()
    for name in m.SCORES:
        thr=m.threshold_from_surrogates(null,split,name)
        assert np.isfinite(thr) and thr>=0
        measured,star=m.audit_surrogate_exceedance(null,split,name,thr)
        assert len(measured)==432
        assert len(star)==36
        assert star.n_surrogates.eq(12).all()
        assert measured.exceeds_threshold.dtype==bool
        assert (
            measured[measured.development_split=="calibration"]
            .exceeds_threshold.mean()<.07
        )
    partial=null[~(
        (null.tic_id==tics[0])&(null.draw==0)&
        (null.block_cadences==24)
    )]
    try:
        m.threshold_from_surrogates(partial,split,"raw_log_power")
    except ValueError as exc:
        assert "Incomplete" in str(exc)
    else:
        raise AssertionError("Incomplete surrogate calibration was accepted")

    print("[PASS] Deterministic disjoint 24/12 development source split")
    print("[PASS] Flux/TGLC jointly resampled within intact time-gap segments")
    print("[PASS] Fixed evidence score, invalid flag rejection and zero cases")
    print("[PASS] Calibrate 95th percentile ONLY on designated surrogate TICs")
    print("[PASS] Independent internal null exceedance tally per source")
    print("[PASS] Missing null replicate rejected; no final holdout access")


if __name__=="__main__":
    main()
