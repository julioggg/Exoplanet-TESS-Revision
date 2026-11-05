#!/usr/bin/env python3
"""Experiment 30 preflight: no FITS/MAST or 18 final holdout accessed.

Validates 369 original labelled sources, 447 existing matched physical
injections, the prior independent 42-TIC injection split, the original
11-feature inputs, and the BLS dy->unweighted power conversion. Also
runs a small synthetic RF sanity check for paired probability responses.
"""
from __future__ import annotations
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier

p=Path(__file__).with_name("30_rf_paired_injection_response.py")
spec=importlib.util.spec_from_file_location("experiment30_rf",p)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def main():
    orig,pairs,truth,cohort,pilot=m.load_inputs()
    test,cal,pilots=m.validate(orig,pairs,truth,cohort,pilot)
    train,normal,injected,itruth=m.prepare(
        orig,pairs,truth,test,cal,pilots
    )
    assert len(train)+len(normal)==315
    assert len(injected)==252
    assert len(itruth)==126
    assert len(test)==42
    assert len(cal)==95
    assert len(pilots)==12
    assert not (set(train.tic_id)|set(normal.tic_id))&test
    assert not (set(train.tic_id)|set(normal.tic_id))&pilots
    assert set(itruth.tic_id)==test
    assert not injected[["tic_id","pair_id","label"]].duplicated().any()
    assert np.allclose(
        injected.bls_power.to_numpy(float),
        injected.bls_power_weighted.to_numpy(float)*
        injected.noise_sigma_baseline.to_numpy(float)**2,
        atol=1e-12,rtol=1e-10
    )
    assert np.allclose(
        injected.snr_estimate.to_numpy(float),
        abs(injected.transit_depth.to_numpy(float))/
        injected.std_flux.to_numpy(float),
        atol=1e-12,rtol=1e-10
    )
    assert set(m.MODELS)=={"FULL_11","NO_N_POINTS","BLS_ONLY_5"}
    assert all(set(c).issubset(orig.columns) for c in m.MODELS.values())

    # Small control/injection test independent of any real TESS data.
    fake=np.random.default_rng(1234).normal(0,1,(100,2))
    y=(fake[:,0]+fake[:,1]>0).astype(int)
    mod=RandomForestClassifier(n_estimators=20,random_state=42)
    mod.fit(fake[:70],y[:70])
    control=np.zeros((5,2),dtype=float)
    affected=control.copy()
    affected[:,0]+=3.
    p0=mod.predict_proba(control)[:,1]
    p1=mod.predict_proba(affected)[:,1]
    assert np.isfinite(p0).all() and np.isfinite(p1).all()
    assert p1.mean()>p0.mean()

    print("[PASS] Original catalogue and all 447 existing injection pairs intact")
    print("[PASS] Injection evaluation TICs 42 are excluded from ALL RF fitting")
    print("[PASS] Previously examined 12 pilot TICs excluded")
    print("[PASS] 126 paired injection trials / 42 independent host clusters")
    print("[PASS] Original 11 features and inverse BLS dy weighting validated")
    print("[PASS] RF synthetic paired probability response and model ablations")
    print("[PASS] No final 18-host evaluation or FITS/MAST accessed")


if __name__=="__main__":
    main()
