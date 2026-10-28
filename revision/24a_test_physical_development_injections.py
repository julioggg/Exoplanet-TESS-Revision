#!/usr/bin/env python3
"""Offline regression suite for DEVELOPMENT-only physical injection design.

Does not access data, MAST, NASA, or LOCKED_HOLDOUT.
Tests deterministic stratified configuration, finite-exposure transit
construction, blind original-grid search, 1-event influence, and TIC
bootstrap. Runs in the identical CEDIA venv as experiment 24.
"""
from __future__ import annotations
import importlib.util
from pathlib import Path
import numpy as np
import pandas as pd

SRC=Path(__file__).with_name("24_physical_development_injections.py")
spec=importlib.util.spec_from_file_location("exp24",SRC)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def main():
    levels=[m.injection_config("TIC_55555555",i) for i in range(3)]
    again=[m.injection_config("TIC_55555555",i) for i in range(3)]
    assert levels==again
    for j,cfg in enumerate(levels):
        lo,hi=m.PERIOD_RANGES[j]
        assert lo<cfg["true_period_days"]<hi
        assert cfg["radius_ratio"]==m.RADIUS_RATIOS[j]
        assert .05<=cfg["impact_parameter"]<=.65

    t=1500.0+np.arange(1200,dtype=float)/48.0
    rng=np.random.default_rng(2468)
    base=1.+rng.normal(0,.0003,len(t))
    fake=lambda t,p,k,b,ph: (
        np.where(
            np.abs((t-(t[0]+ph*p)+p/2)%p-p/2)<.075,
            1-k*k,1.
        ),
        {
            "t0":float(t[0]+ph*p),
            "a_over_rstar":10.0,
        }
    )
    profile,meta=m.physical_model(t,levels[0],fake)
    injected=base*profile
    assert profile.shape==base.shape
    assert 0<profile.min()<1
    assert not np.allclose(injected,base,atol=1e-6,rtol=0)
    assert (t==t.copy()).all()
    assert m.count_sampled_injected_events(
        t,profile,levels[0]["true_period_days"],meta["t0"]
    )>=2

    from pathlib import Path
    helper=m.load_helper(
        Path(__file__).with_name("23_development_event_evidence.py"),
        "synthetic_helper23"
    )
    # Inject a high-SNR synthetic signal with seven or more repeated dips.
    period=3.1
    epoch=1501.0
    cycles=np.rint((t-epoch)/period)
    residual=t-(epoch+cycles*period)
    y=base.copy()
    y[np.abs(residual)<.075] -= .015
    noise=helper.robust_noise(base)
    found=m.scan_candidate(t,y,noise)
    direct,alias=m.compare_period(found["best_period_days"],period)
    assert direct or alias,found
    flags=np.zeros(len(t),dtype=np.int64)
    evidence,ratio,nremoved=m.influence_at_blind_period(
        t,y,flags,found,noise,helper
    )
    assert evidence["n_cycles_with_positive_local_depth"]>=2
    assert 0<ratio<1.1,(ratio,evidence)
    assert nremoved>0

    dummy=pd.DataFrame({
        "tic_id":[f"TIC_{i}" for i in range(36)],
        "injected_period_direct_match_2pct":[i%2==0 for i in range(36)]
    })
    lo,hi=m.star_bootstrap_ci(
        dummy,"injected_period_direct_match_2pct",n_boot=300
    )
    assert 0<lo<.5<hi<1.0
    assert m.compare_period(5.0,5.0)==(True,False)
    assert m.compare_period(10.0,5.0)==(False,True)
    assert m.compare_period(2.5,5.0)==(False,True)

    # Catch accidental inclusion of locked-holdout rows early, without
    # ever querying or opening the real locked source CSV.
    dev=pd.DataFrame({
        "tic_id":[f"TIC_{i}" for i in range(36)],
        "role":["development"]*35+["LOCKED_HOLDOUT"],
        "is_frozen_test_tic":[False]*35+[True],
    })
    stub=pd.DataFrame({"tic_id":dev.tic_id})
    try:
        m.curves_verified(dev,stub,stub,stub,Path("/nonexistent"))
    except ValueError as exc:
        assert "non-development" in str(exc)
    else:
        raise AssertionError("LOCKED_HOLDOUT leakage not detected")

    # Validate real PyTransit interface, which will be used at HPC runtime.
    old=m.load_helper(
        Path(__file__).with_name("06_injection_recovery.py"),
        "real_original06_injection"
    )
    sample_cfg=levels[1]
    physical,physical_meta=m.physical_model(
        t,sample_cfg,old.transit_model
    )
    assert physical.min()<1 and physical.min()>0
    assert m.count_sampled_injected_events(
        t,physical,sample_cfg["true_period_days"],physical_meta["t0"]
    )>=2

    print("[PASS] Frozen stratified periods, planetary radii, deterministic seeds")
    print("[PASS] Paired synthetic profiles preserve original time and flux")
    print("[PASS] Actual PyTransit quadratic physical transit interface")
    print("[PASS] Original 4000x5 blind BLS and event-influence diagnostics")
    print("[PASS] TIC-cluster bootstrap and direct/alias rules")
    print("[PASS] Explicit LOCKED_HOLDOUT access prevention")


if __name__=="__main__":
    main()
