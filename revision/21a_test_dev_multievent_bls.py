#!/usr/bin/env python3
"""Offline tests of experiment-21 blind multi-event BLS; NEVER opens holdout.

Synthetic signal verifies that a nominally stronger isolated-event period
can be refused while the prespecified repeated-event candidate is accepted.
Also checks exact frozen policy file names, TGLC flag veto, and the real
4000-period BLS scan.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path
from types import SimpleNamespace

import numpy as np

CODE=Path(__file__).with_name("21_dev_blind_multievent_bls.py")
spec=importlib.util.spec_from_file_location("exoplanet_dev21",CODE)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def main():
    assert m.file_for(
        "TIC_98765",4,Path("/tmp/engineering_test")
    ).name=="TIC_98765_sector_04_TGLC_APERTURE_1800S_V1.csv"
    rng=np.random.default_rng(20261027)
    t=1500.0+np.arange(1200,dtype=float)/48.0
    period=3.2
    epoch=1501.0
    duration=.15
    cyc=np.rint((t-epoch)/period)
    residual=t-(epoch+cyc*period)
    dip=np.abs(residual)<=duration/2
    y=1.0+rng.normal(0,0.00030,len(t))
    y[dip]-=.02
    flags=np.zeros(len(t),dtype=np.int32)
    good=m.event_support(t,y,flags,period,epoch,duration)
    assert good["n_supported_events"]>=7,good
    assert good["n_positive_events"]>=7,good
    assert good["largest_positive_drop_share"]<.25,good
    assert m.accepted_events(good)

    single=m.event_support(t,y,flags,14.0,1502.0,.10)
    assert single["n_supported_events"]<3,single
    assert not m.accepted_events(single)

    trial=SimpleNamespace(
        period=np.array([14.0,period]),
        transit_time=np.array([1502.0,epoch]),
        duration=np.array([.10,duration]),
        depth=np.array([.2,.02]),
        depth_snr=np.array([300.0,100.0]),
    )
    peak,desc,multi,visits=m.choose_candidates(
        t,y,flags,trial,np.array([200.0,160.0])
    )
    assert peak==0
    assert multi is not None and multi[0]==1
    assert visits==2
    assert desc["n_supported_events"]<3

    # Flag every point in the true transit windows; the same signal must
    # fail the prespecified contamination sensitivity rule.
    marked=flags.copy()
    marked[dip]=1
    contaminated=m.event_support(t,y,marked,period,epoch,duration)
    assert contaminated["fraction_TGLC_flagged_inside"]==1.0
    assert not m.accepted_events(contaminated)

    # Real Astropy BLS smoke test uses the frozen 4000 x 5 grid. It is not
    # fitted on the NASA development or test stars.
    scan,power,noise=m.blind_bls(t,y)
    assert len(power)==4000
    assert noise>0
    original=int(np.nanargmax(power))
    direct,alias=m.match_period(float(scan.period[original]),period)
    assert direct or alias, (
        f"Strong synthetic transit failed BLS: P_BLS={scan.period[original]}"
    )
    assert m.match_period(2.5,5.0)==(False,True)
    assert m.match_period(10.,5.0)==(False,True)
    assert m.match_period(float("nan"),5.0)==(False,False)
    assert m.exact_binomial_ci(0,36)[0]==0.0
    assert m.exact_binomial_ci(36,36)[1]==1.0

    print("[PASS] Exact original TGLC policy-14 cleaned filename")
    print("[PASS] 7+ independent synthetic events, positive depth and coherence")
    print("[PASS] Highest isolated-event peak rejected; repeated peak selected")
    print("[PASS] TGLC-flag veto is label-blind and transparent")
    print("[PASS] Original full-grid BLS remains unchanged")
    print("[PASS] Period aliases and exact binomial bounds")
    print("All development-only multi-event BLS synthetic tests passed.")


if __name__=="__main__":
    main()
