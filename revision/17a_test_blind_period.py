#!/usr/bin/env python3
"""Synthetic offline smoke tests for the FROZEN blind TESS period scan.

No MAST access or reference planet ephemerides are passed into the BLS scan.
A synthetic Gaussian-noise + box-transit light curve must permit approximate
period recovery under the published, unmodified grid and duration choices.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np

SCRIPT=Path(__file__).with_name("17_blind_confirmed_host_period_recovery.py")
spec=importlib.util.spec_from_file_location("period_scan17",SCRIPT)
mod=importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def main():
    true_period=3.4
    rng=np.random.default_rng(20261027)
    t=1500 + np.arange(1200)/48.0
    phase=(t-(1500+1.22)+true_period/2)%true_period-true_period/2
    flux=1.0-0.025*(np.abs(phase)<=0.10/2)
    flux=flux+rng.normal(0,0.0009,len(flux))
    fit,power=mod.blind_scan(t,flux)
    recovered,alias=mod.period_match(
        fit["best_period_days"],true_period
    )
    assert recovered, (
        f"Strong synthetic transit was not recovered: "
        f"{fit['best_period_days']:.6f} vs {true_period:.6f}"
    )
    assert not alias
    assert len(power)==4000
    assert mod.period_match(1.70,3.4)==(False,True)
    assert mod.period_match(6.80,3.4)==(False,True)
    assert mod.period_match(4.4,3.4)==(False,False)
    assert mod.binomial_exact(0,18)[0]==0
    assert mod.binomial_exact(18,18)[1]==1
    assert mod.binomial_exact(9,18)[0]<.5
    assert mod.binomial_exact(9,18)[1]>.5
    assert 0<fit["point_to_point_white_noise_proxy"]<.01
    print("[PASS] Blind BLS recovers strong synthetic P=3.4-day transit")
    print("[PASS] Predeclared direct and 2:1 period-match criteria")
    print("[PASS] Exact binomial intervals including boundary cases")
    print("[PASS] BLS did not receive catalogue ephemerides or periods")
    print("All blind confirmed-host period benchmark tests passed.")


if __name__=="__main__":
    main()
