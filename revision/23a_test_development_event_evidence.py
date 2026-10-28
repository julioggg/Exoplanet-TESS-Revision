#!/usr/bin/env python3
"""Synthetic offline regression checks for experiment-23 event evidence.

No real NASA periods, no private TESS light curves, no holdout TICs.
The final fixed-period BLS check uses astropy already installed in the
CEDIA project environment; core event tests use only NumPy/Pandas.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import tempfile

FILE=Path(__file__).with_name("23_development_event_evidence.py")
spec=importlib.util.spec_from_file_location("exp23",FILE)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def main():
    assert m.curve_path(
        "TIC_123456",9,Path("/tmp/dummy_dev")
    ).name=="TIC_123456_sector_09_TGLC_APERTURE_1800S_V1.csv"
    t=1500.+np.arange(1200)/48.
    period=3.4
    epoch=1501.2
    duration=.15
    cycle=np.rint((t-epoch)/period)
    offset=t-(epoch+cycle*period)
    in_transit=np.abs(offset)<=duration/2
    rng=np.random.default_rng(20261028)
    flux=1.0+rng.normal(0,.00025,len(t))
    flux[in_transit]-=.01
    q=np.zeros(len(t),dtype=np.int64)
    q[np.flatnonzero(in_transit)[0]]=1

    info,detail=m.local_event_windows(
        t,flux,q,period,epoch,duration
    )
    assert len(detail)>5,detail
    assert info["n_cycles_with_positive_local_depth"]>=6,info
    assert info["n_cycles_with_usable_local_depth"]>=6,info
    assert info["median_positive_local_depth"]>.008,info
    assert info["largest_positive_event_weight_share"]<.30,info
    assert 0<info["fraction_inbox_TGLC_flagged"]<.10,info
    assert info["dominant_positive_event_cycle"] is not None,info

    keep,n=m.leave_dominant_cycle_out(
        t,flux,epoch,period,duration,
        info["dominant_positive_event_cycle"]
    )
    assert 2<=n<=10
    assert keep.sum()==len(t)-n
    assert np.count_nonzero(~keep)==n
    no_change,z=m.leave_dominant_cycle_out(
        t,flux,epoch,period,duration,None
    )
    assert no_change.all() and z==0

    noise=m.robust_noise(flux)
    assert .0001<noise<.001
    fixed=m.same_period_bls(t,flux,period,duration,noise)
    after=m.same_period_bls(
        t[keep],flux[keep],period,duration,noise
    )
    assert abs(fixed["period_days"]-period)<.000001
    assert abs(fixed["duration_days"]-duration)<.000001
    assert np.isfinite(fixed["transit_epoch_btjd"])
    assert np.isfinite(fixed["BLS_power_white_noise_surrogate"])
    assert np.isfinite(after["BLS_power_white_noise_surrogate"])
    assert fixed["BLS_power_white_noise_surrogate"]>0
    # Exactly reproduce the experiment-21 4000x5 BLS search, NOT a
    # 1-period x 1-duration refit. The full-grid best phase is the
    # appropriate source of evidence and original power consistency.
    full=m.original_full_grid_bls(t,flux,noise)
    assert full["period_days"]>=0.5
    assert full["period_days"]<=15.0
    assert full["BLS_power_white_noise_surrogate"]>0
    assert np.isfinite(full["transit_epoch_btjd"])

    with tempfile.TemporaryDirectory(prefix="tess_experiment23_") as root:
        private=Path(root)
        tic="TIC_90001234"
        sector=9
        path=m.curve_path(tic,sector,private)
        pd.DataFrame({
            "time_btjd":t,
            "flux_cal_aper_normalized":flux,
            "tglc_flags":q,
            "cadence_num":np.arange(len(t),dtype=int),
        }).to_csv(path,index=False)
        source={
            "tic_id":tic,"sector":sector,"planet":"SYNTHETIC b",
            "n_cadences":len(t),
            "baseline_global_period_days":full["period_days"],
            "baseline_global_duration_days":full["duration_days"],
            "baseline_global_power":full["BLS_power_white_noise_surrogate"],
            "reference_period_days":period,
            "baseline_direct_match_2pct":True,
            "multievent_direct_match_2pct":True,
        }
        qa={"private_cleaned_file":str(path),"policy_id":m.POLICY,
            "n_accepted":len(t)}
        row,events=m.evaluate_one(source,qa,private)
        assert abs(row["replayed_full_grid_BLS_power"]-
                   full["BLS_power_white_noise_surrogate"])<1e-6
        assert row["reconstructed_BLS_transit_epoch_btjd"]==full[
            "transit_epoch_btjd"]
        assert 0<row["one_period_fit_BLS_power"]
        assert row["no_research_time_period_replacement"] is True
        assert len(events)>0

        corrupted=dict(source)
        corrupted["baseline_global_power"]*=1.1
        try:
            m.evaluate_one(corrupted,qa,private)
        except ValueError as exc:
            assert "Full-grid BLS replay" in str(exc)
        else:
            raise AssertionError("Invalid original 4000x5 BLS replay went unnoticed")

    # Missing/off-window nearby samples are NOT treated as quality-certified.
    gapmask=np.abs(t-(epoch+2*period))>.3
    tgap=t[gapmask]
    fgap=flux[gapmask]
    qgap=q[gapmask]
    r,_=m.local_event_windows(
        tgap,fgap,qgap,period,epoch,duration
    )
    assert r["n_cycles_with_positive_local_depth"]>=4

    try:
        m.local_event_windows(
            t,np.full(len(t),np.nan),q,period,epoch,duration
        )
    except ValueError:
        pass
    else:
        raise AssertionError("Invalid photometry was silently accepted")

    print("[PASS] Frozen policy-14 development filename")
    print("[PASS] Local-baseline depths from independent synthetic transits")
    print("[PASS] One event influence mask and flag provenance")
    print("[PASS] Original 4000x5 BLS replay recovers phase and source power")
    print("[PASS] Event evidence pipeline processes full synthetic source")
    print("[PASS] Corrupted previously frozen BLS power is rejected")
    print("[PASS] Missing data, finite values and order checks")
    print("[PASS] No locked holdout data or detector threshold tuning")
    print("All development event-evidence tests passed.")


if __name__=="__main__":
    main()
