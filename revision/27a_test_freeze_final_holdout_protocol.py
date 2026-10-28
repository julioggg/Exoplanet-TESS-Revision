#!/usr/bin/env python3
"""Offline preflight test of experiment-27 freeze BEFORE opening heldout FITS.

All 54 TICs below are SYNTHETIC and intentionally fictitious.
No original or final TESS photometry or MAST source accessed.
"""
from __future__ import annotations
import importlib.util
from pathlib import Path
import numpy as np
import pandas as pd

FILE=Path(__file__).with_name("27_freeze_final_holdout_protocol.py")
spec=importlib.util.spec_from_file_location("exp27_freeze",FILE)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def main():
    def sample(num,prefix,role):
        return pd.DataFrame([
            {
                "tic_id":f"TIC_{prefix*100000+i}",
                "role":role,
                "is_frozen_test_tic":role=="LOCKED_HOLDOUT",
                "selection_rank":i+1,
                "period_bin":["1.5-3.5","3.5-6.5","6.5-9.5"][
                    i%(3)
                ],
                "author":"TGLC",
                "exptime_seconds":1800.0,
                "sector":(i%12)+1,
                "mast_uri":"mast:HLSP/tglc/s0001/dummy",
            } for i in range(num)
        ])
    dev=sample(36,5,"development")
    hold=sample(18,6,"LOCKED_HOLDOUT")
    prior=pd.DataFrame({"tic_id":[
        f"TIC_{800000+i}" for i in range(18)
    ]})
    old=pd.DataFrame({"tic_id":[
        f"TIC_{900000+i}" for i in range(369)
    ]})
    m.validate_selected(dev,hold,old,prior)
    assert dev.period_bin.value_counts().eq(12).all()
    assert hold.period_bin.value_counts().eq(6).all()
    corrupted=hold.copy()
    corrupted.loc[0,"tic_id"]=dev.tic_id.iloc[0]
    try:
        m.validate_selected(dev,corrupted,old,prior)
    except ValueError as e:
        assert "Leakage" in str(e)
    else:
        raise AssertionError("Development-to-holdout overlap accepted")

    changed=hold.copy()
    changed.loc[0,"is_frozen_test_tic"]=False
    try:
        m.validate_selected(dev,changed,old,prior)
    except ValueError as e:
        assert "flag" in str(e)
    else:
        raise AssertionError("Unmarked final star accepted")

    altered=hold.copy()
    altered.loc[0,"exptime_seconds"]=120.
    try:
        m.validate_selected(dev,altered,old,prior)
    except ValueError as e:
        assert "1800" in str(e)
    else:
        raise AssertionError("Changed cadence accepted")

    thresholds=m.fixed_thresholds()
    assert set(thresholds)=={"raw_log_power","multievent_evidence"}
    assert np.isclose(thresholds["raw_log_power"],10.93563724046706,atol=1e-9)
    assert np.isclose(thresholds["multievent_evidence"],5.329805494852721,atol=1e-9)

    # This is PRE-FREEZE metadata inspection only. It is forbidden to
    # open the actual 18 final TGLC FITS at this stage.
    protocol=m.prepare()
    assert protocol["experiment"]==27
    assert len(protocol["locked_holdout_tics"])==18
    assert len(protocol["frozen_sources_hashes"])==6
    assert protocol["final_test_permitted_only_once"]
    assert protocol["no_fits_or_photometric_holdout_access_in_freeze"]
    assert protocol["main_thresholds"]==thresholds
    assert all(len(x)==64 for x in protocol["frozen_sources_hashes"].values())

    print("[PASS] Synthetic 36/18 disjoint TGLC source-role freeze")
    print("[PASS] Any development-to-holdout leakage rejected")
    print("[PASS] Changed TGLC exposure or lock flag rejected")
    print("[PASS] Exactly original experiment-25 95%-quantile thresholds")
    print("[PASS] Fully specified code+data hashes and final endpoints")
    print("[PASS] Metadata-only freeze; NO heldout flux accessed")


if __name__=="__main__":
    main()
