#!/usr/bin/env python3
"""Offline tests for non-retuning BLS peak forensic event support.

One fake BLS period and frozen curve with two/one dips show event counting.
No catalog download, no MAST query, no BLS rescan or target selection.
"""
from __future__ import annotations
import importlib.util
from pathlib import Path

import numpy as np

FILE=Path(__file__).with_name("19_forensic_bls_event_support.py")
spec=importlib.util.spec_from_file_location("frozen_bls19",FILE)
mod=importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def main():
    t=1400.0+np.arange(1200,dtype=float)/48
    period=7.0
    epoch=1401.0
    residual,k=mod.phase_and_cycle(t,period,epoch)
    assert np.abs(residual).max()<=period/2+1e-10
    f=np.ones(len(t),dtype=float)
    box=np.abs(residual)<=.15/2
    f[box]=.99
    flag=np.zeros(len(t),dtype=np.int64)
    flag[np.flatnonzero(box)[0]]=1
    d=mod.count_event_support(t,f,flag,period,epoch,.15)
    assert d["n_bls_cycles_with_at_least_two_cadences"]>=3
    assert d["n_bls_events_with_positive_median_drop"]>=3
    assert d["n_TGLC_flagged_in_bls_windows"]==1
    assert abs(d["observed_bls_window_depth_median"]-.01)<1e-7
    assert 0<d["largest_event_fraction_in_transit_cadences"]<.6

    # A single in-sector event can appear very deep despite no periodic
    # support. It MUST remain descriptive, not a score or confirmation.
    p_single=18.0
    ep_single=1404.0
    r_single,k_single=mod.phase_and_cycle(t,p_single,ep_single)
    box_single=np.abs(r_single)<=.15/2
    f_single=np.ones(len(t))
    f_single[box_single]=.8
    d_single=mod.count_event_support(
        t,f_single,np.zeros(len(t),dtype=np.int64),
        p_single,ep_single,.15
    )
    assert d_single["n_bls_cycles_with_at_least_two_cadences"]<=2
    assert d_single["largest_event_fraction_in_transit_cadences"]>=.5
    assert d_single["observed_bls_window_depth_median"]>.15

    known=mod.catalogue_windows(
        t,period,2457000.0+epoch,3.6,"BJD-TDB"
    )
    assert known["n_frozen_catalogue_midpoints_in_sector"]>=3
    assert known["n_catalogue_windows_with_2plus_cadences"]>=3
    assert known["catalogue_window_time_standard"]=="explicit_BJD-TDB"

    ambiguous=mod.catalogue_windows(
        t,period,2457000.0+epoch,3.6,"BJD"
    )
    assert ambiguous["catalogue_window_time_standard"]==(
        "provisional_unresolved_BJD"
    )
    missing=mod.catalogue_windows(
        t,period,2457000.0+epoch,float("nan"),"BJD-TDB"
    )
    assert np.isnan(missing["n_catalogue_windows_with_2plus_cadences"])
    gaps=mod.gap_proximity(t,period,epoch)
    assert gaps["n_time_gaps_over_3h"]==0
    print("[PASS] Independent sampled transit-cycle support")
    print("[PASS] A deep single-event BLS peak does not imply a detection")
    print("[PASS] Recovered flux-drop and TGLC flagged cadence counts")
    print("[PASS] Catalogue windows keep unresolved BJD explicitly provisional")
    print("[PASS] No new BLS search or tuning performed")


if __name__=="__main__":
    main()
