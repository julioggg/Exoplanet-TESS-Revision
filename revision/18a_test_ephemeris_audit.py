#!/usr/bin/env python3
"""Offline tests of strict NASA BJD-TDB transit timing verification.

No NASA requests, no FITS downloads. Ensures source period and epoch are
FROZEN, BJD with unspecified dynamical time standard is not silently treated
as BJD-TDB, and long-period ephemeris extrapolation cannot be overstated.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path
import numpy as np
import pandas as pd

SRC=Path(__file__).with_name("18_ephemeris_reference_audit.py")
spec=importlib.util.spec_from_file_location("audit18",SRC)
mod=importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def main():
    p=3.0
    t=np.arange(1350.0, 1370.0, 1.0/48.0)
    curve=pd.DataFrame({"time_btjd":t})
    selected={
        "tic_id":"TIC_90000707","planet":"Synthetic confirmed-host b",
        "sector":3,"planet_period_days":p,
        "transit_epoch_jd":2457000.0+1351.2,
        "transit_epoch_time_standard":"BJD-TDB",
        "transit_epoch_error1":.001,
        "transit_epoch_error2":-.001,
        "planet_transit_duration_hours":3.0,
    }
    blind={
        "blind_direct_period_match_2pct":True,
        "best_period_days":p,
        "best_transit_time_btjd":1351.2,
    }
    live={
        "pl_orbper":p,
        "pl_orbpererr1":0.00001,
        "pl_orbpererr2":-0.00001,
        "ttv_flag":0,
    }
    v=mod.audit_one(selected,blind,live,curve)
    assert v["strict_timing_reference_and_uncertainty_usable"] is True
    assert v["blind_period_and_ephemeris_consistent"] is True
    assert v["n_predicted_transit_midpoints_in_sector"]>=5
    assert v["n_midpoints_with_cadence_within_half_duration"]>=5
    assert abs(v["blind_BLS_vs_ephemeris_phase_offset_hours"])<.001
    assert abs(mod.phase_offset(1.0,4.0,3.0))<1e-10

    unknown=dict(selected,transit_epoch_time_standard="BJD")
    v=mod.audit_one(unknown,blind,live,curve)
    assert v["strict_timing_reference_and_uncertainty_usable"] is False
    assert v["blind_period_and_ephemeris_consistent"] is None
    assert "time_system_not_explicit_BJD-TDB" in v[
        "timing_not_scoreable_reasons"
    ]

    ttv=dict(live,ttv_flag=1)
    v=mod.audit_one(selected,blind,ttv,curve)
    assert v["strict_timing_reference_and_uncertainty_usable"] is False
    assert "ttv_or_unknown_ttv_flag" in v["timing_not_scoreable_reasons"]

    period_changed=dict(live,pl_orbper=p+.01)
    v=mod.audit_one(selected,blind,period_changed,curve)
    assert v["strict_timing_reference_and_uncertainty_usable"] is False
    assert "current_NASA_period_differs" in v["timing_not_scoreable_reasons"]

    huge=dict(live,pl_orbpererr1=.2,pl_orbpererr2=-.2)
    v=mod.audit_one(selected,blind,huge,curve)
    assert v["strict_timing_reference_and_uncertainty_usable"] is False
    assert "ephemeris_uncertainty_exceeds" in v[
        "timing_not_scoreable_reasons"
    ]

    missed=dict(blind,best_transit_time_btjd=1352.0)
    v=mod.audit_one(selected,missed,live,curve)
    assert v["strict_timing_reference_and_uncertainty_usable"] is True
    assert v["blind_period_and_ephemeris_consistent"] is False
    print("[PASS] Strict BJD-TDB epoch scoring and frozen period")
    print("[PASS] Unknown time system and TTV prohibit timing claims")
    print("[PASS] Ephemeris changes and large uncertainty flagged")
    print("[PASS] Blind BLS phase mismatch is detected")
    print("All ephemeris provenance tests passed.")


if __name__=="__main__":
    main()
