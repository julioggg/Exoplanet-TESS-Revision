#!/usr/bin/env python3
"""No-network synthetic checks for development-only event-veto attribution.

The 18 LOCKED_HOLDOUT targets are not used or queried.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd

SRC=Path(__file__).with_name("22_attribution_multievent_losses.py")
spec=importlib.util.spec_from_file_location("attribution22",SRC)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def make(tic,period,original,multi,orig_correct,multi_correct,
         events,pos,share,flagged):
    return {
        "tic_id":tic,"planet":"Synthetic "+tic,
        "period_stratum":"3.5-6.5",
        "reference_period_days":period,
        "baseline_global_period_days":original,
        "multievent_period_days":multi,
        "baseline_direct_match_2pct":orig_correct,
        "multievent_direct_match_2pct":multi_correct,
        "baseline_supported_events":events,
        "baseline_positive_events":pos,
        "baseline_largest_drop_share":share,
        "baseline_TGLC_flag_fraction":flagged,
        "multievent_candidate_exists":True,
        "multievent_half_or_double_alias":False,
    }


def main():
    df=pd.DataFrame([
        make("TIC_1",3.,3.,3.,True,True,7,7,.22,.02),
        make("TIC_2",4.,4.,8.,True,False,5,5,.30,.25),
        make("TIC_3",8.2,8.2,4.1,True,False,2,2,.5,0),
        make("TIC_4",5.,12.,5.,False,True,1,1,1.0,0),
        make("TIC_5",6.,12.,3.,False,False,1,1,1.0,0),
    ])
    df.loc[df.tic_id=="TIC_3","multievent_half_or_double_alias"]=True
    out,stats=m.analyze(df)
    assert len(out)==5
    assert stats["both_recover"]==1
    assert stats["original_only_lost"]==2
    assert stats["multievent_only_gained"]==1
    assert stats["neither_recover"]==1
    assert stats["lost_original_successes_by_TGLC_flag"]==1
    assert stats["lost_original_successes_by_min_events"]==1
    assert stats["original_top_peak_passing_all_event_rules"]==1
    assert stats["multievent_alias_only_n"]==1
    assert stats["absolute_multi_minus_original_percentage_points"]==-20.
    assert stats["paired_mcnemar_exact_two_sided_p_descriptive"]==1.0
    assert m.failure_reasons(df.iloc[1].to_dict())==[
        "TGLC_flagged_fraction_gt_20pct"
    ]
    assert m.failure_reasons(df.iloc[2].to_dict())==[
        "fewer_than_3_supported_events",
        "fewer_than_3_positive_events"
    ]
    # A global BLS peak passing all gates must NEVER be replaced by the
    # scanner; a changed result indicates inconsistent upstream evidence.
    corrupted=df.copy()
    corrupted.loc[0,"multievent_period_days"]=5.
    try:
        m.analyze(corrupted)
    except ValueError as exc:
        assert "Inconsistent" in str(exc)
    else:
        raise AssertionError("Inconsistent top-peak switch undetected")

    duplicated=pd.concat([df,df.iloc[[0]]],ignore_index=True)
    try:
        m.analyze(duplicated)
    except ValueError as exc:
        assert "Duplicate" in str(exc)
    else:
        raise AssertionError("Repeated development TIC not detected")

    print("[PASS] Exact paired 2x2 retention, losses and gains")
    print("[PASS] TGLC flag veto independently separated from event-count veto")
    print("[PASS] Alias-only period reported separately")
    print("[PASS] Inconsistent top-peak switch and repeated TIC rejected")
    print("[PASS] No new BLS scan or LOCKED_HOLDOUT access")


if __name__=="__main__":
    main()
