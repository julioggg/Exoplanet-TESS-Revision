#!/usr/bin/env python3
"""Experiment 23: separate blind PERIOD estimation and EVENT evidence.

EXPLORATORY DEVELOPMENT ANALYSIS on the EXACT 36 previously used
development TICs. This program NEVER reads or accesses LOCKED_HOLDOUT.

Principles:
 - Keep original BLS global best periods and frozen 6/36 period agreement.
 - Reconstruct original BLS epoch at the *already chosen* period/duration,
   using only accepted flux/time: no planet ephemeris enters any fit.
 - Measure independent transit-window support and local-baseline depth.
 - Remove the most influential transit window and RE-FIT only at the
   original period/duration, to quantify dependence on one event.
 - Produce diagnostic features, NOT a new optimized confidence threshold.
 - NASA truth is read only for reporting unchanged original period
   concordance. Neither novel period ranking nor planet confirmation here.

The original 18-host results remain separate and previously examined;
the final new 18-TIC holdout stays unopened.

Run:
  python -u revision/23_development_event_evidence.py

Astropy is imported lazily inside the fixed-period BLS function so the
event-statistics regression test can run even without astropy.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
COMPARISON = ROOT / "results/development_multievent_bls/per_host_development_comparison.csv"
PHOTOMETRY = ROOT / "results/new_tess_development_photometry/aperture_curve_metrics.csv"
SELECTED = ROOT / "results/future_tess_cohorts/development_products.csv"
OUT = ROOT / "results/development_period_event_evidence"
PRIVATE = Path.home() / "exoplanet_revision/data/new_tess_development_cleaned"
POLICY = "TGLC_APERTURE_1800S_V1"
# EXACT search grid from experiment 21; the phase discretization can depend
# on the full period/duration search. A single-period BLS fit is not an
# exact reconstruction of the periodogram peak.
ORIGINAL_PERIOD_GRID = np.linspace(0.5, 15.0, 4000)
ORIGINAL_DURATION_GRID = np.array([0.05, 0.10, 0.15, 0.20, 0.25])

MIN_TWO_CADENCE_EVENTS = 2
MIN_LOCAL_FLANK_POINTS = 4
LOCAL_FLANK_OUTER_MULTIPLIER = 5.0


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for piece in iter(lambda: f.read(1024*1024), b""):
            h.update(piece)
    return h.hexdigest()


def curve_path(tic, sector, private):
    return private / f"{tic}_sector_{int(sector):02d}_{POLICY}.csv"


def robust_noise(f):
    diffs = np.diff(np.asarray(f, dtype=float))
    mad = np.median(np.abs(diffs - np.median(diffs)))
    noise = float(1.4826 * mad / np.sqrt(2))
    if not np.isfinite(noise) or noise <= 0:
        raise ValueError("No usable frozen point-to-point noise proxy")
    return noise


def original_full_grid_bls(t, f, noise):
    """Replay EXACT original blind BLS grid to reconstruct the original epoch.

    Unlike a 1x1 BLS search, this uses the original five-duration phase
    resolution. It never reads NASA periods or the holdout and NEVER
    replaces any previously published period-recovery classification.
    """
    from astropy.timeseries import BoxLeastSquares
    model = BoxLeastSquares(
        t, f, dy=np.full(len(f), noise, dtype=float)
    )
    periodogram = model.power(
        ORIGINAL_PERIOD_GRID, ORIGINAL_DURATION_GRID
    )
    powers = np.asarray(periodogram.power, dtype=float)
    if len(powers) != len(ORIGINAL_PERIOD_GRID) or not np.isfinite(powers).any():
        raise ValueError("Original BLS full-grid replay returned invalid values")
    idx = int(np.nanargmax(powers))
    return {
        "period_days": float(periodogram.period[idx]),
        "duration_days": float(periodogram.duration[idx]),
        "transit_epoch_btjd": float(periodogram.transit_time[idx]),
        "BLS_power_white_noise_surrogate": float(powers[idx]),
        "original_best_index": idx,
    }


def same_period_bls(t, f, period, duration, noise):
    """A NEW 1x1 fixed-(P,d) phase refit for sensitivity analysis.

    This is NOT an exact replay of the global periodogram's phase/power.
    Only compare single-period fits to other single-period fits.
    """
    from astropy.timeseries import BoxLeastSquares
    if not (np.isfinite(period) and np.isfinite(duration) and
            0 < duration < period):
        raise ValueError("Invalid frozen BLS candidate")
    bls = BoxLeastSquares(
        t, f, dy=np.full(len(f), noise, dtype=float)
    )
    scan = bls.power(
        np.asarray([period], dtype=float),
        np.asarray([duration], dtype=float),
    )
    return {
        "period_days": float(scan.period[0]),
        "duration_days": float(scan.duration[0]),
        "transit_epoch_btjd": float(scan.transit_time[0]),
        "BLS_power_white_noise_surrogate": float(scan.power[0]),
        "depth_fraction": float(scan.depth[0]),
    }


def local_event_windows(t, flux, flags, period, epoch, duration):
    """Local median depth of each distinct observed BLS transit window.

    Label-blind, and DOES NOT change period/epoch/quality masks. A local
    out-of-transit flank avoids equating full-sector trends with box depth.
    Require at least four local off-transit cadences for a valid local depth.

    No number here has been calibrated to a false-alarm probability.
    """
    t=np.asarray(t,dtype=float)
    flux=np.asarray(flux,dtype=float)
    q=np.asarray(flags,dtype=np.int64)
    if len(t)!=len(flux) or len(t)!=len(q) or len(t)<300:
        raise ValueError("Invalid matched frozen photometry")
    if not (np.isfinite(t).all() and np.isfinite(flux).all()
            and (np.diff(t)>0).all()):
        raise ValueError("Nonfinite or unordered frozen time series")
    if not (np.isfinite(period) and np.isfinite(epoch)
            and np.isfinite(duration) and 0<duration<period):
        raise ValueError("Invalid period, epoch, or transit box duration")

    cycles=np.rint((t-epoch)/period).astype(np.int64)
    midpoints=epoch+cycles*period
    offset=t-midpoints
    in_box=np.abs(offset)<=duration/2
    all_cycle_ids, counts=np.unique(cycles[in_box],return_counts=True)
    local_flank_outer=min(LOCAL_FLANK_OUTER_MULTIPLIER*duration,.35*period)
    rows=[]
    for cycle,count in zip(all_cycle_ids,counts):
        c=int(cycle)
        center=float(epoch+c*period)
        inside=in_box & (cycles==c)
        # Outer flank must be inside half an orbital cycle; this prevents
        # contaminating the baseline with the adjacent expected transit.
        flank=(
            (np.abs(t-center)>=duration)
            & (np.abs(t-center)<=local_flank_outer)
            & (cycles==c)
        ) if local_flank_outer>duration else np.zeros(len(t),bool)
        n_flank=int(flank.sum())
        n_in=int(inside.sum())
        depth=float(np.median(flux[flank])-np.median(flux[inside])
                    ) if n_flank>=MIN_LOCAL_FLANK_POINTS else np.nan
        rows.append({
            "cycle":c,
            "estimated_midpoint_btjd":center,
            "n_in_box":n_in,
            "n_local_out_of_transit":n_flank,
            "local_median_depth_fraction":depth,
            "TGLC_flagged_in_box":int(np.count_nonzero(q[inside])),
            "fraction_TGLC_flagged_in_box":float(np.count_nonzero(q[inside])/n_in),
            "usable_local_depth":bool(n_in>=MIN_TWO_CADENCE_EVENTS
                                      and np.isfinite(depth)),
            "positive_local_depth":bool(n_in>=MIN_TWO_CADENCE_EVENTS
                                        and np.isfinite(depth) and depth>0),
        })
    events=pd.DataFrame(rows,columns=[
        "cycle","estimated_midpoint_btjd","n_in_box",
        "n_local_out_of_transit","local_median_depth_fraction",
        "TGLC_flagged_in_box","fraction_TGLC_flagged_in_box",
        "usable_local_depth","positive_local_depth",
    ])
    eligible=events[events.usable_local_depth] if len(events) else events
    positives=eligible[
        eligible.local_median_depth_fraction>0
    ] if len(eligible) else eligible
    signal_weight=(
        positives.local_median_depth_fraction.to_numpy(float) *
        positives.n_in_box.to_numpy(float)
    ) if len(positives) else np.array([])
    dominant_cycle=(
        int(positives.iloc[int(np.argmax(signal_weight))].cycle)
        if len(signal_weight) else None
    )
    peak_share=(
        float(signal_weight.max()/signal_weight.sum())
        if len(signal_weight) and signal_weight.sum()>0 else np.nan
    )
    return {
        "n_cycles_any_box_sample":int(len(events)),
        "n_cycles_with_2plus_points":int(
            (events.n_in_box>=2).sum()) if len(events) else 0,
        "n_cycles_with_usable_local_depth":int(len(eligible)),
        "n_cycles_with_positive_local_depth":int(len(positives)),
        "positive_local_depth_fraction_of_evaluable_events":(
            float(len(positives)/len(eligible)) if len(eligible) else np.nan
        ),
        "median_positive_local_depth":(
            float(positives.local_median_depth_fraction.median())
            if len(positives) else np.nan
        ),
        "largest_positive_event_weight_share":peak_share,
        "dominant_positive_event_cycle":dominant_cycle,
        "n_inbox_TGLC_flagged":int(events.TGLC_flagged_in_box.sum())
            if len(events) else 0,
        "n_inbox_total":int(events.n_in_box.sum()) if len(events) else 0,
        "fraction_inbox_TGLC_flagged":(
            float(events.TGLC_flagged_in_box.sum()/events.n_in_box.sum())
            if len(events) and events.n_in_box.sum() else np.nan
        ),
    },events


def leave_dominant_cycle_out(t,flux,epoch,period,duration,cycle):
    """Mask only cadences inside one original fitted transit event."""
    if cycle is None:
        return np.ones(len(t),dtype=bool),0
    cyc=np.rint((t-epoch)/period).astype(np.int64)
    in_box=(
        (cyc==int(cycle))
        & (np.abs(t-(epoch+cyc*period))<=duration/2)
    )
    keep=~in_box
    if keep.sum()<300:
        raise ValueError("Dominant transit removal leaves <300 cadences")
    return keep,int(in_box.sum())


def evaluate_one(source,qa,private):
    tic=str(source["tic_id"])
    sector=int(source["sector"])
    path=curve_path(tic,sector,private)
    if path.resolve()!=Path(str(qa["private_cleaned_file"])).resolve():
        raise ValueError(f"Frozen cleaned-file path mismatch for {tic}")
    if str(qa["policy_id"])!=POLICY:
        raise ValueError(f"Unexpected photometric policy for {tic}")
    if not path.is_file():
        raise FileNotFoundError(path)
    frame=pd.read_csv(path)
    need={"time_btjd","flux_cal_aper_normalized","tglc_flags","cadence_num"}
    if not need.issubset(frame):
        raise ValueError("Frozen source missing mandatory data columns")
    if len(frame)!=int(qa["n_accepted"]) or len(frame)!=int(source["n_cadences"]):
        raise ValueError("Development curve length differs from prior analyses")
    t=frame.time_btjd.to_numpy(float)
    f=frame.flux_cal_aper_normalized.to_numpy(float)
    flags=frame.tglc_flags.to_numpy(np.int64)
    period=float(source["baseline_global_period_days"])
    duration=float(source["baseline_global_duration_days"])
    noise=robust_noise(f)
    # Replay the FULL exact original BLS search to recover its phase
    # at the already-published maximum. This changes no period decision.
    original=original_full_grid_bls(t,f,noise)
    original_power=float(source["baseline_global_power"])
    if (
        not np.isclose(original["period_days"],period,rtol=0,atol=1e-8)
        or not np.isclose(original["duration_days"],duration,rtol=0,atol=1e-8)
        or not np.isclose(original["BLS_power_white_noise_surrogate"],
                           original_power,rtol=2e-6,atol=1e-5)
    ):
        raise ValueError(
            f"Full-grid BLS replay does not reproduce frozen original "
            f"period/duration/power for {tic}: replay={original}; "
            f"original_P={period}, original_D={duration}, "
            f"original_power={original_power}"
        )
    # A single-(P,d) BLS call has a different phase-grid construction;
    # use its power ONLY as the denominator of an identical-method
    # single-(P,d) fit after removing the strongest observed event.
    one_grid= same_period_bls(t,f,period,duration,noise)
    epoch=original["transit_epoch_btjd"]
    evidence,events=local_event_windows(
        t,f,flags,period,epoch,duration
    )
    cycle=evidence["dominant_positive_event_cycle"]
    keep,n_removed=leave_dominant_cycle_out(
        t,f,epoch,period,duration,cycle
    )
    counterfactual=(
        same_period_bls(t[keep],f[keep],period,duration,noise)
        if n_removed else None
    )
    cfpower=(
        counterfactual["BLS_power_white_noise_surrogate"]
        if counterfactual else np.nan
    )
    result={
        "tic_id":tic,"planet":str(source["planet"]),
        "sector":sector,
        "catalogue_period_days_reporting_only":float(
            source["reference_period_days"]
        ),
        "frozen_original_BLS_period_days":period,
        "frozen_original_BLS_duration_days":duration,
        "frozen_original_direct_period_match_2pct":bool(
            source["baseline_direct_match_2pct"]
        ),
        "frozen_multievent_direct_period_match_2pct":bool(
            source["multievent_direct_match_2pct"]
        ),
        "source_curve_sha256":sha256(path),
        "reconstructed_BLS_transit_epoch_btjd":epoch,
        "replayed_full_grid_BLS_power":original[
            "BLS_power_white_noise_surrogate"],
        "one_period_fit_BLS_epoch_btjd":one_grid["transit_epoch_btjd"],
        "one_period_fit_BLS_power":one_grid[
            "BLS_power_white_noise_surrogate"],
        "one_period_vs_full_grid_power_ratio":float(
            one_grid["BLS_power_white_noise_surrogate"]/original_power
        ) if original_power>0 else np.nan,
        "n_dominant_transit_box_points_removed":n_removed,
        "fixed_period_BLS_power_after_one_event_removed":cfpower,
        "post_removal_to_original_BLS_power_ratio":(
            # Historical field name retained for report compatibility,
            # but numerator/denominator are BOTH one-period re-fits.
            float(cfpower/one_grid["BLS_power_white_noise_surrogate"])
            if n_removed and one_grid["BLS_power_white_noise_surrogate"]>0
            else np.nan
        ),
        "no_research_time_period_replacement":True,
        **evidence,
    }
    events.insert(0,"tic_id",tic)
    events.insert(1,"sector",sector)
    return result,events


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--private-dir",type=Path,default=PRIVATE)
    ap.add_argument("--outdir",type=Path,default=OUT)
    args=ap.parse_args()
    source=pd.read_csv(COMPARISON,dtype={"tic_id":str})
    qa=pd.read_csv(PHOTOMETRY,dtype={"tic_id":str})
    selected=pd.read_csv(SELECTED,dtype={"tic_id":str})
    if any(len(x)!=36 for x in (source,qa,selected)):
        raise ValueError("Only the complete 36-source development set is allowed")
    if any(df.tic_id.duplicated().any() for df in (source,qa,selected)):
        raise ValueError("Duplicate TIC")
    if set(source.tic_id)!=set(qa.tic_id)!=set(selected.tic_id):
        raise ValueError("Development TIC sets disagree")
    if not selected.role.eq("development").all():
        raise ValueError("Locked holdout or other role among selected sources")
    if not (qa.policy_id==POLICY).all():
        raise ValueError("Non-frozen photometric policy")
    # We never read any heldout file or reference its file path.
    row_by_tic={x["tic_id"]:x for x in qa.to_dict("records")}
    output=[]
    all_events=[]
    failures=[]
    for r in source.to_dict("records"):
        try:
            row,events=evaluate_one(r,row_by_tic[str(r["tic_id"])],
                                    args.private_dir)
            output.append(row)
            all_events.append(events)
            ratio=row["post_removal_to_original_BLS_power_ratio"]
            print(
                f"[{len(output):02d}/36] {r['tic_id']}: "
                f"events local_positive={row['n_cycles_with_positive_local_depth']}, "
                f"single_event_power_ratio={ratio:.3f}"
                if np.isfinite(ratio) else
                f"[{len(output):02d}/36] {r['tic_id']}: "
                f"no dominant positive sampled local event",
                flush=True
            )
        except Exception as err:
            failures.append({"tic_id":str(r["tic_id"]),
                             "error":repr(err)})
            print(f"[FAIL] {r['tic_id']}: {err}",flush=True)
    args.outdir.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(failures,columns=["tic_id","error"]).to_csv(
        args.outdir/"evidence_errors.csv",index=False
    )
    if len(output)!=36:
        raise RuntimeError(
            f"{len(output)}/36 source scans validated. No partial evidence "
            "summary written; inspect evidence_errors.csv."
        )
    df=pd.DataFrame(output)
    events=pd.concat(all_events,ignore_index=True)
    df.to_csv(args.outdir/"per_host_period_and_evidence.csv",index=False)
    events.to_csv(args.outdir/"per_event_local_depths.csv",index=False)

    primary_success=int(df.frozen_original_direct_period_match_2pct.sum())
    multievent_success=int(df.frozen_multievent_direct_period_match_2pct.sum())
    if (primary_success,multievent_success)!=(6,3):
        raise ValueError("Previously published 6/36 and 3/36 changed")
    config={
        "experiment":23,
        "role":"development-only, exploratory evidence audit",
        "original_global_BLS_period_estimate_changed":False,
        "frozen_original_direct_period_concordance":"6/36",
        "frozen_multievent_direct_period_concordance":"3/36",
        "original_BLS_full_grid_replayed_for_epoch_provenance":True,
        "replay_grid_n_periods":len(ORIGINAL_PERIOD_GRID),
        "replay_duration_grid_days":ORIGINAL_DURATION_GRID.tolist(),
        "no_NASA_truth_used_in_replay":True,
        "BLS_evaluation_only_at_frozen_period_and_duration_for_event_influence":True,
        "single_period_fit_not_expected_to_equal_full_grid_peak":True,
        "leave_one_event_power_denominator":"single-period single-duration fit before removal",
        "leave_one_sampled_event_box_out_applied_at_same_period":True,
        "one_event_removal_does_not_claim_independent_validation":True,
        "no_new_detector_thresholds_set":True,
        "locked_holdout_accessed":False,
        "number_sources":36,
        "prior_BLS_comparison_sha256":sha256(COMPARISON),
        "prior_photometry_QA_sha256":sha256(PHOTOMETRY),
        "prior_source_selection_sha256":sha256(SELECTED),
        "original_flux_quality_policy":POLICY,
        "event_local_baseline_flank_min_points":MIN_LOCAL_FLANK_POINTS,
        "local_outer_flank_max_duration_multiplier":
            LOCAL_FLANK_OUTER_MULTIPLIER,
        "local_outer_flank_max_period_fraction":0.35,
        "limitations":[
            "Event-level local median depth is descriptive; not a detection probability.",
            "Refitting BLS at the original period after deleting one event may shift phase.",
            "Positive host status does not guarantee in-sector planet transit.",
            "The full original 4000x5 grid is replayed only to reconstruct the exact original BLS epoch.",
            "1x1 (P,d) BLS refits use potentially different phase binning; their raw power is not expected to reproduce the full grid peak.",
            "The one-event power ratio compares identical 1x1 methods before and after removal; phase may reoptimize.",
            "No signal-free or contaminant-verified negative control cohort is included.",
            "BLS power is not a validated false-alarm probability.",
            "Hyperparameter selection and calibrated null tests still needed.",
            "Final locked holdout untouched; no independent performance claim."
        ],
    }
    (args.outdir/"evidence_protocol.json").write_text(
        json.dumps(config,indent=2),encoding="utf-8"
    )
    report=[
        "# Experiment 23 — Period estimation vs independent-event evidence",
        "",
        "DEVELOPMENT COHORT ONLY: 36 sources, previous 18-host and final locked",
        "18-host cohorts not accessed. No period reranking or signal thresholds.",
        "",
        "## Unchanged frozen development results",
        f"Original blind BLS period matches: {primary_success}/36.",
        f"Previous hard multi-event filter period matches: {multievent_success}/36.",
        "",
        "## Descriptive event-evidence audit",
        "Original 4000x5 BLS grid replayed without NASA truth solely to recover the published peak epoch.",
        "1x1 fixed-period fits are used for the leave-one-event power ratio; they may have different phase binning.",
        f"Sources with 3 or more independent locally positive transit windows: "
        f"{int((df.n_cycles_with_positive_local_depth>=3).sum())}/36.",
        f"Sources with 2 or fewer locally positive transit windows: "
        f"{int((df.n_cycles_with_positive_local_depth<=2).sum())}/36.",
        f"Sources with at least one locally positive event: "
        f"{int((df.n_cycles_with_positive_local_depth>0).sum())}/36.",
        f"Sources where an event removal experiment was feasible: "
        f"{int(df.n_dominant_transit_box_points_removed.gt(0).sum())}/36.",
        "",
        "## Individual source diagnostics",
        df[[
            "tic_id","planet","frozen_original_BLS_period_days",
            "catalogue_period_days_reporting_only",
            "frozen_original_direct_period_match_2pct",
            "n_cycles_with_usable_local_depth",
            "n_cycles_with_positive_local_depth",
            "largest_positive_event_weight_share",
            "fraction_inbox_TGLC_flagged",
            "n_dominant_transit_box_points_removed",
            "post_removal_to_original_BLS_power_ratio"
        ]].round(5).to_string(index=False),
        "", "## Interpretation",
        *["- "+s for s in config["limitations"]],
        "- Do not tune or test on the 18 final locked-holdout sources."
    ]
    (args.outdir/"evidence_report.md").write_text(
        "\n".join(report)+"\n",encoding="utf-8"
    )
    print("\n"+"\n".join(report[:18]),flush=True)


if __name__=="__main__":
    main()
