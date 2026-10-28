#!/usr/bin/env python3
"""Experiment 24: blinded PHYSICAL transit injections into 36 development TICs.

Do not open the 18 final LOCKED_HOLDOUT sources.
Do not reuse the previously inspected 18 pilot-positive TICs.
No catalogue planet period/epoch is ever supplied to BLS.

Each ORIGINAL MAST TGLC 1800s aperture light curve is left unchanged.
Three deterministic injected versions are generated from this same baseline:
  radius ratio k = 0.03, 0.06, 0.12;
  injection period in each of 1.5-3.5, 3.5-6.5, 6.5-9.5 day bins;
  3x3 Latin-square schedule using prior frozen TIC selection rank.
This makes orbital-period stratum and injected radius ratio independently
balanced (12 trials per 3x3 cell), while using only 3 trials per TIC.
Synthetic PHYSICAL profiles use the unchanged PyTransit RoadRunner function
from experiment 06: quadratic limb darkening and finite-exposure sampling,
but an explicitly ASSUMED solar-like star, not observed stellar parameters.
Injection is multiplicative; no new clipping, normalization, or rebinning.

BLS is 4000 x 5 frozen period/duration search. A common white-noise
surrogate derived from each unmodified source is used for each of its three
synthetic pairs. After blind period selection, report:
 - recovery of INJECTED period within 2% (or 1:2 / 2:1 alias),
 - distinct sampled injected events (not an exclusion criterion),
 - local evidence and leave-one-event influence at recovered BLS period,
 - whether the unmodified source's BLS period matches injected period,
 - star-cluster bootstrap intervals for descriptive injected recovery.

Do NOT call these true-negative, real-world FPR, generalizable sensitivity,
or statistical planet-confirmation results. Baseline curves contain unknown
other variability and possibly real transits. This is DEVELOPMENT only.

This does NOT calibrate a detection threshold; separate null-surrogate and
verified astrophysical-control validation will follow.

Run:
 python -u revision/24_physical_development_injections.py
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path

os.environ.setdefault("NUMBA_NUM_THREADS","1")

import numpy as np
import pandas as pd
from astropy.timeseries import BoxLeastSquares
from scipy.stats import beta

ROOT=Path(__file__).resolve().parents[1]
SELECTION=ROOT/"results/future_tess_cohorts/development_products.csv"
ACQUISITION=ROOT/"results/new_tess_development_acquisition/verified_original_fits.csv"
QA=ROOT/"results/new_tess_development_photometry/aperture_curve_metrics.csv"
PRIOR=ROOT/"results/development_multievent_bls/per_host_development_comparison.csv"
PRIVATE=Path.home()/"exoplanet_revision/data/new_tess_development_cleaned"
OUT=ROOT/"results/development_physical_injection"
POLICY="TGLC_APERTURE_1800S_V1"
SEED=20261028
RADIUS_RATIOS=[0.03,0.06,0.12]
PERIOD_RANGES=[(1.5,3.5),(3.5,6.5),(6.5,9.5)]
PGRID=np.linspace(.5,15.0,4000)
DGRID=np.array([.05,.10,.15,.20,.25],dtype=float)
DIRECT_TOL=.02
CI_BOOTSTRAPS=2000


def sha256(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda:stream.read(1024*1024),b""):
            h.update(block)
    return h.hexdigest()


def rng_for(tic,level):
    digest=hashlib.sha256(
        f"EXOPLANET24|{SEED}|{tic}|{level}".encode()
    ).digest()
    return np.random.default_rng(
        int.from_bytes(digest[:8],"little")
    )


def load_helper(path,module_name):
    spec=importlib.util.spec_from_file_location(module_name,path)
    obj=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


def curves_verified(selected,acquired,qa,prior,private):
    tables=(selected,acquired,qa,prior)
    if any(len(d)!=36 for d in tables):
        raise ValueError("Exactly 36 frozen development TICs required")
    if any(d.tic_id.duplicated().any() for d in tables):
        raise ValueError("Repeated TIC in upstream development table")
    if not all(set(d.tic_id)==set(selected.tic_id) for d in tables):
        raise ValueError("Inconsistent development TIC membership")
    if (
        set(selected.selection_rank.astype(int))!=set(range(1,37))
        or selected.selection_rank.duplicated().any()
    ):
        raise ValueError("Frozen development selection_rank 1..36 is required")
    if not selected.role.astype(str).eq("development").all():
        raise ValueError("Attempt to access non-development TIC")
    if selected.is_frozen_test_tic.astype(str).str.lower().eq("true").any():
        raise ValueError("The final locked holdout must NOT be used")
    if not acquired.status.eq("verified_original_fits").all():
        raise ValueError("Non-verified original FITS in manifest")
    if not qa.policy_id.eq(POLICY).all():
        raise ValueError("Policy-14 photometry mismatch")
    lookup={x["tic_id"]:x for x in qa.to_dict("records")}
    for r in selected.itertuples():
        tic=r.tic_id
        expected=private/f"{tic}_sector_{int(r.sector):02d}_{POLICY}.csv"
        if not expected.is_file():
            raise FileNotFoundError(f"Missing frozen development flux: {expected}")
        if Path(str(lookup[tic]["private_cleaned_file"])).resolve()!=expected.resolve():
            raise ValueError(f"Frozen QA path drift for {tic}")
    for attr in ("sector","mast_uri","author","exptime_seconds"):
        merged=selected[["tic_id",attr]].merge(
            acquired[["tic_id",attr]],on="tic_id",validate="one_to_one",
            suffixes=("_selected","_actual")
        )
        x=merged[f"{attr}_selected"]
        y=merged[f"{attr}_actual"]
        if attr=="exptime_seconds":
            ok=np.allclose(x.astype(float),y.astype(float))
        else:
            ok=x.astype(str).equals(y.astype(str))
        if not ok:
            raise ValueError("Frozen source MAST identity mismatch: "+attr)


def read_one(selected,qa,prior,private):
    tic=str(selected["tic_id"])
    path=private/f"{tic}_sector_{int(selected['sector']):02d}_{POLICY}.csv"
    d=pd.read_csv(path)
    req={"time_btjd","flux_cal_aper_normalized","tglc_flags","cadence_num"}
    if not req.issubset(d):
        raise ValueError("Missing frozen photometric source columns")
    if (
        len(d)!=int(qa["n_accepted"])
        or len(d)!=int(prior["n_cadences"])
        or sha256(path)!=str(prior["cleaned_curve_sha256"])
    ):
        raise ValueError("Cleaned development flux provenance changed")
    t=d.time_btjd.to_numpy(float)
    f=d.flux_cal_aper_normalized.to_numpy(float)
    flags=d.tglc_flags.to_numpy(np.int64)
    if (
        len(t)<300 or not np.isfinite(t).all() or not np.isfinite(f).all()
        or (np.diff(t)<=0).any()
    ):
        raise ValueError("Invalid time or flux in original frozen curve")
    return t,f,flags,sha256(path)


def injection_config(tic,level,selection_rank):
    if level not in (0,1,2):
        raise ValueError("Expected one of 3 fixed synthetic planet levels")
    if not 1<=int(selection_rank)<=36:
        raise ValueError("Invalid frozen 36-TIC selection rank")
    rng=rng_for(tic,level)
    lo,hi=PERIOD_RANGES[level]
    radius_idx=(int(selection_rank)-1+level)%3
    period=float(rng.uniform(lo+.05,hi-.05))
    # Choose epoch safely within the observed span, with no knowledge of
    # catalogue ephemerides, source transit labels or BLS outputs.
    phase=float(rng.uniform(.0,1.0))
    impact=float(rng.uniform(.05,.65))
    return {
        "injection_bin":f"{lo:.1f}-{hi:.1f}",
        "injection_level":level,
        "true_period_days":period,
        "radius_ratio":float(RADIUS_RATIOS[radius_idx]),
        "radius_level":int(radius_idx),
        "impact_parameter":impact,
        "phase_0_to_1":phase
    }


def physical_model(t,config,transit_model):
    profile,meta=transit_model(
        t,
        config["true_period_days"],
        config["radius_ratio"],
        config["impact_parameter"],
        config["phase_0_to_1"]
    )
    profile=np.asarray(profile,dtype=float)
    if profile.shape!=t.shape or not np.isfinite(profile).all():
        raise ValueError("Invalid physical transit profile")
    if profile.min()>=1. or profile.min()<=0.:
        raise ValueError("Synthetic transit profile has no valid dips")
    return profile,meta


def count_sampled_injected_events(t,profile,period,t0):
    """Observed distinct injected cycles with >=2 deeply affected cadences.

    Sampled threshold is HALF of the maximum modeled synthetic depth.
    This describes observation coverage, not a detection label.
    """
    depths=1-np.asarray(profile,float)
    high=float(np.max(depths))
    if high<=0:
        return 0
    affected=np.flatnonzero(depths>=.5*high)
    cycles=np.rint((t[affected]-t0)/period).astype(np.int64)
    _,n=np.unique(cycles,return_counts=True)
    return int((n>=2).sum())


def scan_candidate(t,flux,noise):
    bls=BoxLeastSquares(t,flux,dy=np.full(len(t),noise))
    power=bls.power(PGRID,DGRID)
    vals=np.asarray(power.power,float)
    if len(vals)!=4000 or not np.isfinite(vals).any():
        raise ValueError("Invalid frozen injected BLS periodogram")
    i=int(np.nanargmax(vals))
    return {
        "best_period_days":float(power.period[i]),
        "best_duration_days":float(power.duration[i]),
        "best_epoch_btjd":float(power.transit_time[i]),
        "best_BLS_power_white_noise_proxy":float(vals[i]),
        "best_BLS_depth_snr_white_noise_proxy":float(power.depth_snr[i]),
    }


def compare_period(best,true):
    ratio=float(best)/float(true)
    return (
        abs(ratio-1)<=DIRECT_TOL,
        abs(ratio-.5)<=.01 or abs(ratio-2)<=.04,
    )


def influence_at_blind_period(t,flux,flags,candidate,noise,helper):
    p=candidate["best_period_days"]
    d=candidate["best_duration_days"]
    epoch=candidate["best_epoch_btjd"]
    summary,_=helper.local_event_windows(t,flux,flags,p,epoch,d)
    keep,excluded=helper.leave_dominant_cycle_out(
        t,flux,epoch,p,d,summary["dominant_positive_event_cycle"]
    )
    # Both powers from comparable FIXED-period/duration optimizations.
    before=helper.same_period_bls(t,flux,p,d,noise)
    after=(
        helper.same_period_bls(t[keep],flux[keep],p,d,noise)
        if excluded else None
    )
    ratio=(
        float(after["BLS_power_white_noise_surrogate"]/
              before["BLS_power_white_noise_surrogate"])
        if after is not None and before["BLS_power_white_noise_surrogate"]>0
        else np.nan
    )
    return summary,ratio,excluded


def star_bootstrap_ci(frame,col,seed=SEED,n_boot=CI_BOOTSTRAPS):
    """Resample source TIC, not multiple injections as independent targets."""
    if frame.tic_id.duplicated().any() or len(frame)<5:
        raise ValueError("Each bootstrap cell requires >=5 distinct source TICs")
    rng=np.random.default_rng(seed)
    x=frame[col].astype(int).to_numpy()
    n=len(x)
    counts=rng.integers(0,n,size=(n_boot,n))
    draws=x[counts].mean(axis=1)
    return np.quantile(draws,[.025,.975]).tolist()


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--private-dir",type=Path,default=PRIVATE)
    ap.add_argument("--outdir",type=Path,default=OUT)
    args=ap.parse_args()
    selected=pd.read_csv(SELECTION,dtype={"tic_id":str})
    acquired=pd.read_csv(ACQUISITION,dtype={"tic_id":str})
    qa=pd.read_csv(QA,dtype={"tic_id":str})
    prior=pd.read_csv(PRIOR,dtype={"tic_id":str})
    curves_verified(selected,acquired,qa,prior,args.private_dir)
    ref={r["tic_id"]:r for r in qa.to_dict("records")}
    previous={r["tic_id"]:r for r in prior.to_dict("records")}
    original=load_helper(ROOT/"revision/06_injection_recovery.py","exp06_transit")
    helper=load_helper(ROOT/"revision/23_development_event_evidence.py","exp23_event")
    completed=[]
    failed=[]
    for j,star in enumerate(
        selected.sort_values("selection_rank").to_dict("records"),1
    ):
        tic=str(star["tic_id"])
        try:
            t,flux,flags,source_sha=read_one(
                star,ref[tic],previous[tic],args.private_dir
            )
            noise=helper.robust_noise(flux)
            for level in range(3):
                cfg=injection_config(tic,level,star["selection_rank"])
                profile,meta=physical_model(t,cfg,original.transit_model)
                injected=flux*profile
                # Re-use the unmodified star's point-to-point noise estimate
                # for all three injected members; same sampling and flags.
                candidate=scan_candidate(t,injected,noise)
                evidence,ratio,nremoved=influence_at_blind_period(
                    t,injected,flags,candidate,noise,helper
                )
                direct,alias=compare_period(
                    candidate["best_period_days"],
                    cfg["true_period_days"]
                )
                base_direct,base_alias=compare_period(
                    previous[tic]["baseline_global_period_days"],
                    cfg["true_period_days"]
                )
                n_covered=count_sampled_injected_events(
                    t,profile,cfg["true_period_days"],meta["t0"]
                )
                if n_covered<1:
                    raise ValueError(
                        "Unexpected no supported synthetic transit cycle; "
                        "do not resample a more favorable epoch"
                    )
                completed.append({
                    "tic_id":tic,
                    "sector":int(star["sector"]),
                    "planet_host_label_not_used_in_BLS":True,
                    "source_sha256":source_sha,
                    "pair_id":f"{tic}_physical_{level}",
                    "n_accepted_cadences":len(t),
                    "original_noise_white_proxy":noise,
                    "original_source_bls_period_days":float(
                        previous[tic]["baseline_global_period_days"]),
                    "original_source_period_matches_injected":bool(base_direct),
                    "original_source_period_alias_matches_injected":bool(base_alias),
                    **cfg,
                    "physical_t0_btjd":float(meta["t0"]),
                    "physical_a_over_rstar_assumed":float(meta["a_over_rstar"]),
                    "profile_min_flux":float(profile.min()),
                    "n_injected_cycles_2plus_half_depth_cadences":n_covered,
                    **candidate,
                    "injected_period_direct_match_2pct":bool(direct),
                    "injected_period_half_or_double_alias":bool(alias),
                    "blind_BLS_positive_local_events":
                        evidence["n_cycles_with_positive_local_depth"],
                    "blind_BLS_positive_event_dominance":
                        evidence["largest_positive_event_weight_share"],
                    "blind_BLS_TGLC_flag_fraction":
                        evidence["fraction_inbox_TGLC_flagged"],
                    "blind_BLS_power_ratio_after_one_event":
                        ratio,
                    "n_removed_in_event_influence":nremoved,
                })
                print(
                    f"[{j:02d}/36] {tic} injection {level+1}/3 "
                    f"k={cfg['radius_ratio']:.3f} "
                    f"P_inj={cfg['true_period_days']:.4f} "
                    f"P_BLS={candidate['best_period_days']:.4f} "
                    f"direct={bool(direct)} cycles={n_covered}",
                    flush=True
                )
        except Exception as exc:
            failed.append({"tic_id":tic,"error":repr(exc)})
            print(f"[FAIL] {tic} {exc}",flush=True)

    args.outdir.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(failed,columns=["tic_id","error"]).to_csv(
        args.outdir/"physical_injection_errors.csv",index=False
    )
    if len(completed)!=108 or failed:
        raise RuntimeError(
            f"Only {len(completed)}/108 physical injection trials complete. "
            "No success-rate files written; inspect physical_injection_errors.csv"
        )
    frame=pd.DataFrame(completed)
    if frame.pair_id.duplicated().any():
        raise ValueError("Injected synthetic pair IDs not unique")
    frame.to_csv(args.outdir/"per_injection_trial.csv",index=False)
    agg=[]
    for (level,radius_level),sub in frame.groupby(
        ["injection_level","radius_level"],sort=True
    ):
        lo,hi=star_bootstrap_ci(
            sub,"injected_period_direct_match_2pct",
            SEED+int(level)*10+int(radius_level)
        )
        agg.append({
            "injection_level":int(level),
            "radius_level":int(radius_level),
            "radius_ratio":float(sub.radius_ratio.iloc[0]),
            "period_stratum":sub.injection_bin.iloc[0],
            "n_unique_source_TICs":sub.tic_id.nunique(),
            "n_trials":len(sub),
            "n_direct":int(sub.injected_period_direct_match_2pct.sum()),
            "direct_fraction":float(sub.injected_period_direct_match_2pct.mean()),
            "cluster_bootstrap95_low":lo,
            "cluster_bootstrap95_high":hi,
            "n_direct_or_alias":int(
                (sub.injected_period_direct_match_2pct |
                 sub.injected_period_half_or_double_alias).sum()
            ),
            "n_at_least_2_observed_injected_cycles":int(
                (sub.n_injected_cycles_2plus_half_depth_cadences>=2).sum()
            ),
            "n_at_least_3_observed_injected_cycles":int(
                (sub.n_injected_cycles_2plus_half_depth_cadences>=3).sum()
            ),
        })
    aggregate=pd.DataFrame(agg)
    if (
        len(aggregate)!=9 or len(frame)!=108
        or aggregate.n_trials.ne(12).any()
        or frame.tic_id.nunique()!=36
        or frame.groupby("tic_id").size().ne(3).any()
    ):
        raise ValueError("3x3 depth-period Latin-square design drift")
    aggregate.to_csv(args.outdir/"physical_injection_aggregate.csv",index=False)
    marginal=[]
    for dim in ("injection_bin","radius_ratio"):
        for value,sub in frame.groupby(dim,sort=True):
            ci=star_bootstrap_ci(
                sub,"injected_period_direct_match_2pct",
                SEED+100+len(marginal)
            )
            marginal.append({
                "factor":dim,
                "value":str(value),
                "n_stars":sub.tic_id.nunique(),
                "n_trials":len(sub),
                "direct_count":int(
                    sub.injected_period_direct_match_2pct.sum()
                ),
                "direct_fraction":float(
                    sub.injected_period_direct_match_2pct.mean()
                ),
                "tic_bootstrap95_low":float(ci[0]),
                "tic_bootstrap95_high":float(ci[1]),
            })
    pd.DataFrame(marginal).to_csv(
        args.outdir/"physical_injection_marginal.csv",index=False
    )
    config={
        "experiment":24,
        "role":"DEVELOPMENT physical injection-recovery, not final validation",
        "frozen_36_hosts_used":36,
        "locked_18_holdout_opened":False,
        "n_trials":108,
        "injections_per_TIC":3,
        "random_seed":SEED,
        "injection_radius_ratios":RADIUS_RATIOS,
        "balanced_3x3_latin_square_12_per_period_depth_cell":True,
        "depth_assignment":"((frozen selection_rank -1)+ period_level) % 3",
        "injection_period_ranges_days":PERIOD_RANGES,
        "injection_model":"RoadRunner PyTransit quadratic limb darkening",
        "stellar_assumptions":"Mstar=1Msun, Rstar=1Rsun (not observed)",
        "limb_darkening_coefficients":[0.3,0.2],
        "finite_exposure_supersampling":7,
        "preserve_original_cadences_and_TGLC_flags":True,
        "no_sigma_clipping_or_detrending":True,
        "no_new_normalization":True,
        "injections_multiply_original_cleaned_flux":True,
        "search_without_catalogue_or_injected_period":True,
        "BLS_period_grid":[.5,15.0,4000],
        "BLS_durations_days":DGRID.tolist(),
        "reference_period_compared_only_after_search":True,
        "same_noise_surrogate_for_baseline_and_injection":True,
        "no_calibrated_detection_threshold_or_real_false_positive_rate":True,
        "existing_planet_host_signals_and_instrumental_artifacts_remain":True,
        "CI_unit":"source TIC, bootstrap within fixed injection depth level",
        "source_selection_sha256":sha256(SELECTION),
        "source_QA_sha256":sha256(QA),
        "prior_development_comparison_sha256":sha256(PRIOR),
        "limitations":[
            "Synthetic solar-star transit population, not actual planet radii or stellar models.",
            "Period and depth are crossed in a balanced 3x3 design; each TIC supplies only 3/9 cells.",
            "These are the 36 previously inspected development stars.",
            "In-sector original planet signals may compete with injected signals.",
            "Original photometry and time sampling is unchanged; no true negatives.",
            "A recovered injected period is not equivalent to a calibrated planet detection.",
            "Empirical noise-only false-alarm rate is NOT reported.",
            "The 18 final locked holdout sources must remain unopened."
        ],
    }
    (args.outdir/"injection_protocol.json").write_text(
        json.dumps(config,indent=2),encoding="utf-8"
    )
    report=[
        "# Experiment 24: blind PHYSICAL synthetic injection recovery",
        "",
        "DEVELOPMENT ONLY, 36 original-TGLC sources, 108 injected trials.",
        "Do NOT interpret as final test accuracy, actual false-alarm rate, or",
        "population transit-detection completeness.",
        "",
        aggregate.to_string(index=False),
        "", "## Marginal summaries (36 source TICs per factor level)",
        pd.DataFrame(marginal).to_string(index=False),
        "", "## Integrity and interpretation", "",
        *["- "+q for q in config["limitations"]],
        "- All 18 LOCKED holdout stellar light curves remain unaccessed."
    ]
    (args.outdir/"physical_injection_report.md").write_text(
        "\n".join(report)+"\n",encoding="utf-8"
    )
    print("\n"+"\n".join(report),flush=True)


if __name__=="__main__":
    main()
