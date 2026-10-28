#!/usr/bin/env python3
"""Experiment 25: DEVELOPMENT-ONLY surrogate-null calibration of BLS evidence.

36 TESS-confirmed hosts were ALREADY examined in experiments 21-24. Here
we ONLY investigate an exploratory *noise-surrogate* null, NOT verified
astrophysical negative stars, real false-positive rates, or detection FAPs.

For each of the 36 same original-TGLC curves generate 12 moving-block
circular-bootstrap surrogates (4 draws at each 24, 48, 96 original
30-minute cadences) within gaps >3 hours. Jointly resample normalized
flux and its TGLC flag; preserve the original BTJD sampling and gaps.
The operation disrupts some phase coherence but MAY retain actual
transit fragments, red noise and artifacts, so IS NOT a clean no-planet null.

Evaluate the frozen BLS 4000x5 scan on each surrogate and TWO prespecified
scores at the original blind highest peak:
  raw_log_power = log1p(BLS power)
  multievent_evidence = raw_log_power
      * sqrt(min(n_positive_local_events,4)/4)
      * sqrt(1 - fraction_nonzero_TGLC_in_box)
      * sqrt(clamp(remaining_1x1_BLS_power_ratio,0,1))

These are heuristic scores, NOT calibrated detection probabilities.
Use ONLY null surrogates of 24 deterministically selected sources for
a nominal 95th-percentile threshold for each score. Check exceedance
on 12 other ALREADY-EXPLORED development TICs. Contrast with the 108
previously obtained synthetic injection-recovery results on those same
TICs. Do not choose or retune score formulas using the observed results.

No LOCKED_HOLDOUT inputs or access; no new MAST request or source selection.
The final holdout remains unopened until genuinely separate astrophysical
controls, calibration protocol, and detector are frozen.

Run: python -u revision/25_development_surrogate_null_calibration.py
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.timeseries import BoxLeastSquares

ROOT=Path(__file__).resolve().parents[1]
SELECTED=ROOT/"results/future_tess_cohorts/development_products.csv"
PHOTOMETRY=ROOT/"results/new_tess_development_photometry/aperture_curve_metrics.csv"
PERIOD=ROOT/"results/development_period_event_evidence/per_host_period_and_evidence.csv"
INJECTED=ROOT/"results/development_physical_injection/per_injection_trial.csv"
PRIVATE=Path.home()/"exoplanet_revision/data/new_tess_development_cleaned"
OUT=ROOT/"results/development_surrogate_null"
POLICY="TGLC_APERTURE_1800S_V1"
SEED=20261028
N_CALIBRATION=24
N_INTERNAL_AUDIT=12
BLOCK_SIZES=[24,48,96]
DRAW_PER_SIZE=4
N_SURROGATES=12
NULL_QUANTILE=.95
PERIODS=np.linspace(.5,15.,4000)
DURATIONS=np.array([.05,.10,.15,.20,.25],dtype=float)
SCORES=["raw_log_power","multievent_evidence"]
MAX_EVENT_SCORE_COUNT=4


def file_sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):
            h.update(b)
    return h.hexdigest()


def import_helper():
    src=ROOT/"revision/23_development_event_evidence.py"
    spec=importlib.util.spec_from_file_location("experiment23_helper",src)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def assign_sources(tics):
    names=sorted(str(t) for t in tics)
    if len(names)!=36 or len(set(names))!=36:
        raise ValueError("Expected 36 distinct development TIC identifiers")
    # Hash + stable ordering avoids referring to catalogue or any scores.
    digest=lambda x: hashlib.sha256(f"EXP25|{SEED}|{x}".encode()).hexdigest()
    shuffled=sorted(names,key=lambda x:(digest(x),x))
    result={t:("calibration" if i<N_CALIBRATION
               else "internal_development_audit")
            for i,t in enumerate(shuffled)}
    if (sum(v=="calibration" for v in result.values())!=24 or
        sum(v=="internal_development_audit" for v in result.values())!=12):
        raise AssertionError("Source split changed")
    return result


def draw_seed(tic,block_size,draw):
    value=f"EXP25|{SEED}|{tic}|{block_size}|{draw}"
    return int.from_bytes(hashlib.sha256(value.encode()).digest()[:8],"little")


def block_surrogate(t,flux,flags,block_size,rng):
    """Moving-block circular bootstrap WITHIN original uninterrupted segments.

    Time grid and gap boundaries NEVER shift. Flux and flags are
    resampled jointly. Blocks are sampled with replacement: their
    marginal distributions need not match the original exactly.
    """
    t=np.asarray(t,float)
    f=np.asarray(flux,float)
    q=np.asarray(flags,np.int64)
    if not (len(t)==len(f)==len(q) and len(t)>=300):
        raise ValueError("Invalid frozen source lengths")
    if (not np.isfinite(t).all() or not np.isfinite(f).all()
        or not (np.diff(t)>0).all()):
        raise ValueError("Invalid frozen BTJD/flux")
    if block_size<2:
        raise ValueError("Unrealistically short time block")
    edges=np.concatenate(([0],np.flatnonzero(np.diff(t)>3/24)+1,[len(t)]))
    sample=np.empty(len(t),dtype=np.int64)
    for start,end in zip(edges[:-1],edges[1:]):
        n=int(end-start)
        if n<1:
            raise AssertionError("Empty observational segment")
        # Preserve singleton segment without inventing a different value.
        if n==1:
            sample[start]=start
            continue
        placed=0
        while placed<n:
            length=min(int(block_size),n-placed)
            origin=int(rng.integers(0,n))
            sample[start+placed:start+placed+length]=(
                start+(origin+np.arange(length))%n
            )
            placed+=length
    return f[sample],q[sample],sample


def score_from_features(power,n_positive,flag_frac,ratio):
    if not np.isfinite(power) or power<0:
        raise ValueError("Invalid BLS power")
    base=float(np.log1p(power))
    positive=0 if pd.isna(n_positive) else int(n_positive)
    if positive<0:
        raise ValueError("Negative positive-event count")
    flag=1. if pd.isna(flag_frac) else float(flag_frac)
    if not 0<=flag<=1:
        raise ValueError("Invalid original TGLC flag proportion")
    stability=0. if pd.isna(ratio) else float(ratio)
    if not np.isfinite(stability):
        stability=0.
    stability=float(np.clip(stability,0,1))
    weighted=base*math.sqrt(min(positive,MAX_EVENT_SCORE_COUNT)/
                            MAX_EVENT_SCORE_COUNT)
    weighted*=math.sqrt(1-flag)*math.sqrt(stability)
    return {
        "raw_log_power":base,
        "multievent_evidence":float(weighted),
        "positive_local_events":positive,
        "TGLC_nonzero_fraction":flag,
        "same_period_event_removal_power_ratio":float(ratio) if
            not pd.isna(ratio) else np.nan,
    }


def candidate_scores(t,flux,flags,noise,helper):
    model=BoxLeastSquares(t,flux,dy=np.full(len(t),noise,dtype=float))
    scan=model.power(PERIODS,DURATIONS)
    powers=np.asarray(scan.power,dtype=float)
    if not np.isfinite(powers).any() or len(powers)!=4000:
        raise ValueError("Nonfinite BLS surrogate power")
    ix=int(np.nanargmax(powers))
    period=float(scan.period[ix])
    epoch=float(scan.transit_time[ix])
    duration=float(scan.duration[ix])
    evidence,_=helper.local_event_windows(
        t,flux,flags,period,epoch,duration
    )
    keep,n_removed=helper.leave_dominant_cycle_out(
        t,flux,epoch,period,duration,
        evidence["dominant_positive_event_cycle"]
    )
    ratio=np.nan
    if n_removed:
        one=helper.same_period_bls(
            t,flux,period,duration,noise
        )["BLS_power_white_noise_surrogate"]
        after=helper.same_period_bls(
            t[keep],flux[keep],period,duration,noise
        )["BLS_power_white_noise_surrogate"]
        if one>0:
            ratio=float(after/one)
    values=score_from_features(
        float(powers[ix]),
        evidence["n_cycles_with_positive_local_depth"],
        evidence["fraction_inbox_TGLC_flagged"],
        ratio
    )
    return {
        "best_BLS_period_days":period,
        "best_BLS_power_white_noise_proxy":float(powers[ix]),
        "n_influential_event_cadences_removed":n_removed,
        **values
    }


def frozen_scores(period,injected):
    """Read already published scores; no injection/BLS rerun or retuning."""
    raw=[]
    for r in period.to_dict("records"):
        score=score_from_features(
            float(r["replayed_full_grid_BLS_power"]),
            r["n_cycles_with_positive_local_depth"],
            r["fraction_inbox_TGLC_flagged"],
            r["post_removal_to_unmasked_1x1_BLS_power_ratio"]
        )
        raw.append({
            "tic_id":str(r["tic_id"]),
            "trial_kind":"original_confirmed_host",
            "true_injected_period_recovered":np.nan,
            "published_original_period_correct":bool(
                r["frozen_original_direct_period_match_2pct"]
            ),
            **score
        })
    synthetic=[]
    for r in injected.to_dict("records"):
        score=score_from_features(
            float(r["best_BLS_power_white_noise_proxy"]),
            r["blind_BLS_positive_local_events"],
            r["blind_BLS_TGLC_flag_fraction"],
            r["blind_BLS_power_ratio_after_one_event"]
        )
        synthetic.append({
            "tic_id":str(r["tic_id"]),
            "pair_id":str(r["pair_id"]),
            "trial_kind":"injection_on_confirmed_host",
            "injection_bin":str(r["injection_bin"]),
            "radius_ratio":float(r["radius_ratio"]),
            "true_injected_period_recovered":bool(
                r["injected_period_direct_match_2pct"]
            ),
            "published_original_period_correct":np.nan,
            **score
        })
    return pd.DataFrame(raw),pd.DataFrame(synthetic)


def threshold_from_surrogates(nulls,groups,score):
    if score not in SCORES:
        raise ValueError("Unexpected score")
    cal=nulls[nulls.tic_id.map(groups)=="calibration"]
    audit=nulls[nulls.tic_id.map(groups)=="internal_development_audit"]
    if (
        len(cal)!=N_CALIBRATION*N_SURROGATES or
        len(audit)!=N_INTERNAL_AUDIT*N_SURROGATES
        or cal.tic_id.nunique()!=N_CALIBRATION or
        audit.tic_id.nunique()!=N_INTERNAL_AUDIT
    ):
        raise ValueError("Incomplete null calibration/audit split")
    threshold=float(np.quantile(cal[score],NULL_QUANTILE,method="higher"))
    if not np.isfinite(threshold):
        raise ValueError("Nonfinite null-score threshold")
    return threshold


def audit_surrogate_exceedance(nulls,groups,score,threshold):
    measured=nulls.copy()
    measured["development_split"]=measured.tic_id.map(groups)
    measured["exceeds_threshold"]=measured[score]>threshold
    frame=measured.groupby(["development_split","tic_id"]).agg(
        mean_exceedance=("exceeds_threshold","mean"),
        n_surrogates=("exceeds_threshold","size")
    ).reset_index()
    if not frame.n_surrogates.eq(N_SURROGATES).all():
        raise ValueError("Null trials per source changed")
    return measured,frame


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--private-dir",type=Path,default=PRIVATE)
    p.add_argument("--outdir",type=Path,default=OUT)
    args=p.parse_args()
    sel=pd.read_csv(SELECTED,dtype={"tic_id":str})
    qa=pd.read_csv(PHOTOMETRY,dtype={"tic_id":str})
    pre=pd.read_csv(PERIOD,dtype={"tic_id":str})
    inj=pd.read_csv(INJECTED,dtype={"tic_id":str})
    if any(len(d)!=36 for d in (sel,qa,pre)) or len(inj)!=108:
        raise ValueError("Experiment 25 requires complete frozen 36+108 data")
    if any(d.tic_id.duplicated().any() for d in (sel,qa,pre)):
        raise ValueError("Repeated TIC in original development manifests")
    if not all(set(d.tic_id)==set(sel.tic_id) for d in (qa,pre,inj)):
        raise ValueError("Development sources differ from selected 36")
    if not sel.role.eq("development").all():
        raise ValueError("Non-development source included")
    if sel.is_frozen_test_tic.astype(str).str.lower().eq("true").any():
        raise ValueError("LOCKED holdout contamination")
    if not qa.policy_id.eq(POLICY).all():
        raise ValueError("Frozen policy mismatch")
    if (inj.pair_id.duplicated().any() or
        inj.groupby("tic_id").size().ne(3).any()):
        raise ValueError("Prior physical injection schedule changed")
    groups=assign_sources(sel.tic_id)
    ref={r["tic_id"]:r for r in qa.to_dict("records")}
    pub={r["tic_id"]:r for r in pre.to_dict("records")}
    helper=import_helper()
    observed,synthetic=frozen_scores(pre,inj)
    observed["development_split"]=observed.tic_id.map(groups)
    synthetic["development_split"]=synthetic.tic_id.map(groups)
    nulls=[]
    errors=[]
    for i,row in enumerate(sel.sort_values("selection_rank").to_dict("records"),1):
        tic=row["tic_id"]
        path=args.private_dir/f"{tic}_sector_{int(row['sector']):02d}_{POLICY}.csv"
        try:
            if not path.is_file():
                raise FileNotFoundError(path)
            if path.resolve()!=Path(ref[tic]["private_cleaned_file"]).resolve():
                raise ValueError("Changed frozen source path")
            if file_sha(path)!=str(pub[tic]["source_curve_sha256"]):
                raise ValueError("Changed original policy-14 flux hash")
            d=pd.read_csv(path)
            if (len(d)!=int(ref[tic]["n_accepted"]) or
                not {"time_btjd","flux_cal_aper_normalized","tglc_flags"}.issubset(d)):
                raise ValueError("Wrong frozen source photometry")
            t=d.time_btjd.to_numpy(float)
            flux=d.flux_cal_aper_normalized.to_numpy(float)
            flag=d.tglc_flags.to_numpy(np.int64)
            noise=helper.robust_noise(flux)
            for block in BLOCK_SIZES:
                for draw in range(DRAW_PER_SIZE):
                    rng=np.random.default_rng(draw_seed(tic,block,draw))
                    bf,bq,indices=block_surrogate(t,flux,flag,block,rng)
                    stat=candidate_scores(t,bf,bq,noise,helper)
                    nulls.append({
                        "tic_id":tic,
                        "sector":int(row["sector"]),
                        "development_split":groups[tic],
                        "block_cadences":block,
                        "draw":draw,
                        "surrogate_seed":draw_seed(tic,block,draw),
                        "original_cleaned_curve_sha256":file_sha(path),
                        "no_original_period_label_used":True,
                        **stat
                    })
            print(
                f"[{i:02d}/36] {tic} surrogate scans complete "
                f"({N_SURROGATES} x {len(PERIODS)} BLS periods)",
                flush=True
            )
        except Exception as exc:
            errors.append({"tic_id":tic,"error":repr(exc)})
            print(f"[FAIL] {tic}: {exc}",flush=True)
    args.outdir.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(errors,columns=["tic_id","error"]).to_csv(
        args.outdir/"null_errors.csv",index=False
    )
    if errors or len(nulls)!=36*N_SURROGATES:
        raise RuntimeError(
            "Incomplete development-only null experiment: "+
            f"{len(nulls)}/{36*N_SURROGATES}. Do not use partial threshold."
        )
    nullframe=pd.DataFrame(nulls)
    nullframe.to_csv(args.outdir/"per_surrogate_null.csv",index=False)
    thresholds=[]
    measured=[]
    star_checks=[]
    for score in SCORES:
        thr=threshold_from_surrogates(nullframe,groups,score)
        n,by_star=audit_surrogate_exceedance(
            nullframe,groups,score,thr
        )
        n["score_name"]=score
        measured.append(n)
        by_star["score_name"]=score
        star_checks.append(by_star)
        for split in ("calibration","internal_development_audit"):
            subset=n[n.development_split==split]
            originals=observed[observed.development_split==split]
            injected=synthetic[synthetic.development_split==split]
            accepted=originals[score]>thr
            synth_accepted=injected[score]>thr
            recovered=injected.true_injected_period_recovered.astype(bool)
            thresholds.append({
                "score_name":score,
                "source_split":split,
                "threshold_95pct_calibration_surrogate":thr,
                "n_surrogate_null_trials":len(subset),
                "n_surrogate_exceedances":int(subset.exceeds_threshold.sum()),
                "surrogate_exceedance_fraction":float(
                    subset.exceeds_threshold.mean()),
                "n_source_hosts":len(originals),
                "n_source_host_blind_peaks_accepted":int(accepted.sum()),
                "n_source_host_correct_periods_accepted":int(
                    (accepted & originals.published_original_period_correct.astype(bool)).sum()
                ),
                "n_injection_trials":len(injected),
                "n_injections_accepted":int(synth_accepted.sum()),
                "n_injection_periods_correct_and_accepted":int(
                    (synth_accepted & recovered).sum()),
                "n_injection_periods_directly_recovered":int(
                    recovered.sum()
                ),
            })
    summary=pd.DataFrame(thresholds)
    summary.to_csv(args.outdir/"calibration_and_internal_audit.csv",index=False)
    pd.concat(star_checks,ignore_index=True).to_csv(
        args.outdir/"per_star_null_exceedance.csv",index=False)
    observed.to_csv(args.outdir/"prior_host_scores.csv",index=False)
    synthetic.to_csv(args.outdir/"prior_injected_scores.csv",index=False)
    details={
        "experiment":"25",
        "research_stage":"exploratory calibration of surrogate null on ALREADY INSPECTED development hosts",
        "n_tics":36,
        "split_calibration_tics":24,
        "split_internal_audit_tics":12,
        "no_independent_heldout_claim":True,
        "n_surrogates_per_tic":N_SURROGATES,
        "total_surrogates":36*N_SURROGATES,
        "block_cadences":BLOCK_SIZES,
        "draws_per_block":DRAW_PER_SIZE,
        "surrogate_method":"joint flux/TGLC moving-block circular bootstrap within >3h time gaps",
        "time_cadences_and_gap_times_preserved":True,
        "surrogates_are_NOT_verified_signal_free":True,
        "source_planetary_period_not_supplied_to_BLS":True,
        "unchanged_BLS_period_grid":[.5,15.,4000],
        "unchanged_BLS_duration_grid":DURATIONS.tolist(),
        "score_raw":"log1p(frozen-grid BLS max power)",
        "score_evidence":"raw * sqrt(min(positive_events,4)/4) * sqrt(1-flag_fraction) * sqrt(clip(one_event_remaining_power_ratio,0,1))",
        "score_selection_optimized_on_period_labels":False,
        "nominal_calibration_quantile":NULL_QUANTILE,
        "strict_greater_than_threshold_acceptance":True,
        "calibration_null_tic_set":sorted(t for t,r in groups.items() if r=="calibration"),
        "internal_audit_tic_set":sorted(t for t,r in groups.items() if r=="internal_development_audit"),
        "no_real_world_false_positive_rate":True,
        "no_detection_FAP_or_planet_significance_claim":True,
        "final_18_locked_holdout_accessed":False,
        "upstream_sha256":{
            "selection":file_sha(SELECTED),
            "photometry":file_sha(PHOTOMETRY),
            "real_event_evidence":file_sha(PERIOD),
            "synthetic_injection":file_sha(INJECTED)
        },
        "limitations":[
            "No real independently vetted non-transiting negative reference cohort.",
            "Moving-block resampling may retain transit fragments; does not produce proven planet-free nulls.",
            "This is a null-model exceedance analysis, NOT real-world FPR, specificity, or calibrated astrophysical false-alarm probability.",
            "The fixed score was chosen before this run, but the 36 development TICs were previously inspected.",
            "36 source TICs and 12 null replicates each: inference is sensitive to source dependence, bootstrap block length, nonstationarity and gap handling.",
            "The 12-source internal split is NOT a scientifically fresh independent test.",
            "Injection period concordance and acceptance are different metrics.",
            "No previous 18-host or final 18 locked holdout curves were accessed."
        ]
    }
    (args.outdir/"null_protocol.json").write_text(
        json.dumps(details,indent=2),encoding="utf-8"
    )
    report=[
        "# Experiment 25: surrogate-null BLS evidence calibration (DEVELOPMENT)",
        "",
        "Only previously explored 36-TIC development sources; 432 model-dependent",
        "moving-block surrogate nulls. The 18 new final heldout sources are untouched.",
        "THIS DOES NOT ESTIMATE AN ASTROPHYSICAL FALSE-POSITIVE RATE.",
        "",
        "## Fixed null-only thresholds and internal development audit", "",
        summary.round(5).to_string(index=False),
        "",
        "## Method restrictions", "",
        *["- "+s for s in details["limitations"]]
    ]
    (args.outdir/"surrogate_null_report.md").write_text(
        "\n".join(report)+"\n",encoding="utf-8"
    )
    print("\n"+"\n".join(report),flush=True)


if __name__=="__main__":
    main()
