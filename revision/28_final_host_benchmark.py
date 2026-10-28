#!/usr/bin/env python3
"""Experiment 28 — ONE-SHOT locked 18-host TESS evaluation after experiment-27 freeze.

Never execute before the experiment-27 frozen protocol is committed to Git.
Run only AFTER acquiring exactly the frozen 18 TGLC 1800s original FITS
and applying unchanged TGLC_APERTURE_1800S_V1 photometry.

Reads science parameters and frozen thresholds solely from experiment-27
protocol, validates source SHA256 and all code hashes BEFORE reading flux.
The original 4000x5 BLS period search and event-evidence score are reused
from experiments 23/25 UNCHANGED. No parameter fitting, hyperparameter
selection, threshold tuning, negative-star claims, or NASA transit epoch
supplied to the blind detector.

Once all 18 input light curves are verified, results are processed in a
deterministic single pass. The run marker prevents reopening the test for
post hoc retuning. If a genuine code failure occurs, report transparently;
do not bypass the marker without documenting and re-registering a repair.

Host-positive period agreement is NOT transit-detection sensitivity; a
confirmed host does not guarantee observable in-sector transits. This
experiment cannot estimate real false-positive rates, specificity or AUC.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import beta

ROOT=Path(__file__).resolve().parents[1]
FROZEN=ROOT/"results/final_holdout_freeze/frozen_protocol.json"
META=ROOT/"results/future_tess_cohorts/LOCKED_DO_NOT_ANALYZE_holdout_products.csv"
ACQUISITION=ROOT/"results/final_holdout_acquisition/verified_original_fits.csv"
QA=ROOT/"results/final_holdout_photometry/aperture_curve_metrics.csv"
PRIVATE=Path.home()/"exoplanet_revision/data/final_holdout_cleaned"
OUT=ROOT/"results/final_holdout_benchmark"
POLICY="TGLC_APERTURE_1800S_V1"
DIRECT_TOL=.02
SCRIPT_FILES=[
    "revision/28_final_host_benchmark.py",
    "revision/23_development_event_evidence.py",
    "revision/25_development_surrogate_null_calibration.py",
    "revision/16_acquire_frozen_new_hosts.py",
    "revision/14_build_uniform_aperture_curves.py",
    "revision/26_paired_development_score_audit.py",
]
DATA_FILES=[
    "results/future_tess_cohorts/LOCKED_DO_NOT_ANALYZE_holdout_products.csv",
    "results/future_tess_cohorts/development_products.csv",
    "results/development_surrogate_null/calibration_and_internal_audit.csv",
    "results/development_surrogate_null/per_surrogate_null.csv",
    "results/development_surrogate_null/null_protocol.json",
    "results/development_paired_score_audit/paired_audit_report.md",
]


def sha256(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for piece in iter(lambda:f.read(1024*1024),b""):
            h.update(piece)
    return h.hexdigest()


def import_module(name,filename):
    spec=importlib.util.spec_from_file_location(name,ROOT/filename)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def frozen_inputs(frozen,metadata,acquired,qa):
    if frozen.get("experiment")!=27 or not frozen.get("ready_for_first_final_test"):
        raise ValueError("A complete experiment-27 FREEZE is required")
    if set(frozen["frozen_sources_hashes"])!=set(SCRIPT_FILES):
        raise ValueError("Unexpected evaluation/photometry code snapshot")
    if set(frozen["frozen_data_hashes"])!=set(DATA_FILES):
        raise ValueError("Unexpected upstream data snapshot")
    for file,h in {**frozen["frozen_sources_hashes"],
                   **frozen["frozen_data_hashes"]}.items():
        if sha256(ROOT/file)!=h:
            raise ValueError(f"Frozen code/data SHA256 changed: {file}")
    if frozen.get("locked_holdout_metadata_sha256")!=sha256(META):
        raise ValueError("Locked holdout selection metadata changed")
    if frozen.get("n_final_host_tics")!=18:
        raise ValueError("Wrong final cohort size")
    if frozen.get("frozen_photometry_policy")!=POLICY:
        raise ValueError("Changed photometry policy")
    if frozen.get("period_search_grid")!=[.5,15.,4000]:
        raise ValueError("Changed period search")
    if frozen.get("duration_grid_days")!=[.05,.10,.15,.20,.25]:
        raise ValueError("Changed BLS durations")
    if frozen.get("main_thresholds") is None or (
        set(frozen["main_thresholds"])!=
        {"raw_log_power","multievent_evidence"}
    ):
        raise ValueError("Frozen development score thresholds not present")
    for df in (metadata,acquired,qa):
        if len(df)!=18 or df.tic_id.duplicated().any():
            raise ValueError("Exactly 18 distinct locked TICs required")
    if not set(metadata.tic_id)==set(acquired.tic_id)==set(qa.tic_id):
        raise ValueError("Holdout acquisition/photometry membership mismatch")
    if not metadata.role.eq("LOCKED_HOLDOUT").all():
        raise ValueError("Non-holdout source entered final evaluation")
    if not metadata.is_frozen_test_tic.astype(str).str.lower().eq("true").all():
        raise ValueError("Original final-sample lock labels changed")
    if not metadata.selection_rank.astype(int).sort_values().to_list()==list(range(1,19)):
        raise ValueError("Frozen source rank not 1..18")
    if not acquired.status.eq("verified_original_fits").all():
        raise ValueError("Non-verified MAST original FITS")
    if not qa.policy_id.eq(POLICY).all():
        raise ValueError("Changed frozen 1800s original aperture policy")
    for col in ("sector","mast_uri","author","exptime_seconds"):
        m=metadata[["tic_id",col]].merge(
            acquired[["tic_id",col]],on="tic_id",validate="one_to_one",
            suffixes=("_chosen","_actual")
        )
        a=m[f"{col}_chosen"]
        b=m[f"{col}_actual"]
        okay=(np.allclose(a.astype(float),b.astype(float))
              if col=="exptime_seconds" else
              a.astype(str).equals(b.astype(str)))
        if not okay:
            raise ValueError("Frozen product identity mismatch: "+col)
    q=qa[["tic_id","sector","fits_sha256"]].merge(
        acquired[["tic_id","sector","sha256"]],on=["tic_id","sector"],
        validate="one_to_one",suffixes=("_qa","_acquired")
    )
    if len(q)!=18 or not q.fits_sha256.astype(str).equals(q.sha256.astype(str)):
        raise ValueError("FITS content hash mismatch between pipeline stages")
    if set(metadata.tic_id)!=set(frozen["locked_holdout_tics"]):
        raise ValueError("Frozen 18 TIC identities changed")


def exact_interval(k,n,alpha=.05):
    return (
        float(beta.ppf(alpha/2,k,n-k+1)) if k else 0.,
        float(beta.ppf(1-alpha/2,k+1,n-k)) if k<n else 1.
    )


def match_period(blind,truth):
    ratio=float(blind)/float(truth)
    return (
        bool(abs(ratio-1)<=DIRECT_TOL),
        bool(abs(ratio-.5)<=.01 or abs(ratio-2)<=.04),
    )


def blind_one(row,metrics_row,acquired_row,private,helper,score_module):
    """No catalogue period/ephemeris/labels accepted among parameters."""
    tic=str(row["tic_id"])
    p=private/f"{tic}_sector_{int(row['sector']):02d}_{POLICY}.csv"
    if p.resolve()!=Path(str(metrics_row["private_cleaned_file"])).resolve():
        raise ValueError("Frozen photometry private file path mismatched")
    if not p.is_file():
        raise FileNotFoundError(p)
    if str(metrics_row["policy_id"])!=POLICY:
        raise ValueError("Wrong original photometry policy")
    frame=pd.read_csv(p)
    cols={"time_btjd","flux_cal_aper_normalized","tglc_flags",
          "cadence_num"}
    if not cols.issubset(frame):
        raise ValueError("Missing original TGLC aperture flux columns")
    if len(frame)!=int(metrics_row["n_accepted"]) or len(frame)<300:
        raise ValueError("Accepted cadence count is invalid")
    if frame.cadence_num.duplicated().any():
        raise ValueError("Repeated original TESS cadence")
    t=frame.time_btjd.to_numpy(float)
    flux=frame.flux_cal_aper_normalized.to_numpy(float)
    flags=frame.tglc_flags.to_numpy(np.int64)
    if not np.isfinite(t).all() or not np.isfinite(flux).all() or (
        np.diff(t)<=0).any():
        raise ValueError("Nonfinite or unordered source flux/time")
    noise=helper.robust_noise(flux)
    features=score_module.candidate_scores(t,flux,flags,noise,helper)
    return {
        "tic_id":tic,
        "sector":int(row["sector"]),
        "n_accepted_cadences":len(t),
        "original_fits_sha256":str(acquired_row["sha256"]),
        "cleaned_curve_sha256":sha256(p),
        "point_to_point_noise_surrogate":float(noise),
        **features
    }


def summarise(blind,truth,thresholds):
    df=blind.merge(
        truth[["tic_id","planet","planet_period_days","period_bin"]],
        on="tic_id",validate="one_to_one"
    )
    if len(df)!=18:
        raise ValueError("Final held-out blind predictions incomplete")
    matches=[match_period(r.best_BLS_period_days,r.planet_period_days)
             for r in df.itertuples(index=False)]
    df["direct_period_match_2pct"]=[m[0] for m in matches]
    df["half_or_double_period_alias"]=[m[1] for m in matches]
    df["raw_threshold_accepted"]=(
        df.raw_log_power>float(thresholds["raw_log_power"])
    )
    df["evidence_threshold_accepted"]=(
        df.multievent_evidence>
        float(thresholds["multievent_evidence"])
    )
    df["direct_correct_and_raw_accepted"]=(
        df.direct_period_match_2pct & df.raw_threshold_accepted
    )
    df["direct_correct_and_evidence_accepted"]=(
        df.direct_period_match_2pct & df.evidence_threshold_accepted
    )
    d=int(df.direct_period_match_2pct.sum())
    alias=int(df.half_or_double_period_alias.sum())
    d_raw=int(df.direct_correct_and_raw_accepted.sum())
    d_evid=int(df.direct_correct_and_evidence_accepted.sum())
    bybin=df.groupby("period_bin",sort=True).agg(
        n_hosts=("tic_id","size"),
        direct=("direct_period_match_2pct","sum"),
        alias=("half_or_double_period_alias","sum"),
        raw_accepted=("raw_threshold_accepted","sum"),
        evidence_accepted=("evidence_threshold_accepted","sum"),
        correct_and_raw=("direct_correct_and_raw_accepted","sum"),
        correct_and_evidence=("direct_correct_and_evidence_accepted","sum"),
    ).reset_index()
    summary={
        "n_final_host_tics":18,
        "n_direct_blind_period_agreements":d,
        "direct_blind_period_concordance":d/18,
        "direct_binomial95_ci_DESCRIPTIVE_ONLY":exact_interval(d,18),
        "n_period_half_or_double_aliases":alias,
        "n_correct_and_raw_threshold_accepted":d_raw,
        "n_correct_and_evidence_threshold_accepted":d_evid,
        "n_wrong_period_but_raw_threshold_accepted":int(
            ((~df.direct_period_match_2pct)&df.raw_threshold_accepted).sum()
        ),
        "n_wrong_period_but_evidence_threshold_accepted":int(
            ((~df.direct_period_match_2pct)&df.evidence_threshold_accepted).sum()
        ),
        "n_raw_threshold_accepted":int(df.raw_threshold_accepted.sum()),
        "n_evidence_threshold_accepted":int(df.evidence_threshold_accepted.sum()),
        "stratified_results":bybin.to_dict("records"),
        "NO_estimation_of_detection_sensitivity":True,
        "NO_negative_host_control_or_real_FPR":True,
        "NO_AUC_specificity_or_planet_confirmation":True,
        "no_NASA_period_or_ephemeris_passed_to_BLS":True
    }
    return df,bybin,summary


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--private-dir",type=Path,default=PRIVATE)
    ap.add_argument("--outdir",type=Path,default=OUT)
    args=ap.parse_args()
    if not FROZEN.is_file():
        raise FileNotFoundError(
            "Run and GIT-COMMIT experiment 27 freeze BEFORE final test"
        )
    freeze=json.loads(FROZEN.read_text(encoding="utf-8"))
    meta=pd.read_csv(META,dtype={"tic_id":str})
    acquired=pd.read_csv(ACQUISITION,dtype={"tic_id":str})
    qa=pd.read_csv(QA,dtype={"tic_id":str})
    frozen_inputs(freeze,meta,acquired,qa)

    if args.outdir.exists() and (
        (args.outdir/"FINAL_TEST_STARTED.json").exists() or
        (args.outdir/"final_18_host_report.md").exists()
    ):
        raise RuntimeError(
            "Final 18-source benchmark already started; do not rerun "
            "or retune on held-out flux. Document any repair explicitly."
        )

    a={r["tic_id"]:r for r in acquired.to_dict("records")}
    q={r["tic_id"]:r for r in qa.to_dict("records")}
    items=meta.sort_values("selection_rank").to_dict("records")
    blind_items=[{
        "tic_id":r["tic_id"],"sector":int(r["sector"])
    } for r in items]
    for row in blind_items:
        p=args.private_dir/f"{row['tic_id']}_sector_{int(row['sector']):02d}_{POLICY}.csv"
        if not p.is_file():
            raise FileNotFoundError(
                f"Missing final holdout cleaned curve {p}; "
                "acquire and process complete frozen 18 before first evaluation"
            )
        if p.resolve()!=Path(q[row["tic_id"]]["private_cleaned_file"]).resolve():
            raise ValueError("Final cleaned photometry QA provenance mismatch")

    helper=import_module(
        "frozen_exp23","revision/23_development_event_evidence.py"
    )
    score_module=import_module(
        "frozen_exp25","revision/25_development_surrogate_null_calibration.py"
    )
    # Lock BEFORE reading final flux. No optimization permitted thereafter.
    args.outdir.mkdir(exist_ok=True,parents=True)
    (args.outdir/"FINAL_TEST_STARTED.json").write_text(
        json.dumps({
            "stage":"one-shot locked host benchmark",
            "freeze_sha256":sha256(FROZEN),
            "n_sources":18,
            "reruns_for_hyperparameter_selection_PROHIBITED":True,
            "results_not_yet_known_at_lock":True
        },indent=2),encoding="utf-8"
    )
    blind=[]
    failures=[]
    for i,row in enumerate(blind_items,1):
        tic=row["tic_id"]
        try:
            result=blind_one(
                row,q[tic],a[tic],args.private_dir,helper,score_module
            )
            blind.append(result)
            print(
                f"[{i:02d}/18] {tic} blind period="
                f"{result['best_BLS_period_days']:.5f}d; "
                f"evidence={result['multievent_evidence']:.3f}",
                flush=True
            )
        except Exception as exc:
            failures.append({"tic_id":tic,"error":repr(exc)})
            print(f"[FAIL] {tic}: {exc}",flush=True)
    pd.DataFrame(failures,columns=["tic_id","error"]).to_csv(
        args.outdir/"final_errors.csv",index=False
    )
    if len(blind)!=18 or failures:
        raise RuntimeError(
            f"Final one-shot evaluation incomplete: {len(blind)}/18. "
            "Do NOT change the protocol. Save logs and document a repair."
        )
    pred=pd.DataFrame(blind)
    pred.to_csv(args.outdir/"final_blind_predictions.csv",index=False)
    scored,bybin,summary=summarise(
        pred,meta,freeze["main_thresholds"]
    )
    scored.to_csv(args.outdir/"final_per_host_scored.csv",index=False)
    bybin.to_csv(args.outdir/"final_by_period_stratum.csv",index=False)
    outcome={
        "experiment":28,
        "stage":"first and only previously untouched 18-host assessment",
        "frozen_protocol_sha256":sha256(FROZEN),
        "fixed_method_sha256":sha256(
            ROOT/"revision/28_final_host_benchmark.py"
        ),
        "n_sources":18,
        "scores_frozen":freeze["main_thresholds"],
        "results":summary,
        "interpretation_limits":[
            "Host-positive period agreement is not detection sensitivity.",
            "Selected archived sector may lack observable planet transits.",
            "Planet orbital period accuracy does not prove a transit signal.",
            "Cannot measure false positives, specificity, ROC AUC or true FAP.",
            "Uncertain BJD ephemerides are not evaluated here.",
            "18 subjects give wide descriptive uncertainty intervals.",
            "No results may be used to retune this one-shot final protocol."
        ]
    }
    (args.outdir/"final_results_protocol.json").write_text(
        json.dumps(outcome,indent=2),encoding="utf-8"
    )
    text=[
        "# Final blinded independent 18-host BLS period/score assessment",
        "",
        "One previously untouched, metadata-only selected cohort,",
        "with frozen methods and thresholds from experiment 27.",
        "NOT a true sensitivity, FPR, AUC, or astrophysical detection study.",
        "",
        f"Direct period concordance: {summary['n_direct_blind_period_agreements']}/18 "
        f"({100*summary['direct_blind_period_concordance']:.2f}%).",
        f"95% exact binomial interval, descriptive only: "
        f"[{summary['direct_binomial95_ci_DESCRIPTIVE_ONLY'][0]:.4f}, "
        f"{summary['direct_binomial95_ci_DESCRIPTIVE_ONLY'][1]:.4f}].",
        f"Direct period AND raw threshold accepted: "
        f"{summary['n_correct_and_raw_threshold_accepted']}/18.",
        f"Direct period AND multievent evidence accepted: "
        f"{summary['n_correct_and_evidence_threshold_accepted']}/18.",
        f"Wrong period, but raw threshold accepted: "
        f"{summary['n_wrong_period_but_raw_threshold_accepted']}/18.",
        f"Wrong period, but multievent evidence accepted: "
        f"{summary['n_wrong_period_but_evidence_threshold_accepted']}/18.",
        "",
        "## Period-stratum results", "",
        bybin.to_string(index=False),
        "",
        "## Important limitations", "",
        *["- "+v for v in outcome["interpretation_limits"]],
    ]
    (args.outdir/"final_18_host_report.md").write_text(
        "\n".join(text)+"\n",encoding="utf-8"
    )
    print("\n"+"\n".join(text),flush=True)


if __name__=="__main__":
    main()
