#!/usr/bin/env python3
"""Experiment 26: source-clustered PAIRED audit of experiment 25.

NO new data, searches, light curves, MAST, original FITS, injection
optimization, score changes, thresholds, or access to LOCKED_HOLDOUT.

Reuses EXACT frozen 95th-percentile null thresholds from experiment 25.
Cohort: 12 already-examined 'internal_development_audit' TESS hosts,
4 surrogate draws at each block size 24/48/96 (144 surrogate rows)
and 3 physical injections per host (36 source-matched trials).

Two distinct paired endpoints:
  - difference in surrogate-null exceedance fractions, score versus raw;
  - difference in *correct injected-period and threshold-passing* fractions.
Each data row carries both scores (paired). Inference resamples the
INDEPENDENT HOST TIC, keeping all trials of that host together; paired
exact sign-flip across 12 host-level average differences is optional
exploratory reference only (never a confirmatory test on a fresh cohort).

Comparisons are conditioned on a known surrogate model and a development
cohort already inspected repeatedly. Results are NOT astrophysical FPR,
transit specificity, FAP, or final test generalization.

Run:
  python -u revision/26_paired_development_score_audit.py
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/"results/development_surrogate_null"
NULL=BASE/"per_surrogate_null.csv"
INJ=BASE/"prior_injected_scores.csv"
SCORES=BASE/"calibration_and_internal_audit.csv"
SPLIT=BASE/"null_protocol.json"
OUT=ROOT/"results/development_paired_score_audit"
SEED=20261028
BOOTSTRAPS=10000
BLOCKS=(24,48,96)
N_HOSTS=12
N_NULL_PER_HOST=12
N_INJ_PER_HOST=3


def hash_file(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for part in iter(lambda:f.read(1024*1024),b""):
            h.update(part)
    return h.hexdigest()


def source_cluster_bootstrap_ci(a,b,seed,n_boot=BOOTSTRAPS):
    """Paired mean(b-a) and percentile CI, resampling host TIC only."""
    av=np.asarray(a,dtype=float)
    bv=np.asarray(b,dtype=float)
    if len(av)!=len(bv) or len(av)<2 or not (
        np.isfinite(av).all() and np.isfinite(bv).all()
    ):
        raise ValueError("Invalid paired host-level bootstrap inputs")
    delta=bv-av
    rng=np.random.default_rng(seed)
    ids=rng.integers(0,len(delta),size=(n_boot,len(delta)))
    resampled=delta[ids].mean(axis=1)
    return {
        "mean_source_cluster_difference":float(delta.mean()),
        "bootstrap95_difference_low":float(np.quantile(resampled,.025)),
        "bootstrap95_difference_high":float(np.quantile(resampled,.975)),
        "n_distinct_host_clusters":len(delta),
    }


def paired_signflip_exact(a,b):
    """Exploratory random-sign reference on source-level paired means.

    Each of 12 TIC clusters has one delta. The null is symmetric around
    zero; no independence of 144 surrogate trials is assumed.
    """
    a=np.asarray(a,dtype=float)
    b=np.asarray(b,dtype=float)
    if len(a)!=len(b) or len(a)<1 or len(a)>14:
        raise ValueError("Source-cluster signflip supports 1..14 sources")
    d=b-a
    observed=abs(float(np.mean(d)))
    combinations=np.array(
        list(itertools.product([-1.,1.],repeat=len(d))),
        dtype=float,
    )
    null=np.abs(combinations@d/len(d))
    return float(np.mean(null>=observed-1e-12))


def validate_inputs(nul,inj,summary,protocol):
    needs={"tic_id","development_split","block_cadences","draw",
           "raw_log_power","multievent_evidence"}
    if not needs.issubset(nul):
        raise ValueError("Missing original experiment25 surrogate fields")
    if not {"tic_id","development_split","pair_id",
            "true_injected_period_recovered","raw_log_power",
            "multievent_evidence"}.issubset(inj):
        raise ValueError("Incomplete original injection score table")
    n_cal=set(protocol["calibration_null_tic_set"])
    n_audit=set(protocol["internal_audit_tic_set"])
    if len(n_cal)!=24 or len(n_audit)!=12 or n_cal&n_audit:
        raise ValueError("Prior 24/12 source split is invalid")
    if len(nul)!=432 or nul.tic_id.nunique()!=36:
        raise ValueError("Expected all frozen 432 surrogate records")
    if (nul.duplicated(["tic_id","block_cadences","draw"]).any() or
        not nul.groupby("tic_id").size().eq(12).all()):
        raise ValueError("Null source-block-draw keys inconsistent")
    if set(nul.block_cadences.astype(int))!=set(BLOCKS):
        raise ValueError("Surrogate block sizes have changed")
    if not nul.groupby(["tic_id","block_cadences"]).size().eq(4).all():
        raise ValueError("Expected 4 surrogate draws per TIC/block")
    if len(inj)!=108 or inj.pair_id.duplicated().any():
        raise ValueError("Frozen 108 injection trials changed")
    if not inj.groupby("tic_id").size().eq(3).all():
        raise ValueError("Expected exactly 3 injections per TIC")
    for table in (nul,inj):
        if set(table.tic_id)!=n_cal|n_audit:
            raise ValueError("Source TIC set mismatch")
        for r in table[["tic_id","development_split"]].drop_duplicates().itertuples():
            expected=("calibration" if r.tic_id in n_cal
                      else "internal_development_audit")
            if r.development_split!=expected:
                raise ValueError("Changed source cluster split")
    selected=summary.pivot(index="score_name",columns="source_split",
                           values="threshold_95pct_calibration_surrogate")
    if set(selected.index)!={"raw_log_power","multievent_evidence"}:
        raise ValueError("Unexpected score list")
    threshold={}
    for score in selected.index:
        value=float(selected.loc[score,"calibration"])
        if not np.isclose(
            value,
            float(selected.loc[score,"internal_development_audit"]),
            rtol=0,atol=1e-10
        ):
            raise ValueError("Experimental threshold changed across split")
        # Validate frozen threshold against CALIBRATION surrogates only.
        cal=nul[(nul.development_split=="calibration")][score]
        recalculated=float(np.quantile(cal,.95,method="higher"))
        if not np.isclose(recalculated,value,rtol=0,atol=1e-10):
            raise ValueError(
                f"Published null quantile threshold drift for {score}"
            )
        threshold[score]=value
    return threshold,n_audit


def evaluate(nul,inj,thresholds,audit_tics,n_boot=BOOTSTRAPS):
    null=nul[nul.tic_id.isin(audit_tics)].copy()
    synthetic=inj[inj.tic_id.isin(audit_tics)].copy()
    if len(null)!=144 or len(synthetic)!=36:
        raise ValueError("Expected exactly 144 internal nulls and 36 injections")
    for frame in (null,synthetic):
        frame["accept_raw"]=frame.raw_log_power>thresholds["raw_log_power"]
        frame["accept_evidence"]=(
            frame.multievent_evidence>thresholds["multievent_evidence"]
        )
    synthetic["correct"]=synthetic.true_injected_period_recovered.astype(bool)
    synthetic["correct_and_raw"]=synthetic.correct & synthetic.accept_raw
    synthetic["correct_and_evidence"]=(
        synthetic.correct & synthetic.accept_evidence
    )

    # Every surrogate draw is paired, but source TIC is resampling unit.
    null_by_star=null.groupby("tic_id",sort=True).agg(
        raw=("accept_raw","mean"),
        evidence=("accept_evidence","mean"),
        n=("accept_raw","size"),
        raw_n=("accept_raw","sum"),
        evidence_n=("accept_evidence","sum"),
    ).reset_index()
    inj_by_star=synthetic.groupby("tic_id",sort=True).agg(
        correct=("correct","sum"),
        raw=("correct_and_raw","mean"),
        evidence=("correct_and_evidence","mean"),
        n=("correct","size"),
        accepted_raw=("accept_raw","sum"),
        accepted_evidence=("accept_evidence","sum"),
    ).reset_index()
    if (
        not null_by_star.n.eq(12).all() or
        not inj_by_star.n.eq(3).all() or
        not null_by_star.tic_id.equals(inj_by_star.tic_id)
    ):
        raise ValueError("Internal source cluster pairing changed")

    null_diff=source_cluster_bootstrap_ci(
        null_by_star.raw,null_by_star.evidence,SEED,n_boot
    )
    inj_diff=source_cluster_bootstrap_ci(
        inj_by_star.raw,inj_by_star.evidence,SEED+1,n_boot
    )
    null_diff["source_level_signflip_exact_p_exploratory"]=paired_signflip_exact(
        null_by_star.raw,null_by_star.evidence
    )
    inj_diff["source_level_signflip_exact_p_exploratory"]=paired_signflip_exact(
        inj_by_star.raw,inj_by_star.evidence
    )

    null_blocks=null.groupby("block_cadences",sort=True).agg(
        n_trials=("tic_id","size"),
        n_sources=("tic_id","nunique"),
        raw_exceedances=("accept_raw","sum"),
        evidence_exceedances=("accept_evidence","sum"),
    ).reset_index()
    null_blocks["raw_exceedance_fraction"]=(
        null_blocks.raw_exceedances/null_blocks.n_trials
    )
    null_blocks["evidence_exceedance_fraction"]=(
        null_blocks.evidence_exceedances/null_blocks.n_trials
    )
    null_blocks["absolute_evidence_minus_raw"]=(
        null_blocks.evidence_exceedance_fraction-
        null_blocks.raw_exceedance_fraction
    )

    summary={
        "analysis_role":"exploratory paired INTERNAL DEVELOPMENT audit",
        "n_host_TIC_clusters":N_HOSTS,
        "n_surrogate_null_trials":len(null),
        "n_physical_injection_trials":len(synthetic),
        "frozen_raw_log_power_threshold":
            thresholds["raw_log_power"],
        "frozen_evidence_threshold":
            thresholds["multievent_evidence"],
        "surrogate_exceedance_raw":int(null.accept_raw.sum()),
        "surrogate_exceedance_evidence":
            int(null.accept_evidence.sum()),
        "surrogate_exceedance_raw_fraction":float(
            null.accept_raw.mean()
        ),
        "surrogate_exceedance_evidence_fraction":float(
            null.accept_evidence.mean()
        ),
        "surrogate_difference_paired_cluster":null_diff,
        "physical_injections_with_correct_blind_period":int(
            synthetic.correct.sum()
        ),
        "correct_period_and_raw_threshold":int(
            synthetic.correct_and_raw.sum()
        ),
        "correct_period_and_evidence_threshold":int(
            synthetic.correct_and_evidence.sum()
        ),
        "all_injections_above_raw_threshold":
            int(synthetic.accept_raw.sum()),
        "all_injections_above_evidence_threshold":
            int(synthetic.accept_evidence.sum()),
        "correct_period_and_accepted_difference_paired_cluster":
            inj_diff,
        "raw_score_accepted_wrong_injection_periods":int(
            (synthetic.accept_raw & ~synthetic.correct).sum()
        ),
        "evidence_score_accepted_wrong_injection_periods":int(
            (synthetic.accept_evidence & ~synthetic.correct).sum()
        ),
        "surrogate_exceedance_by_block":null_blocks.to_dict("records"),
    }
    return null,synthetic,null_by_star,inj_by_star,null_blocks,summary


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--outdir",type=Path,default=OUT)
    ap.add_argument("--bootstrap",type=int,default=BOOTSTRAPS)
    args=ap.parse_args()
    if not 100<=args.bootstrap<=100000:
        raise ValueError("Bootstrap replicates must be 100..100000")
    nul=pd.read_csv(NULL,dtype={"tic_id":str})
    inj=pd.read_csv(INJ,dtype={"tic_id":str})
    thresholds=pd.read_csv(SCORES)
    protocol=json.loads(SPLIT.read_text(encoding="utf8"))
    frozen,audit=validate_inputs(nul,inj,thresholds,protocol)
    null,synthetic,stars_null,stars_inj,blocks,results=evaluate(
        nul,inj,frozen,audit,n_boot=args.bootstrap
    )
    args.outdir.mkdir(parents=True,exist_ok=True)
    stars_null.to_csv(
        args.outdir/"paired_null_exceedance_by_star.csv",index=False
    )
    stars_inj.to_csv(
        args.outdir/"paired_injection_acceptance_by_star.csv",index=False
    )
    blocks.to_csv(
        args.outdir/"block_duration_sensitivity.csv",index=False
    )
    details={
        "experiment":26,
        "source_experiment":25,
        "source_tics_previously_explored":True,
        "internal_audit_NOT_final_independent_test":True,
        "original_24_TIC_calibration_thresholds_unchanged":True,
        "calibration_scores_not_refitted":True,
        "novel_detector_hyperparameters_not_tuned":True,
        "no_additional_FITS_or_null_simulations":True,
        "locked_18_host_holdout_not_accessed":True,
        "bootstrap_resampling_unit":"TIC, with all 3 injections or 12 null surrogates retained as a block",
        "source_cluster_exact_signflip_test":"exploratory only; 12 host clusters, not 144 independent nulls",
        "results":results,
        "SHA256":{
            "per_surrogate_null":hash_file(NULL),
            "prior_injected_scores":hash_file(INJ),
            "original_thresholds":hash_file(SCORES),
            "original_protocol":hash_file(SPLIT)
        },
        "limitations":[
            "The original stars, detection scores and null model have been examined during development.",
            "Nulls are moving-block resamples of positive hosts; they may contain authentic transit fragments.",
            "Repeated synthetic null draws from each host are dependent; unit of inference is TIC.",
            "No inference of real astrophysical false-positive rate, true-negative specificity, or detection false-alarm probability.",
            "Exact sign-flip p values are descriptive on 12 previously examined source clusters; cannot establish novelty or final validation.",
            "The synthetic injection profiles used assumed solar host properties and may overlap real transit signals.",
            "12 host clusters are insufficient for precise or broadly generalizable calibration.",
            "A truly independent astrophysical-control cohort is still needed.",
        ]
    }
    (args.outdir/"paired_audit_protocol.json").write_text(
        json.dumps(details,indent=2,default=float),encoding="utf8"
    )
    r=[
        "# Experiment 26: paired source-clustered audit of development null scoring",
        "",
        "All results exploratory, from 12 **already studied** development hosts.",
        "Does not establish astrophysical real-world FPR, specificity or FAP.",
        "",
        "## Frozen thresholds (no retuning)",
        f"Log-power: {frozen['raw_log_power']:.6f}",
        f"Multi-event evidence: {frozen['multievent_evidence']:.6f}",
        "",
        "## Original 12-host / 144-surrogate null audit",
        f"Raw null exceedances: {results['surrogate_exceedance_raw']}/144.",
        f"Evidence null exceedances: {results['surrogate_exceedance_evidence']}/144.",
        f"Mean paired source-level difference (evidence minus raw): "
        f"{results['surrogate_difference_paired_cluster']['mean_source_cluster_difference']:.6f}.",
        f"Paired host-bootstrap 95% CI: "
        f"[{results['surrogate_difference_paired_cluster']['bootstrap95_difference_low']:.6f},"
        f" {results['surrogate_difference_paired_cluster']['bootstrap95_difference_high']:.6f}].",
        f"Exploratory host-signflip p: "
        f"{results['surrogate_difference_paired_cluster']['source_level_signflip_exact_p_exploratory']:.6f}.",
        "",
        "## Frozen 36 physical-injection trials on the same 12 hosts",
        f"Correct blind injected periods before threshold: "
        f"{results['physical_injections_with_correct_blind_period']}/36.",
        f"Correct-period AND accepted raw: "
        f"{results['correct_period_and_raw_threshold']}/36.",
        f"Correct-period AND accepted evidence: "
        f"{results['correct_period_and_evidence_threshold']}/36.",
        f"Accepted raw with wrong injected period: "
        f"{results['raw_score_accepted_wrong_injection_periods']}/36.",
        f"Accepted evidence with wrong injected period: "
        f"{results['evidence_score_accepted_wrong_injection_periods']}/36.",
        f"Paired host-level difference (correct AND accepted; evidence minus raw): "
        f"{results['correct_period_and_accepted_difference_paired_cluster']['mean_source_cluster_difference']:.6f}.",
        f"Paired host-bootstrap 95% CI: "
        f"[{results['correct_period_and_accepted_difference_paired_cluster']['bootstrap95_difference_low']:.6f},"
        f" {results['correct_period_and_accepted_difference_paired_cluster']['bootstrap95_difference_high']:.6f}].",
        f"Exploratory host-signflip p: "
        f"{results['correct_period_and_accepted_difference_paired_cluster']['source_level_signflip_exact_p_exploratory']:.6f}.",
        "",
        "## Dependence on surrogate block size", "",
        blocks.to_string(index=False),
        "", "## Limitations", "",
        *["- "+x for x in details["limitations"]],
        "- The 18 final locked host sources remain unexamined.",
    ]
    (args.outdir/"paired_audit_report.md").write_text(
        "\n".join(r)+"\n",encoding="utf8"
    )
    print("\n".join(r),flush=True)


if __name__=="__main__":
    main()
