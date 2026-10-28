#!/usr/bin/env python3
"""Experiment 29: READ-ONLY publication synthesis of frozen exoplanet results.

The final 18-host BLS test is OVER (commit 66446fe). Do NOT rerun it.
This script reads only published CSV/JSON files of experiments 21, 24, 25,
26, 27 and 28. It NEVER touches raw original FITS, cleaned light curves,
the original BLS detector, NASA archive, or holdout acquisition pipelines.

Produces descriptive English-ready methods/results text, transparent
per-source tables and figures. No model fits, threshold selection,
rescanning, injections, null simulations, or post-hoc change to the
frozen 18-host evaluation. All primary interpretations remain fixed.

Usage:
 python -u revision/29_readonly_publication_synthesis.py
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
FINAL=ROOT/"results/final_holdout_benchmark"
FREEZE=ROOT/"results/final_holdout_freeze/frozen_protocol.json"
DEV=ROOT/"results/development_multievent_bls/development_aggregate.csv"
PHYSICAL=ROOT/"results/development_physical_injection/physical_injection_aggregate.csv"
MARGINAL=ROOT/"results/development_physical_injection/physical_injection_marginal.csv"
CALIB=ROOT/"results/development_surrogate_null/calibration_and_internal_audit.csv"
PAIRED=ROOT/"results/development_paired_score_audit/paired_audit_protocol.json"
OUT=ROOT/"results/final_publication_synthesis"
POLICY="TGLC_APERTURE_1800S_V1"


def hash_file(path):
    digest=hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda:stream.read(1024*1024),b""):
            digest.update(block)
    return digest.hexdigest()


def validate(final,blind,freeze,protocol,dev,physical,marginal,calib,paired):
    if len(final)!=18 or len(blind)!=18:
        raise ValueError("Expected immutable 18 final source records")
    if (final.tic_id.duplicated().any() or blind.tic_id.duplicated().any()
        or set(final.tic_id)!=set(blind.tic_id)):
        raise ValueError("Final published blind and scored star sets differ")
    if set(final.tic_id)!=set(freeze["locked_holdout_tics"]):
        raise ValueError("Final 18-star metadata differs from prior freeze")
    if protocol.get("experiment")!=28 or freeze.get("experiment")!=27:
        raise ValueError("Unknown final/frozen protocol")
    if not freeze.get("ready_for_first_final_test"):
        raise ValueError("Final protocol not originally locked")
    if protocol["frozen_protocol_sha256"]!=hash_file(FREEZE):
        raise ValueError("Frozen source code/protocol provenance mismatch")
    if protocol.get("scores_frozen")!=freeze.get("main_thresholds"):
        raise ValueError("Final thresholds differ from frozen development")
    if not (FINAL/"FINAL_TEST_STARTED.json").exists():
        raise ValueError("Original one-shot final-test marker missing")
    marker=json.loads((FINAL/"FINAL_TEST_STARTED.json").read_text())
    if marker["freeze_sha256"]!=hash_file(FREEZE):
        raise ValueError("Original one-shot test marker provenance mismatch")
    d=final.set_index("tic_id").sort_index()
    b=blind.set_index("tic_id").sort_index()
    for field in ("best_BLS_period_days","raw_log_power",
                  "multievent_evidence","cleaned_curve_sha256"):
        if not d[field].astype(str).equals(b[field].astype(str)):
            raise ValueError("Final blind prediction changed: "+field)
    if final.direct_period_match_2pct.sum()!=4:
        raise ValueError("Original 4/18 direct-period outcome changed")
    if final.half_or_double_period_alias.sum()!=0:
        raise ValueError("Original 0/18 alias outcome changed")
    for field,expected in [
        ("raw_threshold_accepted",2),
        ("evidence_threshold_accepted",2),
        ("direct_correct_and_raw_accepted",0),
        ("direct_correct_and_evidence_accepted",2),
    ]:
        if int(final[field].sum())!=expected:
            raise ValueError("Original final result changed: "+field)
    if not final.groupby("period_bin").size().eq(6).all():
        raise ValueError("Frozen final period-stratum quota drift")
    if len(dev)<2 or (
        dev.loc[dev.metric=="global_BLS_direct","success"].iloc[0]!=6
        or dev.loc[dev.metric=="multievent_BLS_direct","success"].iloc[0]!=3
    ):
        raise ValueError("Original 36-host development results changed")
    if (len(physical)!=9 or physical.n_trials.sum()!=108
        or physical.n_direct.sum()!=29):
        raise ValueError("Original 108 synthetic injection metrics changed")
    if len(marginal)!=6 or len(calib)!=4:
        raise ValueError("Development injection/null summary changed")
    if paired.get("source_experiment")!=25 or (
        paired["results"]["n_host_TIC_clusters"]!=12
    ):
        raise ValueError("Frozen internal null/paired experiment provenance differs")
    main=freeze["main_thresholds"]
    for score in main:
        s=calib[calib.score_name==score]
        if len(s)!=2 or not np.allclose(
            s.threshold_95pct_calibration_surrogate,
            main[score],atol=1e-12,rtol=0
        ):
            raise ValueError("Post-freeze threshold drift: "+score)
    # These four known hosts are a result, NOT a newly selected subset.
    if final.tic_id.nunique()!=18 or not (
        np.isfinite(final.best_BLS_period_days.to_numpy(float)).all()
    ):
        raise ValueError("Missing final 18-host period predictions")
    outcomes=protocol["results"]
    if not (
        outcomes["n_direct_blind_period_agreements"]==4
        and outcomes["n_correct_and_evidence_threshold_accepted"]==2
        and outcomes["n_wrong_period_but_evidence_threshold_accepted"]==0
        and outcomes["n_wrong_period_but_raw_threshold_accepted"]==2
    ):
        raise ValueError("Published final protocol outcome changed")


def grouped_success(final):
    records=[]
    for name,col in [
        ("BLS period concordant with known host (2%)","direct_period_match_2pct"),
        ("Correct period AND conventional raw threshold","direct_correct_and_raw_accepted"),
        ("Correct period AND fixed evidence threshold","direct_correct_and_evidence_accepted"),
        ("Wrong period AND conventional raw threshold","raw_threshold_accepted"),
        ("Wrong period AND fixed evidence threshold","evidence_threshold_accepted"),
    ]:
        if name.startswith("Wrong period"):
            count=int((final[col] & ~final.direct_period_match_2pct).sum())
        else:
            count=int(final[col].sum())
        records.append({
            "endpoint":name,
            "n":count,
            "denominator_n_hosts":18,
            "descriptive_fraction":count/18
        })
    return pd.DataFrame(records)


def fig_period_agreement(final,dest):
    f,ax=plt.subplots(figsize=(7.4,5.4),layout="constrained")
    correct=final[final.direct_period_match_2pct]
    wrong=final[~final.direct_period_match_2pct]
    ax.scatter(wrong.planet_period_days,wrong.best_BLS_period_days,
               s=50,marker="o",label=f"Not concordant (n={len(wrong)})")
    ax.scatter(correct.planet_period_days,correct.best_BLS_period_days,
               s=75,marker="D",label=f"Concordant within 2% (n={len(correct)})")
    ax.plot([0,15],[0,15],"k--",lw=1,label="Exact period agreement")
    ax.set(xlim=(0,15),ylim=(0,15),
           xlabel="Catalogue orbital period (days)",
           ylabel="Blind highest-power BLS period (days)",
           title="One-shot independent 18-host period-concordance test")
    ax.legend(fontsize=8,loc="best")
    f.savefig(dest/"figure_1_independent_period_concordance.png",dpi=220)
    f.savefig(dest/"figure_1_independent_period_concordance.pdf")
    plt.close(f)


def fig_correctness_acceptance(final,dest):
    raw_correct=int(final.direct_correct_and_raw_accepted.sum())
    evidence_correct=int(final.direct_correct_and_evidence_accepted.sum())
    raw_wrong=int((final.raw_threshold_accepted &
                   ~final.direct_period_match_2pct).sum())
    evid_wrong=int((final.evidence_threshold_accepted &
                    ~final.direct_period_match_2pct).sum())
    fig,ax=plt.subplots(figsize=(8,5),layout="constrained")
    x=np.arange(2)
    ax.bar(x,[raw_correct,evidence_correct],width=.65,
           label="Concordant periods accepted")
    ax.bar(x,[raw_wrong,evid_wrong],width=.65,
           bottom=[raw_correct,evidence_correct],
           label="Nonconcordant periods accepted")
    ax.set(xticks=x,xticklabels=[
        "Conventional log(BLS power+1)",
        "Frozen multi-event evidence"
    ],ylim=(0,4.5),ylabel="Number of accepted candidate hosts (out of 18)",
        title="Held-out threshold decisions (descriptive, not specificity)")
    ax.legend(loc="upper right")
    for i,(a,b) in enumerate(zip(
        [raw_correct,evidence_correct],[raw_wrong,evid_wrong]
    )):
        ax.text(i,a+b+.10,f"{a} concordant / {b} nonconcordant",
                ha="center",fontsize=9)
    fig.savefig(dest/"figure_2_independent_threshold_decisions.png",dpi=220)
    fig.savefig(dest/"figure_2_independent_threshold_decisions.pdf")
    plt.close(fig)


def fig_injection_grid(physical,dest):
    bins=["1.5-3.5","3.5-6.5","6.5-9.5"]
    radius=[.03,.06,.12]
    cells=np.full((3,3),np.nan)
    counts=np.empty((3,3),dtype=object)
    for i,b in enumerate(bins):
        for j,k in enumerate(radius):
            rows=physical[
                (physical.period_stratum==b) &
                np.isclose(physical.radius_ratio,k)
            ]
            if len(rows)!=1 or rows.n_trials.iloc[0]!=12:
                raise ValueError("Fixed 3x3 injection design differs")
            cells[i,j]=100*float(rows.direct_fraction.iloc[0])
            counts[i,j]=f"{int(rows.n_direct.iloc[0])}/12"
    fig,ax=plt.subplots(figsize=(6.8,5.3),layout="constrained")
    img=ax.imshow(cells,aspect="auto",vmin=0,vmax=100,cmap="viridis")
    ax.set(xticks=range(3),yticks=range(3),
           xticklabels=[f"{x:.2f}" for x in radius],
           yticklabels=bins,
           xlabel="Injected planet/star radius ratio (assumed solar host)",
           ylabel="Injected period range (days)",
           title="Development-only physical injections: recovered period")
    for i in range(3):
        for j in range(3):
            color="white" if cells[i,j]<50 else "black"
            ax.text(j,i,f"{cells[i,j]:.1f}%\n({counts[i,j]})",
                    ha="center",va="center",color=color,fontsize=10)
    cbar=fig.colorbar(img,ax=ax,label="Direct recovered periods (%)")
    fig.savefig(dest/"figure_3_development_injections.png",dpi=220)
    fig.savefig(dest/"figure_3_development_injections.pdf")
    plt.close(fig)


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--outdir",type=Path,default=OUT)
    args=p.parse_args()
    final=pd.read_csv(FINAL/"final_per_host_scored.csv",dtype={"tic_id":str})
    blind=pd.read_csv(FINAL/"final_blind_predictions.csv",dtype={"tic_id":str})
    freeze=json.loads(FREEZE.read_text(encoding="utf-8"))
    protocol=json.loads((FINAL/"final_results_protocol.json").read_text("utf-8"))
    dev=pd.read_csv(DEV)
    physical=pd.read_csv(PHYSICAL)
    marginal=pd.read_csv(MARGINAL)
    calib=pd.read_csv(CALIB)
    paired=json.loads(PAIRED.read_text("utf-8"))
    validate(final,blind,freeze,protocol,dev,physical,marginal,calib,paired)
    output=args.outdir
    charts=output/"figures"
    charts.mkdir(parents=True,exist_ok=True)
    grouped_success(final).to_csv(
        output/"table_final_endpoints.csv",index=False
    )
    final[[
        "tic_id","planet","period_bin","planet_period_days",
        "best_BLS_period_days","direct_period_match_2pct",
        "raw_log_power","multievent_evidence",
        "raw_threshold_accepted","evidence_threshold_accepted",
        "positive_local_events","TGLC_nonzero_fraction",
        "same_period_event_removal_power_ratio"
    ]].to_csv(output/"table_final_18_per_source.csv",index=False)
    physical.to_csv(output/"table_development_3x3_injections.csv",index=False)
    fig_period_agreement(final,charts)
    fig_correctness_acceptance(final,charts)
    fig_injection_grid(physical,charts)
    ci=protocol["results"]["direct_binomial95_ci_DESCRIPTIVE_ONLY"]
    n_raw=int(final.raw_threshold_accepted.sum())
    n_evid=int(final.evidence_threshold_accepted.sum())
    ncorrect=int(final.direct_period_match_2pct.sum())
    accepted_e=final[final.evidence_threshold_accepted]
    accepted_r=final[final.raw_threshold_accepted]
    correct=final[final.direct_period_match_2pct]
    text=[
        "# Read-only, publication-oriented synthesis: TESS transit period validation",
        "",
        "**Status:** FROZEN INDEPENDENT TEST FINISHED. No model refitting, "
        "new test or hyperparameter changes. Read-only reporting on commits "
        "preceding this script.",
        "",
        "## Methods-ready text (English)",
        "",
        "We selected one original 1800-s TGLC product per TIC using a "
        "metadata-only, prespecified protocol, excluding the 369-source "
        "original study and a separately examined 18-host engineering cohort. "
        "Thirty-six confirmed-planet hosts supported detector development; "
        "18 distinct hosts were reserved for a one-shot independent assessment "
        "with six in each of three orbital-period strata (1.5–3.5, "
        "3.5–6.5, and 6.5–9.5 days). All products were subjected to "
        "the identical previously frozen TGLC_APERTURE_1800S_V1 processing "
        "policy. The blind period estimator returned the maximum of an "
        "unchanged 4,000-period (0.5–15 d), five-duration (0.05–0.25 d) "
        "box-least-squares (BLS) search. The independent observational "
        "reference period was read only after blind candidate estimation. "
        "Direct period concordance was defined as a relative discrepancy "
        "within 2%; half/double-period aliases were tabulated separately. "
        "Two confidence criteria were fixed from prior development-only "
        "surrogate-null calibration: log(1+BLS maximum power) and an "
        "event-evidence score combining locally positive events, TGLC "
        "flags, and a same-period leave-one-event influence ratio. "
        "Their thresholds and the evaluation source-code hashes were "
        "Git-committed before any final-holdout flux was accessed. The "
        "final cohort was evaluated exactly once, without retuning.",
        "",
        "## Results-ready text (English)",
        "",
        f"All 18 original TGLC source files were verified, and all 18 "
        f"uniformly processed light curves passed the frozen pipeline. "
        f"In the previously untouched holdout, the global blind BLS "
        f"period agreed with the catalogue period for {ncorrect}/18 hosts "
        f"({100*ncorrect/18:.1f}%; exact binomial 95% interval "
        f"{100*ci[0]:.1f}–{100*ci[1]:.1f}%), with no half/double-period aliases. "
        f"The unmodified BLS-power criterion accepted {n_raw}/18 source "
        f"candidates, neither having a concordant period; the frozen "
        f"multi-event confidence criterion accepted {n_evid}/18, both "
        f"having concordant periods. The remaining two period-concordant "
        f"hosts did not pass either confidence criterion. "
        f"Across the 36 development hosts, direct concordance was 6/36 "
        f"with the original BLS maximum; a previously tested hard "
        f"three-event period-replacement rule recovered only 3/36, "
        f"motivating the separation of period estimation from evidence-based "
        f"candidate acceptance. In 108 development-only physical-injection "
        f"trials, 29/108 injected orbital periods were recovered directly "
        f"by blind BLS; the results are synthetic-injection recovery "
        f"rather than astrophysical detection sensitivity.",
        "",
        "## Limitations-ready text (English)",
        "",
        "The final 18-source sample contains confirmed planet HOSTS, "
        "not independently established in-sector transit detections. "
        "In a given selected sector, the planet may have no observable "
        "transit or its signal may not dominate the light curve. "
        "Therefore catalogue-period concordance should not be interpreted "
        "as transit sensitivity, and catalogue-period discordance does "
        "not automatically demonstrate a false astrophysical detection. "
        "The sample contains no verified transit-free control stars, "
        "so these results cannot quantify astrophysical false-positive "
        "rate, specificity, precision or ROC AUC. The development "
        "null model used block-resampled signals from previously "
        "inspected positive hosts and cannot establish a true no-planet "
        "null. Both accepted candidate groups consist of only two "
        "objects; the observed discrimination is descriptive, not "
        "statistically conclusive or broadly generalizable. "
        "The frozen final cohort must not be used for any subsequent "
        "algorithm or threshold selection.",
        "",
        "## Exact per-host endpoints",
        "",
        "Concordant catalogued-period hosts: "+
        ", ".join(correct.planet.tolist())+".",
        "Frozen event-score accepted hosts: "+
        ", ".join(accepted_e.planet.tolist())+".",
        "Frozen raw-power accepted hosts: "+
        ", ".join(accepted_r.planet.tolist())+".",
        "",
        "## Reproducibility and integrity",
        "",
        f"- Frozen protocol SHA256: {hash_file(FREEZE)}",
        f"- Final results protocol SHA256: "
        f"{hash_file(FINAL/'final_results_protocol.json')}",
        f"- Final per-host scored SHA256: "
        f"{hash_file(FINAL/'final_per_host_scored.csv')}",
        f"- Main thresholds: "
        f"raw={freeze['main_thresholds']['raw_log_power']:.12f}, "
        f"evidence={freeze['main_thresholds']['multievent_evidence']:.12f}",
        "- Read-only figures and tables are DERIVED PUBLICATION MATERIAL, "
        "not a rerun of the locked holdout detector.",
        "- Changes to the frozen code, frozen thresholds or final result "
        "numbers are not permitted.",
    ]
    (output/"publication_results_methods_limitations.md").write_text(
        "\n".join(text)+"\n",encoding="utf-8"
    )
    manifest={
        "study":"Exoplanet TESS frozen first independent 18-host benchmark",
        "experiment":29,
        "type":"READ_ONLY_PUBLICATION_SYNTHESIS",
        "no_BLS_rerun":True,
        "no_fits_or_private_flux_read":True,
        "no_detector_refit":True,
        "final_holdout_results_unchanged":True,
        "frozen_source_code_unchanged":True,
        "original_final_one_shot_commit":"66446fe",
        "input_hashes_sha256":{
            str(p.relative_to(ROOT)):hash_file(p)
            for p in [
                FREEZE,FINAL/"final_per_host_scored.csv",
                FINAL/"final_blind_predictions.csv",
                FINAL/"final_results_protocol.json",DEV,PHYSICAL,
                MARGINAL,CALIB,PAIRED
            ]
        }
    }
    (output/"publication_synthesis_manifest.json").write_text(
        json.dumps(manifest,indent=2)+"\n",encoding="utf-8"
    )
    print("EXP29 PUBLICATION-ONLY SUMMARY: final direct 4/18; "
          "raw 0 concordant of 2 accepted; event 2 concordant of "
          "2 accepted; no rerun or threshold tuning.",flush=True)
    print("Outputs:",output,flush=True)


if __name__=="__main__":
    main()
