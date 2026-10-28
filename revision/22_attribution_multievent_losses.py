#!/usr/bin/env python3
"""Experiment 22: frozen, post-hoc FAILURE ATTRIBUTION for development BLS.

Read ONLY the published experiment-21 per-source CSV; NO TESS FITS, private
curves, new BLS scans, catalogue queries, method changes, or LOCKED_HOLDOUT.

Purpose:
  (1) preserve ORIGINAL paired 6/36 vs 3/36 period-concordance results;
  (2) enumerate each baseline direct-period success vetoed by event rules;
  (3) distinguish alias-only replacements from direct period recovery;
  (4) characterize whether the multi-event veto is best framed as a
      confidence/quality FLAG rather than a replacement period-search rule.

All claims are explicitly development/exploratory. This experiment CANNOT
provide sensitivity, specificity, calibrated FPR, or independent performance.
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
from scipy.stats import binomtest

ROOT=Path(__file__).resolve().parents[1]
INPUT=ROOT/"results/development_multievent_bls/per_host_development_comparison.csv"
EXPERIMENT_RULES=ROOT/"results/development_multievent_bls/experimental_rules.json"
OUT=ROOT/"results/development_bls_failure_attribution"

MIN_EVENTS=3
MIN_POSITIVE=3
MAX_SHARE=.60
MAX_FLAGGED=.20


def sha256(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda:f.read(1024*1024),b""):
            h.update(block)
    return h.hexdigest()


def boolean(x):
    if isinstance(x,(bool,np.bool_)):
        return bool(x)
    if str(x).strip().lower() in ("true","1"):
        return True
    if str(x).strip().lower() in ("false","0"):
        return False
    raise ValueError("Unexpected boolean input: "+str(x))


def failure_reasons(r):
    reasons=[]
    if int(r["baseline_supported_events"])<MIN_EVENTS:
        reasons.append("fewer_than_3_supported_events")
    if int(r["baseline_positive_events"])<MIN_POSITIVE:
        reasons.append("fewer_than_3_positive_events")
    share=float(r["baseline_largest_drop_share"])
    if not np.isfinite(share) or share>MAX_SHARE:
        reasons.append("single_event_dominance_gt_60pct")
    flagged=float(r["baseline_TGLC_flag_fraction"])
    if not np.isfinite(flagged) or flagged>MAX_FLAGGED:
        reasons.append("TGLC_flagged_fraction_gt_20pct")
    return reasons


def analyze(df):
    if not len(df):
        raise ValueError("Empty experiment-21 results")
    if df.tic_id.duplicated().any():
        raise ValueError("Duplicate development TIC")
    requirements={
        "tic_id","planet","period_stratum","reference_period_days",
        "baseline_global_period_days","multievent_period_days",
        "baseline_direct_match_2pct","multievent_direct_match_2pct",
        "baseline_supported_events","baseline_positive_events",
        "baseline_largest_drop_share","baseline_TGLC_flag_fraction",
        "multievent_candidate_exists","multievent_half_or_double_alias",
    }
    if not requirements.issubset(df):
        raise ValueError("Missing experiment-21 columns: "+
                         str(requirements-set(df.columns)))
    rows=[]
    for r in df.to_dict("records"):
        orig=boolean(r["baseline_direct_match_2pct"])
        multi=boolean(r["multievent_direct_match_2pct"])
        found=boolean(r["multievent_candidate_exists"])
        aliases=boolean(r["multievent_half_or_double_alias"])
        reasons=failure_reasons(r)
        # If the original global BLS peak passes the coherence rule,
        # the same ranked period MUST remain selected (no rescanning).
        if not reasons:
            if not found or not np.isclose(
                float(r["baseline_global_period_days"]),
                float(r["multievent_period_days"]),
                rtol=0,atol=1e-10
            ):
                raise ValueError(
                    f"Inconsistent original-vs-multi selection for {r['tic_id']}"
                )
        outcome=(
            "both_correct" if orig and multi else
            "baseline_only" if orig else
            "multievent_only" if multi else
            "neither"
        )
        rows.append({
            "tic_id":str(r["tic_id"]),
            "planet":str(r["planet"]),
            "period_stratum":str(r["period_stratum"]),
            "reference_period_days":float(r["reference_period_days"]),
            "original_period_days":float(r["baseline_global_period_days"]),
            "multievent_period_days":float(r["multievent_period_days"])
                if found else np.nan,
            "original_direct":orig,"multievent_direct":multi,
            "multievent_alias_only":aliases and not multi,
            "multievent_candidate_exists":found,
            "outcome_category":outcome,
            "original_peak_passed_all_event_rules":not reasons,
            "original_peak_veto_reasons":";".join(reasons),
            "original_peak_n_supported_events":int(
                r["baseline_supported_events"]
            ),
            "original_peak_n_positive_events":int(
                r["baseline_positive_events"]
            ),
            "original_peak_max_event_share":float(
                r["baseline_largest_drop_share"]
            ),
            "original_peak_TGLC_flag_fraction":float(
                r["baseline_TGLC_flag_fraction"]
            ),
            "original_peak_event_count_failed":int(
                r["baseline_supported_events"]
            )<MIN_EVENTS,
            "original_peak_quality_flag_failed":float(
                r["baseline_TGLC_flag_fraction"]
            )>MAX_FLAGGED,
        })
    out=pd.DataFrame(rows)
    n=len(out)
    both=int((out.outcome_category=="both_correct").sum())
    lost=int((out.outcome_category=="baseline_only").sum())
    gained=int((out.outcome_category=="multievent_only").sum())
    neither=int((out.outcome_category=="neither").sum())
    assert both+lost+gained+neither==n
    assert both+lost==int(out.original_direct.sum())
    assert both+gained==int(out.multievent_direct.sum())
    discordant=lost+gained
    p=(
        float(binomtest(min(lost,gained),discordant,.5).pvalue)
        if discordant else 1.
    )
    summary={
        "n_stars":n,
        "baseline_direct_success":both+lost,
        "multievent_direct_success":both+gained,
        "absolute_multi_minus_original_percentage_points":
            100*(gained-lost)/n,
        "both_recover":both,
        "original_only_lost":lost,
        "multievent_only_gained":gained,
        "neither_recover":neither,
        "paired_mcnemar_exact_two_sided_p_descriptive":p,
        "original_top_peak_passing_all_event_rules":int(
            out.original_peak_passed_all_event_rules.sum()
        ),
        "original_top_peaks_vetoed":int(
            (~out.original_peak_passed_all_event_rules).sum()
        ),
        "multievent_candidate_found":int(
            out.multievent_candidate_exists.sum()
        ),
        "multievent_alias_only_n":int(out.multievent_alias_only.sum()),
        "lost_original_successes_by_TGLC_flag":int(
            ((out.outcome_category=="baseline_only")&
             out.original_peak_quality_flag_failed).sum()
        ),
        "lost_original_successes_by_min_events":int(
            ((out.outcome_category=="baseline_only")&
             out.original_peak_event_count_failed).sum()
        ),
    }
    return out,summary


def plot(out,folder):
    fig,ax=plt.subplots(figsize=(7.5,5.5),layout="constrained")
    for cat,g in out.groupby("outcome_category"):
        ax.scatter(
            g.reference_period_days,g.original_period_days,
            s=40,label=cat,alpha=.8
        )
    ax.plot([0,15],[0,15],"k--",lw=1,label="identity")
    ax.set(
        xlim=(0,15),ylim=(0,15),
        xlabel="NASA catalogued planet period (days)",
        ylabel="Original blind BLS period (days)",
        title="Development-only outcomes; event rules do not change original result"
    )
    ax.legend(fontsize=8)
    fig.savefig(folder/"original_blind_period_vs_catalogue_by_outcome.png",dpi=150)
    plt.close(fig)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--input",type=Path,default=INPUT)
    parser.add_argument("--outdir",type=Path,default=OUT)
    args=parser.parse_args()
    if not EXPERIMENT_RULES.is_file():
        raise FileNotFoundError("Original experiment-21 rule snapshot absent")
    policy=json.loads(EXPERIMENT_RULES.read_text(encoding="utf-8"))
    if (policy.get("min_distinct_2_cadence_events")!=MIN_EVENTS
        or policy.get("min_positive_depth_events")!=MIN_POSITIVE
        or policy.get("max_one_event_positive_depth_share")!=MAX_SHARE
        or policy.get("max_in_box_nonzero_TGLC_flags_fraction")!=MAX_FLAGGED):
        raise ValueError(
            "Experiment-21 fixed event rules differ from audit: do not "
            "silently change post-hoc thresholds"
        )
    df=pd.read_csv(args.input,dtype={"tic_id":str})
    if len(df)!=36:
        raise ValueError(f"Expected frozen 36 development rows, found {len(df)}")
    out,summary=analyze(df)
    args.outdir.mkdir(parents=True,exist_ok=True)
    figs=args.outdir/"figures"
    figs.mkdir(parents=True,exist_ok=True)
    plot(out,figs)
    out.to_csv(args.outdir/"per_host_failure_attribution.csv",index=False)
    lost=out[out.outcome_category=="baseline_only"]
    lost.to_csv(args.outdir/"lost_correct_periods.csv",index=False)
    cfg={
        "experiment":22,
        "source_experiment":21,
        "source_csv_sha256":sha256(args.input),
        "source_rules_sha256":sha256(EXPERIMENT_RULES),
        "analysis_role":"POST-HOC attribution on known DEVELOPMENT sample",
        "no_model_fitting":True,
        "no_BLS_rerun":True,
        "no_new_thresholds_or_tuning":True,
        "original_metrics_preserved":True,
        "LOCKED_HOLDOUT_data_access":False,
        "important_distinction":(
            "Period estimate from the highest BLS peak and reliability of "
            "its repeated-event evidence are separate scientific outputs."
        ),
        "limits":[
            "Paired exact McNemar p-value on 36 sources is exploratory.",
            "This is not a confirmed transit-detection sensitivity estimate.",
            "A failed 3-event gate may reflect valid two-transit coverage.",
            "Not all planet-host observations show detectable transits.",
            "No negative control stars or calibrated false alarm rates used."
        ],
        "aggregate":summary
    }
    (args.outdir/"audit_config.json").write_text(
        json.dumps(cfg,indent=2),encoding="utf-8"
    )
    report=[
        "# Experiment 22 — Why multi-event BLS lost previously recovered periods",
        "",
        "EXPLORATORY DEVELOPMENT FORENSICS: no model or threshold changes.",
        "",
        f"Original BLS direct recovery: {summary['baseline_direct_success']}/36.",
        f"Multi-event BLS direct recovery: {summary['multievent_direct_success']}/36.",
        f"Net absolute difference: {summary['absolute_multi_minus_original_percentage_points']:.3f} percentage points.",
        f"Paired counts: both={summary['both_recover']}, "
        f"original_only={summary['original_only_lost']}, "
        f"multi_only={summary['multievent_only_gained']}, "
        f"neither={summary['neither_recover']}.",
        f"McNemar exact two-sided p (exploratory, not a confirmatory test): "
        f"{summary['paired_mcnemar_exact_two_sided_p_descriptive']:.5f}.",
        f"Original highest BLS peaks passing every coherence gate: "
        f"{summary['original_top_peak_passing_all_event_rules']}/36.",
        f"Lost directly correct original BLS periods due to TGLC flags: "
        f"{summary['lost_original_successes_by_TGLC_flag']}.",
        f"Lost directly correct original BLS periods due to 3-event gate: "
        f"{summary['lost_original_successes_by_min_events']}.",
        "", "## Previously correct period estimates that were vetoed", "",
        lost[[
            "tic_id","planet","reference_period_days","original_period_days",
            "multievent_period_days","multievent_alias_only",
            "original_peak_veto_reasons",
            "original_peak_n_supported_events",
            "original_peak_TGLC_flag_fraction"
        ]].to_string(index=False) if len(lost) else "(none)",
        "", "## Methodological conclusion", "",
        "Do NOT replace a period estimate solely because its evidence-quality "
        "flag does not pass a fixed multi-event gate.",
        "Instead separate period SEARCH from candidate VETTING / UNCERTAINTY; "
        "these 36 development TICs can guide future methods but not provide "
        "an independent final score.",
        "Do not open the 18 LOCKED holdout TICs before algorithm/protocol freeze."
    ]
    (args.outdir/"attribution_report.md").write_text(
        "\n".join(report)+"\n",encoding="utf-8"
    )
    print("\n".join(report[:17]),flush=True)


if __name__=="__main__":
    main()
