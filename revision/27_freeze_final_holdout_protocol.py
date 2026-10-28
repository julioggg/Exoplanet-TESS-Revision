#!/usr/bin/env python3
"""Experiment 27 — FREEZE the final independent 18-host protocol.

Run and COMMIT this metadata-only protocol BEFORE downloading or opening
any original TGLC FITS or cleaned photometric series of the 18 locked TICs.

This experiment:
  * locks source TICs and selected TGLC 1800s sectors/Mast URIs;
  * hashes the ENTIRE final detector, preprocessing and acquisition scripts;
  * hashes unchanged development score thresholds, null draws and reports;
  * confirms the 36 DEVELOPMENT / 18 LOCKED sets are disjoint;
  * verifies previous 18-host labels and 369 old source TICs remain excluded;
  * freezes all reportable endpoints and scientific limitations.

NO FITS, private photometry, BLS runs, catalogue future updates or holdout
flux is accessed. The final evaluation must run with these exact hashes.
No tuning on final 18 test stars and NO claim of planet-detection FPR,
sensitivity or specificity without separate vetted controls.

Run:
 python -u revision/27_freeze_final_holdout_protocol.py

Then git add results/final_holdout_freeze && git commit && git push.
Only AFTER that is committed, acquire the frozen 18 TGLC originals.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"results/final_holdout_freeze"
PROTOCOL=OUT/"frozen_protocol.json"
COHORT=ROOT/"results/future_tess_cohorts"
SOURCE_LIST=[
    "revision/28_final_host_benchmark.py",
    "revision/23_development_event_evidence.py",
    "revision/25_development_surrogate_null_calibration.py",
    "revision/16_acquire_frozen_new_hosts.py",
    "revision/14_build_uniform_aperture_curves.py",
    "revision/26_paired_development_score_audit.py",
]
DATA_LIST=[
    "results/future_tess_cohorts/LOCKED_DO_NOT_ANALYZE_holdout_products.csv",
    "results/future_tess_cohorts/development_products.csv",
    "results/development_surrogate_null/calibration_and_internal_audit.csv",
    "results/development_surrogate_null/per_surrogate_null.csv",
    "results/development_surrogate_null/null_protocol.json",
    "results/development_paired_score_audit/paired_audit_report.md",
]
OLD=ROOT/"data/dataset/final_dataset.csv"
PREVIOUS=ROOT/"results/new_confirmed_tess_inventory/selected_prospective_hosts.csv"


def hash_file(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):
            h.update(b)
    return h.hexdigest()


def fixed_thresholds():
    source=ROOT/"results/development_surrogate_null/calibration_and_internal_audit.csv"
    nulls=ROOT/"results/development_surrogate_null/per_surrogate_null.csv"
    cfg=ROOT/"results/development_surrogate_null/null_protocol.json"
    s=pd.read_csv(source)
    n=pd.read_csv(nulls,dtype={"tic_id":str})
    p=json.loads(cfg.read_text("utf-8"))
    if (
        len(s)!=4 or len(n)!=432 or n.tic_id.nunique()!=36
        or len(p["calibration_null_tic_set"])!=24
        or len(p["internal_audit_tic_set"])!=12
    ):
        raise ValueError("Wrong development null calibration size")
    thresholds={}
    for score in ("raw_log_power","multievent_evidence"):
        rows=s[s.score_name==score]
        if (len(rows)!=2 or
            set(rows.source_split)!=
            {"calibration","internal_development_audit"}):
            raise ValueError("Incomplete frozen threshold rows for "+score)
        a=float(rows[
            rows.source_split=="calibration"
        ].threshold_95pct_calibration_surrogate.iloc[0])
        b=float(rows[
            rows.source_split=="internal_development_audit"
        ].threshold_95pct_calibration_surrogate.iloc[0])
        if not np.isclose(a,b,atol=1e-12,rtol=0):
            raise ValueError("Threshold differs by audit group")
        cal=n[n.tic_id.isin(p["calibration_null_tic_set"])]
        if (len(cal)!=288 or cal.tic_id.nunique()!=24 or
            not (cal.development_split=="calibration").all()):
            raise ValueError("Null calibration source membership mismatch")
        expected=float(np.quantile(cal[score],.95,method="higher"))
        if not np.isclose(a,expected,atol=1e-12,rtol=0):
            raise ValueError("The original 95%-quantile threshold has changed")
        thresholds[score]=a
    return thresholds


def validate_selected(dev,hold,old,prev):
    if (len(dev)!=36 or len(hold)!=18 or
        dev.tic_id.duplicated().any() or hold.tic_id.duplicated().any()):
        raise ValueError("Wrong fixed development/holdout source counts")
    if not (dev.role=="development").all() or not (
        hold.role=="LOCKED_HOLDOUT").all():
        raise ValueError("Holdout roles changed")
    if not dev.is_frozen_test_tic.astype(str).str.lower().eq("false").all():
        raise ValueError("Development contains frozen final source")
    if not hold.is_frozen_test_tic.astype(str).str.lower().eq("true").all():
        raise ValueError("Frozen final role flag changed")
    if dev.selection_rank.astype(int).sort_values().tolist()!=list(range(1,37)):
        raise ValueError("Development selected ranks drift")
    if hold.selection_rank.astype(int).sort_values().tolist()!=list(range(1,19)):
        raise ValueError("Holdout selected ranks drift")
    tics=set(dev.tic_id)|set(hold.tic_id)
    if (
        len(tics)!=54 or tics&set(old.tic_id) or tics&set(prev.tic_id)
    ):
        raise ValueError("Leakage to previously inspected source TICs")
    if len(set(old.tic_id))!=369 or len(set(prev.tic_id))!=18:
        raise ValueError("Previously excluded source counts changed")
    for tab in (dev,hold):
        if not tab.author.eq("TGLC").all():
            raise ValueError("Only frozen TGLC products permitted")
        if not pd.to_numeric(
            tab.exptime_seconds,errors="coerce"
        ).between(1797,1803).all():
            raise ValueError("Expected original 1800-second TGLC cadence")
        if not tab.sector.astype(int).between(1,26).all():
            raise ValueError("Original primary mission sectors required")
        if not tab.mast_uri.str.startswith("mast:HLSP/tglc/").all():
            raise ValueError("Unknown original MAST product URI")
        if (len(tab.groupby("period_bin"))!=3 or
            tab.groupby("period_bin").size().nunique()!=1):
            raise ValueError("Period strata altered")
    expected_dev=set(dev.groupby("period_bin").size())
    expected_hold=set(hold.groupby("period_bin").size())
    if expected_dev!={12} or expected_hold!={6}:
        raise ValueError("Period-stratum quotas changed")


def prepare():
    paths={p:ROOT/p for p in SOURCE_LIST+DATA_LIST}
    for name,file in paths.items():
        if not file.is_file():
            raise FileNotFoundError("Cannot freeze missing prerequisite: "+name)
    dev=pd.read_csv(COHORT/"development_products.csv",
                    dtype={"tic_id":str})
    hold=pd.read_csv(
        COHORT/"LOCKED_DO_NOT_ANALYZE_holdout_products.csv",
        dtype={"tic_id":str}
    )
    old=pd.read_csv(OLD,dtype={"tic_id":str})
    prev=pd.read_csv(PREVIOUS,dtype={"tic_id":str})
    validate_selected(dev,hold,old,prev)
    thresholds=fixed_thresholds()
    selected_cohort=json.loads(
        (COHORT/"cohort_protocol.json").read_text("utf-8")
    )
    if (
        selected_cohort.get("complete_predeclared_cohorts") is not True
        or selected_cohort.get("target_untouched_holdout_hosts")!=18
        or selected_cohort.get("metadata_only") is not True
    ):
        raise ValueError("Original experiment 20 metadata freeze invalid")
    return {
        "experiment":27,
        "stage":"PRE-FIRST-FINAL-HOLDOUT-FLUX FREEZE",
        "ready_for_first_final_test":True,
        "n_development_tics":36,
        "n_final_host_tics":18,
        "locked_holdout_tics":hold.tic_id.tolist(),
        "locked_holdout_metadata_sha256":hash_file(
            COHORT/"LOCKED_DO_NOT_ANALYZE_holdout_products.csv"
        ),
        "development_metadata_sha256":hash_file(
            COHORT/"development_products.csv"
        ),
        "frozen_sources_hashes":{
            p:hash_file(paths[p]) for p in SOURCE_LIST
        },
        "frozen_data_hashes":{
            p:hash_file(paths[p]) for p in DATA_LIST
        },
        "excluded_old_369_source_tics_sha256":hash_file(OLD),
        "excluded_prior_18_source_tics_sha256":hash_file(PREVIOUS),
        "frozen_photometry_policy":"TGLC_APERTURE_1800S_V1",
        "selected_TGLC_source":"original 1800s primary mission 1-26",
        "period_search_grid":[.5,15.,4000],
        "duration_grid_days":[.05,.10,.15,.20,.25],
        "BLS_noise_method":"1.4826*MAD(diff(flux))/sqrt(2)",
        "blind_candidate_choice":"highest unfiltered BLS power",
        "decision_score_formula":(
            "log1p(BLS_power) * sqrt(min(n_positive_local_events,4)/4)"
            " * sqrt(1-TGLC_nonzero_fraction)"
            " * sqrt(clip(leave_one_event_remaining_1x1_BLS_power_ratio,0,1))"
        ),
        "main_thresholds":thresholds,
        "thresholds_derived_only_from":"experiment-25 original 24-development-host surrogate calibration",
        "threshold_rule":"score STRICTLY GREATER THAN threshold",
        "no_refit_on_holdout":True,
        "no_hyperparameter_tuning_on_holdout":True,
        "final_test_permitted_only_once":True,
        "final_test_primary_descriptive_endpoint":(
            "known-host blind BLS period within 2 percent of NASA orbital period"
        ),
        "secondary_endpoints":[
            "blind period half/double alias",
            "number of BLS periods correctly matched AND passing raw threshold",
            "number of BLS periods correctly matched AND passing evidence threshold",
            "number of wrong BLS periods passing each threshold",
            "predeclared 3-bin orbital period stratification",
            "exact 95 percent binomial interval on 18 host period agreements"
        ],
        "no_true_planet_detection_sensitivity_claim":True,
        "no_real_astrophysical_false_positive_rate_claim":True,
        "no_specificity_AUC_or_FAP_claim":True,
        "host_positive_not_necessarily_in_sector_visible":True,
        "no_precise_ephemeris_metric_for_unresolved_BJD":True,
        "no_fits_or_photometric_holdout_access_in_freeze":True,
        "planned_holdout_acquisition_folder":"results/final_holdout_acquisition",
        "planned_holdout_photometry_folder":"results/final_holdout_photometry",
        "holdout_cleaned_private_path":"~/exoplanet_revision/data/final_holdout_cleaned",
        "planned_first_only_final_output":"results/final_holdout_benchmark",
    }


def main():
    if PROTOCOL.exists():
        raise RuntimeError(
            "Freeze already exists; never overwrite protocol after"
            " possible exposure to final heldout flux."
        )
    # Validate everything and hash the fully written evaluation script
    # before touching any holdout photometry.
    frozen=prepare()
    OUT.mkdir(parents=True,exist_ok=True)
    PROTOCOL.write_text(
        json.dumps(frozen,indent=2,sort_keys=True)+"\n",encoding="utf-8"
    )
    report=[
        "# EXPERIMENT 27: LOCKED, one-shot independent host benchmark",
        "",
        "## Freeze must be GIT-COMMITTED before MAST FITS acquisition",
        "This is a METADATA-ONLY freeze; no final holdout curves inspected.",
        "",
        "Prior known positive host cohort: 18 stars, 4 direct period matches.",
        "Development cohort: 36 stars; frozen independent final cohort: 18.",
        "Final per-bin period quota: 6 hosts × 3 bins.",
        "",
        "## Immutable detector and evaluation",
        "Original TGLC aperture V1; full 4000x5 BLS, highest peak unchanged.",
        "Retain original BLS period; evaluate fixed multi-event confidence.",
        f"Raw BLS log-power threshold: {frozen['main_thresholds']['raw_log_power']:.9f}.",
        f"Multi-event evidence threshold: "
        f"{frozen['main_thresholds']['multievent_evidence']:.9f}.",
        "",
        "## Scientific restrictions",
        "- No post-hoc changes to period search, scoring formula or threshold.",
        "- Report ALL 18 previously selected hosts; never replace nonrecoveries.",
        "- Report direct known-host period agreement, not planet sensitivity.",
        "- No true-negative host group: cannot estimate FPR or specificity.",
        "- No threshold p-value significance claimed; prior work exploratory.",
        "- If data processing fails, document and fix without selecting new TICs.",
        "- DO NOT open a single final host FITS before this commit is pushed.",
        "",
        "## Required order after freeze commit",
        "1. Reuse script 16 to acquire all 18 selected original FITS.",
        "2. Reuse script 14 for unchanged original aperture photometry.",
        "3. Execute script 28 once, after all 18 source checks pass.",
        "4. Preserve independent result even if it is worse than development.",
    ]
    (OUT/"frozen_protocol_report.md").write_text(
        "\n".join(report)+"\n",encoding="utf8"
    )
    print("\n".join(report),flush=True)


if __name__=="__main__":
    main()
