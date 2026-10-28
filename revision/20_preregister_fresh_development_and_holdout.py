#!/usr/bin/env python3
"""Pre-register TWO NEW NASA-confirmed TESS positive-host cohorts.

Experiment 19 exposed failure patterns in an 18-TIC benchmark. Those 18
targets are FOREVER excluded from model optimization and any new independent
performance claim. This inventory draws from the SAME FROZEN NASA snapshot
used for experiment 15 but selects distinct, previously UNUSED source TICs.

Intended use (positive-host recovery only):
   DEVELOPMENT   = 12 TESS confirmed hosts per 1.5-3.5, 3.5-6.5, 6.5-9.5d bin
   LOCKED HOLDOUT=  6 additional hosts per orbital-period bin

Total goal: 36 development + 18 final holdout, no shared TICs. Selection
occurs from archive metadata ONLY, before downloading any new FITS or testing
any new detection algorithm. No labels of non-hosts are invented.

No FITS download; no BLS; no photometric inspection. The LOCKED holdout must
NOT be downloaded, displayed, examined or analyzed during development. Its
performance can be computed ONCE after all model choices are frozen.

Run (metadata only):
   python -u revision/20_preregister_fresh_development_and_holdout.py

Output:
  results/future_tess_cohorts/development_products.csv
  results/future_tess_cohorts/LOCKED_DO_NOT_ANALYZE_holdout_products.csv
  results/future_tess_cohorts/candidate_metadata_audit.csv
  results/future_tess_cohorts/cohort_protocol.json
  results/future_tess_cohorts/cohort_report.md

For methodological transparency, selecting a confirmed planet host is NOT
enough to ensure that the planet transits during its chosen TESS sector.
These data are not true negative detections; specificity is not estimable.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import time
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/"revision/15_new_confirmed_tess_cohort.py"
SNAPSHOT=ROOT/"results/new_confirmed_tess_inventory/archive_snapshot.csv"
PREVIOUS=ROOT/"results/new_confirmed_tess_inventory/selected_prospective_hosts.csv"
OUT=ROOT/"results/future_tess_cohorts"
DEV_PER_BIN=12
HOLDOUT_PER_BIN=6
PERIOD_BINS=["1.5-3.5","3.5-6.5","6.5-9.5"]
MAX_QUERIES=240


def load_discovery_module():
    spec=importlib.util.spec_from_file_location("frozen_tess15",SRC)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sha256(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda:stream.read(1024*1024),b""):
            h.update(chunk)
    return h.hexdigest()


def excluded_tics(discovery):
    """Prior dataset + all previous engineering and experiment-15 targets."""
    old, prior_hash=discovery.previous_tics()
    prev=pd.read_csv(PREVIOUS,dtype={"tic_id":str})
    if len(prev)!=18 or prev.tic_id.duplicated().any():
        raise ValueError("Existing 18-host benchmark not intact")
    frozen18=set(prev.tic_id)
    if not all(t.startswith("TIC_") for t in frozen18):
        raise ValueError("Unexpected frozen test TIC format")
    if old & frozen18:
        raise ValueError("An original TIC was reused in 18-host experiment")
    return old|frozen18,prior_hash,frozen18


def plan(queue, search_result_getter, excluded, query_cap, pause=0.0):
    """Returns selections based only on MAST product availability.

    The first DEV_PER_BIN successfully found distinct TICs within each
    stratum are development; subsequent HOLDOUT_PER_BIN TICs are held out.
    All selection is determined without downloading or viewing any flux.
    """
    counts={b:Counter() for b in PERIOD_BINS}
    chosen=[]
    audit=[]
    complete=lambda:all(
        counts[b]["development"]>=DEV_PER_BIN
        and counts[b]["LOCKED_HOLDOUT"]>=HOLDOUT_PER_BIN
        for b in PERIOD_BINS
    )
    for item in queue.head(query_cap).to_dict("records"):
        if complete():
            break
        tic=item["tic_norm"]
        bin_name=str(item["bin"])
        if tic in excluded:
            raise AssertionError(f"Leakage: previously known {tic}")
        if bin_name not in counts:
            raise ValueError(f"Unexpected orbital-period bin: {bin_name}")
        if counts[bin_name]["LOCKED_HOLDOUT"]>=HOLDOUT_PER_BIN:
            audit.append({
                "queue_rank":int(item["queue_rank"]),"tic_id":tic,
                "period_bin":bin_name,"status":"bin_quota_already_complete",
                "n_products":0,"error":""
            })
            continue
        try:
            products=search_result_getter(item)
            by_sector={}
            for product in products:
                by_sector.setdefault(int(product["sector"]),[]).append(product)
            unique=[
                rows[0] for sector,rows in sorted(by_sector.items())
                if len(rows)==1
            ]
            if not unique:
                status="no_unique_TGLC1800_primary_product"
            else:
                selected=min(unique,key=lambda x:x["sector"])
                purpose=(
                    "development" if counts[bin_name]["development"]<DEV_PER_BIN
                    else "LOCKED_HOLDOUT"
                )
                counts[bin_name][purpose]+=1
                chosen.append({
                    "role":purpose,
                    "selection_rank_within_role":counts[bin_name][purpose],
                    "queue_rank":int(item["queue_rank"]),
                    "tic_id":tic,
                    "hostname":str(item["hostname"]),
                    "planet":str(item["pl_name"]),
                    "planet_period_days":float(item["planet_period_days"]),
                    "period_bin":bin_name,
                    "transit_epoch_jd":item.get("pl_tranmid"),
                    "transit_epoch_time_standard":item.get("pl_tsystemref"),
                    "transit_duration_hours":item.get("pl_trandur"),
                    "transit_epoch_error1":item.get("pl_tranmiderr1"),
                    "transit_epoch_error2":item.get("pl_tranmiderr2"),
                    "planet_radius_ratio":item.get("pl_ratror"),
                    "planet_transit_depth":item.get("pl_trandep"),
                    "prior_benchmark_source":False,
                    "is_frozen_test_tic":purpose=="LOCKED_HOLDOUT",
                    **selected
                })
                status="selected_"+purpose
            audit.append({
                "queue_rank":int(item["queue_rank"]),"tic_id":tic,
                "period_bin":bin_name,"status":status,
                "n_products":len(products),"error":""
            })
        except Exception as exc:
            audit.append({
                "queue_rank":int(item["queue_rank"]),"tic_id":tic,
                "period_bin":bin_name,"status":"MAST_metadata_error",
                "n_products":0,"error":repr(exc)
            })
        if pause>0:
            time.sleep(pause)
    selected=pd.DataFrame(chosen)
    checked=pd.DataFrame(audit)
    if not selected.empty:
        # Each role is a separately frozen cohort; acquisition script 16
        # requires a unique selection_rank in EACH exported CSV. It must
        # not depend on the period-bin-local order.
        selected.insert(
            0,"selection_rank",
            selected.groupby("role",sort=False).cumcount()+1
        )
        if selected.tic_id.duplicated().any():
            raise AssertionError("Duplicate selected TIC")
        if set(selected.tic_id)&excluded:
            raise AssertionError("Leakage into new cohort")
        for role,n in [("development",DEV_PER_BIN),("LOCKED_HOLDOUT",HOLDOUT_PER_BIN)]:
            for b in PERIOD_BINS:
                actual=len(selected[(selected.role==role)&(selected.period_bin==b)])
                if actual>n:
                    raise AssertionError("Period-bin quota overflow")
    return selected,checked,counts,complete()


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--max-queries",type=int,default=MAX_QUERIES)
    parser.add_argument("--pause-seconds",type=float,default=.10)
    parser.add_argument("--outdir",type=Path,default=OUT)
    args=parser.parse_args()
    if args.max_queries<54 or args.max_queries>600:
        raise ValueError("Need query budget >=54 and <=600")
    args.outdir.mkdir(exist_ok=True,parents=True)

    discovery=load_discovery_module()
    excluded,prior_hash,old18=excluded_tics(discovery)
    archive=pd.read_csv(SNAPSHOT,dtype={"tic_id":str})
    queue=discovery.prepare_queue(archive,excluded)
    print(
        f"FROZEN NASA rows={len(archive)}; unused TIC queue={len(queue)}; "
        f"excluded old+previous_18={len(excluded)}",
        flush=True
    )

    import lightkurve as lk
    chosen,audit,counters,complete=plan(
        queue,lambda item:discovery.metadata_rows(item,lk),
        excluded, args.max_queries,args.pause_seconds
    )
    audit.to_csv(args.outdir/"candidate_metadata_audit.csv",index=False)
    dev=chosen[chosen.role=="development"].copy() if len(chosen) else pd.DataFrame()
    hold=chosen[chosen.role=="LOCKED_HOLDOUT"].copy() if len(chosen) else pd.DataFrame()
    # Do not use an incomplete cohort for an independent final evaluation.
    dev.to_csv(args.outdir/"development_products.csv",index=False)
    hold.to_csv(
        args.outdir/"LOCKED_DO_NOT_ANALYZE_holdout_products.csv",index=False
    )
    protocol={
        "study_phase":"preregistered new cohort metadata selection",
        "prior_18_host_result":"4/18 frozen, NEVER retune based on those 18",
        "selection_metadata_source":"same frozen NASA PS snapshot as experiment 15",
        "archive_snapshot_sha256":sha256(SNAPSHOT),
        "original_dataset_sha256":prior_hash,
        "previous_18_selection_sha256":sha256(PREVIOUS),
        "excluded_prior_and_experiment_15_tic_count":len(excluded),
        "prior_18_TICs_EXCLUDED":sorted(old18),
        "period_bins":PERIOD_BINS,
        "development_per_bin":DEV_PER_BIN,
        "holdout_per_bin":HOLDOUT_PER_BIN,
        "target_development_hosts":36,
        "target_untouched_holdout_hosts":18,
        "target_total":54,
        "metadata_only":True,
        "number_of_MAST_queries_recorded":int(
            (audit.status!="bin_quota_already_complete").sum()
        ) if not audit.empty else 0,
        "MAST_query_cap":args.max_queries,
        "complete_predeclared_cohorts":complete,
        "development_count":len(dev),
        "LOCKED_holdout_count":len(hold),
        "selection_tied_to_flux_or_BLS":False,
        "original_TGLC_author":"TGLC",
        "exposure_seconds":1800,
        "primary_mission_sectors":"1-26",
        "sector_choice":"first unique original 1800-second TGLC product",
        "predeclared_photometry":"TGLC_APERTURE_1800S_V1",
        "LOCKED_holdout_data_access_prohibited_during_development":True,
        "no_replacement_after_photometric_quality_or_detection_failure":True,
        "source_positive_label":"NASA-confirmed TRANSITING planet host only",
        "in_sector_planet_transit_required_for_positive_detection_claim":True,
        "no_planet_catalogue_entry_is_not_a_true_negative":True,
        "other_controls_needed_for_specificity":True,
        "planned_training_calibration":[
            "Develop multi-event-support and robust event-coherence algorithms ONLY on new development TICs.",
            "Calibrate detector scoring with null surrogates, physically injected events and separate control cohorts.",
            "Distinguish catalogue host membership from observable in-sector transits.",
            "Freeze all model hyperparameters, outlier and flagging rules before opening LOCKED_HOLDOUT."
        ],
    }
    (args.outdir/"cohort_protocol.json").write_text(
        json.dumps(protocol,indent=2,default=str),encoding="utf-8"
    )
    counts_table=chosen.groupby(["role","period_bin"]).size().to_string() if len(chosen) else "(none)"
    report=[
        "# Prespecified independent TESS model-development and LOCKED holdout cohorts","",
        f"New model-development TICs: {len(dev)}/36",
        f"New locked and untouched final TICs: {len(hold)}/18",
        f"Total distinct newly selected TICs: {len(chosen)}/54",
        f"MAST metadata queries performed: {protocol['number_of_MAST_queries_recorded']}",
        f"MAST metadata errors: {int((audit.status=='MAST_metadata_error').sum()) if len(audit) else 0}",
        f"Frozen design complete: {complete}",
        "", "## Selection by period stratum and role", "",
        counts_table,"",
        "## Integrity and limitations","",
        "- No FITS, light curves, transit metrics, BLS, or model output were inspected.",
        "- The previous 18 hosts and all 369 original study source TICs are excluded.",
        "- DO NOT download or inspect LOCKED_HOLDOUT until model pipeline is frozen.",
        "- Only host-positive labels are known; no true-negative inference is made.",
        "- Sources with missing observed transits must not be silently replaced.",
        "- A future genuinely independent test can only occur ONCE per frozen model.",
    ]
    (args.outdir/"cohort_report.md").write_text(
        "\n".join(report)+"\n",encoding="utf-8"
    )
    print("\n"+"\n".join(report[:13]),flush=True)
    if not complete:
        raise RuntimeError(
            "The prespecified balanced 36/18 independent cohorts are not "
            "complete. Inspect cohort_report.md and metadata audit; do not "
            "analyze the locked holdout or silently change quotas."
        )


if __name__=="__main__":
    main()
