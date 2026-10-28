#!/usr/bin/env python3
"""Experiment 21: development-only, blind BLS versus multi-event evidence.

This is a MODEL-DEVELOPMENT experiment on the 36 *development* TICs frozen
in experiment 20. It NEVER accesses the 18 LOCKED_HOLDOUT sources.

Photometry remains original-MAST TGLC_APERTURE_1800S_V1; the scan is
exactly the prior BLS grid (4000 periods from 0.5 to 15 days, durations
0.05..0.25 days). Every period candidate is fitted WITHOUT NASA period/epoch.
NASA period concordance is checked AFTER all candidate selection.

Two fixed experimental methods:
 (A) original global highest BLS power;
 (B) highest BLS power meeting all:
     - >=3 distinct observed transit windows, with >=2 cadences each;
     - >=3 of those have a positive median flux deficit;
     - no single positive window contributes >60% of summed deficits;
     - <=20% of points inside the model box have nonzero TGLC flags.
The B rule is a development hypothesis, NOT a validated detection rule.
No transit significance/FPR, specificity, AUC or true-positive rate is
claimed. Catalogued positive hosts are NOT a control-negative cohort.
No sigma-clipping, detrending, quality-filter modifications, BLS rescan
after reference-period comparison, source replacement or heldout access.

Example:
 python -u revision/21_dev_blind_multievent_bls.py

Raw/processed curves live outside Git; output aggregate CSV and diagnostic
figures only. Output reports remain labelled DEVELOPMENT/EXPLORATORY.
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
from astropy.timeseries import BoxLeastSquares
from scipy.stats import beta

ROOT=Path(__file__).resolve().parents[1]
SELECTION=ROOT/"results/future_tess_cohorts/development_products.csv"
ACQUISITION=ROOT/"results/new_tess_development_acquisition/verified_original_fits.csv"
PHOTOMETRY=ROOT/"results/new_tess_development_photometry/aperture_curve_metrics.csv"
PRIVATE=Path.home()/"exoplanet_revision/data/new_tess_development_cleaned"
OUT=ROOT/"results/development_multievent_bls"
POLICY="TGLC_APERTURE_1800S_V1"
PERIOD_GRID=np.linspace(.5,15.0,4000)
DURATION_GRID=np.array([.05,.10,.15,.20,.25],dtype=float)
MIN_2PT_EVENTS=3
MIN_POSITIVE_EVENT_DROPS=3
MAX_SINGLE_EVENT_DROP_SHARE=.60
MAX_INBOX_TGLC_FLAG_FRACTION=.20
MAX_RELATIVE_PERIOD_ERROR=.02


def sha256(p):
    h=hashlib.sha256()
    with Path(p).open("rb") as f:
        for chunk in iter(lambda:f.read(1024*1024), b""):
            h.update(chunk)
    return h.hexdigest()


def exact_binomial_ci(k,n):
    if n<1 or not 0<=k<=n:
        raise ValueError("Invalid binomial count")
    return (
        float(beta.ppf(.025,k,n-k+1)) if k else 0.,
        float(beta.ppf(.975,k+1,n-k)) if k<n else 1.,
    )


def file_for(tic, sector, directory):
    return directory/f"{tic}_sector_{int(sector):02d}_{POLICY}.csv"


def verify_sources(selection,acquisition,metrics,private):
    if any(len(x)!=36 for x in (selection,acquisition,metrics)):
        raise ValueError("Expected EXACTLY 36 development source rows in all inputs")
    if not all(set(x.tic_id)==set(selection.tic_id) for x in
               (acquisition,metrics)):
        raise ValueError("Development acquisition and photometry TIC sets differ")
    for df in (selection,acquisition,metrics):
        if df.tic_id.duplicated().any():
            raise ValueError("Duplicate TIC in input")
    if not selection.role.astype(str).eq("development").all():
        raise ValueError("Selection includes something other than development TICs")
    if "is_frozen_test_tic" in selection.columns:
        if selection.is_frozen_test_tic.astype(str).str.lower().eq("true").any():
            raise ValueError("LOCKED_HOLDOUT TIC in development CSV")
    if not metrics.policy_id.eq(POLICY).all():
        raise ValueError("Unexpected photometry policy")
    if not acquisition.status.eq("verified_original_fits").all():
        raise ValueError("Unverified original MAST FITS")
    for field in ("sector","mast_uri","author","exptime_seconds"):
        m=selection[["tic_id",field]].merge(
            acquisition[["tic_id",field]],on="tic_id",
            validate="one_to_one",suffixes=("_sel","_fits"))
        if field=="exptime_seconds":
            if not np.allclose(m.exptime_seconds_sel,m.exptime_seconds_fits):
                raise ValueError("Changed MAST original exposure")
        else:
            if not m[f"{field}_sel"].astype(str).equals(
                m[f"{field}_fits"].astype(str)):
                raise ValueError("Selected MAST product was changed")
    missing=[
        file_for(row.tic_id,row.sector,private)
        for row in selection.itertuples(index=False)
        if not file_for(row.tic_id,row.sector,private).is_file()
    ]
    if missing:
        raise FileNotFoundError(
            f"{len(missing)} frozen development curves absent; "
            f"first missing {missing[0]}. No benchmark is run."
        )


def load_one(row,metrics_row,private):
    p=file_for(row["tic_id"],row["sector"],private)
    expected=Path(str(metrics_row["private_cleaned_file"]))
    if p.resolve()!=expected.resolve():
        raise ValueError(f"Expected policy-14 filename differs for {row['tic_id']}")
    data=pd.read_csv(p)
    required={"time_btjd","cadence_num","flux_cal_aper_normalized",
              "tess_flags","tglc_flags"}
    if not required.issubset(data):
        raise ValueError("Missing frozen aperture photometry columns")
    if len(data)!=int(metrics_row["n_accepted"]):
        raise ValueError("Accepted cadence count changed")
    if data.cadence_num.duplicated().any():
        raise ValueError("Duplicate original TESS cadence identifiers")
    time=data.time_btjd.to_numpy(dtype=float)
    flux=data.flux_cal_aper_normalized.to_numpy(dtype=float)
    flags=data.tglc_flags.to_numpy(dtype=np.int64)
    if (
        len(time)<300 or not np.isfinite(time).all()
        or not np.isfinite(flux).all()
        or not (np.diff(time)>0).all()
    ):
        raise ValueError("Frozen photometry has invalid times or flux")
    return time,flux,flags,sha256(p)


def robust_noise(flux):
    delta=np.diff(np.asarray(flux,float))
    value=1.4826*np.median(np.abs(delta-np.median(delta)))/np.sqrt(2)
    if not np.isfinite(value) or value<=0:
        raise ValueError("No valid white-noise surrogate")
    return float(value)


def blind_bls(time,flux):
    """Important: no NASA orbital truth/ephemerides among function inputs."""
    noise=robust_noise(flux)
    model=BoxLeastSquares(time,flux,dy=np.full(len(time),noise))
    scan=model.power(PERIOD_GRID,DURATION_GRID)
    powers=np.asarray(scan.power,dtype=float)
    if len(powers)!=len(PERIOD_GRID) or not np.isfinite(powers).any():
        raise ValueError("Invalid BLS periodogram")
    return scan,powers,noise


def event_support(time,flux,flags,period,epoch,duration):
    """Descriptive event coherence; not a statistical detection threshold."""
    if not (0<duration<period):
        raise ValueError("Inadmissible BLS duration")
    cyc=np.rint((time-epoch)/period).astype(np.int64)
    residual=time-(epoch+cyc*period)
    box=np.abs(residual)<=duration/2
    n_box=int(box.sum())
    q=np.asarray(flags,dtype=np.int64)
    flag_fraction=float(np.count_nonzero(q[box])/n_box) if n_box else 1.
    outside=flux[~box]
    if len(outside)<10:
        return {
            "n_box_points":n_box,"n_supported_events":0,
            "n_positive_events":0,
            "largest_positive_drop_share":1.,
            "fraction_TGLC_flagged_inside":flag_fraction,
        }
    baseline=float(np.median(outside))
    unique,counts=np.unique(cyc[box],return_counts=True)
    keep=unique[counts>=2]
    drops=[
        float(baseline-np.median(flux[box & (cyc==i)]))
        for i in keep
    ]
    positive=np.asarray([x for x in drops if x>0],dtype=float)
    share=(
        float(positive.max()/positive.sum())
        if len(positive) and positive.sum()>0 else 1.
    )
    return {
        "n_box_points":n_box,
        "n_supported_events":len(keep),
        "n_positive_events":len(positive),
        "largest_positive_drop_share":share,
        "fraction_TGLC_flagged_inside":flag_fraction,
    }


def accepted_events(d):
    return bool(
        d["n_supported_events"]>=MIN_2PT_EVENTS
        and d["n_positive_events"]>=MIN_POSITIVE_EVENT_DROPS
        and d["largest_positive_drop_share"]<=MAX_SINGLE_EVENT_DROP_SHARE
        and d["fraction_TGLC_flagged_inside"]<=MAX_INBOX_TGLC_FLAG_FRACTION
    )


def choose_candidates(time,flux,flags,scan,powers,progress=None):
    """Compare a baseline peak and a filtered peak, WITHOUT NASA truth."""
    valid=np.flatnonzero(np.isfinite(powers))
    if not len(valid):
        raise ValueError("No finite blind BLS powers")
    sorted_ix=valid[np.argsort(-powers[valid],kind="stable")]
    original=int(sorted_ix[0])

    # No top-K truncation. Every BLS grid point is eligible for the
    # independent-event check. The first supported one wins by fixed power.
    constrained=None
    visits=0
    for ix in sorted_ix:
        visits+=1
        if progress is not None and visits%500==0:
            progress(visits,len(sorted_ix))
        d=event_support(
            time,flux,flags,
            float(scan.period[ix]),
            float(scan.transit_time[ix]),
            float(scan.duration[ix])
        )
        if accepted_events(d):
            constrained=(int(ix),d)
            break
    original_support=event_support(
        time,flux,flags,
        float(scan.period[original]),
        float(scan.transit_time[original]),
        float(scan.duration[original])
    )
    return original,original_support,constrained,visits


def result_at(scan,powers,index):
    if index is None:
        return {
            "period_days":np.nan,
            "duration_days":np.nan,
            "transit_epoch_btjd":np.nan,
            "bls_power_white_noise_proxy":np.nan,
            "depth_fraction":np.nan,
            "depth_snr_white_noise_proxy":np.nan,
        }
    return {
        "period_days":float(scan.period[index]),
        "duration_days":float(scan.duration[index]),
        "transit_epoch_btjd":float(scan.transit_time[index]),
        "bls_power_white_noise_proxy":float(powers[index]),
        "depth_fraction":float(scan.depth[index]),
        "depth_snr_white_noise_proxy":float(scan.depth_snr[index]),
    }


def match_period(blind,true):
    if not np.isfinite(blind):
        return False,False
    ratio=float(blind)/float(true)
    return (bool(abs(ratio-1)<=.02),
            bool(abs(ratio-.5)<=.01 or abs(ratio-2)<=.04))


def make_figure(tic,t,f,scan,powers,global_ix,constrained_ix,out):
    fig,axs=plt.subplots(2,1,figsize=(9.5,7),layout="constrained")
    axs[0].plot(t-t[0],f,".",ms=2)
    axs[0].set(title=f"{tic} | exploratory 36-TIC development cohort",
               xlabel="Days since first accepted cadence",
               ylabel="Original normalized aperture flux")
    axs[1].plot(PERIOD_GRID,powers,lw=.6)
    axs[1].axvline(float(scan.period[global_ix]),lw=1,ls=":",
                   color="black",label="Global BLS peak")
    if constrained_ix is not None:
        axs[1].axvline(float(scan.period[constrained_ix]),lw=1,
                       ls="--",label="Multi-event eligible peak")
    axs[1].set(xlabel="Blind period (days)",
               ylabel="BLS power (uncalibrated)")
    axs[1].legend(fontsize=8)
    fig.savefig(out/f"{tic}_development_peak_comparison.png",dpi=125)
    plt.close(fig)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--private-dir",type=Path,default=PRIVATE)
    ap.add_argument("--outdir",type=Path,default=OUT)
    args=ap.parse_args()
    selected=pd.read_csv(SELECTION,dtype={"tic_id":str})
    acquired=pd.read_csv(ACQUISITION,dtype={"tic_id":str})
    metrics=pd.read_csv(PHOTOMETRY,dtype={"tic_id":str})
    verify_sources(selected,acquired,metrics,args.private_dir)
    # This script intentionally DOES NOT open or inspect any locked holdout file.
    m={row["tic_id"]:row for row in metrics.to_dict("records")}
    a={row["tic_id"]:row for row in acquired.to_dict("records")}
    results=[]
    errors=[]
    figs=args.outdir/"figures"
    figs.mkdir(parents=True,exist_ok=True)
    for item in selected.sort_values("selection_rank").to_dict("records"):
        tic=str(item["tic_id"])
        try:
            print(f"[START {tic}] Original-TGLC blind development scan",flush=True)
            t,flux,flags,source_sha=load_one(item,m[tic],args.private_dir)
            # The blind scan and candidate selection run BEFORE consulting
            # item["planet_period_days"].
            scan,powers,noise=blind_bls(t,flux)
            global_idx,global_ev,constrained,visits=choose_candidates(
                t,flux,flags,scan,powers,
                progress=lambda n,total:print(
                    f"[SCAN {tic}] inspected {n}/{total} BLS grid candidates",
                    flush=True
                )
            )
            constrained_idx=constrained[0] if constrained is not None else None
            constrained_ev=constrained[1] if constrained is not None else None
            best=result_at(scan,powers,global_idx)
            eligible=result_at(scan,powers,constrained_idx)
            # The reference catalogued period is read ONLY after BLS output.
            planet_p=float(item["planet_period_days"])
            base_ok,base_alias=match_period(best["period_days"],planet_p)
            multi_ok,multi_alias=match_period(eligible["period_days"],planet_p)
            result={
                "tic_id":tic,
                "planet":str(item["planet"]),
                "sector":int(item["sector"]),
                "period_stratum":str(item["period_bin"]),
                "reference_period_days":planet_p,
                "n_cadences":len(t),
                "cleaned_curve_sha256":source_sha,
                "original_fits_sha256":str(a[tic]["sha256"]),
                "point_to_point_white_noise_proxy":noise,
                "baseline_global_period_days":best["period_days"],
                "baseline_global_duration_days":best["duration_days"],
                "baseline_global_power":best["bls_power_white_noise_proxy"],
                "baseline_global_depth_snr_uncalibrated":
                    best["depth_snr_white_noise_proxy"],
                "baseline_direct_match_2pct":base_ok,
                "baseline_half_or_double_alias":base_alias,
                "baseline_supported_events":
                    global_ev["n_supported_events"],
                "baseline_positive_events":
                    global_ev["n_positive_events"],
                "baseline_largest_drop_share":
                    global_ev["largest_positive_drop_share"],
                "baseline_TGLC_flag_fraction":
                    global_ev["fraction_TGLC_flagged_inside"],
                "multievent_period_days":eligible["period_days"],
                "multievent_duration_days":eligible["duration_days"],
                "multievent_power":eligible["bls_power_white_noise_proxy"],
                "multievent_depth_snr_uncalibrated":
                    eligible["depth_snr_white_noise_proxy"],
                "multievent_candidate_exists":constrained is not None,
                "multievent_direct_match_2pct":multi_ok,
                "multievent_half_or_double_alias":multi_alias,
                "multievent_supported_events":(
                    constrained_ev["n_supported_events"]
                    if constrained_ev else 0
                ),
                "multievent_positive_events":(
                    constrained_ev["n_positive_events"]
                    if constrained_ev else 0
                ),
                "multievent_largest_drop_share":(
                    constrained_ev["largest_positive_drop_share"]
                    if constrained_ev else np.nan
                ),
                "multievent_TGLC_flag_fraction":(
                    constrained_ev["fraction_TGLC_flagged_inside"]
                    if constrained_ev else np.nan
                ),
                "grid_periods_inspected_for_multievent":visits,
            }
            results.append(result)
            make_figure(tic,t,flux,scan,powers,global_idx,
                        constrained_idx,figs)
            print(
                f"[{len(results):02d}/36] {tic} "
                f"global={best['period_days']:.4f}d "
                f"multi={eligible['period_days']:.4f}d "
                f"NASA={planet_p:.4f}d "
                f"global_match={base_ok} multi_match={multi_ok} "
                f"tested={visits}",
                flush=True,
            )
        except Exception as exc:
            errors.append({"tic_id":tic,"error":repr(exc)})
            print(f"[FAIL] {tic}: {exc}",flush=True)

    args.outdir.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(errors,columns=["tic_id","error"]).to_csv(
        args.outdir/"development_errors.csv",index=False
    )
    if len(results)!=36:
        raise RuntimeError(
            f"{len(results)}/36 benchmark scans succeeded. "
            "No partial success-rate report will be written. "
            "Inspect development_errors.csv and repair source issues."
        )
    df=pd.DataFrame(results)
    df.to_csv(args.outdir/"per_host_development_comparison.csv",index=False)
    summaries=[]
    for name,col in [
        ("global_BLS_direct", "baseline_direct_match_2pct"),
        ("multievent_BLS_direct", "multievent_direct_match_2pct"),
        ("global_BLS_direct_or_alias", "baseline_direct_match_2pct"),
    ]:
        if name.endswith("_or_alias"):
            success=df.baseline_direct_match_2pct|df.baseline_half_or_double_alias
        else:
            success=df[col]
        k=int(success.sum())
        low,high=exact_binomial_ci(k,len(df))
        summaries.append({
            "metric":name,"success":k,"n_hosts":len(df),
            "fraction":k/len(df),
            "exact_binomial_ci95_low":low,
            "exact_binomial_ci95_high":high,
        })
    scores=pd.DataFrame(summaries)
    scores.to_csv(args.outdir/"development_aggregate.csv",index=False)
    bybin=df.groupby("period_stratum",sort=True).agg(
        n_hosts=("tic_id","size"),
        global_direct=("baseline_direct_match_2pct","sum"),
        multievent_direct=("multievent_direct_match_2pct","sum"),
        multievent_found=("multievent_candidate_exists","sum"),
    ).reset_index()
    bybin.to_csv(args.outdir/"development_by_period_stratum.csv",index=False)

    rules={
        "experiment":"21",
        "role":"MODEL DEVELOPMENT ONLY; NOT AN INDEPENDENT FINAL TEST",
        "sources":36,
        "locked_holdout_FITS_opened":False,
        "previous_18_host_TICs_used_to_fit_detector":False,
        "selection_sha256":sha256(SELECTION),
        "acquisition_manifest_sha256":sha256(ACQUISITION),
        "photometry_summary_sha256":sha256(PHOTOMETRY),
        "policy":POLICY,
        "BLS_period_grid":[.5,15.,4000],
        "BLS_durations_days":DURATION_GRID.tolist(),
        "min_distinct_2_cadence_events":MIN_2PT_EVENTS,
        "min_positive_depth_events":MIN_POSITIVE_EVENT_DROPS,
        "max_one_event_positive_depth_share":MAX_SINGLE_EVENT_DROP_SHARE,
        "max_in_box_nonzero_TGLC_flags_fraction":
            MAX_INBOX_TGLC_FLAG_FRACTION,
        "rank":"descending original BLS power, first passing event support",
        "inspects_all_grid_periods_without_NASA_ephemeris":True,
        "planet_period_is_used_only_for_posthoc_matching":True,
        "no_statistical_detection_threshold_calibrated":True,
        "no_negative_or_false_alarm_rates":True,
        "known_limitations":[
            "Development sample is not an independent final evaluation.",
            "Transit host membership does not prove observed sector contains planet transit.",
            "Positive event depths computed against global out-of-box median may be distorted by trends.",
            "Event coherence and TGLC flags are heuristic; threshold tuning requires separate controls.",
            "BLS SNR assumes a white-noise surrogate, not calibrated physical significance.",
            "Positive-host period agreement may occur by chance, and alias handling is descriptive.",
            "No causal proof of planet signal, false-positive calibration or planet classification."
        ],
    }
    (args.outdir/"experimental_rules.json").write_text(
        json.dumps(rules,indent=2),encoding="utf-8"
    )
    report=[
        "# Experiment 21: blind BLS multi-event heuristic, DEVELOPMENT ONLY","",
        "All 36 previously selected and uniformly processed development TICs.",
        "No final locked-holdout curves accessed. Period truth withheld from scans.","",
        "## Blind-period concordance with NASA (descriptive, NOT sensitivity)", "",
        scores.round(6).to_string(index=False), "",
        "## Orbital-period stratification", "",
        bybin.to_string(index=False), "",
        "## Exploratory diagnostics", "",
        f"Global BLS peaks >=12d: {int((df.baseline_global_period_days>=12).sum())}/36",
        f"Multi-event candidates found: {int(df.multievent_candidate_exists.sum())}/36",
        f"Median grid ranks examined: {df.grid_periods_inspected_for_multievent.median():.0f}",
        "", "## Warning", "",
        *["- "+x for x in rules["known_limitations"]]
    ]
    (args.outdir/"development_report.md").write_text(
        "\n".join(report)+"\n",encoding="utf-8"
    )
    print("\n"+"\n".join(report[:18]),flush=True)


if __name__=="__main__":
    main()
