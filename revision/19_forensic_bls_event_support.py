#!/usr/bin/env python3
"""POST-HOC forensic audit of FROZEN experiment-17 BLS period failures.

DO NOT use this script to change the 4/18 primary period recovery rate.
It explains why the largest BLS periodogram peak can be misleading:
  - count independent BLS transit cycles actually sampled;
  - detect whether a single transit window dominates the in-transit deficit;
  - measure how close BLS midpoints fall to observed >3-hour data gaps;
  - inventory TGLC flags inside blind BLS transit windows;
  - count observed catalogue-predicted transit windows WITHOUT fitting them;
  - report implausibly deep BLS box fits and 12-15d peak concentrations.

No BLS rerun, detrending, parameter retuning, model training or per-target
exclusions. The selection and blind scan are those frozen in experiments
15-17. Catalogue timing with ambiguous "BJD" is PROVISIONAL ONLY.

Inputs:
  results/new_confirmed_tess_blind_period/blind_period_per_host.csv
  results/new_confirmed_tess_inventory/selected_prospective_hosts.csv
  results/new_confirmed_tess_photometry/aperture_curve_metrics.csv
  private original-policy-14 cleaned aperture curves.

Run:
  python -u revision/19_forensic_bls_event_support.py
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
BLIND = ROOT/"results/new_confirmed_tess_blind_period/blind_period_per_host.csv"
SELECTED = ROOT/"results/new_confirmed_tess_inventory/selected_prospective_hosts.csv"
METRICS = ROOT/"results/new_confirmed_tess_photometry/aperture_curve_metrics.csv"
OUT = ROOT/"results/new_confirmed_tess_bls_forensics"
PRIVATE = Path.home()/"exoplanet_revision/data/new_confirmed_tess_cleaned"
POLICY="TGLC_APERTURE_1800S_V1"
BTJD_ZERO=2457000.0


def file_hash(p):
    h=hashlib.sha256()
    with Path(p).open("rb") as fp:
        for chunk in iter(lambda: fp.read(1024*1024), b""):
            h.update(chunk)
    return h.hexdigest()


def numeric(v):
    try:
        val=float(v)
        return val if math.isfinite(val) else float("nan")
    except (TypeError, ValueError):
        return float("nan")


def phase_and_cycle(t, period, epoch):
    """For each cadence, return signed time to nearest cycle and its number."""
    t=np.asarray(t,dtype=float)
    if not (period>0 and np.isfinite(period) and np.isfinite(epoch)):
        raise ValueError("Invalid positive period / epoch")
    cycle=np.rint((t-epoch)/period).astype(np.int64)
    residual=t-(epoch+cycle*period)
    return residual,cycle


def count_event_support(t, flux, flag, period, epoch, duration):
    """Descriptive transit-window support with no refitting.

    An 'event with >=2 cadences' counts *distinct orbital cycles*, not
    cadences, and is not automatically an astrophysical transit.
    """
    t=np.asarray(t,dtype=float)
    f=np.asarray(flux,dtype=float)
    q=np.asarray(flag,dtype=np.int64)
    if not (len(t)==len(f)==len(q)) or len(t)<300:
        raise ValueError("Missing/invalid frozen per-star curve")
    if not (np.isfinite(t).all() and np.isfinite(f).all()):
        raise ValueError("Nonfinite time or photometry")
    if not (0<duration<period):
        raise ValueError("Unphysical box duration/period")
    residual, cycle=phase_and_cycle(t,period,epoch)
    inside=np.abs(residual)<=duration/2
    count=int(inside.sum())
    unique, nums=np.unique(cycle[inside],return_counts=True)
    events=[int(n) for n in nums]
    independent_with_two=int(np.count_nonzero(nums>=2))
    out=f[~inside]
    event_drops=[]
    if len(out)>5:
        outside_median=float(np.median(out))
        for c in unique:
            local=f[inside & (cycle==c)]
            event_drops.append(float(outside_median-np.median(local)))
    qinside=q[inside]
    # Original source BLS SNR uses point-to-point white-noise surrogate.
    # These event supports are independent descriptive diagnostics only.
    return {
        "n_in_bls_transit_cadences":count,
        "n_bls_cycles_with_at_least_one_cadence":len(events),
        "n_bls_cycles_with_at_least_two_cadences":independent_with_two,
        "n_bls_event_cadences_min":min(events) if events else 0,
        "n_bls_event_cadences_max":max(events) if events else 0,
        "largest_event_fraction_in_transit_cadences":(
            max(events)/count if count else np.nan
        ),
        "n_TGLC_flagged_in_bls_windows":int(np.count_nonzero(qinside)),
        "fraction_TGLC_flagged_in_bls_windows":(
            float(np.count_nonzero(qinside)/count) if count else np.nan
        ),
        "observed_bls_window_depth_median":(
            float(np.median(out)-np.median(f[inside]))
            if count and len(out)>5 else np.nan
        ),
        "n_bls_events_with_positive_median_drop":int(
            np.count_nonzero(np.asarray(event_drops)>0)
        ),
        "smallest_event_median_drop":(
            min(event_drops) if event_drops else np.nan
        ),
        "largest_event_median_drop":(
            max(event_drops) if event_drops else np.nan
        ),
    }


def gap_proximity(t,period,epoch):
    """How close any predicted BLS box CENTER is to an actual >3h gap EDGE."""
    gaps=np.flatnonzero(np.diff(t)>3/24)
    if not len(gaps):
        return {"n_time_gaps_over_3h":0,
                "nearest_bls_midpoint_to_gap_edge_hours":np.nan}
    left=t[gaps]
    right=t[gaps+1]
    n0=int(np.floor((t[0]-epoch)/period))-1
    n1=int(np.ceil((t[-1]-epoch)/period))+1
    centers=epoch+np.arange(n0,n1+1)*period
    centers=centers[(centers>=t[0])&(centers<=t[-1])]
    if not len(centers):
        nearest=np.nan
    else:
        edges=np.concatenate([left,right])
        nearest=float(24*np.min(np.abs(centers[:,None]-edges[None,:])))
    return {
        "n_time_gaps_over_3h":int(len(gaps)),
        "nearest_bls_midpoint_to_gap_edge_hours":nearest,
    }


def catalogue_windows(t,period,epoch_jd,duration_hours,ts):
    """Count cadence support at frozen catalogue epochs; NEVER use for BLS.

    Missing duration or epoch => missing support, NOT zero. For ambiguous
    time standard "BJD", report provisional counts, NOT phase evidence.
    """
    if not np.isfinite(epoch_jd) or not np.isfinite(period) or period<=0:
        return {"n_frozen_catalogue_midpoints_in_sector":np.nan,
                "n_catalogue_windows_with_2plus_cadences":np.nan,
                "catalogue_window_time_standard":"unusable"}
    t0=epoch_jd-BTJD_ZERO
    a=int(math.ceil((t[0]-t0)/period))
    b=int(math.floor((t[-1]-t0)/period))
    if b-a>1000:
        raise ValueError("Too many predicted catalogue cycles")
    predicted=max(0,b-a+1)
    strict=str(ts).strip().upper().replace("_","-")=="BJD-TDB"
    time_status="explicit_BJD-TDB" if strict else "provisional_unresolved_BJD"
    if not np.isfinite(duration_hours) or duration_hours<=0:
        return {"n_frozen_catalogue_midpoints_in_sector":predicted,
                "n_catalogue_windows_with_2plus_cadences":np.nan,
                "catalogue_window_time_standard":time_status}
    residual,cycle=phase_and_cycle(t,period,t0)
    incl=(np.abs(residual)<=duration_hours/(24*2))
    within=(cycle>=a)&(cycle<=b)&incl
    _,nums=np.unique(cycle[within],return_counts=True)
    return {"n_frozen_catalogue_midpoints_in_sector":predicted,
            "n_catalogue_windows_with_2plus_cadences":int((nums>=2).sum()),
            "catalogue_window_time_standard":time_status}


def one_row(blind,selected,source,private):
    tic=str(blind["tic_id"])
    sector=int(blind["sector"])
    path=private/f"{tic}_sector_{sector:02d}_{POLICY}.csv"
    if not path.is_file():
        raise FileNotFoundError(f"Missing frozen curve: {path}")
    expected=Path(str(source["private_cleaned_file"]))
    if expected.resolve()!=path.resolve():
        raise ValueError(f"Photometry manifest filename mismatch for {tic}")
    if str(source["policy_id"])!=POLICY:
        raise ValueError("Frozen photometric policy mismatch")
    curve=pd.read_csv(path)
    need={"time_btjd","flux_cal_aper_normalized","tglc_flags","cadence_num"}
    if not need.issubset(curve):
        raise ValueError(f"Missing original policy curve columns: {need-set(curve)}")
    if int(source["n_accepted"])!=len(curve) or int(blind["n_points"])!=len(curve):
        raise ValueError("Accepted cadence count changed after experiment 17")
    t=curve.time_btjd.to_numpy(float)
    f=curve.flux_cal_aper_normalized.to_numpy(float)
    q=curve.tglc_flags.to_numpy(np.int64)
    if not np.all(np.diff(t)>0):
        raise ValueError("Frozen curve times not increasing")
    observed=count_event_support(
        t,f,q,float(blind["best_period_days"]),
        float(blind["best_transit_time_btjd"]),
        float(blind["best_duration_days"])
    )
    gaps=gap_proximity(
        t,float(blind["best_period_days"]),
        float(blind["best_transit_time_btjd"])
    )
    cat=catalogue_windows(
        t,float(selected["planet_period_days"]),
        numeric(selected.get("transit_epoch_jd")),
        numeric(selected.get("planet_transit_duration_hours")),
        selected.get("transit_epoch_time_standard","")
    )
    return {
        "tic_id":tic, "planet":str(blind["planet"]),
        "sector":sector, "period_stratum":str(blind["period_bin"]),
        "blind_period_days":float(blind["best_period_days"]),
        "frozen_nasa_planet_period_days":float(selected["planet_period_days"]),
        "direct_period_match_frozen":bool(blind["blind_direct_period_match_2pct"]),
        "blind_alias_match_frozen":bool(blind["blind_2to1_alias_match"]),
        "blind_period_fraction_of_time_baseline":(
            float(blind["best_period_days"])/float(blind["baseline_days"])
        ),
        "best_BLS_duration_hours":24*float(blind["best_duration_days"]),
        "best_BLS_depth_fraction":float(blind["best_depth_fraction"]),
        "best_BLS_SNR_uncalibrated":float(
            blind["best_depth_snr_white_noise_proxy"]
        ),
        "n_original_policy_cadences":len(t),
        "fraction_TGLC_flagged_all_accepted":float(
            np.count_nonzero(q)/len(q)
        ),
        "frozen_cleaned_curve_sha256":file_hash(path),
        **observed,**gaps,**cat,
    }


def plot_summary(frame,where):
    x=frame.frozen_nasa_planet_period_days.to_numpy(float)
    y=frame.blind_period_days.to_numpy(float)
    fig,ax=plt.subplots(figsize=(8,6.5),layout="constrained")
    for ok,part in frame.groupby("direct_period_match_frozen"):
        ax.scatter(part.frozen_nasa_planet_period_days,
                   part.blind_period_days,label=(
                       "Blind period matches NASA" if ok else
                       "Blind period does not match NASA"
                   ),s=38)
    ax.plot([0,15],[0,15],"k:",lw=1,label="Identity")
    ax.set(xlim=(0,15),ylim=(0,15),
           xlabel="Frozen NASA planet period (days)",
           ylabel="Blind BLS strongest period (days)",
           title="Post hoc BLS peak forensics (18 untouched hosts)")
    ax.legend(fontsize=8)
    fig.savefig(where/"blind_peak_vs_planet_period.png",dpi=150)
    plt.close(fig)

    fig,ax=plt.subplots(figsize=(8,5.5),layout="constrained")
    ax.scatter(frame.blind_period_days,
               frame.n_bls_cycles_with_at_least_two_cadences,
               s=42)
    ax.set(xlabel="BLS strongest period (days)",
           ylabel="Distinct supported box cycles (>=2 cadences)",
           title="How many independent transit-like events support BLS?")
    fig.savefig(where/"blind_event_support_vs_period.png",dpi=150)
    plt.close(fig)


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--private-dir",type=Path,default=PRIVATE)
    p.add_argument("--outdir",type=Path,default=OUT)
    args=p.parse_args()
    blind=pd.read_csv(BLIND,dtype={"tic_id":str})
    selected=pd.read_csv(SELECTED,dtype={"tic_id":str})
    sources=pd.read_csv(METRICS,dtype={"tic_id":str})
    if (
        len(blind)!=18 or len(selected)!=18 or len(sources)!=18
        or blind.tic_id.duplicated().any()
        or selected.tic_id.duplicated().any()
        or sources.tic_id.duplicated().any()
        or set(blind.tic_id)!=set(selected.tic_id)
        or set(blind.tic_id)!=set(sources.tic_id)
    ):
        raise ValueError("Original frozen 18-host cohort is incomplete or inconsistent")
    sel={r.tic_id:r._asdict() for r in selected.itertuples(index=False)}
    src={r.tic_id:r._asdict() for r in sources.itertuples(index=False)}
    data=[]
    for b in blind.sort_values("selection_rank").to_dict("records"):
        tic=b["tic_id"]
        r=one_row(b,sel[tic],src[tic],args.private_dir)
        data.append(r)
        print(
            f"[{len(data):02d}/18] {tic}: P_BLS={r['blind_period_days']:.4f}, "
            f"events_2+={r['n_bls_cycles_with_at_least_two_cadences']}, "
            f"largest_event_fraction="
            f"{r['largest_event_fraction_in_transit_cadences']:.2f}",
            flush=True
        )
    frame=pd.DataFrame(data)
    args.outdir.mkdir(parents=True,exist_ok=True)
    figures=args.outdir/"figures"
    figures.mkdir(parents=True,exist_ok=True)
    plot_summary(frame,figures)
    frame.to_csv(args.outdir/"per_host_blind_peak_support.csv",index=False)
    # Report DESCRIPTIVE counts only; do not create a new 'success rate'
    # based on post hoc event criteria or reclassify any 17 outcome.
    n_long=int(frame.blind_period_days.between(12,15,inclusive="both").sum())
    n_sparse=int((frame.n_bls_cycles_with_at_least_two_cadences<2).sum())
    n_deep=int((frame.best_BLS_depth_fraction>.05).sum())
    config={
        "design":"Exploratory, POST-HOC forensic audit of FROZEN BLS outputs",
        "not_a_new_detector_or_retuned_pipeline":True,
        "original_blind_direct_period_recoveries":int(
            frame.direct_period_match_frozen.sum()
        ),
        "original_blind_n_hosts":18,
        "original_scores_are_unchanged":True,
        "bls_period_range_high_days":[12,15],
        "n_BLS_peaks_12_to_15_days":n_long,
        "n_BLS_peaks_with_less_than_2_sampled_events":n_sparse,
        "n_box_depths_greater_than_5_percent":n_deep,
        "catalogue_BJD_without_TDB_unresolved":True,
        "science_goal":"diagnose long-period/gap/single-event artifacts",
        "frozen_upstream_blind_csv_sha256":file_hash(BLIND),
        "frozen_upstream_selection_sha256":file_hash(SELECTED),
        "photometric_policy":POLICY,
        "must_not_train_or_tune_using_these_18":True,
        "notes":[
            "Large BLS SNR is white-noise proxy, not calibrated significance.",
            "Single supported BLS transit event is not periodic confirmation.",
            "A BLS box depth exceeding 5% is a descriptive flag, not exclusion.",
            "BJD unqualified time system is not a valid precise timing claim.",
            "All 18 hosts remain in the frozen experiment-17 denominator.",
            "Future model improvements require separate new development TICs."
        ]
    }
    (args.outdir/"forensic_config.json").write_text(
        json.dumps(config,indent=2),encoding="utf-8"
    )
    columns=[
        "tic_id","planet","period_stratum","blind_period_days",
        "frozen_nasa_planet_period_days","direct_period_match_frozen",
        "best_BLS_depth_fraction","best_BLS_SNR_uncalibrated",
        "n_bls_cycles_with_at_least_two_cadences",
        "largest_event_fraction_in_transit_cadences",
        "nearest_bls_midpoint_to_gap_edge_hours",
        "fraction_TGLC_flagged_in_bls_windows",
        "n_frozen_catalogue_midpoints_in_sector",
        "n_catalogue_windows_with_2plus_cadences",
        "catalogue_window_time_standard"
    ]
    report=[
        "# Post-hoc support for BLIND BLS transit windows", "",
        "This is a diagnostic, NOT a revised detection benchmark.",
        "Original experiment-17 direct period concordance retained: "
            f"{int(frame.direct_period_match_frozen.sum())}/18.",
        f"Strongest BLS peaks at 12-15 days: {n_long}/18.",
        f"Peaks supported by fewer than two independently sampled events: {n_sparse}/18.",
        f"Peak box depths >5% (descriptive, not an exclusion): {n_deep}/18.",
        "", "## Frozen per-target results plus diagnostics", "",
        frame[columns].round(5).to_string(index=False),
        "", "## Strict interpretation", "",
        *["- "+x for x in config["notes"]],
    ]
    (args.outdir/"forensic_report.md").write_text(
        "\n".join(report)+"\n",encoding="utf-8"
    )
    print("\n"+"\n".join(report[:8]),flush=True)


if __name__=="__main__":
    main()
