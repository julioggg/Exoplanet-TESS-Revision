#!/usr/bin/env python3
"""LOCKED blind BLS period-recovery benchmark: 18 NEW NASA-confirmed TESS hosts.

Before this new cohort was downloaded, a previous injection experiment fixed:
  period grid 0.5..15 days, 4000 points;
  allowed transit durations 0.05, 0.10, 0.15, 0.20, 0.25 days;
  period agreement |P_BLS/P_archive - 1| <= 0.02;
  harmonics 0.5x and 2x evaluated separately.

Data preparation is frozen as TGLC_APERTURE_1800S_V1. BLS receives only
time and normalized calibrated aperture flux, never true period/epoch/planet.
Period matching occurs strictly AFTER scanning each source.

Descriptive cohort sensitivity ONLY. A confirmed planet HOST may not yield
recoverable in-sector transits. A period match ALONE is NOT exoplanet
confirmation and may happen by chance. No negatives here => NO specificity,
false positive rates, AUROC or calibrated detection threshold.

Inputs:
  results/new_confirmed_tess_inventory/selected_prospective_hosts.csv
  results/new_confirmed_tess_acquisition/verified_original_fits.csv
  results/new_confirmed_tess_photometry/aperture_curve_metrics.csv
  ~/exoplanet_revision/data/new_confirmed_tess_cleaned/*.csv
Output:
  results/new_confirmed_tess_blind_period/

Use:
  python -u revision/17_blind_confirmed_host_period_recovery.py
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
INVENTORY=ROOT/"results/new_confirmed_tess_inventory/selected_prospective_hosts.csv"
PROVENANCE=ROOT/"results/new_confirmed_tess_acquisition/verified_original_fits.csv"
METRICS=ROOT/"results/new_confirmed_tess_photometry/aperture_curve_metrics.csv"
PRIVATE=Path.home()/"exoplanet_revision/data/new_confirmed_tess_cleaned"
OUT=ROOT/"results/new_confirmed_tess_blind_period"
PERIODS=np.linspace(0.5,15.0,4000)
DURATIONS=np.array([0.05,0.10,0.15,0.20,0.25],dtype=float)
DIRECT_TOL=0.02
ALIAS_TOL_HALF=0.01
ALIAS_TOL_TWICE=0.04
MIN_POINTS=300


def sha256(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):
            h.update(b)
    return h.hexdigest()


def mad_sigma(x):
    arr=np.asarray(x,dtype=float)
    return float(1.4826*np.median(np.abs(arr-np.median(arr))))


def binomial_exact(k,n,alpha=.05):
    if not 0<=k<=n or n<=0:
        raise ValueError("Invalid binomial parameters")
    lo=float(beta.ppf(alpha/2,k,n-k+1)) if k>0 else 0.0
    hi=float(beta.ppf(1-alpha/2,k+1,n-k)) if k<n else 1.0
    return lo,hi


def period_match(best,true):
    r=float(best)/float(true)
    direct=abs(r-1.0)<=DIRECT_TOL
    alias=(abs(r-.5)<=ALIAS_TOL_HALF or abs(r-2.0)<=ALIAS_TOL_TWICE)
    return bool(direct),bool(alias)


def blind_scan(t,flux):
    """No catalogue ephemeris/period in signature or BLS search."""
    t=np.asarray(t,dtype=float)
    f=np.asarray(flux,dtype=float)
    if len(t)<MIN_POINTS:
        raise ValueError("Insufficient accepted cadences")
    if not np.isfinite(t).all() or not np.isfinite(f).all():
        raise ValueError("Nonfinite input to BLS")
    if not np.all(np.diff(t)>0):
        raise ValueError("BTJD timestamps are not strictly increasing")
    noise=mad_sigma(np.diff(f))/np.sqrt(2)
    if not np.isfinite(noise) or noise<=0:
        raise ValueError(
            "Point-to-point MAD noise is zero or nonfinite: no artificial "
            "error floor applied to an originally frozen light curve"
        )
    bls=BoxLeastSquares(t,f,dy=np.full(len(t),noise))
    scan=bls.power(PERIODS,DURATIONS)
    power=np.asarray(scan.power,dtype=float)
    if not np.isfinite(power).any():
        raise ValueError("BLS scan contains no finite powers")
    ix=int(np.nanargmax(power))
    return {
        "best_period_days":float(scan.period[ix]),
        "best_duration_days":float(scan.duration[ix]),
        "best_transit_time_btjd":float(scan.transit_time[ix]),
        "best_depth_fraction":float(scan.depth[ix]),
        "best_depth_snr_white_noise_proxy":float(scan.depth_snr[ix]),
        "best_bls_power_white_noise_proxy":float(power[ix]),
        "point_to_point_white_noise_proxy":float(noise),
        "n_points":int(len(t)),
        "baseline_days":float(t[-1]-t[0]),
        "median_cadence_min":float(np.median(np.diff(t))*1440),
    },power


def load_curve(tic,sector,private):
    path=private/f"{tic}_sector_{int(sector):02d}.csv"
    if not path.is_file():
        raise FileNotFoundError(
            f"Original policy-14 cleaned file missing: {path}"
        )
    data=pd.read_csv(path)
    required={
        "time_btjd","flux_cal_aper_normalized",
        "tess_flags","tglc_flags","cadence_num"
    }
    if not required.issubset(data):
        raise ValueError(
            f"Wrong uniform-aperture CSV layout for {tic}: "
            +str(required-set(data.columns))
        )
    if data.cadence_num.duplicated().any():
        raise ValueError("Repeated cadence identifiers")
    if not np.isfinite(data[[
        "time_btjd","flux_cal_aper_normalized"
    ]].to_numpy(float)).all():
        raise ValueError("Nonfinite cleaned curve")
    return data


def verify_cohort(sel,manifest,metrics):
    if any(x.empty for x in [sel,manifest,metrics]):
        raise ValueError("Empty source manifest")
    for name,df in [
        ("frozen selection",sel),("acquisition",manifest),("processing",metrics)
    ]:
        if df.tic_id.duplicated().any():
            raise ValueError(f"Duplicate TIC: {name}")
    if len(sel)!=18:
        raise ValueError(
            "Frozen benchmark must contain EXACTLY 18 selected new TICs"
        )
    if set(sel.tic_id)!=set(manifest.tic_id)!=set(metrics.tic_id):
        raise ValueError("Cohort not fully acquired and uniformly processed")
    if not (metrics.policy_id=="TGLC_APERTURE_1800S_V1").all():
        raise ValueError("Incorrect predeclared photometric policy")
    if not (manifest.status=="verified_original_fits").all():
        raise ValueError("Unverified MAST original FITS")
    for field in ["sector","author","exptime_seconds","mast_uri"]:
        if field in sel and field in manifest:
            compared=sel[["tic_id",field]].merge(
                manifest[["tic_id",field]],on="tic_id",
                suffixes=("_selected","_download"),validate="one_to_one"
            )
            a=compared[field+"_selected"]
            b=compared[field+"_download"]
            if field=="exptime_seconds":
                if not np.allclose(a.to_numpy(float),b.to_numpy(float)):
                    raise ValueError(f"Frozen {field} changed")
            elif not a.astype(str).equals(b.astype(str)):
                raise ValueError(f"Frozen {field} changed")


def plot_one(tic,sector,data,fit,power,figures):
    t=data.time_btjd.to_numpy(float)
    y=data.flux_cal_aper_normalized.to_numpy(float)
    per=fit["best_period_days"]
    epoch=fit["best_transit_time_btjd"]
    phase=((t-epoch+.5*per)%per)/per-.5
    fig,axs=plt.subplots(3,1,figsize=(10.5,8.5),layout="constrained")
    axs[0].plot(t-t[0],y,".",markersize=2.2)
    axs[0].set(
        title=f"{tic}, TESS sector {sector} | blind BLS",
        ylabel="Normalized aperture flux",
        xlabel="Days since first accepted exposure"
    )
    axs[1].plot(PERIODS,power,lw=.8)
    axs[1].axvline(per,ls="--",lw=1,color="black")
    axs[1].set(xlabel="Blind searched period (days)",
               ylabel="BLS power (white-noise proxy)")
    axs[2].plot(phase,y,".",markersize=2.4)
    axs[2].set(xlabel="Phase at BLIND best period (cycles)",
               ylabel="Normalized aperture flux")
    fig.savefig(
        figures/f"{tic}_sector_{int(sector):02d}_blind_bls.png",dpi=140
    )
    plt.close(fig)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--private-dir",type=Path,default=PRIVATE)
    ap.add_argument("--outdir",type=Path,default=OUT)
    args=ap.parse_args()
    args.outdir.mkdir(parents=True,exist_ok=True)
    figures=args.outdir/"figures"
    figures.mkdir(exist_ok=True)

    selection=pd.read_csv(INVENTORY,dtype={"tic_id":str})
    manifest=pd.read_csv(PROVENANCE,dtype={"tic_id":str})
    prior_metrics=pd.read_csv(METRICS,dtype={"tic_id":str})
    verify_cohort(selection,manifest,prior_metrics)
    selection=selection.sort_values("selection_rank")
    per_target=[]
    failures=[]
    for item in selection.to_dict("records"):
        tic=str(item["tic_id"])
        sector=int(item["sector"])
        try:
            curve=load_curve(tic,sector,args.private_dir)
            # SCIENTIFIC SEPARATION: scan first, compare to NASA second.
            fit,power=blind_scan(
                curve.time_btjd.to_numpy(float),
                curve.flux_cal_aper_normalized.to_numpy(float)
            )
            measured_n=len(curve)
            row_qa=prior_metrics.loc[prior_metrics.tic_id==tic].iloc[0]
            if measured_n!=int(row_qa.n_accepted):
                raise ValueError("Inconsistent policy 14 accepted-cadence count")
            if abs(fit["baseline_days"]-float(row_qa.baseline_days))>1e-6:
                raise ValueError("Different time baseline from frozen curve QA")
            reference=float(item["planet_period_days"])
            if not 1.5<=reference<=9.5:
                raise ValueError("NASA reference period outside prespecified window")
            match,alias=period_match(fit["best_period_days"],reference)
            result={
                "selection_rank":int(item["selection_rank"]),
                "tic_id":tic,"planet":str(item["planet"]),
                "sector":sector,
                "period_truth_days":reference,
                "period_bin":str(item["period_bin"]),
                "blind_direct_period_match_2pct":match,
                "blind_2to1_alias_match":alias,
                "blind_any_direct_or_2to1_alias":match or alias,
                "relative_period_error":float(
                    (fit["best_period_days"]-reference)/reference
                ),
                "known_transit_epoch_jd":item.get("transit_epoch_jd"),
                "known_epoch_time_standard":item.get(
                    "transit_epoch_time_standard"
                ),
                "n_TGLC_flagged_kept":int(
                    np.count_nonzero(curve.tglc_flags.to_numpy())
                ),
                "raw_FITS_sha256":str(
                    manifest.loc[manifest.tic_id==tic,"sha256"].iloc[0]
                ),
                **fit
            }
            per_target.append(result)
            plot_one(tic,sector,curve,fit,power,figures)
            print(
                f"[{len(per_target):02d}/18] {tic} {item['planet']}: "
                f"P_BLS={fit['best_period_days']:.5f}d, "
                f"P_NASA={reference:.5f}d, "
                f"direct={match}, alias={alias}",
                flush=True,
            )
        except Exception as exc:
            failures.append({
                "tic_id":tic,"sector":sector,"error":repr(exc)
            })
            print(f"[FAIL] {tic}: {exc}",flush=True)

    df=pd.DataFrame(per_target)
    df.to_csv(args.outdir/"blind_period_per_host.csv",index=False)
    pd.DataFrame(failures,columns=[
        "tic_id","sector","error"
    ]).to_csv(args.outdir/"benchmark_errors.csv",index=False)
    tables=[]
    if len(df):
        for metric in [
            "blind_direct_period_match_2pct",
            "blind_2to1_alias_match",
            "blind_any_direct_or_2to1_alias"
        ]:
            k=int(df[metric].astype(bool).sum())
            n=len(df)
            lo,hi=binomial_exact(k,n)
            tables.append({
                "metric":metric,"n_success":k,"n_hosts_evaluated":n,
                "rate":float(k/n),
                "exact_binomial_ci95_low":lo,
                "exact_binomial_ci95_high":hi
            })
    results=pd.DataFrame(tables)
    results.to_csv(args.outdir/"blind_period_aggregate.csv",index=False)

    by_period=[]
    if len(df):
        for name,g in df.groupby("period_bin",sort=True):
            by_period.append({
                "period_bin":name,"n_hosts":len(g),
                "direct_period_matches":int(g.blind_direct_period_match_2pct.sum()),
                "direct_match_rate":float(
                    g.blind_direct_period_match_2pct.mean()
                ),
                "alias_matches":int(g.blind_2to1_alias_match.sum())
            })
    pd.DataFrame(by_period).to_csv(
        args.outdir/"period_stratum_breakdown.csv",index=False
    )

    protocol={
        "design":"FROZEN blind BLS period recovery on 18 new confirmed hosts",
        "source_cohort":"experiment 15, exactly 18 new distinct TICs",
        "frozen_selection_sha256":sha256(INVENTORY),
        "verified_fits_manifest_sha256":sha256(PROVENANCE),
        "uniform_qa_csv_sha256":sha256(METRICS),
        "photometric_policy":"TGLC_APERTURE_1800S_V1",
        "BLS_input_truth_or_epoch":False,
        "blind_period_grid":[float(PERIODS[0]),float(PERIODS[-1]),len(PERIODS)],
        "blind_duration_days":DURATIONS.tolist(),
        "BLS_noise":"MAD(diff(normalized flux))*1.4826/sqrt(2), a white-noise surrogate",
        "primary_metric":"direct period within 2% of NASA default solution",
        "additional_metric":"half- and double-period alias, separately",
        "phase_timing_scored":False,
        "classification_AUC_or_specificity_scored":False,
        "no_BLS_threshold_applied":True,
        "selection_after_inspecting_period_results":False,
        "n_source_TICs_expected":18,
        "n_source_TICs_scanned":len(df),
        "n_source_TICs_failed":len(failures),
        "cautions":[
            "This is positive-host-only period concordance, not planet validation.",
            "A close orbital-period match may arise from correlated noise or chance.",
            "An in-sector transit is not established by a NASA planet-host flag.",
            "No statistically calibrated false alarm rate can be derived here.",
            "The NASA archive snapshot was frozen before downloading these curves.",
            "Timing reference and propagated ephemeris uncertainties not yet audited.",
            "New independent samples required for future detector tuning."
        ]
    }
    (args.outdir/"benchmark_protocol.json").write_text(
        json.dumps(protocol,indent=2),encoding="utf-8"
    )
    report=[
        "# Blind BLS period recovery on previously unused NASA-confirmed TESS hosts","",
        "## Frozen design",
        "Photometry: original TGLC calibrated aperture, fixed 1800-second exposure.",
        "BLS scanning was completed WITHOUT using the NASA planet period or epoch.",
        f"New unique source stars: {len(df)}/18; failures: {len(failures)}",
        "", "## Cohort-level period concordance (Clopper-Pearson 95% CI)","",
        results.round(6).to_string(index=False) if len(results) else "(none)",
        "", "## By catalogue orbital period stratum", "",
        pd.DataFrame(by_period).to_string(index=False) if by_period else "(none)",
        "", "## Per-source BLIND searches", "",
        df[[
            "tic_id","planet","sector","period_truth_days",
            "best_period_days","relative_period_error",
            "blind_direct_period_match_2pct","blind_2to1_alias_match",
            "best_depth_snr_white_noise_proxy"
        ]].round(6).to_string(index=False) if len(df) else "(none)",
        "", "## Interpretation limitations","",
        *["- "+s for s in protocol["cautions"]]
    ]
    (args.outdir/"blind_period_report.md").write_text(
        "\n".join(report)+"\n",encoding="utf-8"
    )
    print("\n"+"\n".join(report[:14]),flush=True)
    if failures:
        raise RuntimeError(
            f"Blind period benchmark had {len(failures)} errors; "
            "read benchmark_errors.csv. Do NOT replace source TICs."
        )


if __name__=="__main__":
    main()
