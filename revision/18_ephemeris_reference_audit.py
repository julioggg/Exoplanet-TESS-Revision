#!/usr/bin/env python3
"""Post-scan ephemeris PROVENANCE and time-reference audit for 18 TESS hosts.

This script never feeds catalogue periods/epochs into blind BLS. It reads
the complete, frozen outputs of revision/17, and then obtains published
period uncertainties / transit timing variation flags from NASA PS.

IMPORTANT limitations:
- The period and transit epochs used for comparison are frozen in the
  experiment-15 source selection, NOT replaced by a new PS solution.
- New NASA query is supplementary metadata; a changed default solution is
  flagged, not silently substituted.
- Strict phase comparison is permissible ONLY for explicitly 'BJD-TDB'
  timestamps, complete P/T0 uncertainties, TTV flag=0, a plausible
  ephemeris extrapolation, and a reported transit duration.
- Sigma_pred = sqrt(sigma_T0^2 + N_cycles^2*sigma_P^2) omits covariance
  and potential timing-systematics; it is a rough lower-bound proxy, not
  a fully propagated authoritative transit ephemeris confidence interval.
- Phase concordance from a BLS peak is descriptive, NOT planet validation.

Run (after experiment 17):
  python -u revision/18_ephemeris_reference_audit.py
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import math
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
SELECTED=ROOT/"results/new_confirmed_tess_inventory/selected_prospective_hosts.csv"
BLIND=ROOT/"results/new_confirmed_tess_blind_period/blind_period_per_host.csv"
PRIVATE=Path.home()/"exoplanet_revision/data/new_confirmed_tess_cleaned"
OUT=ROOT/"results/new_confirmed_tess_ephemeris_audit"
TAP="https://exoplanetarchive.ipac.caltech.edu/TAP/sync"
TESS_BTJD_ZERO=2457000.0
LIVE_FIELDS=[
    "pl_name","pl_orbper","pl_orbpererr1","pl_orbpererr2",
    "pl_tranmid","pl_tranmiderr1","pl_tranmiderr2","pl_tsystemref",
    "pl_trandur","ttv_flag","default_flag"
]


def hash_bytes(value):
    return hashlib.sha256(value).hexdigest()


def live_tap(names):
    allowed = [
        str(name).replace("'","''") for name in sorted(set(names))
    ]
    quoted=",".join("'" +v+"'" for v in allowed)
    query=(
        "select "+",".join(LIVE_FIELDS)+" from ps where "
        "default_flag=1 and pl_name in ("+quoted+")"
    )
    args=urllib.parse.urlencode({"query":query,"format":"csv"})
    endpoint=TAP+"?"+args
    request=urllib.request.Request(
        endpoint,headers={"User-Agent":"TESS-Ephemeris-Provenance-Audit/1.0"}
    )
    with urllib.request.urlopen(request,timeout=90) as response:
        content=response.read(4*1024*1024+1)
    if len(content)>4*1024*1024:
        raise ValueError("Unexpectedly large NASA PS ephemeris response")
    frame=pd.read_csv(io.BytesIO(content))
    if frame.empty or not set(LIVE_FIELDS).issubset(frame):
        raise ValueError(
            "NASA PS query returned no rows or incorrect schema: "
            +str(frame.columns.tolist())
        )
    if frame.pl_name.duplicated().any():
        raise ValueError("Duplicated NASA default solution for a planet")
    return frame,endpoint,hash_bytes(content)


def errorbar(row,pos_name,neg_name):
    a=pd.to_numeric(row.get(pos_name),errors="coerce")
    b=pd.to_numeric(row.get(neg_name),errors="coerce")
    if not np.isfinite(a) or not np.isfinite(b):
        return float("nan")
    return float(max(abs(a),abs(b)))


def phase_offset(a,b,p):
    """Signed residual modulo true period, in days."""
    return float((a-b+.5*p)%p-.5*p)


def audit_one(selected,blind,live,curve):
    period=float(selected["planet_period_days"])
    if not np.isfinite(period) or period<=0:
        raise ValueError("Invalid frozen planet period")
    jd_epoch=float(selected["transit_epoch_jd"])
    ts=str(selected.get("transit_epoch_time_standard","")).strip().upper()
    t=np.asarray(curve.time_btjd,dtype=float)
    if len(t)<300 or not np.isfinite(t).all() or not (np.diff(t)>0).all():
        raise ValueError("Invalid frozen aperture time series")
    obs_first=float(t[0]+TESS_BTJD_ZERO)
    obs_last=float(t[-1]+TESS_BTJD_ZERO)
    obs_mid=float(np.median(t)+TESS_BTJD_ZERO)
    n=int(round((obs_mid-jd_epoch)/period)) if np.isfinite(jd_epoch) else 0
    closest=jd_epoch+n*period if np.isfinite(jd_epoch) else np.nan

    # Epoch error was frozen in the original, verified NASA catalogue snapshot.
    sigma_epoch=errorbar(
        selected,"transit_epoch_error1","transit_epoch_error2"
    )
    sigma_period=errorbar(live,"pl_orbpererr1","pl_orbpererr2")
    sigma_pred=(
        float(np.hypot(sigma_epoch,abs(n)*sigma_period))
        if np.isfinite(sigma_epoch) and np.isfinite(sigma_period)
        else np.nan
    )
    period_current=pd.to_numeric(live.get("pl_orbper"),errors="coerce")
    delta_p=(
        abs(float(period_current)-period)/period
        if np.isfinite(period_current) else np.nan
    )
    mismatch=not np.isfinite(delta_p) or delta_p>5e-5
    ttv=pd.to_numeric(live.get("ttv_flag"),errors="coerce")
    duration_h=pd.to_numeric(
        selected.get("planet_transit_duration_hours"),errors="coerce"
    )
    duration_d=float(duration_h)/24 if np.isfinite(duration_h) and duration_h>0 else np.nan
    n_transits=0
    n_sampled=0
    if np.isfinite(jd_epoch):
        first=math.ceil((obs_first-jd_epoch)/period)
        last=math.floor((obs_last-jd_epoch)/period)
        n_transits=max(0,last-first+1)
        if n_transits>1000:
            raise ValueError("Implausibly many predicted transit epochs")
        if n_transits and np.isfinite(duration_d):
            for cycle in range(first,last+1):
                tm=jd_epoch+cycle*period-TESS_BTJD_ZERO
                if np.any(np.abs(t-tm)<=duration_d/2):
                    n_sampled+=1
    is_tdb=ts.replace("_","-")=="BJD-TDB"
    timing_scoreable=bool(
        is_tdb and not mismatch and np.isfinite(jd_epoch)
        and np.isfinite(sigma_pred) and np.isfinite(duration_d)
        and sigma_pred<=duration_d/2 and ttv==0
        and n_transits>0 and n_sampled>0
    )
    # The blind BLS best transit epoch is a model output, not a NASA input.
    best_bls_time=float(blind["best_transit_time_btjd"])+TESS_BTJD_ZERO
    offset=(
        phase_offset(best_bls_time,closest,period)
        if np.isfinite(jd_epoch) else np.nan
    )
    period_ok=str(blind["blind_direct_period_match_2pct"]).lower()=="true"
    match_timing=(
        bool(period_ok and abs(offset)<=duration_d/2+3*sigma_pred)
        if timing_scoreable else None
    )
    failure_reasons=[]
    if not is_tdb:
        failure_reasons.append("time_system_not_explicit_BJD-TDB")
    if mismatch:
        failure_reasons.append("current_NASA_period_differs_from_frozen_period")
    if not np.isfinite(jd_epoch):
        failure_reasons.append("no_frozen_epoch")
    if not np.isfinite(sigma_epoch):
        failure_reasons.append("no_frozen_epoch_uncertainty")
    if not np.isfinite(sigma_period):
        failure_reasons.append("no_period_uncertainty")
    if not np.isfinite(duration_d):
        failure_reasons.append("no_transit_duration")
    if ttv!=0:
        failure_reasons.append("ttv_or_unknown_ttv_flag")
    if n_transits==0:
        failure_reasons.append("no_predicted_midpoint_in_selected_sector")
    if n_sampled==0:
        failure_reasons.append("no_cadences_close_to_predicted_midpoint")
    if np.isfinite(sigma_pred) and np.isfinite(duration_d):
        if sigma_pred>duration_d/2:
            failure_reasons.append("ephemeris_uncertainty_exceeds_half_duration")
    return {
        "tic_id":str(selected["tic_id"]),
        "planet":str(selected["planet"]),
        "sector":int(selected["sector"]),
        "known_period_frozen_days":period,
        "known_period_current_days":float(period_current)
            if np.isfinite(period_current) else np.nan,
        "relative_archive_period_change":delta_p,
        "frozen_epoch_jd":jd_epoch,
        "frozen_time_reference":ts,
        "n_cycles_to_observed_median":n,
        "sigma_epoch_frozen_days":sigma_epoch,
        "sigma_period_live_days":sigma_period,
        "sigma_propagated_days_without_covariance":sigma_pred,
        "transit_duration_hours":float(duration_h)
            if np.isfinite(duration_h) else np.nan,
        "ttv_flag_live":ttv,
        "n_predicted_transit_midpoints_in_sector":n_transits,
        "n_midpoints_with_cadence_within_half_duration":n_sampled,
        "blind_BLS_period_days":float(blind["best_period_days"]),
        "blind_direct_period_match":bool(period_ok),
        "blind_BLS_epoch_BTJD":float(blind["best_transit_time_btjd"]),
        "blind_BLS_vs_ephemeris_phase_offset_hours":offset*24,
        "strict_timing_reference_and_uncertainty_usable":timing_scoreable,
        "blind_period_and_ephemeris_consistent":match_timing,
        "timing_not_scoreable_reasons":";".join(failure_reasons)
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--private-dir",type=Path,default=PRIVATE)
    ap.add_argument("--outdir",type=Path,default=OUT)
    ap.add_argument("--ps-live-csv",type=Path,default=None,
                    help="Use a saved NASA PS ephemeris metadata snapshot")
    args=ap.parse_args()
    args.outdir.mkdir(parents=True,exist_ok=True)

    selected=pd.read_csv(SELECTED,dtype={"tic_id":str})
    blind=pd.read_csv(BLIND,dtype={"tic_id":str})
    if (
        len(selected)!=18 or len(blind)!=18
        or set(selected.tic_id)!=set(blind.tic_id)
        or selected.tic_id.duplicated().any()
        or blind.tic_id.duplicated().any()
    ):
        raise ValueError("Incomplete or inconsistent frozen 18-TIC benchmark")
    names=selected.planet.astype(str).tolist()
    if args.ps_live_csv:
        live=pd.read_csv(args.ps_live_csv)
        nasa_url="offline_supplied_snapshot"
        nasa_hash=hash_bytes(args.ps_live_csv.read_bytes())
    else:
        live,nasa_url,nasa_hash=live_tap(names)
    live.to_csv(args.outdir/"nasa_ps_ephemeris_metadata_snapshot.csv",
                index=False)
    if live.pl_name.duplicated().any():
        raise ValueError("NASA PS returned duplicate planets")

    rows,failures=[],[]
    for entry in selected.to_dict("records"):
        name=str(entry["planet"])
        tic=str(entry["tic_id"])
        try:
            bl=blind.loc[blind.tic_id==tic].iloc[0].to_dict()
            ref=live.loc[live.pl_name==name]
            if len(ref)!=1:
                raise ValueError(
                    f"No unique NASA PS updated ephemeris for {name}"
                )
            path=args.private_dir/f"{tic}_sector_{int(entry['sector']):02d}.csv"
            if not path.is_file():
                raise FileNotFoundError(path)
            curve=pd.read_csv(path)
            result=audit_one(entry,bl,ref.iloc[0].to_dict(),curve)
            rows.append(result)
            print(
                f"[{len(rows):02d}/18] {tic}: "
                f"cycles={result['n_cycles_to_observed_median']}, "
                f"sigma_pred={result['sigma_propagated_days_without_covariance']:.5g} d, "
                f"phase_scoreable={result['strict_timing_reference_and_uncertainty_usable']}",
                flush=True
            )
        except Exception as exc:
            failures.append({"tic_id":tic,"planet":name,"error":repr(exc)})
            print(f"[FAIL] {tic}: {exc}",flush=True)
    df=pd.DataFrame(rows)
    df.to_csv(args.outdir/"per_host_ephemeris_audit.csv",index=False)
    pd.DataFrame(failures,columns=[
        "tic_id","planet","error"
    ]).to_csv(args.outdir/"ephemeris_errors.csv",index=False)
    n_eligible=(
        int(df.strict_timing_reference_and_uncertainty_usable.sum())
        if not df.empty else 0
    )
    n_joint=(
        int(df.blind_period_and_ephemeris_consistent.fillna(False).sum())
        if not df.empty else 0
    )
    cfg={
        "snapshot_source":"NASA Exoplanet Archive PS current default rows",
        "metadata_query_url":nasa_url,
        "nasa_ephemeris_snapshot_sha256":nasa_hash,
        "frozen_catalogue_period_and_epoch_NOT_updated":True,
        "BJD_to_BTJD_offset":TESS_BTJD_ZERO,
        "phase_scoring_restricted_to_explicit_BJD-TDB":True,
        "sigma_prediction_ignores_T0_P_covariance_and_TTV":True,
        "timing_score_rule":"only direct BLS period and abs(phase offset)<=duration/2+3*sigma",
        "timing_uncertainty_eligibility":"sigma_pred <= duration/2; no TTV; unchanged period",
        "in_sector_predicted_transits_not_necessarily_observed":True,
        "selection_or_BLS_tuning_performed":False,
        "n_host_sources_expected":18,
        "n_host_sources_audited":len(rows),
        "n_host_sources_failed":len(failures),
        "n_with_strictly_scoreable_timing":n_eligible,
        "n_with_strict_scoreable_BLS_period_and_epoch_match":n_joint,
        "not_a_formal_detector_true_positive_rate":True
    }
    (args.outdir/"ephemeris_config.json").write_text(
        json.dumps(cfg,indent=2),encoding="utf-8"
    )
    report=[
        "# NASA confirmed-host ephemeris provenance / prediction audit", "",
        f"Sources audited: {len(df)}/18; failures: {len(failures)}",
        f"Timing scoreable (strict time/uncertainty criteria): {n_eligible}/18",
        f"Blind BLS period and timing jointly consistent: {n_joint}/{n_eligible} scoreable",
        "", "## Source-by-source diagnostics","",
        df.round(6).to_string(index=False) if not df.empty else "(none)",
        "", "## Interpretation limits", "",
        "- NASA default periods and epochs from cohort selection remain frozen.",
        "- Live NASA PS metadata only contributes period uncertainty and TTV flags.",
        "- BJD without an explicit TDB standard is not treated as BJD-TDB.",
        "- Propagated timing errors assume independent T0 and period, omit covariance.",
        "- Epoch propagation over hundreds of cycles may be inconclusive.",
        "- A favorable BLS match is NOT proof of a transit of the catalogued planet.",
        "- This is an exploratory post-scan timing audit, not calibrated planet detection."
    ]
    (args.outdir/"ephemeris_report.md").write_text(
        "\n".join(report)+"\n",encoding="utf-8"
    )
    print("\n"+"\n".join(report[:6]),flush=True)
    if failures:
        raise RuntimeError(
            f"{len(failures)} ephemeris-audit failures: "
            "inspect ephemeris_errors.csv, do not substitute targets"
        )


if __name__=="__main__":
    main()
