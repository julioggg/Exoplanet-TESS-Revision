#!/usr/bin/env python3
"""Prospective inventory: NEW NASA-confirmed transiting TESS planet hosts.

Uses NASA Exoplanet Archive PS/default solution, independently of the old
369 source TICs. Selection and source filtering are specified BEFORE any
light-curve extraction, classification or detection. It checks MAST TGLC
1800-second ORIGINAL product availability in primary sectors 1-26, but DOES
NOT download photometry.

This creates POSITIVE *planet-host* candidates only. A confirmed planet host
does not guarantee a transit in any particular TESS sector. Absence of a known
planet does NOT establish a negative/nontransiting class. Later blind period
recovery on new light curves is a different endpoint from planet validation.

Run:
  python -u revision/15_new_confirmed_tess_cohort.py \
    --target-stars 18 --max-mast-queries 70

Files in results/new_confirmed_tess_inventory:
  archive_snapshot.csv                  NASA TAP rows (new+old, for audit)
  candidate_queue.csv                   fixed-seed prespecified search order
  available_original_mast_products.csv  exact TIC+TGLC+1800+sector results
  selected_prospective_hosts.csv        18 or fewer NEW distinct source TICs
  mast_query_audit.csv                   all searched/rejected candidates
  discovery_report.md, discovery_config.json

Credentials: none. Network: NASA TAP + MAST metadata only. No FITS downloads.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import re
import time
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/new_confirmed_tess_inventory"
PRIOR = ROOT / "data/dataset/final_dataset.csv"
SEED = 20261027
RANGES = [1.5, 3.5, 6.5, 9.5]
PERIOD_LABELS = ("1.5-3.5", "3.5-6.5", "6.5-9.5")
BASE_URL = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"
QUERY = """select pl_name,hostname,tic_id,pl_orbper,pl_tranmid,
 pl_tranmiderr1,pl_tranmiderr2,pl_tsystemref,pl_ratror,pl_trandep,
 pl_trandur,disc_facility,discoverymethod,tran_flag,default_flag,
 pl_controv_flag from ps
 where default_flag=1 and tran_flag=1 and tic_id is not null
 and pl_orbper>=1.5 and pl_orbper<=9.5
 and disc_facility='Transiting Exoplanet Survey Satellite (TESS)'"""


def normalized_tic(value):
    if value is None or pd.isna(value):
        return ""
    text = str(value).strip()
    if re.fullmatch(r"\d+(?:\.0+)?", text):
        numeric = text.split(".")[0]
        return "TIC_" + numeric
    digits = re.fullmatch(r"(?i)TIC[\s_-]*(\d+)", text)
    return "TIC_" + digits.group(1) if digits else ""


def download_archive_csv():
    params = urllib.parse.urlencode({"query": QUERY, "format": "csv"})
    url = BASE_URL + "?" + params
    req = urllib.request.Request(
        url, headers={"User-Agent":"TESS-Original-Photometry-Replication/1.0"}
    )
    with urllib.request.urlopen(req, timeout=90) as resp:
        payload = resp.read(20*1024*1024 + 1)
    if len(payload) > 20*1024*1024:
        raise ValueError("Unexpectedly large archive TAP result")
    frame = pd.read_csv(io.BytesIO(payload), dtype={"tic_id": str})
    expected = {
        "pl_name", "hostname", "tic_id", "pl_orbper",
        "pl_tranmid", "disc_facility", "default_flag", "tran_flag"
    }
    if frame.empty or not expected.issubset(frame):
        raise ValueError(
            f"NASA TAP query returned unexpected schema or no rows: "
            f"columns={list(frame)}; first bytes={payload[:120]!r}"
        )
    return frame, url, hashlib.sha256(payload).hexdigest()


def previous_tics():
    prior = pd.read_csv(PRIOR, dtype={"tic_id":str})
    if "tic_id" not in prior:
        raise ValueError("Previous dataset lacks tic_id; cannot exclude")
    tics = {normalized_tic(x) for x in prior.tic_id}
    # Also exclude the four original-MAST engineering targets even if a
    # later dataset refactor no longer includes every old positive.
    engineering_manifest = (
        ROOT / "results/mast_original_pilot/verified_original_fits.csv"
    )
    if engineering_manifest.is_file():
        used = pd.read_csv(engineering_manifest, dtype={"tic_id":str})
        if "tic_id" in used:
            tics.update(normalized_tic(v) for v in used.tic_id)
    tics.discard("")
    if len(tics) < 300:
        raise ValueError(
            f"Only {len(tics)} previous target TICs found, expected >=300; "
            "refuse potential reuse/leakage"
        )
    return tics, hashlib.sha256(PRIOR.read_bytes()).hexdigest()


def prepare_queue(source, excluded):
    x = source.copy()
    x["tic_norm"] = x["tic_id"].map(normalized_tic)
    x["planet_period_days"] = pd.to_numeric(x.pl_orbper,errors="coerce")
    x = x.loc[
        x.tic_norm.ne("")
        & x.planet_period_days.between(1.5,9.5,inclusive="both")
        & pd.to_numeric(x.default_flag,errors="coerce").eq(1)
        & pd.to_numeric(x.tran_flag,errors="coerce").eq(1)
        & ~x.tic_norm.isin(excluded)
        & (pd.to_numeric(x.pl_controv_flag,errors="coerce").fillna(0) == 0)
    ].copy()
    x["bin"] = pd.cut(
        x.planet_period_days, RANGES,
        labels=PERIOD_LABELS, include_lowest=True
    ).astype(str)
    # A given TIC may host multiple planets. A single fixed representative
    # planet/ephemeris is chosen deterministically by planet name, not SNR.
    x = x.sort_values(["tic_norm", "pl_name"], kind="mergesort")
    duplicates = x.groupby("tic_norm").size()
    x["known_planets_matching_query_on_host"] = x.tic_norm.map(duplicates)
    x = x.drop_duplicates("tic_norm", keep="first")
    rng = np.random.default_rng(SEED)
    batches = []
    for lab in PERIOD_LABELS:
        group = x[x.bin == lab].copy()
        if len(group):
            ranks = rng.permutation(len(group))
            group["random_order_within_bin"] = ranks
            batches.append(
                group.sort_values("random_order_within_bin").reset_index(drop=True)
            )
    if not batches:
        raise ValueError("No unused NASA-confirmed TESS transit hosts to query")
    # Interleave orbital-period strata so low-period planets do not dominate
    # solely due to source catalogue ordering.
    queue = []
    maxlen = max(len(part) for part in batches)
    for idx in range(maxlen):
        for part in batches:
            if idx < len(part):
                queue.append(part.iloc[idx].to_dict())
    q = pd.DataFrame(queue)
    q.insert(0, "queue_rank", np.arange(1,len(q)+1))
    if q.tic_norm.duplicated().any():
        raise AssertionError("Repeated TIC in NASA host queue")
    if q.tic_norm.isin(excluded).any():
        raise AssertionError("Previously-used TIC in prospective queue")
    return q


def parse_sector(mission):
    match = re.search(r"(?i)\bsector\s*0*(\d+)\b", str(mission))
    return int(match.group(1)) if match else None


def metadata_rows(host, lk):
    tic = host["tic_norm"]
    result = lk.search_lightcurve(
        "TIC " + tic[4:],
        mission="TESS", author="TGLC", exptime=1800
    )
    rows = []
    for ix in range(len(result)):
        if normalized_tic(result.target_name[ix]) != tic:
            continue
        sector = parse_sector(result.mission[ix])
        if sector is None or not 1 <= sector <= 26:
            continue
        author = str(result.author[ix]).strip()
        exposure = result.exptime[ix]
        exp = float(exposure.value if hasattr(exposure,"value") else exposure)
        if author.upper() != "TGLC" or abs(exp - 1800) > 3:
            continue
        table = result.table
        namecols = set(table.colnames)
        def find_field(keys):
            for key in keys:
                if key in namecols:
                    value = table[key][ix]
                    if not np.ma.is_masked(value):
                        return str(value).strip()
            return ""
        uri = find_field(["dataURI","dataUri","data_uri"])
        filename = find_field(["productFilename","productfilename","filename"])
        obsid = find_field(["obsid","obs_id","obsID"])
        if not uri.startswith("mast:"):
            raise RuntimeError(
                f"Available TGLC metadata for {tic} missing original MAST URI"
            )
        rows.append({
            "tic_id": tic,
            "sector": sector,
            "author": author,
            "exptime_seconds": exp,
            "mast_uri": uri,
            "product_filename": filename,
            "mast_obsid": obsid,
            "search_result_index": int(ix),
        })
    # Avoid arbitrary selection if archive returns duplicate product views.
    rows.sort(key=lambda x:(x["sector"],x["mast_uri"]))
    unique = {}
    for row in rows:
        key = (row["sector"],row["mast_uri"])
        unique[key] = row
    return list(unique.values())


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--target-stars",type=int,default=18)
    p.add_argument("--max-mast-queries",type=int,default=70)
    p.add_argument("--pause-seconds",type=float,default=0.20)
    p.add_argument("--archive-csv",type=Path,default=None,
                   help="Use existing NASA snapshot instead of a new TAP query")
    p.add_argument("--outdir",type=Path,default=OUT)
    args = p.parse_args()
    if not 1 <= args.target_stars <= 100:
        raise ValueError("target-stars must be between 1 and 100")
    if args.max_mast_queries < args.target_stars:
        raise ValueError("max-mast-queries must be >= target-stars")
    args.outdir.mkdir(parents=True,exist_ok=True)
    excluded,old_hash = previous_tics()
    if args.archive_csv:
        snapshot = pd.read_csv(args.archive_csv,dtype={"tic_id":str})
        archive_url = "provided_archive_snapshot"
        archive_sha = hashlib.sha256(args.archive_csv.read_bytes()).hexdigest()
    else:
        snapshot,archive_url,archive_sha = download_archive_csv()
    snapshot.to_csv(args.outdir/"archive_snapshot.csv",index=False)
    queue = prepare_queue(snapshot,excluded)
    queue.to_csv(args.outdir/"candidate_queue.csv",index=False)
    print(
        f"NASA hosts in deterministic queue: {len(queue)}; "
        f"excluded previous source TICs: {len(excluded)}",
        flush=True
    )

    import lightkurve as lk
    selected = []
    inventory = []
    queries = []
    for candidate in queue.head(args.max_mast_queries).to_dict("records"):
        tic = candidate["tic_norm"]
        rank = int(candidate["queue_rank"])
        try:
            found = metadata_rows(candidate,lk)
            unique_sectors = sorted({x["sector"] for x in found})
            # Prefer a single predeclared sector: earliest primary-mission
            # 30-min sector with an unambiguous original MAST product.
            usable = []
            for sector in unique_sectors:
                products = [x for x in found if x["sector"]==sector]
                if len(products)==1:
                    usable.append(products[0])
            for prod in found:
                inventory.append({
                    "queue_rank":rank,"planet":candidate["pl_name"],**prod
                })
            if usable:
                chosen = min(usable,key=lambda x:x["sector"])
                selected.append({
                    "selection_rank":len(selected)+1,
                    "queue_rank":rank,
                    "tic_id":tic,
                    "hostname":str(candidate["hostname"]),
                    "planet":str(candidate["pl_name"]),
                    "planet_period_days":float(candidate["planet_period_days"]),
                    "period_bin":str(candidate["bin"]),
                    "n_known_plans_on_host_in_query":int(
                        candidate["known_planets_matching_query_on_host"]
                    ),
                    "transit_epoch_jd":candidate.get("pl_tranmid"),
                    "transit_epoch_time_standard":candidate.get("pl_tsystemref"),
                    "transit_epoch_error1":candidate.get("pl_tranmiderr1"),
                    "transit_epoch_error2":candidate.get("pl_tranmiderr2"),
                    "planet_transit_depth_percent":candidate.get("pl_trandep"),
                    "planet_transit_duration_hours":candidate.get("pl_trandur"),
                    "planet_to_star_radius_ratio":candidate.get("pl_ratror"),
                    "old_dataset_source":False,
                    **chosen,
                })
                status="selected"
            else:
                status="no_unique_primary_mission_1800s_TGLC_product"
            queries.append({
                "queue_rank":rank,"tic_id":tic,"period_bin":candidate["bin"],
                "n_products":len(found),"n_unique_sectors":len(unique_sectors),
                "status":status,"error":"",
            })
            print(f"[{len(queries):03d}] {tic} {candidate['pl_name']}: "
                  f"{status}; accepted={len(selected)}",flush=True)
        except Exception as exc:
            queries.append({
                "queue_rank":rank,"tic_id":tic,"period_bin":candidate["bin"],
                "n_products":0,"n_unique_sectors":0,
                "status":"search_error","error":repr(exc),
            })
            print(f"[{len(queries):03d}] {tic}: SEARCH ERROR {exc}",
                  flush=True)
        if len(selected)>=args.target_stars:
            break
        if args.pause_seconds:
            time.sleep(args.pause_seconds)

    sel=pd.DataFrame(selected)
    inv=pd.DataFrame(inventory)
    audit=pd.DataFrame(queries)
    sel.to_csv(args.outdir/"selected_prospective_hosts.csv",index=False)
    inv.to_csv(args.outdir/"available_original_mast_products.csv",index=False)
    audit.to_csv(args.outdir/"mast_query_audit.csv",index=False)
    if len(sel):
        assert not sel.tic_id.isin(excluded).any()
        assert not sel.tic_id.duplicated().any()
        assert sel.sector.between(1,26).all()
        assert (sel.exptime_seconds-1800).abs().le(3).all()

    cfg={
        "selection_type":"NEW independently inventoried NASA-confirmed transiting hosts",
        "source_archive":"NASA Exoplanet Archive ps default_flag=1 and tran_flag=1",
        "archive_tap_url":archive_url,
        "archive_download_sha256":archive_sha,
        "old_dataset_sha256":old_hash,
        "old_dataset_tics_excluded":len(excluded),
        "candidate_planet_source":"NASA archive default PS row, not the old model label",
        "candidate_queue_seed":SEED,
        "queue_interleaving":"period strata 1.5-3.5;3.5-6.5;6.5-9.5 days",
        "target_stars_requested":args.target_stars,
        "query_cap":args.max_mast_queries,
        "mast_queries_performed":len(queries),
        "hosts_available":len(selected),
        "mast_search_errors":int((audit.status=="search_error").sum()),
        "product_author":"TGLC",
        "product_exptime_seconds":1800,
        "product_sector":"1 through 26",
        "product_selection":"earliest unambiguous available sector",
        "quality_pipeline_frozen":"TGLC_APERTURE_1800S_V1",
        "no_FITS_downloads":True,
        "no_performance_evaluation":True,
        "limitations":[
            "Confirmed host status does NOT guarantee a transit in chosen sector.",
            "TIC identification and NASA default solutions must be audited.",
            "Ephemeris time reference may require correction/verification.",
            "TGLC product metadata matching is necessary but not sufficient.",
            "This is positive-host inventory only; no transit-free negatives.",
            "No future model fitting should use independently frozen test TICs."
        ],
    }
    (args.outdir/"discovery_config.json").write_text(
        json.dumps(cfg,indent=2,default=str),encoding="utf-8"
    )
    counts=sel.groupby("period_bin").size() if len(sel) else pd.Series(dtype=int)
    report=[
        "# NASA confirmed TESS planet-host prospective metadata inventory","",
        f"NASA archive rows returned: {len(snapshot)}",
        f"Prior source TICs excluded: {len(excluded)}",
        f"Unique new TIC candidates queued: {len(queue)}",
        f"MAST metadata queries: {len(queries)}",
        f"Original TGLC 1800-s products identified for NEW TICs: {len(selected)}",
        f"MAST errors: {cfg['mast_search_errors']}","",
        "## Chosen new host sample by orbital period", "",
        counts.to_string() if len(sel) else "(none)", "",
        "## Selected source/planet observations", "",
        sel[[
            "tic_id","planet","planet_period_days","sector","mast_uri"
        ]].to_string(index=False) if len(sel) else "(none)", "",
        "## Interpretation safeguards", "",
        *["- "+x for x in cfg["limitations"]],"",
        "NO original FITS or light curves downloaded in this step."
    ]
    (args.outdir/"discovery_report.md").write_text(
        "\n".join(report)+"\n", encoding="utf-8"
    )
    print("\n"+"\n".join(report[:14]),flush=True)
    if not len(sel):
        raise RuntimeError(
            "No new original TGLC hosts available. Examine search errors, "
            "candidate queue and MAST archive metadata."
        )


if __name__=="__main__":
    main()
