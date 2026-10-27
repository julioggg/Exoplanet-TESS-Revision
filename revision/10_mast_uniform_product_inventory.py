#!/usr/bin/env python3
"""MAST feasibility inventory for homogeneous TESS products.

Technological feasibility study ONLY: metadata search, NO downloads.
Select matched positive-negative TIC pairs from the already examined sector
matching study; search each TIC for Lightkurve/MAST products in its recorded
sector. Determine whether both members share exactly the SAME:
  - data product author / reduction pipeline;
  - observing sector;
  - nominal exposure seconds.

This pilot is for acquiring a NEW study sample later. It is NOT an independent
heldout analysis, not label verification, and not photometric validation.

Usage:
  python -u revision/10_mast_uniform_product_inventory.py --pairs 8
"""
from __future__ import annotations

import argparse
import json
import re
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PAIRS = ROOT / "results/sector_matched/sector_matched_pairs.csv"
OUTPUT = ROOT / "results/mast_inventory"


def tic_digits(tic):
    d = re.sub(r"[^0-9]", "", str(tic))
    if not d:
        raise ValueError(f"Invalid TIC: {tic}")
    return d


def scalar(value):
    if hasattr(value, "value"):
        return value.value
    return value


def do_search(tic, sector, lk):
    query = f"TIC {tic_digits(tic)}"
    result = lk.search_lightcurve(
        query, mission="TESS", sector=int(sector)
    )
    n = len(result)
    rows = []
    if n == 0:
        return rows
    authors = list(result.author)
    missions = list(result.mission)
    expir = scalar(result.exptime)
    exposures = list(np.atleast_1d(expir))
    targets = list(result.target_name)
    if len(authors) != n or len(missions) != n or len(exposures) != n:
        raise ValueError(f"MAST metadata length mismatch for TIC {tic}")
    for i in range(n):
        exposure = float(scalar(exposures[i]))
        author = str(authors[i]).strip()
        target_name = str(targets[i]).strip() if len(targets) == n else ""
        target_digits = tic_digits(target_name) if re.search(
            r"\d", target_name
        ) else ""
        strict_target_match = target_digits == tic_digits(tic)
        rows.append({
            "tic_id": str(tic),
            "sector_query": int(sector),
            "author": author,
            "exptime_seconds": exposure,
            "exptime_rounded_seconds": int(round(exposure)),
            "mission": str(missions[i]),
            "target_name_from_mast": target_name,
            "exact_tic_match": bool(strict_target_match),
            "data_key": f"{author}|{int(round(exposure))}",
            "product_index_in_search": i,
        })
    return rows


def select_pairs(all_pairs, limit, tight_only):
    if tight_only:
        all_pairs = all_pairs[
            all_pairs["passes_calipers"].astype(str).str.lower() == "true"
        ].copy()
    if not len(all_pairs):
        raise ValueError("No eligible matched pairs")
    # Across sectors rather than repeatedly querying many TICs in one sector.
    selected = (
        all_pairs.sort_values(
            ["match_cost_unpenalized", "sector", "pair_id"]
        ).drop_duplicates(subset=["sector"], keep="first")
        .head(limit).copy()
    )
    if len(selected) < limit:
        rest = all_pairs[~all_pairs["pair_id"].isin(selected["pair_id"])]
        selected = pd.concat([selected, rest.head(limit - len(selected))])
    return selected.reset_index(drop=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", type=int, default=8)
    ap.add_argument("--all-sector-pairs", action="store_true")
    ap.add_argument("--pause-seconds", type=float, default=0.25)
    ap.add_argument("--outdir", type=Path, default=OUTPUT)
    args = ap.parse_args()
    if args.pairs < 1 or args.pause_seconds < 0:
        raise ValueError("Invalid pair count or pause")
    args.outdir.mkdir(parents=True, exist_ok=True)

    import lightkurve as lk

    candidate = pd.read_csv(PAIRS)
    selected = select_pairs(
        candidate, args.pairs, not args.all_sector_pairs
    )
    selected.to_csv(args.outdir / "queried_matched_pairs.csv", index=False)

    records = []
    issues = []
    pair_results = []
    cache = {}

    for i, pair in selected.iterrows():
        sector = int(pair.sector)
        print(
            f"[{i+1:02d}/{len(selected):02d}] sector {sector}, "
            f"{pair.negative_tic} vs {pair.positive_tic}",
            flush=True,
        )
        per = {}
        for label, tic in [
            (0, str(pair.negative_tic)), (1, str(pair.positive_tic))
        ]:
            key = (tic, sector)
            if key not in cache:
                try:
                    cache[key] = do_search(tic, sector, lk)
                except Exception as err:
                    issues.append({
                        "tic_id": tic, "sector": sector,
                        "pair_id": pair.pair_id,
                        "error": repr(err)
                    })
                    cache[key] = None
                if args.pause_seconds:
                    time.sleep(args.pause_seconds)
            found = cache[key]
            per[label] = [] if found is None else found
            if found is not None:
                for item in found:
                    records.append({
                        "pair_id": pair.pair_id, "label": label, **item
                    })
        neg = {
            v["data_key"] for v in per[0] if v["exact_tic_match"]
        }
        pos = {
            v["data_key"] for v in per[1] if v["exact_tic_match"]
        }
        common = sorted(neg & pos)
        pair_results.append({
            "pair_id": str(pair.pair_id),
            "recorded_sector": sector,
            "negative_tic": str(pair.negative_tic),
            "positive_tic": str(pair.positive_tic),
            "n_negative_exact_products": sum(
                int(x["exact_tic_match"]) for x in per[0]
            ),
            "n_positive_exact_products": sum(
                int(x["exact_tic_match"]) for x in per[1]
            ),
            "n_common_author_exptime_modes": len(common),
            "common_author_exptime_modes": ";".join(common),
            "query_success_both": not any(
                (tic, sector) in cache and cache[(tic, sector)] is None
                for tic in [str(pair.negative_tic), str(pair.positive_tic)]
            ),
        })

    pd.DataFrame(records).to_csv(
        args.outdir / "mast_product_inventory.csv", index=False
    )
    summary = pd.DataFrame(pair_results)
    summary.to_csv(
        args.outdir / "same_sector_same_pipeline_pairs.csv", index=False
    )
    pd.DataFrame(issues, columns=[
        "tic_id", "sector", "pair_id", "error"
    ]).to_csv(args.outdir / "mast_search_errors.csv", index=False)

    ok_pairs = int((summary.n_common_author_exptime_modes > 0).sum())
    config = {
        "queried_pairs": len(selected),
        "unique_queried_sectors": selected.sector.nunique(),
        "pairs_with_matching_author_exptime_sector": ok_pairs,
        "queries_with_errors": len(issues),
        "no_science_fits_downloaded": True,
        "exact_target_match_required": True,
        "notes": [
            "This is only a metadata availability survey of previously used stars.",
            "Matching author and exposure does not prove identical photometric quality.",
            "A source outside its recorded exported sector may still have other products.",
            "MAST network errors are not proof of missing products.",
            "Final independent validation will require newly chosen target TICs.",
        ],
    }
    (args.outdir / "survey_config.json").write_text(
        json.dumps(config, indent=2), encoding="utf-8"
    )
    report = [
        "# MAST availability survey for homogeneous TESS products", "",
        f"Matched pairs queried: {len(selected)}",
        f"Unique sectors queried: {selected.sector.nunique()}",
        f"Pairs with both TICs available at identical author/exptime: {ok_pairs}",
        f"MAST search errors: {len(issues)}", "",
        "## Pair-level product matching", "",
        summary.to_string(index=False), "",
        "A matched author/exposure sector only establishes basic product compatibility.",
        "This script did not download original FITS or validate astrophysical labels.",
        "The next data collection must use independent source stars.",
    ]
    (args.outdir / "mast_survey_report.md").write_text(
        "\n".join(report) + "\n", encoding="utf-8"
    )
    print("\n" + "\n".join(report[:10]), flush=True)


if __name__ == "__main__":
    main()
