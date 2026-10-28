#!/usr/bin/env python3
"""Acquire EXACTLY the frozen new NASA-confirmed TESS-host MAST products.

Input is experiment-15 selected_prospective_hosts.csv (18 NEW TICs, ONE
preselected original TGLC 1800-second sector and URI each). This operation
does not update the sample, substitute failed hosts, use transit labels to
choose photometry, or tune a detection threshold.

Downloads ORIGINAL archived dataURI via astroquery.mast.Observations.
Original FITS are stored OUTSIDE Git; only sha256-backed provenance,
acquisition failures and HDU inventories are published in results/.
The output verified_original_fits.csv is directly compatible with
revision/14_build_uniform_aperture_curves.py --manifest.

Run:
  python -u revision/16_acquire_frozen_new_hosts.py
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
from astroquery.mast import Observations

ROOT = Path(__file__).resolve().parents[1]
SELECTED = (
    ROOT / "results/new_confirmed_tess_inventory/"
    "selected_prospective_hosts.csv"
)
OLD_DATASET = ROOT / "data/dataset/final_dataset.csv"
ENGINEERING = ROOT / "results/mast_original_pilot/verified_original_fits.csv"
OUT = ROOT / "results/new_confirmed_tess_acquisition"
PRIVATE = Path.home() / "exoplanet_revision/data/new_confirmed_tess_original"
MIN_ROWS = 300


def canonical_tic(tic):
    m = re.fullmatch(r"TIC_(\d+)", str(tic).strip())
    if m is None:
        raise ValueError(f"Noncanonical TIC identifier: {tic!r}")
    return m.group(1)


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as file:
        for block in iter(lambda: file.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_name_for(row):
    """Filename is derived from dataURI hash, never from a remote path."""
    tic = canonical_tic(row["tic_id"])
    sector = int(row["sector"])
    uri = str(row["mast_uri"])
    if not uri.startswith("mast:HLSP/tglc/"):
        raise ValueError("Expected frozen original MAST TGLC HLSP dataURI")
    m = re.search(r"/s(\d{4})/", uri)
    if m is None or int(m.group(1)) != sector:
        raise ValueError("MAST URI does not encode the selected TESS sector")
    extension = ".fits.gz" if uri.lower().endswith(".fits.gz") else ".fits"
    if not uri.lower().endswith((".fits", ".fits.gz")):
        raise ValueError("Selection URI is not a FITS product")
    short = hashlib.sha256(uri.encode("utf-8")).hexdigest()[:12]
    return f"TIC_{tic}_s{sector:04d}_TGLC1800_{short}{extension}"


def check_selection(selected):
    required = {
        "tic_id", "sector", "selection_rank", "mast_uri", "author",
        "exptime_seconds", "planet", "planet_period_days",
    }
    if selected.empty or not required.issubset(selected):
        raise ValueError("Missing or empty frozen selection")
    if selected.tic_id.duplicated().any():
        raise ValueError("Frozen selection has repeated TIC IDs")
    if selected.selection_rank.duplicated().any():
        raise ValueError("Frozen selection has repeated ranks")
    if not selected.tic_id.map(lambda x: bool(re.fullmatch(r"TIC_\d+", str(x)))).all():
        raise ValueError("Noncanonical frozen TIC identifier")
    if not selected.sector.between(1, 26).all():
        raise ValueError("Only primary-mission TGLC sectors 1-26 allowed")
    if not selected.author.astype(str).str.upper().eq("TGLC").all():
        raise ValueError("Product author changed after freezing")
    if not (pd.to_numeric(selected.exptime_seconds)-1800).abs().le(3).all():
        raise ValueError("Expected 1800-second products")
    for entry in selected.to_dict("records"):
        file_name_for(entry)

    prior = pd.read_csv(OLD_DATASET, dtype={"tic_id":str})
    if "tic_id" not in prior:
        raise ValueError("Original full dataset lacks tic_id")
    old = set(prior.tic_id)
    if ENGINEERING.is_file():
        pilot = pd.read_csv(ENGINEERING, dtype={"tic_id":str})
        old.update(pilot.tic_id)
    if set(selected.tic_id) & old:
        raise ValueError(
            "Leakage: new frozen sample intersects previously inspected stars"
        )


def inspect_original(path, row):
    """Validate expected TGLC provenance and raw science table, no tuning."""
    if path.stat().st_size < 2048:
        raise ValueError("MAST original FITS is unexpectedly small")
    with fits.open(path, memmap=True, checksum=False) as hdul:
        if len(hdul) < 2:
            raise ValueError("No science HDU in original FITS")
        telescope = str(hdul[0].header.get("TELESCOP", "")).strip()
        origin = str(hdul[0].header.get("ORIGIN", "")).strip()
        if telescope and telescope.upper() != "TESS":
            raise ValueError("Original FITS telescope is not TESS")
        if origin and "TGLC" not in origin.upper():
            raise ValueError(
                "Unexpected ORIGIN in nominal TGLC original FITS: "+origin
            )
        target = canonical_tic(row["tic_id"])
        detected = []
        science_rows = 0
        column_names = []
        for idx, hd in enumerate(hdul):
            for key in ("TICID", "TIC_ID"):
                if key in hd.header:
                    found = str(int(hd.header[key]))
                    if found != target:
                        raise ValueError(
                            f"Header {key}={found} mismatches selected TIC_{target}"
                        )
            if "SECTOR" in hd.header:
                if int(hd.header["SECTOR"]) != int(row["sector"]):
                    raise ValueError("Header sector mismatches frozen selection")
            c = getattr(hd, "columns", None)
            names = list(c.names or []) if c is not None else []
            lowered = {str(n).lower():str(n) for n in names}
            if {"time", "cal_aper_flux", "tess_flags",
                "cadence_num"}.issubset(lowered):
                n = len(hd.data) if hd.data is not None else 0
                detected.append(idx)
                science_rows = n
                column_names = names
        if len(detected) != 1:
            raise ValueError(
                f"Expected exactly one TGLC calibrated-aperture science "
                f"extension, found {len(detected)}"
            )
        if science_rows < MIN_ROWS:
            raise ValueError(
                f"Too few original time-series cadences: {science_rows}"
            )
        meta = {
            "n_science_rows":science_rows,
            "science_hdu_index":detected[0],
            "original_fits_origin":origin,
            "original_fits_telescope":telescope,
            "original_hdu_count":len(hdul),
            "time_column_name":next(
                col for col in column_names if col.lower()=="time"
            ),
            "flux_column_name":next(
                col for col in column_names if col.lower()=="cal_aper_flux"
            ),
        }
        columns = [
            {
                "tic_id":str(row["tic_id"]),
                "sector":int(row["sector"]),
                "hdu":detected[0],
                "column_name":str(col),
            } for col in column_names
        ]
    return meta, columns


def download_once(uri, path):
    """Do not silently replace already obtained files."""
    if path.exists():
        return "previously_downloaded"
    path.parent.mkdir(parents=True, exist_ok=True)
    # No alternate catalogue/product selection or mirror fallback.
    response = Observations.download_file(
        uri, local_path=str(path), cache=True, verbose=False
    )
    if not response or str(response[0]) not in {"COMPLETE", "SKIPPED"}:
        raise RuntimeError(
            f"MAST download returned non-success status: {response!r}"
        )
    if not path.exists():
        raise FileNotFoundError("MAST claims success but file is absent")
    return "new_download"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--selected",type=Path,default=SELECTED)
    ap.add_argument("--outdir",type=Path,default=OUT)
    ap.add_argument("--private-dir",type=Path,default=PRIVATE)
    ap.add_argument("--limit",type=int,default=0,
                    help="0 = entire frozen cohort. Lower values are partial engineering runs.")
    args = ap.parse_args()
    if args.limit < 0:
        raise ValueError("--limit must be >=0")
    if args.private_dir.resolve().is_relative_to(ROOT.resolve()):
        raise ValueError("Original FITS must remain outside the Git repo")
    args.outdir.mkdir(parents=True,exist_ok=True)
    args.private_dir.mkdir(parents=True,exist_ok=True)

    selected = pd.read_csv(
        args.selected, dtype={"tic_id":str,"planet":str}
    ).sort_values("selection_rank").reset_index(drop=True)
    check_selection(selected)
    original_count = len(selected)
    frozen_sha = sha256(args.selected)
    if args.limit:
        selected=selected.head(args.limit)

    results, diagnostics, errors = [], [], []
    for i, row in enumerate(selected.to_dict("records"),1):
        tic = str(row["tic_id"])
        sector = int(row["sector"])
        try:
            dst = args.private_dir / f"sector_{sector:02d}" / file_name_for(row)
            status = download_once(str(row["mast_uri"]),dst)
            meta, columns = inspect_original(dst,row)
            digest = sha256(dst)
            diagnostics.extend(columns)
            result = {
                **row,
                **meta,
                "local_path":str(dst),
                "file_bytes":dst.stat().st_size,
                "sha256":digest,
                "source_selection_sha256":frozen_sha,
                "status":"verified_original_fits",
                "acquisition_action":status,
            }
            results.append(result)
            print(
                f"[{i:02d}/{len(selected):02d}] OK {tic} "
                f"sector {sector}: {meta['n_science_rows']} raw cadences; "
                f"{dst.stat().st_size/1024:.1f} KiB",
                flush=True,
            )
        except Exception as exc:
            errors.append({
                "selection_rank":int(row["selection_rank"]),
                "tic_id":tic,
                "sector":sector,
                "mast_uri":str(row["mast_uri"]),
                "error":repr(exc),
            })
            print(
                f"[{i:02d}/{len(selected):02d}] FAIL {tic}: {exc}",
                flush=True,
            )

    manifest = pd.DataFrame(results)
    manifest.to_csv(args.outdir/"verified_original_fits.csv",index=False)
    pd.DataFrame(diagnostics,columns=[
        "tic_id","sector","hdu","column_name"
    ]).to_csv(args.outdir/"science_column_inventory.csv",index=False)
    pd.DataFrame(errors,columns=[
        "selection_rank","tic_id","sector","mast_uri","error"
    ]).to_csv(args.outdir/"acquisition_errors.csv",index=False)
    config = {
        "source_selection_csv":str(args.selected),
        "source_selection_sha256":frozen_sha,
        "frozen_cohort_size":original_count,
        "frozen_source_tics":selected.tic_id.tolist() if not args.limit else (
            "partial run: inspect input frozen selection for all source TICs"
        ),
        "sources_attempted":len(selected),
        "sources_verified":len(results),
        "source_failures":len(errors),
        "no_failed_source_replacement":True,
        "only_preselected_mast_uris":True,
        "no_lightcurve_reexport":True,
        "product":"original TGLC 1800-second primary-mission FITS",
        "next_policy":"TGLC_APERTURE_1800S_V1",
        "science_fits_outside_git":str(args.private_dir),
        "notes":[
            "Some confirmed hosts may have no transit observable in chosen sector.",
            "A confirmed planet does not imply a model-detected transit.",
            "No TIC is substituted for a failed acquisition or poor quality.",
            "This acquisition step cannot validate exoplanet classifications.",
        ],
    }
    (args.outdir/"acquisition_config.json").write_text(
        json.dumps(config,indent=2,default=str),encoding="utf-8"
    )
    report = [
        "# Frozen NEW TESS confirmed-host original-MAST FITS acquisition","",
        f"Preselected unique TICs: {original_count}",
        f"TICs attempted in this run: {len(selected)}",
        f"Original FITS verified: {len(manifest)}",
        f"Acquisition / integrity failures: {len(errors)}",
        f"Original FITS storage (outside Git): {args.private_dir}",
        f"Policy to apply unchanged: TGLC_APERTURE_1800S_V1","",
        "Source selection is FROZEN; failed acquisitions are NOT replaced.",
        "No planet-classification or transit-detection scores were computed.",
        "Selected objects are catalogued planet HOSTS, not verified "
        "in-sector recovered transits.",
    ]
    (args.outdir/"acquisition_report.md").write_text(
        "\n".join(report)+"\n",encoding="utf-8"
    )
    print("\n"+"\n".join(report),flush=True)
    if errors:
        raise RuntimeError(
            f"{len(errors)} acquisitions failed; read acquisition_errors.csv. "
            "Do not replace failed targets with a different star."
        )


if __name__=="__main__":
    main()
