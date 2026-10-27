#!/usr/bin/env python3
"""Download a SMALL pilot of ORIGINAL, provenance-verified MAST TESS FITS.

This is an engineering feasibility check using previously examined TICs, NOT
new independent validation. No trained detector, labels, BLS or numerical
results are updated.

Default: up to two same-sector positive/negative pairs with TGLC 1800-s
products identified in the existing 8-pair MAST survey. Files are downloaded
verbatim using the MAST dataURI (not LightCurve.to_fits()). Exact TIC, author,
sector, exposure and header metadata are audited. Original FITS stay OUTSIDE
the git repository; only compact manifests and column inventories are saved.

Example:
    python -u revision/11_download_original_mast_pilot.py --pairs 2 \
      --author TGLC --exptime 1800
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
PAIR_INPUT = ROOT / "results/mast_inventory/same_sector_same_pipeline_pairs.csv"
RESULTS_DIR = ROOT / "results/mast_original_pilot"
STORAGE_DIR = Path.home() / "exoplanet_revision/data/mast_original_pilot"


def digits(value):
    match = re.findall(r"\d+", str(value))
    return "".join(match) if match else ""


def stringify(v):
    try:
        if np.ma.is_masked(v):
            return ""
    except Exception:
        pass
    return str(v)


def sha256_file(path):
    sha = hashlib.sha256()
    with path.open("rb") as fp:
        for chunk in iter(lambda: fp.read(1024 * 1024), b""):
            sha.update(chunk)
    return sha.hexdigest()


def get_row_value(table, ix, choices):
    for name in choices:
        if name in table.colnames:
            val = table[name][ix]
            if np.ma.is_masked(val):
                continue
            value = stringify(val).strip()
            if value and value.lower() not in {"none", "nan", "null"}:
                return value
    return ""


def find_selected_product(tic, sector, author, exptime, lk):
    """
    Require unique, EXACT TIC target match, sector, author and integration time.
    Do NOT silently select search result[0] from an unfiltered table.
    """
    result = lk.search_lightcurve(
        f"TIC {digits(tic)}",
        mission="TESS",
        sector=int(sector),
        author=author,
        exptime=int(exptime),
    )
    chosen = []
    for ix in range(len(result)):
        target = stringify(result.target_name[ix])
        if digits(target) != digits(tic):
            continue
        auth = stringify(result.author[ix]).strip()
        exposure = float(result.exptime[ix].value if hasattr(
            result.exptime[ix], "value") else result.exptime[ix])
        mission = stringify(result.mission[ix])
        if auth.upper() != author.upper():
            continue
        if abs(exposure - exptime) > 3.0:
            continue
        # Query includes sector, but check returned sector encoded in mission.
        if not re.search(r"(?i)sector\s*0*" + str(int(sector)) + r"\b", mission):
            continue
        tab = result.table
        uri = get_row_value(
            tab, ix, ("dataURI", "dataUri", "data_uri", "dataURL")
        )
        filename = get_row_value(
            tab, ix, ("productFilename", "productfilename", "filename")
        )
        obsid = get_row_value(tab, ix, ("obsid", "obs_id", "obsID"))
        if not uri.startswith("mast:"):
            raise RuntimeError(
                f"Matched science product for {tic} {sector} "
                "but original dataURI is not available in SearchResult.table. "
                "Please provide its columns: " + str(tab.colnames)
            )
        chosen.append({
            "tic_id": tic, "sector": int(sector), "author": author,
            "exptime_seconds": exposure, "mast_uri": uri,
            "product_filename": filename, "mast_obsid": obsid,
            "search_index": ix, "mission": mission,
            "query_result_count": len(result),
        })
    if len(chosen) != 1:
        raise RuntimeError(
            f"Expected exactly 1 original FITS product for TIC={tic} "
            f"sector={sector}, author={author}, exptime={exptime}, "
            f"found {len(chosen)}. Do not arbitrarily choose."
        )
    return chosen[0]


def inspect_fits(path, entry):
    with fits.open(path, memmap=True, checksum=False) as hdul:
        if len(hdul) < 2:
            raise ValueError("Original FITS has no science extension")
        primary = hdul[0].header
        columns = []
        n_science = 0
        for hduix, hdu in enumerate(hdul):
            if getattr(hdu, "columns", None) is not None:
                names = list(hdu.columns.names or [])
                columns.extend({
                    "tic_id": entry["tic_id"], "sector": entry["sector"],
                    "author": entry["author"], "hdu": hduix,
                    "extension": str(hdu.header.get("EXTNAME", "")),
                    "column": name,
                    "unit": str(hdu.columns[name].unit or "")
                } for name in names)
                if "TIME" in names:
                    n_science = len(hdu.data) if hdu.data is not None else 0
        hdr = {
            "telescop": str(primary.get("TELESCOP", "")),
            "origin": str(primary.get("ORIGIN", "")),
            "creator": str(primary.get("CREATOR", "")),
            "procver": str(primary.get("PROCVER", "")),
            "sector_in_header": str(primary.get("SECTOR", "")),
            "ticid_in_header": str(primary.get("TICID", "")),
            "n_science_rows": int(n_science),
            "n_hdus": len(hdul),
        }
        # If provided, verify a header TIC matches expected target.
        for hdu in hdul:
            for keyword in ["TICID"]:
                hv = hdu.header.get(keyword, None)
                if hv is not None and digits(hv) and digits(hv) != digits(
                    entry["tic_id"]
                ):
                    raise ValueError(
                        "Source TIC conflicts with a science-product header"
                    )
        if n_science < 300:
            raise ValueError(f"Insufficient science rows: {n_science}")
        if hdr["telescop"] and hdr["telescop"].upper() != "TESS":
            raise ValueError("Unexpected telescope keyword in original FITS")
    return hdr, columns


def choose_pairs(input_table, author, exptime, count):
    mode = f"{author}|{int(exptime)}"
    eligible = input_table[
        input_table["common_author_exptime_modes"].fillna("").map(
            lambda s: mode in s.split(";")
        )
    ].copy()
    if not len(eligible):
        raise RuntimeError(
            f"No known same-author same-exptime pairs for {mode}"
        )
    eligible = eligible.sort_values(
        ["recorded_sector", "pair_id"]
    ).head(count)
    return eligible


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pairs", type=int, default=2)
    parser.add_argument("--author", default="TGLC")
    parser.add_argument("--exptime", type=int, default=1800)
    parser.add_argument("--store", type=Path, default=STORAGE_DIR)
    parser.add_argument("--outdir", type=Path, default=RESULTS_DIR)
    args = parser.parse_args()
    if args.pairs < 1 or args.pairs > 8:
        raise ValueError("Only 1-8 paired products allowed for engineering pilot")
    if args.exptime < 1:
        raise ValueError("Exposure time must be positive")
    args.outdir.mkdir(parents=True, exist_ok=True)
    args.store.mkdir(parents=True, exist_ok=True)
    import lightkurve as lk
    selected = choose_pairs(
        pd.read_csv(PAIR_INPUT), args.author, args.exptime, args.pairs
    )
    print(f"ORIGINAL MAST PILOT: {len(selected)} matched pairs "
          f"using {args.author} at {args.exptime}s", flush=True)

    rows, cols, errors = [], [], []
    for ix, pair in selected.iterrows():
        for label, tic in [
            (0, str(pair["negative_tic"])),
            (1, str(pair["positive_tic"]))
        ]:
            sector = int(pair["recorded_sector"])
            base = {
                "pair_id": str(pair["pair_id"]),
                "label_from_old_catalogue_not_verified": label,
                "tic_id": tic, "sector": sector,
            }
            try:
                entry = find_selected_product(
                    tic, sector, args.author, args.exptime, lk
                )
                # Original archive filename is preserved in a separate column.
                # The local name cannot collide across TIC/sector/product mode.
                outfolder = args.store / f"sector_{sector:02d}"
                outfolder.mkdir(parents=True, exist_ok=True)
                # Retain .fits.gz when the archived science product is gzip
                # encoded. Never trust its original basename as a local path.
                suffix = (
                    ".fits.gz" if entry["product_filename"].lower().endswith(".gz")
                    else ".fits"
                )
                path = outfolder / (
                    f"{tic}_{args.author}_{args.exptime}_ORIGINAL{suffix}"
                )
                if not path.is_file():
                    result = Observations.download_file(
                        entry["mast_uri"], local_path=str(path),
                        cache=True, verbose=False
                    )
                    status = str(result[0]) if result else "UNKNOWN"
                    if status not in {"COMPLETE", "SKIPPED"}:
                        raise RuntimeError(
                            f"MAST original download returned {result}"
                        )
                if not path.is_file() or path.stat().st_size < 2048:
                    raise RuntimeError("Original FITS missing or too small")
                header, column_rows = inspect_fits(path, entry)
                cols.extend(column_rows)
                rows.append({
                    **base, **entry, **header,
                    "local_path": str(path),
                    "file_bytes": path.stat().st_size,
                    "sha256": sha256_file(path),
                    "status": "verified_original_fits",
                })
                print(f"[OK] TIC {tic} sector {sector}: "
                      f"{path.stat().st_size / 1024:.0f} KiB; "
                      f"rows={header['n_science_rows']}", flush=True)
            except Exception as exc:
                errors.append({
                    **base, "error": repr(exc),
                })
                print(f"[FAIL] TIC {tic} sector {sector}: {exc}",
                      flush=True)

    manifest = pd.DataFrame(rows)
    manifest.to_csv(args.outdir / "verified_original_fits.csv", index=False)
    pd.DataFrame(cols, columns=[
        "tic_id", "sector", "author", "hdu", "extension", "column", "unit"
    ]).to_csv(args.outdir / "original_fits_columns.csv", index=False)
    pd.DataFrame(errors, columns=[
        "pair_id", "label_from_old_catalogue_not_verified",
        "tic_id", "sector", "error"
    ]).to_csv(args.outdir / "mast_original_errors.csv", index=False)
    selected.to_csv(args.outdir / "requested_pairs.csv", index=False)

    grouped = (
        manifest.groupby("pair_id")["label_from_old_catalogue_not_verified"]
        .agg(["count", "sum"]).reset_index()
        if len(manifest) else pd.DataFrame()
    )
    n_complete = int(
        ((grouped["count"] == 2) & (grouped["sum"] == 1)).sum()
    ) if len(grouped) else 0
    cfg = {
        "goal": "original FITS acquisition feasibility, not model validation",
        "n_pairs_requested": len(selected),
        "n_pairs_both_verified": n_complete,
        "n_fits_verified": len(manifest),
        "n_errors": len(errors),
        "author": args.author, "exptime_seconds": args.exptime,
        "files_live_outside_repository": str(args.store),
        "download_used_original_dataURI": True,
        "no_lightkurve_to_fits": True,
        "no_classification_training": True,
        "limitation": (
            "Prior labels and TICs were used only for engineering. "
            "Need new independently drawn and verified positive/negative "
            "labels before prospective science validation."
        ),
    }
    (args.outdir / "acquisition_config.json").write_text(
        json.dumps(cfg, indent=2), encoding="utf-8"
    )
    lines = [
        "# MAST original-product acquisition pilot", "",
        f"Matched pairs requested: {len(selected)}",
        f"Original FITS verified: {len(manifest)}",
        f"Complete same-sector/author/exptime pairs: {n_complete}",
        f"Failures: {len(errors)}",
        f"Product mode: {args.author}, {args.exptime} s",
        f"Raw storage outside git: {args.store}", "",
        "See original_fits_columns.csv before selecting a photometric flux column.",
        "A local SHA256 content hash and original MAST URI give provenance; "
        "no science validation is implied.",
    ]
    (args.outdir / "acquisition_report.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )
    print("\n" + "\n".join(lines), flush=True)


if __name__ == "__main__":
    main()
