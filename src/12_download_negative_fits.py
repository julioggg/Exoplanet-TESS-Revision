import pandas as pd
import lightkurve as lk

from pathlib import Path

print("=" * 60)
print("MASS DOWNLOAD OF NEGATIVE LIGHT CURVES")
print("=" * 60)

NEGATIVES = pd.read_csv(
    "data/catalog/negative_targets_valid.csv"
)

OUTDIR = Path("data/raw/negatives")
OUTDIR.mkdir(parents=True, exist_ok=True)

downloaded = 0
skipped = 0

for i, tic in enumerate(NEGATIVES["tic_id"], start=1):

    tic = str(tic)

    outfile = (
        OUTDIR /
        f"{tic.replace(' ', '_')}.fits"
    )

    if outfile.exists():

        skipped += 1

        print(
            f"\n[{i}/{len(NEGATIVES)}] {tic}"
        )
        print("  Already downloaded")

        continue

    print(
        f"\n[{i}/{len(NEGATIVES)}] {tic}"
    )

    try:

        search = lk.search_lightcurve(
            tic,
            mission="TESS"
        )

        if len(search) == 0:

            print("  No observations")
            continue

        collection = search.download_all()

        if collection is None:

            print("  Download failed")
            continue

        stitched = collection.stitch()

        stitched.to_fits(
            path=str(outfile),
            overwrite=True
        )

        print(
            f"  Saved ({len(stitched.time)} cadences)"
        )

        downloaded += 1

    except Exception as e:

        print("  ERROR:", e)

print()
print("=" * 60)
print("SUMMARY")
print("=" * 60)

print("Downloaded :", downloaded)
print("Skipped    :", skipped)
print("Total      :", len(NEGATIVES))