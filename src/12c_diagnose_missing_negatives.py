import pandas as pd
import lightkurve as lk

print("=" * 60)
print("DIAGNOSING MISSING NEGATIVE TARGETS")
print("=" * 60)

missing = pd.read_csv(
    "data/catalog/missing_negative_targets.csv"
)

results = []

for i, tic in enumerate(missing["tic_id"], start=1):

    tic = str(tic)
    message = ""
    print(f"\n[{i}/{len(missing)}] {tic}")

    try:

        search = lk.search_lightcurve(
            tic,
            mission="TESS"
        )

        n_results = len(search)

        print(f"  Search results: {n_results}")

        if n_results == 0:

            status = "NO_SEARCH_RESULTS"

        else:

            try:

                lc = search.download_all()

                if lc is None:

                    status = "DOWNLOAD_RETURNED_NONE"

                else:

                    status = "DOWNLOAD_OK"

            except Exception as e:

                status = type(e).__name__

                message = str(e)

                print("  ERROR:", message)

        print(f"  Status: {status}")

    except Exception as e:

        status = f"SEARCH_ERROR: {type(e).__name__}"

        print(f"  Status: {status}")

    results.append({
    "tic_id": tic,
    "status": status,
    "message": message
    })

results = pd.DataFrame(results)

outfile = (
    "data/catalog/missing_negative_diagnosis.csv"
)

results.to_csv(
    outfile,
    index=False
)

print()
print("=" * 60)
print("SUMMARY")
print("=" * 60)

print(results["status"].value_counts())

print()
print("Saved:")
print(outfile)