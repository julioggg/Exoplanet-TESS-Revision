import pandas as pd
import lightkurve as lk

targets = pd.read_csv("data/catalog/targets.csv")

N = 20

print("=" * 60)
print("ESTIMACION DE SECTORES TESS")
print("=" * 60)

total_sectors = 0

for i, tic in enumerate(targets["tic_id"].head(N), start=1):

    tic_clean = str(tic).replace("TIC", "").strip()

    try:

        search = lk.search_lightcurve(
            f"TIC {tic_clean}",
            mission="TESS"
        )

        n = len(search)

        total_sectors += n

        print(
            f"[{i:02d}/{N}] "
            f"TIC {tic_clean} -> {n} sectores"
        )

    except Exception as e:

        print(
            f"[{i:02d}/{N}] "
            f"TIC {tic_clean} -> ERROR {e}"
        )

print("\nTotal sectores:", total_sectors)
print("Promedio:", total_sectors / N)