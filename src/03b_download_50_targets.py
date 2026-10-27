import pandas as pd
from pathlib import Path
import lightkurve as lk
import time

N_TARGETS = 50

targets = pd.read_csv(
    "data/catalog/targets.csv"
)

base_dir = Path(
    "data/raw/tess"
)

base_dir.mkdir(
    parents=True,
    exist_ok=True
)

log = []

print("=" * 60)
print("DESCARGA MASIVA TESS")
print("=" * 60)

for i, tic in enumerate(
    targets["tic_id"].head(N_TARGETS),
    start=1
):

    tic_clean = (
        str(tic)
        .replace("TIC", "")
        .strip()
    )

    print(
        f"\n[{i}/{N_TARGETS}] TIC {tic_clean}"
    )

    try:

        search = lk.search_lightcurve(
            f"TIC {tic_clean}",
            mission="TESS"
        )

        n_sectors = len(search)

        print(
            f"  Sectores encontrados: {n_sectors}"
        )

        if n_sectors == 0:

            log.append(
                [tic_clean, 0, "NO_DATA"]
            )

            continue

        target_dir = (
            base_dir /
            f"TIC_{tic_clean}"
        )

        target_dir.mkdir(
            parents=True,
            exist_ok=True
        )

        products = search.download_all()

        saved = 0

        if products is not None:

            for j, lc in enumerate(products):

                outfile = (
                    target_dir /
                    f"sector_{j+1:02d}.fits"
                )

                lc.to_fits(
                    path=str(outfile),
                    overwrite=True
                )

                saved += 1

        log.append(
            [
                tic_clean,
                saved,
                "OK"
            ]
        )

        print(
            f"  Guardados: {saved}"
        )

    except Exception as e:

        print(
            f"  ERROR: {e}"
        )

        log.append(
            [
                tic_clean,
                0,
                str(e)
            ]
        )

    time.sleep(0.5)

log_df = pd.DataFrame(
    log,
    columns=[
        "tic_id",
        "n_files",
        "status"
    ]
)

log_df.to_csv(
    "data/raw/download_log.csv",
    index=False
)

print("\nProceso completado.")