import numpy as np
import pandas as pd

from pathlib import Path

PROCESSED_DIR = Path("data/processed")
FEATURES_FILE = Path("data/features/bls_features.csv")

OUT_DIR = Path("data/dataset")
OUT_DIR.mkdir(parents=True, exist_ok=True)

bls = pd.read_csv(FEATURES_FILE)

rows = []

print("=" * 60)
print("CONSTRUCCION DEL DATASET")
print("=" * 60)

for _, row in bls.iterrows():

    tic = row["tic_id"]

    csv_file = (
        PROCESSED_DIR /
        f"{tic}_processed.csv"
    )

    if not csv_file.exists():
        continue

    df = pd.read_csv(csv_file)

    flux = df["FLUX_NORM"].values

    mean_flux = np.mean(flux)
    median_flux = np.median(flux)
    std_flux = np.std(flux)

    mad_flux = np.median(
        np.abs(
            flux - median_flux
        )
    )

    amplitude = (
        np.max(flux)
        - np.min(flux)
    )

    rows.append({

        "tic_id": tic,

        "mean_flux": mean_flux,
        "median_flux": median_flux,
        "std_flux": std_flux,
        "mad_flux": mad_flux,
        "amplitude": amplitude,

        "best_period": row["best_period"],
        "transit_duration": row["transit_duration"],
        "transit_depth": row["transit_depth"],
        "bls_power": row["bls_power"],
        "snr_estimate": row["snr_estimate"],

        "n_points": row["n_points"],

        "label": 1
    })

dataset = pd.DataFrame(rows)

outfile = (
    OUT_DIR /
    "master_dataset.csv"
)

dataset.to_csv(
    outfile,
    index=False
)

print("\nDataset generado:")
print(dataset.head())

print(
    f"\nGuardado: {outfile}"
)