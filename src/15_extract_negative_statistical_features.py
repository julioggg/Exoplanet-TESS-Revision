import numpy as np
import pandas as pd

from pathlib import Path

INPUT_DIR = Path("data/processed/negatives")
OUTPUT_FILE = Path(
    "data/features/negatives/negative_statistical_features.csv"
)

records = []

print("=" * 60)
print("FEATURES ESTADISTICAS NEGATIVOS")
print("=" * 60)

for file in INPUT_DIR.glob("*_processed.csv"):

    df = pd.read_csv(file)

    flux = df["FLUX_NORM"]

    records.append({
        "tic_id": file.stem.replace("_processed", ""),
        "mean_flux": flux.mean(),
        "median_flux": flux.median(),
        "std_flux": flux.std(),
        "mad_flux": np.median(
            np.abs(flux - np.median(flux))
        ),
        "amplitude": flux.max() - flux.min(),
        "n_points": len(flux)
    })

result = pd.DataFrame(records)

print(result.head())

OUTPUT_FILE.parent.mkdir(
    parents=True,
    exist_ok=True
)

result.to_csv(
    OUTPUT_FILE,
    index=False
)

print("\nGuardado:")
print(OUTPUT_FILE)