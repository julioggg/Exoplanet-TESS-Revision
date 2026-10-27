import numpy as np
import pandas as pd

from pathlib import Path
from astropy.timeseries import BoxLeastSquares

INPUT_DIR = Path("data/processed/negatives")
OUT_DIR = Path("data/features/negatives")

OUT_DIR.mkdir(parents=True, exist_ok=True)

results = []

print("=" * 60)
print("BLS NEGATIVOS")
print("=" * 60)

for file in INPUT_DIR.glob("*_processed.csv"):

    print(f"\nProcesando {file.name}")

    try:

        df = pd.read_csv(file)

        time = df["TIME"].values
        flux = df["FLUX_NORM"].values

        bls = BoxLeastSquares(time, flux)

        periods = np.linspace(0.5, 30, 500)

        durations = np.array([
            0.05,
            0.10,
            0.20,
            0.25
        ])

        power = bls.power(
            periods,
            durations
        )

        idx = np.argmax(power.power)

        best_period = power.period[idx]
        best_duration = power.duration[idx]
        best_depth = power.depth[idx]
        best_power = power.power[idx]

        snr = (
            abs(best_depth)
            / np.std(flux)
        )

        tic_id = file.stem.replace(
            "_processed",
            ""
        )

        results.append({
            "tic_id": tic_id,
            "best_period": best_period,
            "transit_duration": best_duration,
            "transit_depth": best_depth,
            "bls_power": best_power,
            "snr_estimate": snr,
            "n_points": len(df),
            "label": 0
        })

        print(
            f"Periodo={best_period:.4f} d | "
            f"Power={best_power:.4f}"
        )

    except Exception as e:

        print("ERROR:", e)

results = pd.DataFrame(results)

outfile = (
    OUT_DIR /
    "negative_bls_features.csv"
)

results.to_csv(
    outfile,
    index=False
)

print("\n")
print("=" * 60)
print("RESULTADOS")
print("=" * 60)

print(results)

print("\nGuardado:", outfile)