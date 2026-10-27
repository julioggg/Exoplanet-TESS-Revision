import numpy as np
import pandas as pd

from pathlib import Path
from astropy.timeseries import BoxLeastSquares

DATA_DIR = Path("data/processed/positives")
OUT_DIR = Path("data/features/positives")

OUT_DIR.mkdir(parents=True, exist_ok=True)

files = sorted(DATA_DIR.glob("*_processed.csv"))

results = []

print("=" * 60)
print("BLS MASIVO - POSITIVOS")
print("=" * 60)

for i, file in enumerate(files, start=1):

    print(f"\n[{i}/{len(files)}] {file.name}")

    try:

        df = pd.read_csv(file)

        time = df["TIME"].values
        flux = df["FLUX_NORM"].values

        bls = BoxLeastSquares(time, flux)

        periods = np.linspace(
            0.5,
            30,
            5000
        )

        durations = np.array([
            0.05,
            0.10,
            0.15,
            0.20,
            0.25
        ])

        result = bls.power(
            periods,
            durations
        )

        idx = np.argmax(result.power)

        best_period = result.period[idx]
        best_duration = result.duration[idx]
        best_power = result.power[idx]
        best_depth = result.depth[idx]

        noise = np.std(flux)

        if noise == 0:
            noise = 1e-8

        snr = np.abs(best_depth) / noise
        
        results.append({
            "tic_id": file.name.replace(
                "_processed.csv",
                ""
            ),
            "best_period": best_period,
            "transit_duration": best_duration,
            "transit_depth": best_depth,
            "bls_power": best_power,
            "snr_estimate": snr,
            "n_points": len(df),
            "label": 1
        })
        pd.DataFrame(results).to_csv(
            OUT_DIR / "positive_bls_features.csv",
            index=False
        )
        print(
            f"Periodo={best_period:.4f} d "
            f"| Power={best_power:.4f}"
        )

    except Exception as e:

        print("ERROR:", e)

out = pd.DataFrame(results)

outfile = (
    OUT_DIR /
    "positive_bls_features.csv"
)

out.to_csv(outfile, index=False)

print("\n")
print(out.head())

print("\nGuardado:")
print(outfile)