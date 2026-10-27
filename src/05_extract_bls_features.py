import numpy as np
import pandas as pd

from pathlib import Path
from astropy.timeseries import BoxLeastSquares

INPUT_DIR = Path("data/processed")

OUTPUT_DIR = Path("data/features")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

results = []

files = list(INPUT_DIR.glob("*_processed.csv"))

print("=" * 60)
print("EXTRACCION DE CARACTERISTICAS BLS")
print("=" * 60)

for file in files:

    print(f"\nProcesando {file.name}")

    try:

        df = pd.read_csv(file)

        time = df["TIME"].values
        flux = df["FLUX_NORM"].values

        mask = np.isfinite(time) & np.isfinite(flux)

        time = time[mask]
        flux = flux[mask]

        if len(time) < 1000:

            print("Muy pocos puntos")
            continue

        model = BoxLeastSquares(
            time,
            flux
        )

        periods = np.linspace(
             1.0,
             30.0,
             5000
        )

        durations = np.array([
             0.05,
             0.10,
             0.15,
             0.20,
             0.25
        ])

        power = model.power(
            periods,
            durations
        )

        idx = np.argmax(
            power.power
        )

        best_period = float(
            power.period[idx]
        )

        best_duration = float(
            power.duration[idx]
        )

        best_power = float(
            power.power[idx]
        )

        best_depth = float(
            power.depth[idx]
        )

        flux_std = np.std(flux)

        snr = (
            abs(best_depth)
            / flux_std
            if flux_std > 0
            else 0
        )

        tic_id = (
            file.stem
            .replace("_processed", "")
        )

        results.append([
            tic_id,
            best_period,
            best_duration,
            best_depth,
            best_power,
            snr,
            len(flux)
        ])

        print(
            f"Periodo={best_period:.4f} d | "
            f"Power={best_power:.4f}"
        )

    except Exception as e:

        print(
            f"ERROR: {e}"
        )

features = pd.DataFrame(
    results,
    columns=[
        "tic_id",
        "best_period",
        "transit_duration",
        "transit_depth",
        "bls_power",
        "snr_estimate",
        "n_points"
    ]
)

outfile = (
    OUTPUT_DIR /
    "bls_features.csv"
)

features.to_csv(
    outfile,
    index=False
)

print("\n")
print("=" * 60)
print("RESULTADOS")
print("=" * 60)

print(features.head())

print(
    f"\nGuardado: {outfile}"
)