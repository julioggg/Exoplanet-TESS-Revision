import numpy as np
import pandas as pd

from pathlib import Path
from astropy.io import fits
from astropy.stats import sigma_clip

RAW_DIR = Path("data/raw/negatives")
OUT_DIR = Path("data/processed/negatives")

OUT_DIR.mkdir(parents=True, exist_ok=True)

files = list(RAW_DIR.glob("*.fits"))

print("=" * 60)
print("PREPROCESAMIENTO NEGATIVOS")
print("=" * 60)

for file in files:

    print(f"\nProcesando {file.name}")

    try:

        hdul = fits.open(file)

        data = hdul[1].data

        time = np.array(data["TIME"], dtype=np.float64)
        flux = np.array(data["FLUX"], dtype=np.float64)
        flux_err = np.array(data["FLUX_ERR"], dtype=np.float64)
        quality = np.array(data["QUALITY"], dtype=np.int64)

        df = pd.DataFrame({
            "TIME": time,
            "FLUX": flux,
            "FLUX_ERR": flux_err,
            "QUALITY": quality
        })

        df = df[df["QUALITY"] == 0]
        df = df.dropna()

        median_flux = np.median(df["FLUX"])

        df["FLUX_NORM"] = df["FLUX"] / median_flux

        clipped = sigma_clip(
            df["FLUX_NORM"],
            sigma=5,
            maxiters=3
        )

        df = df[~clipped.mask]

        outfile = (
            OUT_DIR /
            f"{file.stem}_processed.csv"
        )

        df.to_csv(outfile, index=False)

        print("Guardado:", outfile)
        print("Observaciones:", len(df))

    except Exception as e:

        print("ERROR:", e)

print("\nProceso completado.")