import numpy as np
import pandas as pd

from pathlib import Path
from astropy.io import fits
from astropy.stats import sigma_clip

RAW_DIR = Path("data/raw/positives")
OUT_DIR = Path("data/processed/positives")

OUT_DIR.mkdir(parents=True, exist_ok=True)

files = sorted(RAW_DIR.glob("*.fits"))

print("=" * 60)
print("PREPROCESAMIENTO MASIVO POSITIVOS")
print("=" * 60)

total = len(files)

for i, file in enumerate(files, start=1):

    print(f"\n[{i}/{total}] {file.name}")

    try:

        with fits.open(file) as hdul:          # ✅ cierra el archivo automáticamente
            data = hdul[1].data

            time     = np.array(data["TIME"]).astype(np.float64)
            flux     = np.array(data["FLUX"]).astype(np.float64)
            flux_err = np.array(data["FLUX_ERR"]).astype(np.float64)
            quality  = np.array(data["QUALITY"]).astype(np.int64)

        df = pd.DataFrame({
            "TIME":     time,
            "FLUX":     flux,
            "FLUX_ERR": flux_err,
            "QUALITY":  quality
        })

        df = df[df["QUALITY"] == 0]
        df = df.dropna()

        if len(df) < 500:                      # ✅ continue DENTRO del if
            print("  Muy pocas observaciones, omitido")
            continue

        median_flux = np.median(df["FLUX"])

        df["FLUX_NORM"] = df["FLUX"] / median_flux

        clipped = sigma_clip(
            df["FLUX_NORM"],
            sigma=5,
            maxiters=3
        )

        df = df[~clipped.mask]

        outfile = OUT_DIR / f"{file.stem}_processed.csv"

        df.to_csv(outfile, index=False)

        print(f"  Observaciones finales: {len(df)}")

    except Exception as e:

        print("  ERROR:", e)

print("\nProceso completado.")