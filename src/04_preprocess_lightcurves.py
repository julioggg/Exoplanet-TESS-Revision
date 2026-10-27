import numpy as np
import pandas as pd

from pathlib import Path
from astropy.io import fits
from astropy.stats import sigma_clip

RAW_DIR = Path("data/raw/fits")
OUT_DIR = Path("data/processed")

OUT_DIR.mkdir(parents=True, exist_ok=True)

files = list(RAW_DIR.glob("*.fits"))

print("=" * 60)
print("PREPROCESAMIENTO DE CURVAS DE LUZ")
print("=" * 60)

for file in files:

    print(f"\nProcesando {file.name}")

    try:

        hdul = fits.open(file)
        
        data = hdul[1].data

        
        time = np.asarray(data["TIME"]).astype(np.float64)

        
        flux = np.asarray(data["FLUX"]).astype(np.float64)

        
        flux_err = np.asarray(data["FLUX_ERR"]).astype(np.float64)

        
        quality = np.asarray(data["QUALITY"]).astype(np.int64)

        

        df = pd.DataFrame({
            "TIME": time,
            "FLUX": flux,
            "FLUX_ERR": flux_err,
            "QUALITY": quality
        })

        # ------------------------------------------------
        # 1. Filtrar QUALITY
        # ------------------------------------------------

        df = df[df["QUALITY"] == 0]

        # ------------------------------------------------
        # 2. Eliminar NaN
        # ------------------------------------------------

        df = df.dropna()

        # ------------------------------------------------
        # 3. Normalización
        # ------------------------------------------------

        median_flux = np.median(df["FLUX"])

        df["FLUX_NORM"] = df["FLUX"] / median_flux

        # ------------------------------------------------
        # 4. Sigma clipping
        # ------------------------------------------------

        clipped = sigma_clip(
            df["FLUX_NORM"],
            sigma=5,
            maxiters=3
        )

        df = df[~clipped.mask]

        # ------------------------------------------------
        # 5. Guardar
        # ------------------------------------------------

        outfile = OUT_DIR / f"{file.stem}_processed.csv"

        df.to_csv(outfile, index=False)

        print(f"Guardado: {outfile}")
        print(f"Observaciones finales: {len(df)}")

    except Exception as e:

        print(f"ERROR: {e}")

print("\nProceso completado.")