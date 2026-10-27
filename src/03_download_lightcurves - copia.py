import pandas as pd
from pathlib import Path
import lightkurve as lk

# =====================================================
# Cargar lista de objetivos
# =====================================================

targets = pd.read_csv("data/catalog/targets.csv")

# =====================================================
# Carpeta de salida
# =====================================================

outdir = Path("data/raw/fits")
outdir.mkdir(parents=True, exist_ok=True)

# =====================================================
# Descarga de prueba
# =====================================================

N_TEST = 5

print("=" * 60)
print("DESCARGA DE CURVAS DE LUZ TESS")
print("=" * 60)

for i, tic in enumerate(targets["tic_id"].head(N_TEST), start=1):

    try:

        tic_clean = str(tic).replace("TIC", "").strip()

        print(f"\n[{i}/{N_TEST}] Buscando TIC {tic_clean}")

        search = lk.search_lightcurve(
            f"TIC {tic_clean}",
            mission="TESS"
        )

        if len(search) == 0:

            print("No se encontraron curvas.")

            continue

        print(f"Encontradas {len(search)} curvas")

        lc = search.download()

        if lc is None:

            print("Descarga vacía")

            continue

        outfile = outdir / f"TIC_{tic_clean}.fits"

        lc.to_fits(
            path=str(outfile),
            overwrite=True
        )

        print(f"Guardado: {outfile}")

    except Exception as e:

        print(f"ERROR TIC {tic}: {e}")

print("\nProceso finalizado.")