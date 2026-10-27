import pandas as pd

from pathlib import Path

print("=" * 60)
print("BUSCANDO FITS NEGATIVOS FALTANTES")
print("=" * 60)

# ----------------------------
# Archivos esperados
# ----------------------------

expected = pd.read_csv(
    "data/catalog/negative_targets_valid.csv"
)

# ----------------------------
# Archivos descargados
# ----------------------------

fits_dir = Path("data/raw/negatives")

downloaded = set()

for file in fits_dir.glob("*.fits"):

    downloaded.add(file.stem.replace("_", " "))

# ----------------------------
# Comparación
# ----------------------------

missing = []

for tic in expected["tic_id"]:

    tic = str(tic)

    if tic not in downloaded:

        missing.append(tic)

# ----------------------------
# Guardar resultado
# ----------------------------

missing_df = pd.DataFrame({
    "tic_id": missing,
    "label": 0
})

outfile = (
    Path("data/catalog")
    / "missing_negative_targets.csv"
)

missing_df.to_csv(
    outfile,
    index=False
)

# ----------------------------
# Resumen
# ----------------------------

print()

print("Negativos esperados :", len(expected))
print("FITS descargados    :", len(downloaded))
print("Faltantes           :", len(missing))

print()

print("Guardado en:")
print(outfile)