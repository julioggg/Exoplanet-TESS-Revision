import pandas as pd
from pathlib import Path

OUTDIR = Path("data/catalog")
OUTDIR.mkdir(parents=True, exist_ok=True)

print("="*60)
print("DESCARGANDO PLANETAS CONFIRMADOS")
print("="*60)

query = (
    "https://exoplanetarchive.ipac.caltech.edu/TAP/sync?"
    "query=select+pl_name,hostname+from+pscomppars"
    "&format=csv"
)

df = pd.read_csv(query)

outfile = OUTDIR / "confirmed_planets_v2.csv"

df.to_csv(outfile, index=False)

print(f"Archivo guardado: {outfile}")
print(f"Registros: {len(df)}")