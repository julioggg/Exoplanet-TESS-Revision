import pandas as pd
from pathlib import Path

OUTDIR = Path("data/catalog")
OUTDIR.mkdir(parents=True, exist_ok=True)

print("=" * 60)
print("DESCARGANDO PLANETAS CONFIRMADOS")
print("=" * 60)

confirmed_url = (
    "https://exoplanetarchive.ipac.caltech.edu/TAP/sync?"
    "query=select+pl_name,hostname,tic_id,disc_facility+from+pscomppars"
    "&format=csv"
)

confirmed = pd.read_csv(confirmed_url)

confirmed.to_csv(
    OUTDIR / "confirmed_planets.csv",
    index=False
)

print(f"Planetas descargados: {len(confirmed)}")
print("Proceso completado.")
