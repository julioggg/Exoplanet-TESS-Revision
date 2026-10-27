import pandas as pd
import numpy as np

from pathlib import Path

print("=" * 60)
print("GENERANDO OBJETIVOS NEGATIVOS")
print("=" * 60)

TARGET_NEGATIVES = 2000

positive = pd.read_csv(
    "data/catalog/targets.csv"
)

positive_ids = set(
    positive["tic_id"].astype(str)
)

print("Positivos:", len(positive_ids))

np.random.seed(42)

negative_ids = set()

while len(negative_ids) < TARGET_NEGATIVES:

    tic = f"TIC {np.random.randint(10000000, 999999999)}"

    if tic not in positive_ids:
        negative_ids.add(tic)

negative = pd.DataFrame({
    "tic_id": sorted(negative_ids),
    "label": 0
})

outdir = Path("data/catalog")
outdir.mkdir(parents=True, exist_ok=True)

outfile = outdir / "negative_targets.csv"

negative.to_csv(
    outfile,
    index=False
)

print()
print("Negativos generados:", len(negative))
print("Guardado:", outfile)