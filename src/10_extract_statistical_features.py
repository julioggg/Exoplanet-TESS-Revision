import numpy as np
import pandas as pd

from pathlib import Path

DATA_DIR = Path(
    "data/processed/positives"
)

OUT_DIR = Path(
    "data/features/positives"
)

OUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)

files = sorted(
    DATA_DIR.glob("*_processed.csv")
)

rows = []

print("=" * 60)
print("FEATURES ESTADISTICAS")
print("=" * 60)

for file in files:

    df = pd.read_csv(file)

    flux = df["FLUX_NORM"].values

    rows.append({

        "tic_id":
            file.name.replace(
                "_processed.csv",
                ""
            ),

        "mean_flux":
            np.mean(flux),

        "median_flux":
            np.median(flux),

        "std_flux":
            np.std(flux),

        "mad_flux":
            np.median(
                np.abs(
                    flux -
                    np.median(flux)
                )
            ),

        "amplitude":
            np.max(flux) -
            np.min(flux),

        "n_points":
            len(flux)
    })

out = pd.DataFrame(rows)

outfile = (
    OUT_DIR /
    "positive_statistical_features.csv"
)

out.to_csv(
    outfile,
    index=False
)

print(out.head())
print("\nGuardado:", outfile)