import pandas as pd

print("=" * 60)
print("LIMPIEZA DE OUTLIERS")
print("=" * 60)

df = pd.read_csv(
    "data/dataset/final_dataset.csv"
)

neg = df[df["label"] == 0].copy()
pos = df[df["label"] == 1].copy()

features = [
    "amplitude",
    "std_flux",
    "transit_depth",
    "bls_power"
]

print("\nNegativos iniciales:", len(neg))

for col in features:

    q1 = neg[col].quantile(0.25)
    q3 = neg[col].quantile(0.75)

    iqr = q3 - q1

    lower = q1 - 1.5 * iqr
    upper = q3 + 1.5 * iqr

    before = len(neg)

    neg = neg[
        (neg[col] >= lower) &
        (neg[col] <= upper)
    ]

    removed = before - len(neg)

    print(
        f"{col}: eliminados {removed}"
    )

clean = pd.concat(
    [pos, neg],
    ignore_index=True
)

print("\nNegativos finales:", len(neg))
print("Positivos:", len(pos))
print("Dataset final:", len(clean))

outfile = (
    "data/dataset/final_dataset_clean.csv"
)

clean.to_csv(
    outfile,
    index=False
)

print("\nGuardado:")
print(outfile)