import pandas as pd

df = pd.read_csv(
    "data/dataset/final_dataset.csv"
)

neg = df[df["label"] == 0]

print("=" * 60)
print("OUTLIERS NEGATIVOS")
print("=" * 60)

print("\nAmplitude:")
print(
    neg[
        ["tic_id", "amplitude"]
    ].sort_values(
        "amplitude",
        ascending=False
    ).head(10)
)

print("\nBLS Power:")
print(
    neg[
        ["tic_id", "bls_power"]
    ].sort_values(
        "bls_power",
        ascending=False
    ).head(10)
)

print("\nTransit Depth:")
print(
    neg[
        ["tic_id", "transit_depth"]
    ].sort_values(
        "transit_depth",
        ascending=False
    ).head(10)
)

print("\nSTD Flux:")
print(
    neg[
        ["tic_id", "std_flux"]
    ].sort_values(
        "std_flux",
        ascending=False
    ).head(10)
)