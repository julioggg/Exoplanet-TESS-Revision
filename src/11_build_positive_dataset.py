import pandas as pd
from pathlib import Path

print("=" * 60)
print("DATASET POSITIVO")
print("=" * 60)

stats = pd.read_csv(
    "data/features/positives/positive_statistical_features.csv"
)

bls = pd.read_csv(
    "data/features/positives/positive_bls_features.csv"
)

dataset = stats.merge(
    bls,
    on=["tic_id", "n_points"],
    how="inner"
)

dataset.to_csv(
    "data/dataset/positive_dataset.csv",
    index=False
)

print(dataset.head())
print()
print("Shape:", dataset.shape)
print()
print("Guardado:")
print("data/dataset/positive_dataset.csv")