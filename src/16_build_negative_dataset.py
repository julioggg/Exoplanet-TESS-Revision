import pandas as pd

STAT_FILE = (
    "data/features/negatives/"
    "negative_statistical_features.csv"
)

BLS_FILE = (
    "data/features/negatives/"
    "negative_bls_features.csv"
)

OUTPUT_FILE = (
    "data/dataset/negative_dataset.csv"
)

print("=" * 60)
print("DATASET NEGATIVO")
print("=" * 60)

stat_df = pd.read_csv(STAT_FILE)

bls_df = pd.read_csv(BLS_FILE)

dataset = pd.merge(
    stat_df,
    bls_df[
        [
            "tic_id",
            "best_period",
            "transit_duration",
            "transit_depth",
            "bls_power",
            "snr_estimate",
            "label"
        ]
    ],
    on="tic_id"
)

print(dataset.head())

print("\nShape:", dataset.shape)

dataset.to_csv(
    OUTPUT_FILE,
    index=False
)

print("\nGuardado:")
print(OUTPUT_FILE)