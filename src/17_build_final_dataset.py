import pandas as pd

POSITIVE_FILE = "data/dataset/positive_dataset.csv"
NEGATIVE_FILE = "data/dataset/negative_dataset.csv"

OUTPUT_FILE = "data/dataset/final_dataset.csv"

print("=" * 60)
print("DATASET FINAL")
print("=" * 60)

positive_df = pd.read_csv(POSITIVE_FILE)
negative_df = pd.read_csv(NEGATIVE_FILE)

dataset = pd.concat(
    [positive_df, negative_df],
    ignore_index=True
)
print("\nDataset information")
print(dataset.info())

print("\nMissing values")
print(dataset.isna().sum())

print("\nDuplicated rows")
print(dataset.duplicated().sum())

print(dataset.head())

print("\nDistribución de clases:")
print(dataset["label"].value_counts())

print("\nShape:")
print(dataset.shape)

dataset.to_csv(
    OUTPUT_FILE,
    index=False
)

print("\nGuardado:")
print(OUTPUT_FILE)