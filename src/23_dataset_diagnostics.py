import pandas as pd

DATASET = "data/dataset/final_dataset.csv"

print("=" * 60)
print("DIAGNOSTICO DEL DATASET")
print("=" * 60)

df = pd.read_csv(DATASET)

print("\nShape:")
print(df.shape)

print("\nDistribucion de clases:")
print(df["label"].value_counts())

print("\nPromedios por clase:")
print(
    df.groupby("label")
      .mean(numeric_only=True)
      .T
)

print("\nDesviaciones estandar por clase:")
print(
    df.groupby("label")
      .std(numeric_only=True)
      .T
)