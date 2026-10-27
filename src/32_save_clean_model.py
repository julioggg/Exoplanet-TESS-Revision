import pandas as pd
import joblib

from sklearn.ensemble import RandomForestClassifier

print("=" * 60)
print("GUARDANDO MODELO FINAL LIMPIO")
print("=" * 60)

df = pd.read_csv(
    "data/dataset/final_dataset_clean.csv"
)

X = df.drop(
    columns=[
        "tic_id",
        "label",
        "n_points"
    ],
    errors="ignore"
)

y = df["label"]

model = RandomForestClassifier(
    n_estimators=200,
    random_state=42
)

model.fit(X, y)

outfile = (
    "models/random_forest_clean.pkl"
)

joblib.dump(
    model,
    outfile
)

print("\nModelo guardado:")
print(outfile)