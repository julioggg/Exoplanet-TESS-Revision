import pandas as pd
import joblib

from pathlib import Path
from sklearn.ensemble import RandomForestClassifier

DATASET = "data/dataset/final_dataset.csv"

print("=" * 60)
print("GUARDANDO MODELO FINAL")
print("=" * 60)

df = pd.read_csv(DATASET)

X = df.drop(
    columns=[
        "tic_id",
        "label",
        "n_points"
    ]
)

y = df["label"]

model = RandomForestClassifier(
    n_estimators=200,
    random_state=42
)

model.fit(X, y)

Path("models").mkdir(
    exist_ok=True
)

joblib.dump(
    model,
    "models/random_forest_exoplanet.pkl"
)

print()
print("Modelo guardado:")
print("models/random_forest_exoplanet.pkl")