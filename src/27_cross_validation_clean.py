# src/27_cross_validation_clean.py

import pandas as pd

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import cross_val_score

DATASET = "data/dataset/final_dataset_clean.csv"

print("=" * 60)
print("VALIDACION CRUZADA - DATASET LIMPIO")
print("=" * 60)

df = pd.read_csv(DATASET)

X = df.drop(
    columns=[
        "tic_id",
        "label",
        "n_points"
    ],
    errors="ignore"
)

y = df["label"]

print("\nShape X:", X.shape)
print("Shape y:", y.shape)

model = RandomForestClassifier(
    n_estimators=200,
    random_state=42
)

cv = StratifiedKFold(
    n_splits=5,
    shuffle=True,
    random_state=42
)

scores = cross_val_score(
    model,
    X,
    y,
    cv=cv,
    scoring="f1"
)

print("\nF1 por fold:")
print(scores)

print("\nF1 promedio:")
print(round(scores.mean(), 4))

print("\nDesviacion estandar:")
print(round(scores.std(), 4))