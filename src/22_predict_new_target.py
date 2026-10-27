import pandas as pd
import joblib

MODEL = "models/random_forest_exoplanet.pkl"

print("=" * 60)
print("PREDICCION DE CANDIDATOS")
print("=" * 60)

model = joblib.load(MODEL)

df = pd.read_csv(
    "data/dataset/final_dataset.csv"
)

sample = df.iloc[[0]]

tic_id = sample["tic_id"].values[0]

X = sample.drop(
    columns=[
        "tic_id",
        "label",
        "n_points"
    ]
)

prob = model.predict_proba(X)[0][1]

pred = model.predict(X)[0]

print()
print("TIC:", tic_id)

print(
    "Probabilidad exoplaneta:",
    round(prob, 4)
)

print(
    "Prediccion:",
    "POSITIVO" if pred == 1 else "NEGATIVO"
)