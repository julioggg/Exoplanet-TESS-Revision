import pandas as pd

from sklearn.ensemble import RandomForestClassifier

print("=" * 60)
print("TABLA DE IMPORTANCIA DE VARIABLES")
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

importance = pd.DataFrame({
    "Feature": X.columns,
    "Importance": model.feature_importances_
})

importance = importance.sort_values(
    "Importance",
    ascending=False
)

print(importance)

importance.to_csv(
    "results/feature_importance_table.csv",
    index=False
)

print("\nGuardado:")
print("results/feature_importance_table.csv")