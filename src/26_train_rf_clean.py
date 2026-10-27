import pandas as pd

from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier

from sklearn.metrics import (
    accuracy_score,
    precision_score,
    recall_score,
    f1_score,
    confusion_matrix
)

print("=" * 60)
print("RANDOM FOREST - DATASET LIMPIO")
print("=" * 60)

df = pd.read_csv(
    "data/dataset/final_dataset_clean.csv"
)

# Eliminamos tic_id, label y n_points
X = df.drop(
    columns=[
        "tic_id",
        "label",
        "n_points"
    ]
)

y = df["label"]

print("\nShape X:", X.shape)
print("Shape y:", y.shape)

X_train, X_test, y_train, y_test = train_test_split(
    X,
    y,
    test_size=0.30,
    random_state=42,
    stratify=y
)

print("\nTrain:", len(X_train))
print("Test:", len(X_test))

model = RandomForestClassifier(
    n_estimators=200,
    random_state=42
)

model.fit(X_train, y_train)

y_pred = model.predict(X_test)

print("\n" + "=" * 60)
print("METRICAS")
print("=" * 60)

print(
    "Accuracy:",
    round(
        accuracy_score(
            y_test,
            y_pred
        ),
        4
    )
)

print(
    "Precision:",
    round(
        precision_score(
            y_test,
            y_pred,
            zero_division=0
        ),
        4
    )
)

print(
    "Recall:",
    round(
        recall_score(
            y_test,
            y_pred,
            zero_division=0
        ),
        4
    )
)

print(
    "F1:",
    round(
        f1_score(
            y_test,
            y_pred,
            zero_division=0
        ),
        4
    )
)

print("\nMatriz de confusion:")

print(
    confusion_matrix(
        y_test,
        y_pred
    )
)

print("\n" + "=" * 60)
print("IMPORTANCIA DE VARIABLES")
print("=" * 60)

importance = pd.DataFrame({
    "feature": X.columns,
    "importance": model.feature_importances_
})

importance = importance.sort_values(
    "importance",
    ascending=False
)

print(importance)