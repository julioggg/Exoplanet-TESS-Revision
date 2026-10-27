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
from sklearn.metrics import classification_report

DATASET = "data/dataset/final_dataset.csv"

print("=" * 60)
print("RANDOM FOREST - DETECCION DE EXOPLANETAS")
print("=" * 60)

df = pd.read_csv(DATASET)

X = df.drop(
    columns=[
        "tic_id",
        "label"
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
    n_estimators=500,
    class_weight="balanced",
    random_state=42,
    n_jobs=-1
)

model.fit(X_train, y_train)

y_pred = model.predict(X_test)
y_prob = model.predict_proba(X_test)[:, 1]

print("\n" + "=" * 60)
print("CLASSIFICATION REPORT")
print("=" * 60)

print(
    classification_report(
        y_test,
        y_pred,
        digits=4
    )
)
results = pd.DataFrame({

    "True": y_test.values,

    "Predicted": y_pred,

    "Probability": y_prob

})

results.to_csv(
    "results/test_predictions.csv",
    index=False
)

print("\nPredictions saved:")
print("results/test_predictions.csv")
print("\n" + "=" * 60)
print("METRICAS")
print("=" * 60)

print(
    "Accuracy:",
    round(
        accuracy_score(y_test, y_pred),
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
cm = confusion_matrix(
    y_test,
    y_pred
)

print(cm)

pd.DataFrame(
    cm,
    index=["True Negative", "True Positive"],
    columns=["Predicted Negative", "Predicted Positive"]
).to_csv(
    "results/confusion_matrix.csv"
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
importance.to_csv(
    "results/feature_importance.csv",
    index=False
)

print("\nFeature importance saved:")
print("results/feature_importance.csv")