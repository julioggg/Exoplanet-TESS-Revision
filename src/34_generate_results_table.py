import pandas as pd

print("=" * 60)
print("TABLA DE RESULTADOS")
print("=" * 60)

results = pd.DataFrame({

    "Modelo": [
        "RF completo",
        "RF limpio",
        "RF reducido"
    ],

    "Accuracy": [
        0.8421,
        0.9286,
        0.8571
    ],

    "Precision": [
        0.7143,
        1.0000,
        0.8333
    ],

    "Recall": [
        0.8333,
        0.8333,
        0.8333
    ],

    "F1": [
        0.7692,
        0.9091,
        0.8333
    ],

    "CV_F1": [
        0.8044,
        0.8833,
        0.9029
    ]
})

print(results)

results.to_csv(
    "results/model_comparison.csv",
    index=False
)

print("\nGuardado:")
print("results/model_comparison.csv")