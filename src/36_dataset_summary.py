import pandas as pd

print("=" * 60)
print("RESUMEN DEL ESTUDIO")
print("=" * 60)

summary = pd.DataFrame({

    "Metric": [
        "Positivos",
        "Negativos",
        "Total muestras",
        "Variables finales",
        "Accuracy",
        "Precision",
        "Recall",
        "F1",
        "CV_F1",
        "CV_STD"
    ],

    "Value": [
        20,
        24,
        44,
        10,
        0.9286,
        1.0000,
        0.8333,
        0.9091,
        0.8833,
        0.1453
    ]
})

print(summary)

summary.to_csv(
    "results/study_summary.csv",
    index=False
)

print("\nGuardado:")
print("results/study_summary.csv")