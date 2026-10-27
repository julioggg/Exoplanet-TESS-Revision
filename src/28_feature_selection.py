import pandas as pd

DATASET = "data/dataset/final_dataset_clean.csv"

print("=" * 60)
print("SELECCION DE VARIABLES")
print("=" * 60)

df = pd.read_csv(DATASET)

features = [
    "mean_flux",
    "median_flux",
    "std_flux",
    "mad_flux",
    "amplitude",
    "best_period",
    "transit_duration",
    "transit_depth",
    "bls_power",
    "snr_estimate"
]

corr = df[features].corr()

print("\nCorrelacion absoluta:")

print(
    corr.abs().round(3)
)