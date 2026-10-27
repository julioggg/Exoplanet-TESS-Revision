import pandas as pd

print("=" * 60)
print("REDUCCION DE VARIABLES")
print("=" * 60)

df = pd.read_csv(
    "data/dataset/final_dataset_clean.csv"
)

drop_cols = [
    "mad_flux",
    "amplitude",
    "transit_depth"
]

df_reduced = df.drop(
    columns=drop_cols
)

print("\nVariables eliminadas:")

for col in drop_cols:
    print("-", col)

print("\nShape original:")
print(df.shape)

print("\nShape reducido:")
print(df_reduced.shape)

outfile = (
    "data/dataset/final_dataset_reduced.csv"
)

df_reduced.to_csv(
    outfile,
    index=False
)

print("\nGuardado:")
print(outfile)