import pandas as pd
import lightkurve as lk

from pathlib import Path

print("=" * 60)
print("VALIDANDO NEGATIVOS - MODO RECUPERABLE")
print("=" * 60)

INPUT_FILE = "data/catalog/negative_targets.csv"

OUTPUT_FILE = (
    Path("data/catalog")
    / "negative_targets_valid.csv"
)

neg = pd.read_csv(INPUT_FILE)

# -----------------------------
# cargar progreso previo
# -----------------------------

if OUTPUT_FILE.exists():

    old = pd.read_csv(OUTPUT_FILE)

    valid = set(
        old["tic_id"].astype(str)
    )

    print(
        f"Ya existen {len(valid)} "
        f"negativos validados"
    )

else:

    valid = set()

# -----------------------------
# validacion
# -----------------------------

for i, tic in enumerate(
    neg["tic_id"],
    start=1
):

    tic = str(tic)

    if tic in valid:

        continue

    print(
        f"[{i}/{len(neg)}] {tic}"
    )

    try:

        search = lk.search_lightcurve(
            tic,
            mission="TESS"
        )

        if len(search) > 0:

            print(
                "  ✓ Observaciones encontradas"
            )

            valid.add(tic)

            pd.DataFrame({
                "tic_id": sorted(valid),
                "label": 0
            }).to_csv(
                OUTPUT_FILE,
                index=False
            )

        else:

            print(
                "  ✗ Sin observaciones"
            )

    except Exception as e:

        print("  ERROR:", e)

# -----------------------------
# resumen
# -----------------------------

print()
print("=" * 60)
print("RESUMEN")
print("=" * 60)

print(
    "Negativos válidos:",
    len(valid)
)

print(
    "Guardado:",
    OUTPUT_FILE
)