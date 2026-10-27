import pandas as pd
import lightkurve as lk
from pathlib import Path

TARGET_POSITIVES = 200

print("=" * 60)
print(f"DESCARGA DE {TARGET_POSITIVES} POSITIVOS")
print("=" * 60)

targets = pd.read_csv("data/catalog/targets.csv")

outdir = Path("data/raw/positives")
outdir.mkdir(parents=True, exist_ok=True)

success = 0
skipped = 0

for tic in targets["tic_id"]:                          # ✅ for correcto

    if success >= TARGET_POSITIVES:                    # ✅ indentado dentro del for
        break

    outfile = outdir / f"{tic.replace(' ', '_')}.fits" # ✅ dentro del for

    if outfile.exists():                               # ✅ dentro del for
        skipped += 1
        print(f"\n[SKIP] {tic}")
        continue

    print(f"\n[{success + 1}/{TARGET_POSITIVES}] {tic}")

    try:                                               # ✅ try dentro del for
        search = lk.search_lightcurve(tic, mission="TESS")

        if len(search) == 0:
            print("  Sin observaciones")
            continue

        lc = search[0].download()

        if lc is None:
            print("  Descarga nula")
            continue

        lc.to_fits(path=str(outfile), overwrite=True)
        print("  Guardado")
        success += 1

    except Exception as e:
        print("  ERROR:", e)

# ✅ Resumen fuera del for (sin indentación)
print()
print("=" * 60)
print("RESUMEN")
print("=" * 60)
print("Positivos descargados:", success)
print("Archivos existentes:  ", skipped)
print("Total objetivo:       ", TARGET_POSITIVES)