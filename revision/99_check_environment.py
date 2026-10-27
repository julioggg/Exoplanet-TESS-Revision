#!/usr/bin/env python3
from __future__ import annotations

import importlib
import platform
import sys

PACKAGES = [
    ("numpy", "numpy"),
    ("pandas", "pandas"),
    ("scipy", "scipy"),
    ("sklearn", "scikit-learn"),
    ("astropy", "astropy"),
    ("pytransit", "PyTransit"),
    ("lightkurve", "lightkurve"),
    ("astroquery", "astroquery"),
    ("matplotlib", "matplotlib"),
    ("shap", "shap"),
    ("joblib", "joblib"),
    ("yaml", "PyYAML"),
]

print("=" * 72)
print("EXOPLANET REVISION ENVIRONMENT CHECK")
print("=" * 72)
print("Python:", sys.version.replace("\n", " "))
print("Executable:", sys.executable)
print("Platform:", platform.platform())
print()

missing = []
for module_name, package_name in PACKAGES:
    try:
        mod = importlib.import_module(module_name)
        version = getattr(mod, "__version__", "unknown")
        print(f"[OK]      {package_name:<16} {version}")
    except Exception as exc:
        missing.append(package_name)
        print(f"[MISSING] {package_name:<16} {exc}")

print()
if missing:
    print("Missing packages:", ", ".join(missing))
    print("Environment is NOT ready.")
    raise SystemExit(1)

print("Environment is ready.")
print("=" * 72)
