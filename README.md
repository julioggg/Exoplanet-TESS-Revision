# Exoplanet TESS Revision

Reproducible revision workspace for the manuscript on automated exoplanet transit detection in TESS light curves.

## Workflow

- `main`: stable baseline and accepted revision checkpoints.
- `reviewer-revision`: active reviewer-response development.
- CEDIA HPC: raw/processed TESS data, FITS files, synthetic injections, and computational runs.
- GitHub: source code, small tabular datasets, configurations, figures, compact results, and revision documentation.

Large raw/processed astronomical files and trained binary models are intentionally excluded from version control.

## Revision priorities

1. Audit labels, TIC-level grouping, observational coverage, and feature definitions.
2. Reproduce the submitted Random Forest baseline.
3. Quantify the effect of `n_points` with controlled ablations.
4. Add grouped/repeated cross-validation and bootstrap uncertainty.
5. Validate positive transit visibility and known-period consistency.
6. Add synthetic transit injection/recovery and false-positive challenge experiments.
7. Recompute feature importance and SHAP explanations.
8. Regenerate manuscript tables/figures and reviewer-response evidence.
