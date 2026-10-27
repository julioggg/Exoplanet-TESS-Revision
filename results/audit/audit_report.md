# Original dataset and pipeline audit

This report describes the submitted baseline. It does not modify the original dataset.

## Dataset

- Rows: 369
- Unique TIC IDs: 369
- Duplicate TIC IDs: 0
- Positive rows: 220
- Negative rows: 149

- n_points negative: median=1029.0, mean=4781.4, IQR=[935.0, 1196.0], range=[651, 91141]
- n_points positive: median=16853.0, mean=21417.1, IQR=[14446.2, 18420.5], range=[695, 115324]

## Source-code asymmetries detected

- Positive download uses only first search result: True
- Negative download uses download_all(): True
- Negative download stitches products: True
- Positive preprocessing enforces >=500 points: True
- Negative preprocessing enforces same threshold: False
- Positive BLS grid: {'period_min': 0.5, 'period_max': 30.0, 'period_grid_points': 5000, 'durations': [0.05, 0.1, 0.15, 0.2, 0.25]}
- Negative BLS grid: {'period_min': 0.5, 'period_max': 30.0, 'period_grid_points': 500, 'durations': [0.05, 0.1, 0.2, 0.25]}
- BLS grids identical: False
- Positive std uses NumPy ddof=0: True
- Negative std uses pandas ddof=1: True
- Existing no-n_points comparison is controlled: False

## Median-flux diagnostic

- Overall median(median_flux): 1
- SD(median_flux): 0.0507857980843
- Fraction within 1e-3 of unity: 0.9783

## Processed-curve observational coverage

- Processed curves audited: 369
- Curves linked to final dataset: 369
- n_points_processed:
  - negative: median=1029, IQR=[935, 1196], n=149
  - positive: median=16853, IQR=[14446.2, 18420.5], n=220
- baseline_days:
  - negative: median=26.6879, IQR=[24.7921, 27.8958], n=149
  - positive: median=26.4669, IQR=[24.4378, 26.8502], n=220
- median_cadence_minutes:
  - negative: median=30, IQR=[29.9996, 30.0003], n=149
  - positive: median=2.00001, IQR=[1.99998, 2.00002], n=220

## Interpretation

The submitted dataset must be treated as potentially confounded by observational coverage until a symmetric positive/negative light-curve construction and feature extraction pipeline is rerun. RF impurity importance and SHAP agreement alone cannot rule out this confounding.
The existing script that removes n_points is not a controlled ablation because its Random Forest hyperparameters differ from the submitted model.
All revised ablations should therefore use identical splits and identical classifier hyperparameters, changing only the feature set.

