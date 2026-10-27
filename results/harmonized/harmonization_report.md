# Cadence-matched symmetric reconstruction

Both classes were rebinned to 30-minute cadence, limited to the first 27 days, and processed with identical clipping, statistical definitions and BLS grids.

- Targets discovered: 369
- Successful: 369
- Failed: 0
- Positive successful: 220
- Negative successful: 149

## n_points_after

- negative: median=987, IQR=[882, 1036], range=[584, 1184], n=149
- positive: median=1100, IQR=[978.75, 1225], range=[477, 1269], n=220

## baseline_days_after

- negative: median=26.0002, IQR=[24.4327, 26.9371], range=[20.5629, 27], n=149
- positive: median=26.4529, IQR=[24.4288, 26.8326], range=[20.2481, 26.9812], n=220

## median_cadence_minutes_after

- negative: median=30.0002, IQR=[29.9998, 30.0004], range=[29.9992, 30.0028], n=149
- positive: median=30.0001, IQR=[29.9997, 30.0004], range=[29.9992, 30.003], n=220

The next step is to statistically test the post-harmonization observational distributions and then train controlled ablations on harmonized features.
