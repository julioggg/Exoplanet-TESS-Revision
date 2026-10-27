# TESS paired transit injection/recovery (full)

Stars requested: 149
Successful injection pairs: 447
Failures: 0

## Held-out recovery metrics

- status: completed
- n_original_stars_requested: 149
- n_unique_stars_usable: 149
- n_pairs: 447
- n_star_train: 104
- n_star_test: 45
- n_test_injections: 135
- n_test_control_stars: 45
- threshold_calibrated_on_training_controls_95pct: 42.03722946539674
- heldout_control_false_positive_rate: 0.08888888888888889
- heldout_injected_detection_rate: 0.1259259259259259
- heldout_direct_period_recovery_fraction: 0.37777777777777777
- heldout_period_or_2to1_alias_recovery_fraction: 0.4074074074074074
- heldout_detection_and_direct_period_recovery_fraction: 0.037037037037037035
- note: Detection uses training-calibrated 95th percentile of control BLS depth_snr. Period recovery requires relative error <=2%, separately allowing 2:1 aliases. These are descriptive point estimates; stars (not injections) are the independent units.

## Protocol notes

A synthetic transit was multiplied into each original Lightkurve-
exported raw negative-catalogue FITS before common preprocessing.
The SAME base-derived clipping mask, normalization and time bins
were applied to injected/control members of each pair.
The split is by TIC; no source star crosses train/test.
The detection threshold is calibrated on training controls only.
Recovered period and thresholded detection are both reported.
Only pilot feasibility should be inferred from small pilot results.
