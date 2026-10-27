# TESS paired transit injection/recovery (pilot)

Stars requested: 12
Successful injection pairs: 12
Failures: 0

## Held-out recovery metrics

- status: completed
- n_original_stars_requested: 12
- n_unique_stars_usable: 12
- n_pairs: 12
- n_star_train: 8
- n_star_test: 4
- n_test_injections: 4
- n_test_control_stars: 4
- threshold_calibrated_on_training_controls_95pct: 31.21266135536353
- heldout_control_false_positive_rate: 0.25
- heldout_injected_detection_rate: 0.25
- heldout_direct_period_recovery_fraction: 0.25
- heldout_period_or_2to1_alias_recovery_fraction: 0.25
- heldout_detection_and_direct_period_recovery_fraction: 0.0
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
