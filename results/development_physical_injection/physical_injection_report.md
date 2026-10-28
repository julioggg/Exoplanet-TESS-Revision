# Experiment 24: blind PHYSICAL synthetic injection recovery

DEVELOPMENT ONLY, 36 original-TGLC sources, 108 injected trials.
Do NOT interpret as final test accuracy, actual false-alarm rate, or
population transit-detection completeness.

 injection_level  radius_level  radius_ratio period_stratum  n_unique_source_TICs  n_trials  n_direct  direct_fraction  cluster_bootstrap95_low  cluster_bootstrap95_high  n_direct_or_alias  n_at_least_2_observed_injected_cycles  n_at_least_3_observed_injected_cycles
               0             0          0.03        1.5-3.5                    12        12         2         0.166667                 0.000000                  0.416667                  2                                     12                                     12
               0             1          0.06        1.5-3.5                    12        12         2         0.166667                 0.000000                  0.416667                  2                                     12                                     12
               0             2          0.12        1.5-3.5                    12        12         5         0.416667                 0.166667                  0.666667                  5                                     12                                     12
               1             0          0.03        3.5-6.5                    12        12         0         0.000000                 0.000000                  0.000000                  0                                     12                                     12
               1             1          0.06        3.5-6.5                    12        12         4         0.333333                 0.083333                  0.583333                  4                                     12                                     12
               1             2          0.12        3.5-6.5                    12        12         6         0.500000                 0.250000                  0.833333                  6                                     12                                     12
               2             0          0.03        6.5-9.5                    12        12         0         0.000000                 0.000000                  0.000000                  1                                     12                                      9
               2             1          0.06        6.5-9.5                    12        12         2         0.166667                 0.000000                  0.416667                  3                                     12                                      7
               2             2          0.12        6.5-9.5                    12        12         8         0.666667                 0.416667                  0.916667                  8                                     12                                     11

## Marginal summaries (36 source TICs per factor level)
       factor   value  n_stars  n_trials  direct_count  direct_fraction  tic_bootstrap95_low  tic_bootstrap95_high
injection_bin 1.5-3.5       36        36             9         0.250000             0.111111              0.388889
injection_bin 3.5-6.5       36        36            10         0.277778             0.138889              0.416667
injection_bin 6.5-9.5       36        36            10         0.277778             0.138889              0.416667
 radius_ratio    0.03       36        36             2         0.055556             0.000000              0.138889
 radius_ratio    0.06       36        36             8         0.222222             0.083333              0.361111
 radius_ratio    0.12       36        36            19         0.527778             0.361111              0.694444

## Integrity and interpretation

- Synthetic solar-star transit population, not actual planet radii or stellar models.
- Period and depth are crossed in a balanced 3x3 design; each TIC supplies only 3/9 cells.
- These are the 36 previously inspected development stars.
- In-sector original planet signals may compete with injected signals.
- Injected trials with no in-sector sampled transit remain in denominator.
- Original photometry and time sampling is unchanged; no true negatives.
- A recovered injected period is not equivalent to a calibrated planet detection.
- Empirical noise-only false-alarm rate is NOT reported.
- The 18 final locked holdout sources must remain unopened.
- All 18 LOCKED holdout stellar light curves remain unaccessed.
