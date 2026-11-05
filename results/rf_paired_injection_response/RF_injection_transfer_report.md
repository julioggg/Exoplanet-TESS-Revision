# Experiment 30: Random Forest response to physical transit injections

EXPLORATORY complementary injection-transfer analysis. This is not a rerun of the 108 TGLC-host injections or the locked 18-host final evaluation.

Original source catalogue: 369 unique TICs; 42 injection test TICs and 12 pilot TICs excluded before RF model fitting.
Remaining RF source pool: 315; training 220; normal test 95.
Independent paired transfer test: 42 stars, 126 synthetic control/injected pairs.

## Ordinary catalogued-label test (same trained models)
      model  n_features  train_n_source_TICs  normal_test_n_source_TICs  normal_catalogue_AUC  normal_catalogue_average_precision  normal_catalogue_F1  normal_catalogue_balanced_accuracy  normal_TN  normal_FP  normal_FN  normal_TP
    FULL_11          11                  220                         95               0.98276                             0.99229              0.96350                             0.91379         24          5          0         66
NO_N_POINTS          10                  220                         95               0.97910                             0.99111              0.96241                             0.93312         26          3          2         64
 BLS_ONLY_5           5                  220                         95               0.96813                             0.98648              0.94737                             0.90831         25          4          3         63

## Paired physical-injection effect (42 source TIC bootstrap clusters)
      model  n_injection_pairs  n_distinct_evaluation_stars  mean_delta_RF_probability  mean_delta_cluster_CI95_low  mean_delta_cluster_CI95_high  median_pair_delta  positive_delta_fraction  positive_delta_cluster_CI95_low  positive_delta_cluster_CI95_high  control_RF_positive_fraction  injected_RF_positive_fraction  negative_to_positive_flip_fraction  negative_to_positive_flip_CI95_low  negative_to_positive_flip_CI95_high  positive_to_negative_flips  direct_BLS_recovery_fraction_in_same_pairs                                                          note
    FULL_11                126                           42                   -0.00430                     -0.02998                       0.01211              0.000                  0.43651                          0.33333                           0.53968                       0.02381                        0.00000                             0.00000                             0.00000                              0.00000                           3                                     0.31746 RF scores are catalogue-host outputs, not transit probability
NO_N_POINTS                126                           42                    0.02671                     -0.01773                       0.06391              0.005                  0.55556                          0.43651                           0.66667                       0.09524                        0.08730                             0.01587                             0.00000                              0.03968                           3                                     0.31746 RF scores are catalogue-host outputs, not transit probability
 BLS_ONLY_5                126                           42                    0.11497                      0.06233                       0.16721              0.002                  0.53968                          0.42857                           0.64306                       0.16667                        0.34921                             0.21429                             0.13492                              0.29365                           4                                     0.31746 RF scores are catalogue-host outputs, not transit probability

## By injected radius-ratio bin
      model  radius_bin  n_pairs  n_stars  mean_delta_probability  fraction_positive_delta  control_RF_positive  injected_RF_positive  flip_negative_to_positive  direct_BLS_period_recovered
 BLS_ONLY_5 0.035–0.070       43       30                -0.01200                  0.25581              0.23256               0.20930                          1                            3
 BLS_ONLY_5 0.070–0.110       44       29                 0.09686                  0.54545              0.18182               0.34091                          8                           12
 BLS_ONLY_5 0.110–0.150       39       25                 0.27538                  0.84615              0.07692               0.51282                         18                           25
    FULL_11 0.035–0.070       43       30                -0.01474                  0.34884              0.04651               0.00000                          0                            3
    FULL_11 0.070–0.110       44       29                -0.00368                  0.36364              0.00000               0.00000                          0                           12
    FULL_11 0.110–0.150       39       25                 0.00651                  0.61538              0.02564               0.00000                          0                           25
NO_N_POINTS 0.035–0.070       43       30                -0.01526                  0.51163              0.11628               0.06977                          0                            3
NO_N_POINTS 0.070–0.110       44       29                -0.00150                  0.47727              0.11364               0.11364                          0                           12
NO_N_POINTS 0.110–0.150       39       25                 0.10482                  0.69231              0.05128               0.07692                          2                           25

## Interpretation
The full 11-feature RF obtained a catalogue-label AUC of 0.983 on its separate normal heldout test, but on 126 paired injected-minus-original samples its mean score shift was -0.0043 (source bootstrap 95% interval -0.0300 to +0.0121). A difference in these endpoints is not a direct difference of RF AUC.

- The paired injection sources come from the catalogue-negative group, not verified transit-free stars.
- The RF was originally built to discriminate confirmed-host catalogue labels, not injection presence.
- Injection BLS preprocessing and trial-period grid differ from those used to build original RF features; power is rescaled for unit comparability, not complete pipeline equivalence.
- The original RF sampling/cadence shortcut can persist even after physical injections.
- RF probabilities are not calibrated astrophysical transit-detection probabilities.
- Normal-test AUC is for catalogue labels among other stars; it is not an injection-detection AUC.
- Only 42 source clusters, not 126 independent stars, define the injection transfer population.
- This is a separate experiment from the 108 TGLC injections in manuscript section 4.7 and from the 18-host frozen final test.
- All results postdate the 18-host frozen benchmark and cannot be used to change its algorithm.
