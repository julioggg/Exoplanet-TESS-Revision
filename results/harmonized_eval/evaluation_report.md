# Harmonized TESS: controlled RF ablation and repeated CV

Each TIC appears once. Same 70/30 holdout and RF hyperparameters.
Repeated 5-fold CV is restricted to the 70% training set.

HOLDOUT METRICS
                   model  n_features  roc_auc  average_precision  accuracy  balanced_accuracy  precision  recall  specificity      f1  tn  fp  fn  tp
                FULL_RAW          11  0.93805            0.96589   0.86486            0.86515    0.90476 0.86364      0.86667 0.88372  39   6   9  57
         NO_N_POINTS_RAW          10  0.82660            0.89351   0.71171            0.70808    0.77419 0.72727      0.68889 0.75000  31  14  18  48
           N_POINTS_ONLY           1  0.75354            0.85734   0.69369            0.69293    0.76667 0.69697      0.68889 0.73016  31  14  20  46
            BLS_ONLY_RAW           5  0.77929            0.85361   0.70270            0.69343    0.75385 0.74242      0.64444 0.74809  29  16  17  49
       FULL_NOISE_SCALED          10  0.93384            0.96231   0.86486            0.86515    0.90476 0.86364      0.86667 0.88372  39   6   9  57
NO_N_POINTS_NOISE_SCALED           9  0.79461            0.86871   0.72072            0.72273    0.79661 0.71212      0.73333 0.75200  33  12  19  47
   BLS_ONLY_NOISE_SCALED           5  0.77778            0.85049   0.72973            0.72677    0.79032 0.74242      0.71111 0.76562  32  13  17  49
      OBSERVATIONAL_ONLY           3  0.92879            0.95210   0.81081            0.79848    0.82609 0.86364      0.73333 0.84444  33  12   9  57
        STATISTICAL_ONLY           4  0.77475            0.86056   0.70270            0.70758    0.78947 0.68182      0.73333 0.73171  33  12  21  45
 TRANSIT_PARAMETERS_ONLY           3  0.74024            0.81337   0.68468            0.67475    0.73846 0.72727      0.62222 0.73282  28  17  18  48

POST-HARMONIZATION OBSERVATIONAL BALANCE
                     feature  negative_median  positive_median  cliffs_delta_pos_minus_neg  ks_statistic  ks_pvalue  mannwhitney_pvalue
                    n_points       987.000000      1100.000000                    0.471660      0.464918   0.000000            0.000000
         baseline_days_after        26.000224        26.452863                    0.032093      0.158786   0.019664            0.601178
median_cadence_minutes_after        30.000164        30.000085                   -0.048566      0.072026   0.710122            0.428786

REPEATED CV FOLD SUMMARY (NOT CONFIDENCE INTERVALS)
                   model            metric  count    mean     std     min     max
   BLS_ONLY_NOISE_SCALED average_precision     15 0.85604 0.03328 0.79569 0.91448
   BLS_ONLY_NOISE_SCALED balanced_accuracy     15 0.68471 0.06508 0.59370 0.81774
   BLS_ONLY_NOISE_SCALED                f1     15 0.76265 0.04228 0.70769 0.87879
   BLS_ONLY_NOISE_SCALED           roc_auc     15 0.78027 0.05004 0.70276 0.87823
            BLS_ONLY_RAW average_precision     15 0.84912 0.04242 0.76516 0.91674
            BLS_ONLY_RAW balanced_accuracy     15 0.66822 0.06588 0.52919 0.80323
            BLS_ONLY_RAW                f1     15 0.73289 0.05225 0.64516 0.83333
            BLS_ONLY_RAW           roc_auc     15 0.76024 0.06004 0.62826 0.84355
       FULL_NOISE_SCALED average_precision     15 0.94903 0.02447 0.90936 0.98161
       FULL_NOISE_SCALED balanced_accuracy     15 0.84044 0.04222 0.77381 0.91167
       FULL_NOISE_SCALED                f1     15 0.86920 0.03332 0.81967 0.91525
       FULL_NOISE_SCALED           roc_auc     15 0.92287 0.03074 0.86587 0.96774
                FULL_RAW average_precision     15 0.94463 0.03230 0.86843 0.97964
                FULL_RAW balanced_accuracy     15 0.83666 0.04020 0.75269 0.89435
                FULL_RAW                f1     15 0.86591 0.03080 0.80702 0.90909
                FULL_RAW           roc_auc     15 0.91928 0.03150 0.85794 0.96621
NO_N_POINTS_NOISE_SCALED average_precision     15 0.88451 0.03433 0.83089 0.94416
NO_N_POINTS_NOISE_SCALED balanced_accuracy     15 0.72463 0.05759 0.60138 0.81644
NO_N_POINTS_NOISE_SCALED                f1     15 0.78009 0.04438 0.71429 0.85714
NO_N_POINTS_NOISE_SCALED           roc_auc     15 0.82362 0.04828 0.74040 0.89597
         NO_N_POINTS_RAW average_precision     15 0.87893 0.03861 0.80995 0.93556
         NO_N_POINTS_RAW balanced_accuracy     15 0.71816 0.05817 0.60829 0.78710
         NO_N_POINTS_RAW                f1     15 0.77973 0.04507 0.66667 0.83077
         NO_N_POINTS_RAW           roc_auc     15 0.81457 0.05364 0.70584 0.88226
           N_POINTS_ONLY average_precision     15 0.83907 0.04130 0.76291 0.90079
           N_POINTS_ONLY balanced_accuracy     15 0.61439 0.05746 0.52619 0.75192
           N_POINTS_ONLY                f1     15 0.69291 0.05642 0.63333 0.79412
           N_POINTS_ONLY           roc_auc     15 0.72794 0.06367 0.61825 0.82796
      OBSERVATIONAL_ONLY average_precision     15 0.91728 0.02462 0.87140 0.95361
      OBSERVATIONAL_ONLY balanced_accuracy     15 0.78138 0.06137 0.67823 0.88018
      OBSERVATIONAL_ONLY                f1     15 0.82278 0.04964 0.75000 0.90909
      OBSERVATIONAL_ONLY           roc_auc     15 0.86783 0.04289 0.79758 0.93548
        STATISTICAL_ONLY average_precision     15 0.86453 0.03962 0.78938 0.93075
        STATISTICAL_ONLY balanced_accuracy     15 0.72349 0.06036 0.62519 0.82857
        STATISTICAL_ONLY                f1     15 0.77632 0.05006 0.69091 0.85714
        STATISTICAL_ONLY           roc_auc     15 0.79951 0.05342 0.70968 0.88326
 TRANSIT_PARAMETERS_ONLY average_precision     15 0.85403 0.03115 0.79899 0.92566
 TRANSIT_PARAMETERS_ONLY balanced_accuracy     15 0.67955 0.05888 0.55376 0.77097
 TRANSIT_PARAMETERS_ONLY                f1     15 0.75166 0.03366 0.67742 0.79310
 TRANSIT_PARAMETERS_ONLY           roc_auc     15 0.76598 0.04081 0.69892 0.85968

CAUTIONS
Class-dependent preprocessing upstream of the harmonized curves remains.
Negative classes denote catalogue status, not proven absence of transits.
BLS noise-scaled variables are uncalibrated constant-noise proxies.
BLS searches periods up to 30 days even though coverage is <=27 days.
