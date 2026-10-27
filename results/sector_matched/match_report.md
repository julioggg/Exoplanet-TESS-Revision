# Within-sector and observational-covariate matched diagnosis

Original targets: 369
Sector-overlap count: 24
Maximum pairs by recorded sector: 64
Pairs matched within sector: 64
Pairs passing all strict calipers: 34

## Covariate balance

                    design             feature  negative_median  positive_median  ks_statistic  ks_pvalue  mean_difference_positive_minus_negative
               ALL_TARGETS            n_points       987.000000      1100.000000      0.464918   0.000000                               125.955400
               ALL_TARGETS baseline_days_after        26.000224        26.452863      0.158786   0.019664                                 0.037046
               ALL_TARGETS  occupancy_fraction         0.773637         0.916726      0.562752   0.000000                                 0.099422
       SECTOR_ONLY_MATCHED            n_points       921.500000      1005.000000      0.281250   0.012310                                89.093750
       SECTOR_ONLY_MATCHED baseline_days_after        25.022547        24.432621      0.203125   0.142948                                -0.326715
       SECTOR_ONLY_MATCHED  occupancy_fraction         0.768059         0.885714      0.406250   0.000041                                 0.084392
SECTOR_AND_CALIPER_MATCHED            n_points      1004.000000      1013.500000      0.235294   0.306758                                38.529412
SECTOR_AND_CALIPER_MATCHED baseline_days_after        24.854482        24.436419      0.235294   0.306758                                -0.131026
SECTOR_AND_CALIPER_MATCHED  occupancy_fraction         0.844479         0.881202      0.352941   0.028185                                 0.037659

## Pair-grouped cross-validation summaries

                    design              model            metric  count    mean     std     min     max
SECTOR_AND_CALIPER_MATCHED           BLS_ONLY average_precision     15 0.84920 0.08237 0.74235 0.92007
SECTOR_AND_CALIPER_MATCHED           BLS_ONLY balanced_accuracy     15 0.67381 0.11809 0.50000 0.78571
SECTOR_AND_CALIPER_MATCHED           BLS_ONLY                f1     15 0.66872 0.13113 0.46154 0.80000
SECTOR_AND_CALIPER_MATCHED           BLS_ONLY           roc_auc     15 0.77937 0.14193 0.61111 0.89796
SECTOR_AND_CALIPER_MATCHED               FULL average_precision     15 0.90448 0.06089 0.82576 1.00000
SECTOR_AND_CALIPER_MATCHED               FULL balanced_accuracy     15 0.82143 0.04676 0.75000 0.85714
SECTOR_AND_CALIPER_MATCHED               FULL                f1     15 0.81022 0.04688 0.72727 0.85714
SECTOR_AND_CALIPER_MATCHED               FULL           roc_auc     15 0.83946 0.12037 0.66667 1.00000
SECTOR_AND_CALIPER_MATCHED        NO_N_POINTS average_precision     15 0.90838 0.06500 0.82576 1.00000
SECTOR_AND_CALIPER_MATCHED        NO_N_POINTS balanced_accuracy     15 0.80714 0.04436 0.75000 0.85714
SECTOR_AND_CALIPER_MATCHED        NO_N_POINTS                f1     15 0.79740 0.04761 0.72727 0.85714
SECTOR_AND_CALIPER_MATCHED        NO_N_POINTS           roc_auc     15 0.84150 0.13274 0.66667 1.00000
SECTOR_AND_CALIPER_MATCHED OBSERVATIONAL_ONLY average_precision     15 0.90786 0.02409 0.87143 0.94048
SECTOR_AND_CALIPER_MATCHED OBSERVATIONAL_ONLY balanced_accuracy     15 0.79048 0.09784 0.66667 0.92857
SECTOR_AND_CALIPER_MATCHED OBSERVATIONAL_ONLY                f1     15 0.79084 0.08213 0.71429 0.92308
SECTOR_AND_CALIPER_MATCHED OBSERVATIONAL_ONLY           roc_auc     15 0.85646 0.04850 0.77551 0.89796
       SECTOR_ONLY_MATCHED           BLS_ONLY average_precision     15 0.87214 0.04148 0.81975 0.94373
       SECTOR_ONLY_MATCHED           BLS_ONLY balanced_accuracy     15 0.78846 0.03561 0.75000 0.84615
       SECTOR_ONLY_MATCHED           BLS_ONLY                f1     15 0.77088 0.04408 0.70000 0.83333
       SECTOR_ONLY_MATCHED           BLS_ONLY           roc_auc     15 0.82072 0.05612 0.76331 0.91716
       SECTOR_ONLY_MATCHED               FULL average_precision     15 0.91772 0.02081 0.89358 0.95330
       SECTOR_ONLY_MATCHED               FULL balanced_accuracy     15 0.80321 0.06321 0.70833 0.88462
       SECTOR_ONLY_MATCHED               FULL                f1     15 0.78876 0.08369 0.66667 0.88889
       SECTOR_ONLY_MATCHED               FULL           roc_auc     15 0.89730 0.03405 0.84024 0.93491
       SECTOR_ONLY_MATCHED        NO_N_POINTS average_precision     15 0.88791 0.02892 0.86024 0.94091
       SECTOR_ONLY_MATCHED        NO_N_POINTS balanced_accuracy     15 0.78141 0.05862 0.73077 0.88462
       SECTOR_ONLY_MATCHED        NO_N_POINTS                f1     15 0.76950 0.06593 0.69565 0.88000
       SECTOR_ONLY_MATCHED        NO_N_POINTS           roc_auc     15 0.84881 0.04091 0.79290 0.91420
       SECTOR_ONLY_MATCHED OBSERVATIONAL_ONLY average_precision     15 0.93265 0.02989 0.87960 0.96478
       SECTOR_ONLY_MATCHED OBSERVATIONAL_ONLY balanced_accuracy     15 0.80449 0.03622 0.76923 0.84615
       SECTOR_ONLY_MATCHED OBSERVATIONAL_ONLY                f1     15 0.80403 0.04119 0.76190 0.85714
       SECTOR_ONLY_MATCHED OBSERVATIONAL_ONLY           roc_auc     15 0.90472 0.04781 0.81944 0.94675

## Limitations

Matching is only a sensitivity check on a restricted support region.
FITS-recorded sector may not describe all observations in stitched curves.
No confirmed planet is not the same as a verified non-transiting target.
Rebuild from homogeneous source products and run injection/recovery next.
