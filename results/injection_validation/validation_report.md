# Expanded paired transit injection/recovery validation

Full cohort with successful injections: 149 distinct source stars
Pilot stars excluded: 12
Independent calibration stars: 95
Independent evaluation stars: 42
Evaluated injected trials: 126

## Primary, calibrated 5%-nominal-FPR results (star-cluster CIs)

 nominal_control_fpr  calibrated_snr_threshold                            metric  point_estimate  ci95_low  ci95_high  n_independent_test_stars bootstrap_unit  n_cluster_bootstrap
                0.05                  36.35483                       control_fpr         0.09524   0.02381    0.19048                        42            TIC                 3000
                0.05                  36.35483                injected_detection         0.15873   0.06349    0.26190                        42            TIC                 3000
                0.05                  36.35483            direct_period_recovery         0.31746   0.21429    0.42063                        42            TIC                 3000
                0.05                  36.35483          period_or_alias_recovery         0.34127   0.23810    0.44444                        42            TIC                 3000
                0.05                  36.35483 joint_detection_and_direct_period         0.06349   0.01587    0.11905                        42            TIC                 3000
                0.05                  36.35483       paired_snr_improvement_rate         0.48413   0.38095    0.58730                        42            TIC                 3000

## Sensitivity by fixed target-FPR levels (not test-tuned)

 nominal_control_fpr  calibrated_snr_threshold                            metric  point_estimate  ci95_low  ci95_high  n_independent_test_stars bootstrap_unit  n_cluster_bootstrap
                0.01                  99.14690                       control_fpr         0.00000   0.00000    0.00000                        42            TIC                 3000
                0.01                  99.14690                injected_detection         0.00794   0.00000    0.02381                        42            TIC                 3000
                0.01                  99.14690            direct_period_recovery         0.31746   0.22222    0.42063                        42            TIC                 3000
                0.01                  99.14690          period_or_alias_recovery         0.34127   0.24603    0.44444                        42            TIC                 3000
                0.01                  99.14690 joint_detection_and_direct_period         0.00794   0.00000    0.02381                        42            TIC                 3000
                0.01                  99.14690       paired_snr_improvement_rate         0.48413   0.38889    0.58730                        42            TIC                 3000
                0.05                  36.35483                       control_fpr         0.09524   0.02381    0.19048                        42            TIC                 3000
                0.05                  36.35483                injected_detection         0.15873   0.06349    0.26190                        42            TIC                 3000
                0.05                  36.35483            direct_period_recovery         0.31746   0.21429    0.42063                        42            TIC                 3000
                0.05                  36.35483          period_or_alias_recovery         0.34127   0.23810    0.44444                        42            TIC                 3000
                0.05                  36.35483 joint_detection_and_direct_period         0.06349   0.01587    0.11905                        42            TIC                 3000
                0.05                  36.35483       paired_snr_improvement_rate         0.48413   0.38095    0.58730                        42            TIC                 3000
                0.10                  27.12734                       control_fpr         0.11905   0.02381    0.21429                        42            TIC                 3000
                0.10                  27.12734                injected_detection         0.21429   0.11111    0.33333                        42            TIC                 3000
                0.10                  27.12734            direct_period_recovery         0.31746   0.21429    0.42857                        42            TIC                 3000
                0.10                  27.12734          period_or_alias_recovery         0.34127   0.24603    0.45238                        42            TIC                 3000
                0.10                  27.12734 joint_detection_and_direct_period         0.09524   0.03175    0.17460                        42            TIC                 3000
                0.10                  27.12734       paired_snr_improvement_rate         0.48413   0.38095    0.58730                        42            TIC                 3000

## Completeness by injection parameters (descriptive)

   dimension           bin  n_injections  n_unique_source_stars  detection_rate  direct_period_recovery  joint_detection_period_rate
period_group  (1.499, 3.5]            32                     25         0.21875                 0.40625                      0.15625
period_group    (3.5, 6.5]            44                     33         0.06818                 0.29545                      0.00000
period_group    (6.5, 9.5]            50                     33         0.20000                 0.28000                      0.06000
radius_group (0.034, 0.07]            43                     30         0.11628                 0.06977                      0.00000
radius_group  (0.07, 0.11]            44                     29         0.11364                 0.27273                      0.02273
radius_group  (0.11, 0.15]            39                     25         0.25641                 0.64103                      0.17949

## Interpretation safeguards

- Pilot TICs were inspected during development, excluded from primary analysis.
- Catalogue negatives are not proven non-transiting stars.
- Threshold learned using CALIBRATION control TICs only.
- CIs reflect sampling of EVALUATION stars with fixed threshold.
- Per-source replicated injections are not independent.
- Assumed solar stellar parameters; injection completeness is conditional.
- Small-bin summaries are descriptive, not confidence intervals.
- This evaluates BLS detectability and period recovery, not an RF exoplanet detector.
