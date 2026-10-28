# Blind BLS period recovery on previously unused NASA-confirmed TESS hosts

## Frozen design
Photometry: original TGLC calibrated aperture, fixed 1800-second exposure.
BLS scanning was completed WITHOUT using the NASA planet period or epoch.
New unique source stars: 18/18; failures: 0

## Cohort-level period concordance (Clopper-Pearson 95% CI)

                        metric  n_success  n_hosts_evaluated     rate  exact_binomial_ci95_low  exact_binomial_ci95_high
blind_direct_period_match_2pct          4                 18 0.222222                 0.064092                  0.476373
        blind_2to1_alias_match          1                 18 0.055556                 0.001406                  0.272944
blind_any_direct_or_2to1_alias          5                 18 0.277778                 0.096949                  0.534802

## By catalogue orbital period stratum

period_bin  n_hosts  direct_period_matches  direct_match_rate  alias_matches
   1.5-3.5        6                      4           0.666667              0
   3.5-6.5        5                      0           0.000000              0
   6.5-9.5        7                      0           0.000000              1

## Per-source BLIND searches

       tic_id          planet  sector  period_truth_days  best_period_days  relative_period_error  blind_direct_period_match_2pct  blind_2to1_alias_match  best_depth_snr_white_noise_proxy
TIC_176956893      TOI-2184 b       1           6.906830         13.836084               1.003247                           False                    True                        221.190133
TIC_142105158      TOI-2200 b       1           2.936028         13.825206               3.708812                           False                   False                        128.930892
 TIC_50618703       TOI-544 b       6           1.548352          1.547887              -0.000300                            True                   False                         13.063361
 TIC_96925186      TOI-2518 b       6           6.897974          3.005501              -0.564292                           False                   False                          6.899097
TIC_262530407       GJ 3090 b       2           2.853102          2.414479              -0.153736                           False                   False                         21.709606
TIC_209459275       TOI-559 b       4           6.983910         13.480745               0.930258                           False                   False                        164.484424
TIC_382602147      TOI-2384 b       1           2.135703          2.135284              -0.000196                            True                   False                         75.833909
TIC_200805484      TOI-2654 b       3           2.967604          2.969242               0.000552                            True                   False                         29.303313
 TIC_94986319       TOI-421 b       5           5.197576         13.078270               1.516225                           False                   False                        103.844598
 TIC_44792534       TOI-954 b       4           3.684973         12.324081               2.344416                           False                   False                        211.357802
TIC_260004324      LHS 1815 b       1           3.814334         13.930358               2.652107                           False                   False                        432.742275
TIC_121490076    HIP 113103 b       1           7.610303         13.828832               0.817120                           False                   False                        403.163838
TIC_409594381      TOI-4214 b       7           3.491388          3.491373              -0.000004                            True                   False                         46.161011
 TIC_77552382  TIC 77552382 b       8           3.510259         12.273318               2.496414                           False                   False                         49.280079
TIC_166527623     HIP 67522 b      11           6.959473         14.390848               1.067807                           False                   False                        299.785401
  TIC_9344899   TIC 9344899 b       5           8.582157         13.078270               0.523891                           False                   False                         84.812698
TIC_117642575 TIC 117642575 b       2           4.944591         14.735309               1.980087                           False                   False                         16.028752
TIC_460205581       TOI-837 b      10           8.324911          2.998250              -0.639846                           False                   False                         69.663952

## Interpretation limitations

- This is positive-host-only period concordance, not planet validation.
- A close orbital-period match may arise from correlated noise or chance.
- An in-sector transit is not established by a NASA planet-host flag.
- No statistically calibrated false alarm rate can be derived here.
- The NASA archive snapshot was frozen before downloading these curves.
- Timing reference and propagated ephemeris uncertainties not yet audited.
- New independent samples required for future detector tuning.
