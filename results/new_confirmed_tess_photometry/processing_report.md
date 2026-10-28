# TGLC original uniform calibrated-aperture photometry pilot

Policy: TGLC_APERTURE_1800S_V1
Processed: 18/18
Errors: 0
Each target uses the identical flux selection and quality rules.
No transit features were extracted; no labels were used.

## Per-star flux QA

       tic_id  sector pair_id  n_raw  n_good_TESS_default  n_positive_finite_aperture  n_accepted  fraction_accepted  n_TGLC_flag_nonzero_raw  n_TGLC_flag_nonzero_accepted  n_TGLC_strict_zero_if_sensitivity  median_cadence_min  baseline_days  n_gaps_over_3h  flux_median_raw  flux_mad_sigma_relative  flux_p2p_sigma_ppm  normalized_flux_p01  normalized_flux_p99              policy_id
TIC_176956893       1           1282                 1267                        1282        1267           0.988300                        0                             0                               1267           29.999674      27.833034               1         1.000007                 0.000776          681.759045             0.973307             1.002094 TGLC_APERTURE_1800S_V1
TIC_142105158       1           1282                 1267                        1282        1267           0.988300                        0                             0                               1267           29.999678      27.833027               1         0.999985                 0.001203         1189.518901             0.975292             1.002933 TGLC_APERTURE_1800S_V1
 TIC_50618703       6            993                  981                         993         981           0.987915                       36                            36                                945           29.999581      21.749678               1         1.000005                 0.000309          259.694360             0.999121             1.000745 TGLC_APERTURE_1800S_V1
 TIC_96925186       6            993                  981                         992         980           0.986908                       20                            18                                962           30.000282      21.750190               1         0.999995                 0.000389          392.419901             0.998680             1.001861 TGLC_APERTURE_1800S_V1
TIC_262530407       2           1245                 1196                        1245        1196           0.960643                      236                           227                                969           30.000341      27.354417               1         1.000023                 0.000340          284.996308             0.998292             1.000875 TGLC_APERTURE_1800S_V1
TIC_209459275       4           1060                 1027                        1060        1027           0.968868                       93                            88                                939           29.999618      25.895507               2         0.999980                 0.000495          463.662049             0.983912             1.011855 TGLC_APERTURE_1800S_V1
TIC_382602147       1           1282                 1267                        1282        1267           0.988300                       25                            25                               1242           30.000292      27.833611               1         0.999825                 0.002732         2570.600027             0.959325             1.007364 TGLC_APERTURE_1800S_V1
TIC_200805484       3           1288                 1077                        1286        1077           0.836180                      208                           167                                910           30.000314      27.312765               5         1.000098                 0.001809         1727.344062             0.987981             1.006262 TGLC_APERTURE_1800S_V1
 TIC_94986319       5           1196                 1176                        1194        1176           0.983278                      126                           114                               1062           30.000292      26.146081               1         0.999995                 0.000327          291.815236             0.996818             1.000868 TGLC_APERTURE_1800S_V1
 TIC_44792534       4           1060                 1027                        1060        1027           0.968868                       95                            90                                937           30.000253      25.896070               2         0.999959                 0.000609          515.288253             0.991794             1.006135 TGLC_APERTURE_1800S_V1
TIC_260004324       1           1282                 1267                        1282        1267           0.988300                        0                             0                               1267           30.000304      27.833622               1         0.999994                 0.000599          505.584604             0.962411             1.003964 TGLC_APERTURE_1800S_V1
TIC_121490076       1           1282                 1267                        1282        1267           0.988300                        0                             0                               1267           30.000382      27.833594               1         1.000030                 0.000488          338.481343             0.955310             1.004381 TGLC_APERTURE_1800S_V1
TIC_409594381       7           1093                 1086                        1093        1086           0.993596                       34                            27                               1059           29.999679      24.416403               1         0.999995                 0.000557          513.332644             0.995234             1.002585 TGLC_APERTURE_1800S_V1
 TIC_77552382       8            968                  924                         950         924           0.954545                      100                            80                                844           30.000452      24.583633               2         1.000032                 0.002051         1904.540743             0.964237             1.024993 TGLC_APERTURE_1800S_V1
TIC_166527623      11           1248                 1037                        1078        1037           0.830929                      348                           314                                723           29.999559      27.062249               4         1.000045                 0.007637         1189.697018             0.980793             1.015818 TGLC_APERTURE_1800S_V1
  TIC_9344899       5           1196                 1176                        1193        1176           0.983278                       69                            63                               1113           30.000292      26.146081               1         0.999951                 0.001848         1860.826604             0.985977             1.004717 TGLC_APERTURE_1800S_V1
TIC_117642575       2           1245                 1196                        1245        1196           0.960643                      311                           294                                902           30.000380      27.354411               1         0.999993                 0.000490          402.418592             0.998147             1.002375 TGLC_APERTURE_1800S_V1
TIC_460205581      10           1211                 1150                        1211        1150           0.949628                      360                           299                                851           30.000303      26.208634               7         1.000037                 0.001478          432.393517             0.995951             1.004232 TGLC_APERTURE_1800S_V1

## Shared sector cadence numbers (audit ONLY)

(none)

## Interpretation

Previous nonfinite PSF values do not reject a valid aperture target.
No per-target fallback is applied; aperture is mandatory for ALL.
Barycentric timestamps are NEVER replaced with a matched source's times.
TGLC_flags remain recorded to support future sensitivity checks.
The four development stars must not serve as a new heldout test set.
