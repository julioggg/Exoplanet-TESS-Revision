# TGLC original uniform calibrated-aperture photometry pilot

Policy: TGLC_APERTURE_1800S_V1
Processed: 18/18
Errors: 0
Each target uses the identical flux selection and quality rules.
No transit features were extracted; no labels were used.

## Per-star flux QA

       tic_id  sector pair_id  n_raw  n_good_TESS_default  n_positive_finite_aperture  n_accepted  fraction_accepted  n_TGLC_flag_nonzero_raw  n_TGLC_flag_nonzero_accepted  n_TGLC_strict_zero_if_sensitivity  median_cadence_min  baseline_days  n_gaps_over_3h  flux_median_raw  flux_mad_sigma_relative  flux_p2p_sigma_ppm  normalized_flux_p01  normalized_flux_p99              policy_id
TIC_235078021       5           1196                 1176                        1196        1176           0.983278                       69                            60                               1116           30.000321      26.146101               1         1.000006                 0.000812          789.954904             0.995948             1.002266 TGLC_APERTURE_1800S_V1
 TIC_98796344       4           1060                 1027                        1060        1027           0.968868                       95                            90                                937           29.999599      25.895500               2         0.999998                 0.001464          765.778818             0.982024             1.006991 TGLC_APERTURE_1800S_V1
TIC_305048087       2           1245                 1196                        1245        1196           0.960643                      306                           294                                902           29.999772      27.353820               1         1.000121                 0.002872         2943.153153             0.992424             1.007158 TGLC_APERTURE_1800S_V1
TIC_408416746       8            968                  954                         962         954           0.985537                      125                           117                                837           30.000395      24.583631               2         1.000087                 0.000563          468.111167             0.869361             1.018693 TGLC_APERTURE_1800S_V1
TIC_255685030       2           1245                 1196                        1245        1196           0.960643                      152                           140                               1056           30.000311      27.354455               2         1.000054                 0.000593          459.191779             0.996893             1.004036 TGLC_APERTURE_1800S_V1
 TIC_39414571       6            993                  981                         993         981           0.987915                       63                            62                                919           29.999581      21.749678               1         0.999986                 0.000743          731.638541             0.989653             1.001525 TGLC_APERTURE_1800S_V1
TIC_451645081      10           1211                 1111                        1140        1111           0.917424                      365                           335                                776           30.000285      25.396124               4         1.000038                 0.000518          410.179968             0.987565             1.002022 TGLC_APERTURE_1800S_V1
TIC_118327550       2           1245                 1196                        1245        1196           0.960643                      311                           294                                902           30.000380      27.354411               1         0.999999                 0.000447          433.097458             0.998356             1.001047 TGLC_APERTURE_1800S_V1
TIC_144193715       1           1282                 1267                        1282        1267           0.988300                       12                            12                               1255           29.999752      27.832999               1         0.999979                 0.000541          478.784665             0.949817             1.004629 TGLC_APERTURE_1800S_V1
TIC_359271092       9           1154                 1135                        1153        1135           0.983536                      278                           260                                875           30.000338      25.125281               2         1.000002                 0.000268          226.983295             0.998747             1.000720 TGLC_APERTURE_1800S_V1
TIC_453211454       7           1093                 1086                        1093        1086           0.993596                       30                            29                               1057           29.999674      24.416400               1         1.000014                 0.000201          180.772099             0.999265             1.000930 TGLC_APERTURE_1800S_V1
TIC_317548889       6            993                  981                         992         980           0.986908                       24                            23                                957           29.999593      21.749688               1         1.000000                 0.000139          127.660776             0.999511             1.000303 TGLC_APERTURE_1800S_V1
TIC_231077395       1           1282                 1267                        1282        1267           0.988300                        0                             0                               1267           30.000317      27.833615               1         0.999987                 0.000497          468.251829             0.976284             1.002237 TGLC_APERTURE_1800S_V1
TIC_259377017       3           1288                 1077                        1282        1077           0.836180                      211                           169                                908           30.000314      27.312775               5         0.999996                 0.000575          565.914894             0.995343             1.001507 TGLC_APERTURE_1800S_V1
TIC_317483660       6            993                  981                         992         980           0.986908                       23                            22                                958           29.999593      21.749688               1         0.999954                 0.000679          692.136597             0.987159             1.001494 TGLC_APERTURE_1800S_V1
 TIC_29960110       4           1060                 1016                        1054        1016           0.958491                       77                            67                                949           29.999567      25.895491               3         1.000093                 0.001126          852.992441             0.940249             1.024761 TGLC_APERTURE_1800S_V1
TIC_176314383       5           1196                 1169                        1182        1169           0.977425                      139                           137                               1032           30.000288      25.979409               1         1.000066                 0.001182         1114.724634             0.996888             1.025201 TGLC_APERTURE_1800S_V1
 TIC_54962195       9           1154                 1125                        1132        1125           0.974870                      197                           196                                929           29.999667      24.833093               1         0.999992                 0.001069          841.892815             0.989612             1.004735 TGLC_APERTURE_1800S_V1

## Shared sector cadence numbers (audit ONLY)

(none)

## Interpretation

Previous nonfinite PSF values do not reject a valid aperture target.
No per-target fallback is applied; aperture is mandatory for ALL.
Barycentric timestamps are NEVER replaced with a matched source's times.
TGLC_flags remain recorded to support future sensitivity checks.
Sources inspected for method development must never be reused as an untouched final test set.
