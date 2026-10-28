# TGLC original uniform calibrated-aperture photometry pilot

Policy: TGLC_APERTURE_1800S_V1
Processed: 4/4
Errors: 0
Each target uses the identical flux selection and quality rules.
No transit features were extracted; no labels were used.

## Per-star flux QA

       tic_id  sector pair_id  n_raw  n_good_TESS_default  n_positive_finite_aperture  n_accepted  fraction_accepted  n_TGLC_flag_nonzero_raw  n_TGLC_flag_nonzero_accepted  n_TGLC_strict_zero_if_sensitivity  median_cadence_min  baseline_days  n_gaps_over_3h  flux_median_raw  flux_mad_sigma_relative  flux_p2p_sigma_ppm  normalized_flux_p01  normalized_flux_p99              policy_id
TIC_144285389       2   P0002   1245                 1196                        1245        1196           0.960643                      297                           285                                911           29.999745      27.353828               1         1.000017                 0.001197         1150.739664             0.996209             1.002800 TGLC_APERTURE_1800S_V1
TIC_183532609       2   P0002   1245                 1196                        1245        1196           0.960643                      309                           293                                903           30.000380      27.354411               1         1.000021                 0.000355          286.871523             0.984024             1.001264 TGLC_APERTURE_1800S_V1
TIC_144415503       7   P0013   1093                 1086                        1093        1086           0.993596                       37                            30                               1056           30.000341      24.416937               1         1.000087                 0.003170         2972.181536             0.992339             1.017787 TGLC_APERTURE_1800S_V1
 TIC_22598327       7   P0013   1093                 1086                        1093        1086           0.993596                       30                            28                               1058           29.999677      24.416399               1         0.999998                 0.000197          148.531678             0.999358             1.000485 TGLC_APERTURE_1800S_V1

## Shared sector cadence numbers (audit ONLY)

pair_id  sector         tic_a         tic_b  n_valid_a  n_valid_b  n_joint_cadence_nums  jaccard_cadence_nums  median_abs_BTJD_difference_seconds  same_source_timestamps_forced
  P0002       2 TIC_144285389 TIC_183532609       1196       1196                  1196                   1.0                            13.26237                          False
  P0013       7 TIC_144415503  TIC_22598327       1086       1086                  1086                   1.0                            80.45887                          False

## Interpretation

Previous nonfinite PSF values do not reject a valid aperture target.
No per-target fallback is applied; aperture is mandatory for ALL.
Barycentric timestamps are NEVER replaced with a matched source's times.
TGLC_flags remain recorded to support future sensitivity checks.
The four development stars must not serve as a new heldout test set.
