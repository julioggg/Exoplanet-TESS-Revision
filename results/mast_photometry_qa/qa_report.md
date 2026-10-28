# Original MAST TGLC photometric QA (label-blind)

FITS audited: 3/4
Failed: 1
Predeclared primary: cal_psf_flux
Secondary sensitivity: cal_aper_flux

## Shared sampling and flag breakdown

       tic_id pair_id  sector  n_total  n_finite_time  n_nonfinite_cal_psf_flux  n_nonfinite_cal_aper_flux  n_nonpositive_cal_psf_flux  n_nonpositive_cal_aper_flux  n_nonzero_TESS_flags  n_failed_TESS_default  n_nonzero_TGLC_flags  n_common_quality_kept  n_strict_TGLC_zero_kept  fraction_kept               tess_flags_unique tglc_flags_unique
TIC_183532609   P0002       2     1245           1245                         0                          0                           0                            0                    49                     49                   309                   1196                      903        0.96064 0:1196;32:1;128:33;160:5;164:10       0:936;1:309
TIC_144415503   P0013       7     1093           1093                         0                          0                           0                            0                     7                      7                    37                   1086                     1056        0.99360               0:1086;32:1;160:6       0:1056;1:37
 TIC_22598327   P0013       7     1093           1093                         0                          0                           0                            0                     7                      7                    30                   1086                     1058        0.99360               0:1086;32:1;160:6       0:1063;1:30

## Per-curve photometric scatter (no clipping/detrending)

 n_valid_cadences  normalized_mean  normalized_median  normalized_std  normalized_robust_sigma  normalized_p2p_sigma  normalized_p2p_ppm  normalized_amplitude_p99_p01  raw_median_flux  min_normalized_flux  max_normalized_flux  baseline_days  median_cadence_minutes  n_gaps_over_1h  n_gaps_over_3h        tic_id pair_id  sector               method  n_valid_common  fraction_kept  tess_default_mask  tglc_flag_zero_required  no_sigma_clipping  no_detrending  n_tglc_strict_valid  strict_tglc_robust_sigma
             1196         0.999722                1.0        0.002184                 0.000933              0.000746          746.429900                      0.016281         0.999968             0.984399             1.006096      27.354411               30.000380              24               1 TIC_183532609   P0002       2          primary_psf            1196       0.960643               True                    False               True           True                  903                  0.000848
             1196         0.999680                1.0        0.002116                 0.000355              0.000287          286.871523                      0.017241         1.000021             0.982616             1.002278      27.354411               30.000380              24               1 TIC_183532609   P0002       2 sensitivity_aperture            1196       0.960643               True                    False               True           True                  903                  0.000317
             1086         1.000003                1.0        0.003856                 0.003630              0.003427         3427.168518                      0.018693         0.999906             0.975267             1.012439      24.416937               30.000341               6               1 TIC_144415503   P0013       7          primary_psf            1086       0.993596               True                    False               True           True                 1056                  0.003520
             1086         1.000271                1.0        0.004074                 0.003170              0.002972         2972.181536                      0.025448         1.000087             0.988672             1.028514      24.416937               30.000341               6               1 TIC_144415503   P0013       7 sensitivity_aperture            1086       0.993596               True                    False               True           True                 1056                  0.003148
             1086         0.999804                1.0        0.001330                 0.000880              0.000789          788.785149                      0.008598         0.999984             0.990963             1.002579      24.416399               29.999677               1               1  TIC_22598327   P0013       7          primary_psf            1086       0.993596               True                    False               True           True                 1058                  0.000864
             1086         0.999983                1.0        0.000224                 0.000197              0.000149          148.531678                      0.001128         0.999998             0.999047             1.000605      24.416399               29.999677               1               1  TIC_22598327   P0013       7 sensitivity_aperture            1086       0.993596               True                    False               True           True                 1058                  0.000194

## Pair timestamp support overlap

pair_id  sector         tic_a        tic_b  camera_a  camera_b  ccd_a  ccd_b  same_camera_and_ccd  n_a  n_b  n_matched_half_minute_bins  jaccard_half_minute_bins
  P0013       7 TIC_144415503 TIC_22598327         2         3      2      2                False 1086 1086                           0                       0.0

## Cautions

- Four previously inspected targets: engineering feasibility, not independent validation.
- Different camera/CCD combinations may cause remaining differences even within sector.
- Primary-mission calibrated TGLC flux can be corrupted for faint/contaminated sources.
- No transit-injection/recovery or false-positive calibration is attempted here.
- Point-to-point scatter is a descriptive proxy, not a rigorous CDPP calculation.
- Rejecting nonpositive calibrated flux may bias rare, strongly variable targets.
- Manual exclusion of TGLC-flagged cadences would require separate validation.
