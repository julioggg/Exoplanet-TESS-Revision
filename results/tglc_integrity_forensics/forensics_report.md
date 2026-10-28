# Original TGLC channel and cadence forensics

Previously obtained FITS rechecked by SHA256.

## Which original channels reject each star?

       tic_id  sector pair_id  n_raw  n_time_finite  n_tess_default_good  n_time_and_tess_good  n_tglc_flags_nonzero  n_corrected_joint_positive  n_all_four_joint_positive  n_corrected_psf_only_good  n_corrected_aperture_only_good  n_unprocessed_psf_good  n_unprocessed_aperture_good
TIC_144285389       2   P0002   1245           1245                 1196                  1196                   297                           0                          0                          0                            1196                       0                         1196
TIC_183532609       2   P0002   1245           1245                 1196                  1196                   309                        1196                       1196                          0                               0                    1196                         1196
TIC_144415503       7   P0013   1093           1093                 1086                  1086                    37                        1086                       1086                          0                               0                    1086                         1086
 TIC_22598327       7   P0013   1093           1093                 1086                  1086                    30                        1086                       1086                          0                               0                    1086                         1086

## Quality-good and signed-flux counts

       tic_id  sector pair_id       channel  n_finite  n_nonfinite  n_positive  n_zero  n_negative  n_quality_good_positive  fraction_quality_good_positive  n_quality_good_negative  n_quality_good_zero  median_finite_flux  median_positive_good_flux  finite_flux_q01  finite_flux_q99  normalized_robust_sigma_if_eligible  eligible_by_channel_300_points
TIC_144285389       2   P0002  cal_psf_flux         0         1245           0       0           0                        0                             0.0                        0                    0                 NaN                        NaN              NaN              NaN                                  NaN                           False
TIC_144285389       2   P0002 cal_aper_flux      1245            0        1245       0           0                     1196                             1.0                        0                    0             0.99995                    1.00002          0.97504          1.00281                              0.00120                            True
TIC_144285389       2   P0002      psf_flux         0         1245           0       0           0                        0                             0.0                        0                    0                 NaN                        NaN              NaN              NaN                                  NaN                           False
TIC_144285389       2   P0002 aperture_flux      1245            0        1245       0           0                     1196                             1.0                        0                    0          2539.58398                 2539.89233       2472.51976       2551.47281                              0.00277                            True
TIC_183532609       2   P0002  cal_psf_flux      1245            0        1245       0           0                     1196                             1.0                        0                    0             0.99994                    0.99997          0.97670          1.00359                              0.00093                            True
TIC_183532609       2   P0002 cal_aper_flux      1245            0        1245       0           0                     1196                             1.0                        0                    0             1.00001                    1.00002          0.97651          1.00128                              0.00036                            True
TIC_183532609       2   P0002      psf_flux      1245            0        1245       0           0                     1196                             1.0                        0                    0         32547.99023                32552.39746      31756.55977      32739.76875                              0.00440                            True
TIC_183532609       2   P0002 aperture_flux      1245            0        1245       0           0                     1196                             1.0                        0                    0         32554.58594                32558.94922      31677.45039      32804.43812                              0.00395                            True
TIC_144415503       7   P0013  cal_psf_flux      1093            0        1093       0           0                     1086                             1.0                        0                    0             0.99991                    0.99991          0.98865          1.00905                              0.00363                            True
TIC_144415503       7   P0013 cal_aper_flux      1093            0        1093       0           0                     1086                             1.0                        0                    0             1.00006                    1.00009          0.99162          1.01850                              0.00317                            True
TIC_144415503       7   P0013      psf_flux      1093            0        1093       0           0                     1086                             1.0                        0                    0           857.66473                  857.66342        846.30240        865.50718                              0.00404                            True
TIC_144415503       7   P0013 aperture_flux      1093            0        1093       0           0                     1086                             1.0                        0                    0           857.84467                  857.86990        844.84108        873.37521                              0.00672                            True
 TIC_22598327       7   P0013  cal_psf_flux      1093            0        1093       0           0                     1086                             1.0                        0                    0             0.99999                    0.99998          0.99279          1.00207                              0.00088                            True
 TIC_22598327       7   P0013 cal_aper_flux      1093            0        1093       0           0                     1086                             1.0                        0                    0             1.00000                    1.00000          0.99919          1.00048                              0.00020                            True
 TIC_22598327       7   P0013      psf_flux      1093            0        1093       0           0                     1086                             1.0                        0                    0        159331.40625               159331.67188     155654.22125     160923.21188                              0.00719                            True
 TIC_22598327       7   P0013 aperture_flux      1093            0        1093       0           0                     1086                             1.0                        0                    0        159331.21875               159331.88281     158926.91688     159562.19125                              0.00068                            True

## Matching actual observing cadence numbers

pair_id  sector         tic_a         tic_b  n_quality_good_a  n_quality_good_b  n_common_cadence_nums  cadence_num_jaccard  n_common_all_four_corrected_valid  median_absolute_time_offset_seconds  median_signed_time_offset_seconds  range_time_offset_seconds                                  timing_comparison_is
  P0002       2 TIC_144285389 TIC_183532609              1196              1196                   1196                  1.0                                  0                              13.2624                             1.0363                    50.3034 same TESS cadence identifier; BTJD differs per target
  P0013       7 TIC_144415503  TIC_22598327              1086              1086                   1086                  1.0                               1086                              80.4589                            80.4589                    46.4794 same TESS cadence identifier; BTJD differs per target

## Design constraints

- No previous flux-selection or clipping rule was changed.
- Independent per-column validity is inspected BEFORE any new policy.
- A TGLC calibrated flux may be corrupt even if a raw column is finite.
- Matched cadence_num is not automatically identical local barycentric time.
- No astrophysical labels were used; these old TICs are engineering samples.
- Next sample must be newly drawn with an a priori quality policy.
