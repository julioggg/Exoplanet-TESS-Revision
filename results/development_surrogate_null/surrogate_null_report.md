# Experiment 25: surrogate-null BLS evidence calibration (DEVELOPMENT)

Only previously explored 36-TIC development sources; 432 model-dependent
moving-block surrogate nulls. The 18 new final heldout sources are untouched.
THIS DOES NOT ESTIMATE AN ASTROPHYSICAL FALSE-POSITIVE RATE.

## Fixed null-only thresholds and internal development audit

         score_name               source_split  threshold_95pct_calibration_surrogate  n_surrogate_null_trials  n_surrogate_exceedances  surrogate_exceedance_fraction  n_source_hosts  n_source_host_blind_peaks_accepted  n_source_host_correct_periods_accepted  n_injection_trials  n_injections_accepted  n_injection_periods_correct_and_accepted  n_injection_periods_directly_recovered
      raw_log_power                calibration                               10.93564                      288                       14                        0.04861              24                                   1                                       0                  72                      7                                         3                                      20
      raw_log_power internal_development_audit                               10.93564                      144                       31                        0.21528              12                                   3                                       0                  36                      9                                         0                                       9
multievent_evidence                calibration                                5.32981                      288                       14                        0.04861              24                                   2                                       2                  72                     21                                        16                                      20
multievent_evidence internal_development_audit                                5.32981                      144                        6                        0.04167              12                                   2                                       1                  36                     10                                         7                                       9

## Method restrictions

- No real independently vetted non-transiting negative reference cohort.
- Moving-block resampling may retain transit fragments; does not produce proven planet-free nulls.
- This is a null-model exceedance analysis, NOT real-world FPR, specificity, or calibrated astrophysical false-alarm probability.
- The fixed score was chosen before this run, but the 36 development TICs were previously inspected.
- 36 source TICs and 12 null replicates each: inference is sensitive to source dependence, bootstrap block length, nonstationarity and gap handling.
- The 12-source internal split is NOT a scientifically fresh independent test.
- Injection period concordance and acceptance are different metrics.
- No previous 18-host or final 18 locked holdout curves were accessed.
