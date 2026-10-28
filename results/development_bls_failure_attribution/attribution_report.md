# Experiment 22 — Why multi-event BLS lost previously recovered periods

EXPLORATORY DEVELOPMENT FORENSICS: no model or threshold changes.

Original BLS direct recovery: 6/36.
Multi-event BLS direct recovery: 3/36.
Net absolute difference: -8.333 percentage points.
Paired counts: both=3, original_only=3, multi_only=0, neither=30.
McNemar exact two-sided p (exploratory, not a confirmatory test): 0.25000.
Original highest BLS peaks passing every coherence gate: 4/36.
Lost directly correct original BLS periods due to TGLC flags: 2.
Lost directly correct original BLS periods due to 3-event gate: 1.

## Previously correct period estimates that were vetoed

       tic_id     planet  reference_period_days  original_period_days  multievent_period_days  multievent_alias_only                                 original_peak_veto_reasons  original_peak_n_supported_events  original_peak_TGLC_flag_fraction
 TIC_36452991 TOI-2969 b               1.823715              1.823456                3.647287                   True                             TGLC_flagged_fraction_gt_20pct                                13                          0.250000
TIC_237922465 TOI-2416 b               8.275479              8.270318                4.140410                   True fewer_than_3_supported_events;fewer_than_3_positive_events                                 2                          0.000000
TIC_280655495 TOI-3364 b               5.876892              5.877219                1.468117                  False                             TGLC_flagged_fraction_gt_20pct                                 4                          0.241379

## Methodological conclusion

Do NOT replace a period estimate solely because its evidence-quality flag does not pass a fixed multi-event gate.
Instead separate period SEARCH from candidate VETTING / UNCERTAINTY; these 36 development TICs can guide future methods but not provide an independent final score.
Do not open the 18 LOCKED holdout TICs before algorithm/protocol freeze.
