# Experiment 21: blind BLS multi-event heuristic, DEVELOPMENT ONLY

All 36 previously selected and uniformly processed development TICs.
No final locked-holdout curves accessed. Period truth withheld from scans.

## Blind-period concordance with NASA (descriptive, NOT sensitivity)

                    metric  success  n_hosts  fraction  exact_binomial_ci95_low  exact_binomial_ci95_high
         global_BLS_direct        6       36  0.166667                 0.063720                  0.328116
     multievent_BLS_direct        3       36  0.083333                 0.017526                  0.224690
global_BLS_direct_or_alias        6       36  0.166667                 0.063720                  0.328116

## Orbital-period stratification

period_stratum  n_hosts  global_direct  multievent_direct  multievent_found
       1.5-3.5       12              3                  2                12
       3.5-6.5       12              2                  1                10
       6.5-9.5       12              1                  0                10

## Exploratory diagnostics

Global BLS peaks >=12d: 25/36
Multi-event candidates found: 32/36
Median grid ranks examined: 2596

## Warning

- Development sample is not an independent final evaluation.
- Transit host membership does not prove observed sector contains planet transit.
- Positive event depths computed against global out-of-box median may be distorted by trends.
- Event coherence and TGLC flags are heuristic; threshold tuning requires separate controls.
- BLS SNR assumes a white-noise surrogate, not calibrated physical significance.
- Positive-host period agreement may occur by chance, and alias handling is descriptive.
- No causal proof of planet signal, false-positive calibration or planet classification.
