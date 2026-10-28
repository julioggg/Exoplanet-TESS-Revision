# Experiment 26: paired source-clustered audit of development null scoring

All results exploratory, from 12 **already studied** development hosts.
Does not establish astrophysical real-world FPR, specificity or FAP.

## Frozen thresholds (no retuning)
Log-power: 10.935637
Multi-event evidence: 5.329805

## Original 12-host / 144-surrogate null audit
Raw null exceedances: 31/144.
Evidence null exceedances: 6/144.
Mean paired source-level difference (evidence minus raw): -0.173611.
Paired host-bootstrap 95% CI: [-0.347222, -0.027778].
Exploratory host-signflip p: 0.093750.

## Frozen 36 physical-injection trials on the same 12 hosts
Correct blind injected periods before threshold: 9/36.
Correct-period AND accepted raw: 0/36.
Correct-period AND accepted evidence: 7/36.
Accepted raw with wrong injected period: 9/36.
Accepted evidence with wrong injected period: 3/36.
Paired host-level difference (correct AND accepted; evidence minus raw): 0.194444.
Paired host-bootstrap 95% CI: [0.055556, 0.333333].
Exploratory host-signflip p: 0.062500.

## Dependence on surrogate block size

 block_cadences  n_trials  n_sources  raw_exceedances  evidence_exceedances  raw_exceedance_fraction  evidence_exceedance_fraction  absolute_evidence_minus_raw
             24        48         12               11                     1                 0.229167                      0.020833                    -0.208333
             48        48         12                8                     1                 0.166667                      0.020833                    -0.145833
             96        48         12               12                     4                 0.250000                      0.083333                    -0.166667

## Limitations

- The original stars, detection scores and null model have been examined during development.
- Nulls are moving-block resamples of positive hosts; they may contain authentic transit fragments.
- Repeated synthetic null draws from each host are dependent; unit of inference is TIC.
- No inference of real astrophysical false-positive rate, true-negative specificity, or detection false-alarm probability.
- Exact sign-flip p values are descriptive on 12 previously examined source clusters; cannot establish novelty or final validation.
- The synthetic injection profiles used assumed solar host properties and may overlap real transit signals.
- 12 host clusters are insufficient for precise or broadly generalizable calibration.
- A truly independent astrophysical-control cohort is still needed.
- The 18 final locked host sources remain unexamined.
