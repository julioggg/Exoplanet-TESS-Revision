# Prespecified independent TESS model-development and LOCKED holdout cohorts

New model-development TICs: 36/36
New locked and untouched final TICs: 18/18
Total distinct newly selected TICs: 54/54
MAST metadata queries performed: 153
MAST metadata errors: 0
Frozen design complete: True

## Selection by period stratum and role

role            period_bin
LOCKED_HOLDOUT  1.5-3.5        6
                3.5-6.5        6
                6.5-9.5        6
development     1.5-3.5       12
                3.5-6.5       12
                6.5-9.5       12

## Integrity and limitations

- No FITS, light curves, transit metrics, BLS, or model output were inspected.
- The previous 18 hosts and all 369 original study source TICs are excluded.
- DO NOT download or inspect LOCKED_HOLDOUT until model pipeline is frozen.
- Only host-positive labels are known; no true-negative inference is made.
- Sources with missing observed transits must not be silently replaced.
- A future genuinely independent test can only occur ONCE per frozen model.
