# EXPERIMENT 27: LOCKED, one-shot independent host benchmark

## Freeze must be GIT-COMMITTED before MAST FITS acquisition
This is a METADATA-ONLY freeze; no final holdout curves inspected.

Prior known positive host cohort: 18 stars, 4 direct period matches.
Development cohort: 36 stars; frozen independent final cohort: 18.
Final per-bin period quota: 6 hosts × 3 bins.

## Immutable detector and evaluation
Original TGLC aperture V1; full 4000x5 BLS, highest peak unchanged.
Retain original BLS period; evaluate fixed multi-event confidence.
Raw BLS log-power threshold: 10.935637240.
Multi-event evidence threshold: 5.329805495.

## Scientific restrictions
- No post-hoc changes to period search, scoring formula or threshold.
- Report ALL 18 previously selected hosts; never replace nonrecoveries.
- Report direct known-host period agreement, not planet sensitivity.
- No true-negative host group: cannot estimate FPR or specificity.
- No threshold p-value significance claimed; prior work exploratory.
- If data processing fails, document and fix without selecting new TICs.
- DO NOT open a single final host FITS before this commit is pushed.

## Required order after freeze commit
1. Reuse script 16 to acquire all 18 selected original FITS.
2. Reuse script 14 for unchanged original aperture photometry.
3. Execute script 28 once, after all 18 source checks pass.
4. Preserve independent result even if it is worse than development.
