# Read-only, publication-oriented synthesis: TESS transit period validation

**Status:** FROZEN INDEPENDENT TEST FINISHED. No model refitting, new test or hyperparameter changes. Read-only reporting on commits preceding this script.

## Methods-ready text (English)

We selected one original 1800-s TGLC product per TIC using a metadata-only, prespecified protocol, excluding the 369-source original study and a separately examined 18-host engineering cohort. Thirty-six confirmed-planet hosts supported detector development; 18 distinct hosts were reserved for a one-shot independent assessment with six in each of three orbital-period strata (1.5–3.5, 3.5–6.5, and 6.5–9.5 days). All products were subjected to the identical previously frozen TGLC_APERTURE_1800S_V1 processing policy. The blind period estimator returned the maximum of an unchanged 4,000-period (0.5–15 d), five-duration (0.05–0.25 d) box-least-squares (BLS) search. The independent observational reference period was read only after blind candidate estimation. Direct period concordance was defined as a relative discrepancy within 2%; half/double-period aliases were tabulated separately. Two confidence criteria were fixed from prior development-only surrogate-null calibration: log(1+BLS maximum power) and an event-evidence score combining locally positive events, TGLC flags, and a same-period leave-one-event influence ratio. Their thresholds and the evaluation source-code hashes were Git-committed before any final-holdout flux was accessed. The final cohort was evaluated exactly once, without retuning.

## Results-ready text (English)

All 18 original TGLC source files were verified, and all 18 uniformly processed light curves passed the frozen pipeline. In the previously untouched holdout, the global blind BLS period agreed with the catalogue period for 4/18 hosts (22.2%; exact binomial 95% interval 6.4–47.6%), with no half/double-period aliases. The unmodified BLS-power criterion accepted 2/18 source candidates, neither having a concordant period; the frozen multi-event confidence criterion accepted 2/18, both having concordant periods. The remaining two period-concordant hosts did not pass either confidence criterion. Across the 36 development hosts, direct concordance was 6/36 with the original BLS maximum; a previously tested hard three-event period-replacement rule recovered only 3/36, motivating the separation of period estimation from evidence-based candidate acceptance. In 108 development-only physical-injection trials, 29/108 injected orbital periods were recovered directly by blind BLS; the results are synthetic-injection recovery rather than astrophysical detection sensitivity.

## Limitations-ready text (English)

The final 18-source sample contains confirmed planet HOSTS, not independently established in-sector transit detections. In a given selected sector, the planet may have no observable transit or its signal may not dominate the light curve. Therefore catalogue-period concordance should not be interpreted as transit sensitivity, and catalogue-period discordance does not automatically demonstrate a false astrophysical detection. The sample contains no verified transit-free control stars, so these results cannot quantify astrophysical false-positive rate, specificity, precision or ROC AUC. The development null model used block-resampled signals from previously inspected positive hosts and cannot establish a true no-planet null. Both accepted candidate groups consist of only two objects; the observed discrimination is descriptive, not statistically conclusive or broadly generalizable. The frozen final cohort must not be used for any subsequent algorithm or threshold selection.

## Exact per-host endpoints

Concordant catalogued-period hosts: TOI-237 b, TOI-2364 b, HD 63935 b, TOI-2346 b.
Frozen event-score accepted hosts: TOI-2364 b, TOI-2346 b.
Frozen raw-power accepted hosts: TOI-4324 b, TOI-4308 b.

## Reproducibility and integrity

- Frozen protocol SHA256: 805a6e5261176a2c392cefa218172112a48333f3c7f43a2029d3134c3bf1c55e
- Final results protocol SHA256: bccd48f3ce10985f7c1c428963df02ab990e3fcebf17fdb1f56ceb5ba2f22ff8
- Final per-host scored SHA256: 7b816c3c896f7e7fbbc945134babd9e219b2173f8cf35d41bffc834d91c3aa1a
- Main thresholds: raw=10.935637240467, evidence=5.329805494853
- Read-only figures and tables are DERIVED PUBLICATION MATERIAL, not a rerun of the locked holdout detector.
- Changes to the frozen code, frozen thresholds or final result numbers are not permitted.
