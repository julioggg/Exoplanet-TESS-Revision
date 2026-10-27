# Controlled ablation on the submitted feature table

All models use the exact same train/test stars and identical Random Forest hyperparameters (500 trees, balanced class weights, random_state=42). Only the feature set changes.

**Important:** this experiment quantifies baseline feature dependence but does not remove the known class-dependent cadence/download/BLS-grid confounding.

## Hold-out metrics

| model | accuracy | precision | recall | specificity | f1 | roc_auc | average_precision | tn | fp | fn | tp |\n| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |\n| FULL | 0.927928 | 0.926471 | 0.954545 | 0.888889 | 0.940299 | 0.979966 | 0.989552 | 40 | 5 | 3 | 63 |\n| NO_N_POINTS | 0.945946 | 0.954545 | 0.954545 | 0.933333 | 0.954545 | 0.976768 | 0.987829 | 42 | 3 | 3 | 63 |\n| N_POINTS_ONLY | 0.918919 | 0.938462 | 0.924242 | 0.911111 | 0.931298 | 0.944613 | 0.937417 | 41 | 4 | 5 | 61 |\n| BLS_ONLY | 0.900901 | 0.898551 | 0.939394 | 0.844444 | 0.918519 | 0.919697 | 0.909532 | 38 | 7 | 4 | 62 |\n| TRANSIT_PARAMETERS_ONLY | 0.738739 | 0.768116 | 0.80303 | 0.644444 | 0.785185 | 0.789562 | 0.853062 | 29 | 16 | 13 | 53 |\n| STATISTICAL_ONLY | 0.828829 | 0.821918 | 0.909091 | 0.711111 | 0.863309 | 0.901178 | 0.931928 | 32 | 13 | 6 | 60 |\n| VARIABILITY_ONLY | 0.765766 | 0.822581 | 0.772727 | 0.755556 | 0.796875 | 0.855387 | 0.910901 | 34 | 11 | 15 | 51 |\n| DETECTION_STRENGTH_ONLY | 0.72973 | 0.736842 | 0.848485 | 0.555556 | 0.788732 | 0.822727 | 0.883815 | 25 | 20 | 10 | 56 |

## Key reviewer comparisons

- FULL ROC-AUC=0.9800, F1=0.9403.
- NO_N_POINTS ROC-AUC=0.9768, F1=0.9545.
- N_POINTS_ONLY ROC-AUC=0.9446, F1=0.9313.
- BLS_ONLY ROC-AUC=0.9197, F1=0.9185.

Interpretation must be deferred until the symmetric reconstruction is run, because n_points is strongly coupled to class-dependent cadence and the BLS features were computed using different positive/negative grids.

