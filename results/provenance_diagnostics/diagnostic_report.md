# TESS observational shortcut diagnostics

Repeated the experiment 03 holdout; only observational variables
are included. A true 30-min cadence, rounded to one minute,
eliminates sub-second timestamp fluctuations that could serve
as product or preprocessing fingerprints.

                   model  n_features  roc_auc  average_precision      f1  balanced_accuracy  tn  fp  fn  tp
         OBS_ALL_PRECISE           3  0.92997            0.95305 0.85075            0.80960  34  11   9  57
     COUNT_PLUS_BASELINE           2  0.91364            0.94201 0.84444            0.79848  33  12   9  57
      COUNT_PLUS_CADENCE           2  0.79444            0.88063 0.74603            0.71162  32  13  19  47
   BASELINE_PLUS_CADENCE           2  0.75337            0.79263 0.73529            0.65657  25  20  16  50
              COUNT_ONLY           1  0.75337            0.84554 0.74016            0.70051  31  14  19  47
           BASELINE_ONLY           1  0.73013            0.77249 0.76471            0.69394  27  18  14  52
    CADENCE_ONLY_PRECISE           1  0.46313            0.56313 0.51240            0.46818  21  24  35  31
OBS_CADENCE_ROUNDED_1MIN           3  0.91061            0.94043 0.85294            0.80606  33  12   8  58
 OBS_CADENCE_FIXED_30MIN           3  0.91061            0.94043 0.85294            0.80606  33  12   8  58
          OCCUPANCY_ONLY           1  0.85320            0.88289 0.77778            0.74899  34  11  17  49
COUNT_BASELINE_OCCUPANCY           3  0.91869            0.95231 0.83212            0.77626  31  14   9  57

FITS headers audited: 369

Do not use these models for exoplanet detection or planet validation.
