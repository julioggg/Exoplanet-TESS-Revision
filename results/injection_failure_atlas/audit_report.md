# TESS injection failure atlas and uncertainty audit

This is a post-hoc descriptive analysis of the FIXED expanded sample.
No thresholds or model parameters have been changed.

## Two-stage bootstrap: threshold AND heldout TIC uncertainty

                                metric  point_estimate  two_stage_ci95_low  two_stage_ci95_high  bootstrap_median  n_calibration_stars  n_evaluation_stars  bootstrap_unit  n_bootstrap
                             threshold        36.35483            25.60857             73.09868          36.35483                   95                  42 TIC both stages         3000
                   false_positive_rate         0.09524             0.00000              0.19048           0.07143                   95                  42 TIC both stages         3000
               injected_detection_rate         0.15873             0.00000              0.30159           0.13492                   95                  42 TIC both stages         3000
joint_detection_and_direct_period_rate         0.06349             0.00000              0.14286           0.04762                   95                  42 TIC both stages         3000

## Baseline SNR tail: highest-scoring heldout control stars

       tic_id  n_injections  control_snr  control_false_alarm  injected_snr_max  injected_snr_mean  mean_delta_snr  fraction_snr_improved  fraction_period_recovered  detection_rate  joint_detection_period  median_nominal_depth  median_depth_to_robust_sigma  median_noise_sigma
TIC_158899044             3     57.13182                 True          57.69389           57.29205         0.16023                0.33333                    0.00000             1.0                     0.0               0.00795                       1.44780             0.00549
TIC_318338452             3     56.51878                 True          56.51565           56.50929        -0.00949                0.00000                    0.00000             1.0                     0.0               0.00233                       0.08196             0.02837
TIC_258001134             3     54.15797                 True          58.58920           55.58181         1.42384                0.33333                    0.00000             1.0                     0.0               0.01451                       4.23955             0.00342
TIC_403866308             3     51.64080                 True          51.63971           51.63183        -0.00897                0.00000                    0.00000             1.0                     0.0               0.00592                       0.00516             1.14820
TIC_459348749             3     29.24363                False          29.23615           29.22836        -0.01527                0.00000                    0.00000             0.0                     0.0               0.00757                       0.10254             0.07381
TIC_181856912             3     26.83958                False          26.79344           26.77118        -0.06841                0.00000                    0.00000             0.0                     0.0               0.00347                       1.13499             0.00306
 TIC_26701250             3     25.41543                False          25.95247           25.72570         0.31027                0.66667                    0.00000             0.0                     0.0               0.00571                       0.66399             0.00860
TIC_211418785             3     19.57736                False          19.65151           19.44331        -0.13405                0.33333                    0.00000             0.0                     0.0               0.00669                       3.64191             0.00184
TIC_325109026             3     17.97626                False          28.94798           22.70645         4.73019                0.66667                    0.66667             0.0                     0.0               0.02135                       6.80644             0.00314
TIC_253794766             3     17.67797                False          20.10280           18.91577         1.23780                0.66667                    0.33333             0.0                     0.0               0.01105                       1.94347             0.00569
TIC_457306533             3     15.52495                False          20.36550           17.08908         1.56413                0.33333                    0.33333             0.0                     0.0               0.00712                       2.24429             0.00317
TIC_925879373             3     15.40701                False          15.37036           15.30410        -0.10291                0.00000                    0.33333             0.0                     0.0               0.00423                       1.22351             0.00346
 TIC_64496748             3     13.89856                False          13.77903           13.74982        -0.14873                0.00000                    0.00000             0.0                     0.0               0.00488                       1.05775             0.00462
TIC_307904524             3     12.93377                False          15.25977           14.02332         1.08955                0.66667                    0.00000             0.0                     0.0               0.01265                       1.98180             0.00638
TIC_908231255             3     11.53604                False          17.83846           13.93477         2.39872                0.66667                    0.33333             0.0                     0.0               0.01428                       2.43291             0.00587

## Interpretive limits

Nominal 5% calibration does not guarantee observed 5% false positives.
Correlated stellar and instrumental noise invalidates simple white-noise significance interpretations.
The baseline source stars are catalogue negatives, not confirmed transit-free objects.
All exploratory findings here must be confirmed using new, independent targets.
Never tune using these heldout stars and then reuse the same data as test.
