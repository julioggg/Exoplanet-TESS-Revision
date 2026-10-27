# Injection diagnostics: pilot

Pairs audited: 12
Failures: 0
Direct blind period recovery: 6/12
Oracle period power increased after injection: 7/12

Important: 'oracle' means true period was supplied. This is only a diagnostic of potential detectability at a known ephemeris.
High BLS depth_snr in control stars may indicate stellar variability, eclipses, systematics or correlated photometric noise; it is not proof that the catalogue-negative star contains a planet.

## Per-injection results

       tic_id split  period_true  period_blind_injected  blind_direct_period_recovery  radius_ratio  blind_control_depth_snr  blind_injected_depth_snr  oracle_delta_power  oracle_delta_depth_snr
 TIC_10070640  test       6.5125                 7.7772                         False        0.0401                   8.4566                    8.4311             -0.1373                 -0.0288
TIC_103144811 train       1.7991                14.7353                         False        0.0827                   6.4758                    6.4226             -0.2545                 -0.1055
TIC_105057280 train       3.7270                 3.7234                          True        0.1153                  42.1950                   70.0928           2357.3813                 51.7140
TIC_107636744 train       1.6486                 0.7393                         False        0.0689                   8.9986                    8.6335              7.8178                  1.4431
TIC_109204497 train       3.7392                 3.7416                          True        0.0603                   4.3312                    6.0066             10.5829                  2.2937
 TIC_11079600 train       5.7107                 5.7104                          True        0.0563                   7.9606                   25.6003            308.9108                 18.8267
TIC_114800566 train       3.1365                 3.1360                          True        0.1437                  10.8169                   36.5976            655.4723                 28.9090
TIC_126595132 train       9.2777                 0.9351                         False        0.0383                   5.0296                    5.1119             -0.0674                 -0.0206
TIC_127105071  test       7.9441                 7.9440                          True        0.0875                  12.0050                   19.9666            186.9170                 14.4820
TIC_128613134  test       5.7585                12.6867                         False        0.0648                 104.7012                  104.6343             -7.3924                 -0.1205
TIC_128975349  test       7.7823                14.1733                         False        0.0549                  25.1825                   25.1709             -0.1603                 -0.0107
TIC_131848736 train       5.6198                 5.6162                          True        0.1482                   9.7775                   16.9463            128.8153                 11.0894

