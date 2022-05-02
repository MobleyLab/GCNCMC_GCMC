Each system was equilibrated with GCMC moves by three stages: uvt1, npt, uvt2 (details in the paper). The input coordinate structure and scripts of performing these simulations for each system can be found in [`gcmc/`](gcmc).

The production phase of GCMC and GCNCMC simulations was performed using the same equilibrated structure. The scripts of performing the initial production run and extended production run can be found in "prod1" and "prod2" folders for each system in [`gcmc/`](gcmc) and [`gcncmc/`](gcncmc).

The atoms and radius used to define the region of interest for water enhanced sampling can be found in the scripts to run production simulations for each system.

The complete list of simulated systems with their PDB IDs: HSP90 (2xab, 2xjg, 3rlp, 3rlr, 3rlq), HIV1-protease (1ec1, 1ec0, 1ebw, 1eby), Trypsin (1c5t, 1gi1, 1f0u, 1o2j), Factor Xa (1ezq, 1lpg, 1lpz, 1f0s), thrombin (2zff), BTK (4zlz), PTP1B (2qbs), TAF1 (5i1q).
