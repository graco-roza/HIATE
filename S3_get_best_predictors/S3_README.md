# S3_get_best_predictors

This step identifies the best environmental predictors for beta-diversity patterns.

## Contents

- `Scripts/get_best_predictors.R`: main analysis script.
- `Scripts/run_getpred.sh`: optional SLURM submission script.
- `Functions/functions_get_best_predictors.R`: helper functions.
- `best_predictors/`: selected predictors (coefficient- and R2-based outputs).
- `relationship_direction/`: direction summaries for predictor-response relationships.

## Inputs

- `S2_get_beta_diversity/betadiv_output/*_beta_Output.rds`
- `S1_Preprocessing/Processed/*_clean.rds`

## Reproducibility notes

- The script can run dataset-by-dataset in HPC arrays or locally.
- Keep naming consistency between S1 and S2 outputs to avoid dropped datasets.
