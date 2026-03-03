# S4_run_BBGDM

This step runs Bayesian bootstrap GDM models and generates Figure 1 assets.

## Contents

- `scripts/main_script_S4.R`: main BBGDM fitting workflow.
- `scripts/plot_figure_1.R`: plotting workflow for Figure 1.
- `scripts/run_final_bbgdm.sh`: optional HPC launcher.
- `functions/helper_functions_S4.R`: helper functions.
- `bbgdm_output/`: model outputs.

## Inputs

- `S2_get_beta_diversity/betadiv_output/*_beta_Output.rds`
- `S3_get_best_predictors/best_predictors/coeff/*`

## Reproducibility notes

- Run `main_script_S4.R` before plotting scripts.
- Keep predictor labels aligned with S3 outputs.
