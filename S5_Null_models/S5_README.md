# S5_Null_models

This step generates null-model expectations and standardized effect sizes (SES) used in synthesis analyses.

## Contents

- `scripts/1-run_null_analysis.R`: runs null BBGDM analyses.
- `scripts/2-null_model_results.R`: summarizes null-model outputs.
- `scripts/3-Null_model_processing.r`: post-processing and export utilities.
- `scripts/check_completeness.R`: sanity-check helper.
- `functions/functions_run_null_bbgdm.R`: null-model core functions.
- `functions/functions_prepare_null_results.R`: aggregation helpers.
- `null_output/`: null-model output files.

## Inputs

- S4 outputs and associated predictor/metadata structures.

## Reproducibility notes

- Run scripts in numeric order.
- Preserve random seeds and model settings if extending null draws.
