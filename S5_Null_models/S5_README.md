# S5_Null_models

This folder is part of the workflow for analyzing taxonomic and functional beta diversity. It focuses on generating null model outputs and synthesizing results for comparison with observed data.

## Folder Structure

```plaintext
S5_Null_models
├── functions
│   ├── functions_prepare_null_results.R   # Functions for preparing null model results.
│   └── functions_run_null_bbgdm.R         # Functions for running null models and BBGDM analysis.
├── null_output
│   ├── <dataset>_null_bbgdm.rds           # Null model outputs for various datasets.
│   └── coeff                              # Subfolder containing coefficient-based outputs.
├── scripts
│   ├── 1-run_null_analysis.R              # Script to run null model analyses.
│   └── 2-null_model_results.R             # Script to compile and analyze null model results.
```

## Overview

### Objective

This step aims to:

1. Generate null models to compare taxonomic and functional beta diversity metrics with observed data.
2. Assess the effects of predictors on beta diversity through the null models.
3. Compile and synthesize results for further interpretation.

### Key Components

- **Null Model Analysis**:
  - Performed using `1-run_null_analysis.R`.
  - Generates null model outputs stored in `null_output/`.

- **Synthesis of Results**:
  - Done via `2-null_model_results.R`.
  - Compiles null model results and compares them with observed data.

### Outputs

- **Null Model Results**:
  - Saved in `null_output/` as `.rds` files.
  - Files are named `<dataset>_null_bbgdm.rds`.

- **Compiled Results**:
  - Results are combined and saved for synthesis in subsequent steps.

## Notes

- **Execution**:
  - Designed for High-Performance Computing (HPC) clusters but can also run on local machines.
  - Ensure required R packages are installed: `tidyverse`, `hypervolume`, `doSNOW`, `glue`, and `readxl`.

- **Data Dependency**:
  - Requires observed model results from `S4_run_BBGDM/bbgdm_output`.
  - Outputs from this step feed into `S6_Synthesis_model`.

- **Custom Functions**:
  - Supporting functions are located in the `functions/` folder for modular use.
