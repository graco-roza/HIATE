# S3_get_best_predictors

This folder is part of the workflow for analyzing taxonomic and functional beta diversity. It focuses on identifying the best environmental predictors for beta diversity using precomputed metrics and candidate predictor variables.

## Folder Structure

``` plaintext
S3_get_best_predictors
├── Functions
│   └── functions_get_best_predictors.R   # Contains auxiliary functions for GDM modeling and data preparation.
├── Scripts
│   ├── get_best_predictors.R            # Main script for identifying the best predictors.
│   └── run_getpred.sh                   # Shell script for submitting batch jobs to HPC.
├── best_predictors
│   ├── coeff                            # Output folder for best predictors based on coefficients.
│   └── r2                               # Output folder for best predictors based on R² values.
└── relationship_direction
    ├── coeff                            # Direction of relationships based on coefficients.
    └── r2                               # Direction of relationships based on R² values.
```

## Overview

### Objective

This module identifies the best environmental predictors for taxonomic and functional beta diversity by:

-   Analyzing precomputed beta diversity metrics.
-   Exploring combinations of predictor variables.
-   Selecting predictors based on effect size (coefficients) and explained variance (R²).
-   Storing results, including predictor relationships and effect directions.

### Key Scripts

1.  **`get_best_predictors.R`**
    -   Main script to process datasets, compute the best predictors, and store the results.
    -   Designed to be run on HPC systems or locally.
2.  **`functions_get_best_predictors.R`**
    -   Provides custom functions for GDM modeling, Bayesian Bootstrap GDM (BBGDM), and data handling.
3.  **`run_getpred.sh`**
    -   Shell script for submitting the main script (`get_best_predictors.R`) to an HPC cluster.

### Outputs

-   **Best Predictors**:\
    Saved in the `best_predictors` folder, grouped by coefficient (`coeff`) and R² (`r2`).

-   **Relationship Directions**:\
    Stored in `relationship_direction` folders, indicating the direction of taxonomic and functional relationships for each predictor.

## Input Data

-   **Beta Diversity Metrics**:\
    Located in `S2_get_beta_diversity/betadiv_output/` (e.g., `*_beta_Output.rds` files).

-   **Predictor Variables**:\
    Processed data from `S1_Preprocessing/Processed/` (e.g., `*_clean.rds` files).

## Notes

-   **Execution on HPC**:\
    The `run_getpred.sh` script facilitates running the workflow in a high-performance computing environment. Modify it to suit your specific job scheduler (e.g., Slurm).

-   **Local Execution**:\
    For local testing, set `ii <- 1` in the R script to process a single dataset.

-   **Dependencies**:\
    Ensure the following R packages are installed: `tidyverse`, `gdm`, `magrittr`, `surveillance`, `glue`, and `plotrix`.
