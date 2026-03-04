# S1_Preprocessing

This step cleans each raw community/trait dataset and appends site-level predictors used throughout the pipeline.

## Contents

- `Functions/auxiliary_functions.R`: helper functions used by preprocessing scripts.
- `Scripts/pre_processing.R`: main preprocessing script (supports local execution or SLURM array execution).
- `Scripts/Figure_S1_and_S4.R`: generates supplementary plotting objects using processed data.
- `Miscellaneous/dataset_info_all.xlsx`: dataset-level metadata used during cleaning.
- `Raw_data/*.xlsx`: raw input datasets.
- `Processed/*_clean.rds`: output produced by `pre_processing.R`.
- `run_step1.sh`: optional HPC submission wrapper.

## Reproducibility notes

1. Run `Scripts/pre_processing.R` first to generate all `*_clean.rds` files.
2. Keep `dataset_info_all.xlsx` synchronized with `Raw_data/` file names.
3. `Figure_S1_and_S4.R` reads `S6_Synthesis_model/data/synthesis_data.xlsx`; run S6 compilation beforehand if you regenerate synthesis inputs.
