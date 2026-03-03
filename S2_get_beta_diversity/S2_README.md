# S2_get_beta_diversity

This step computes taxonomic and functional beta diversity from preprocessed datasets.

## Scripts

1. `Scripts/1-get_functional_dissimilarities.R`
   - Builds functional dissimilarity matrices (saved to `betadiv_input/`).
2. `Scripts/2-get_beta_diversity.R`
   - Computes beta-diversity outputs (saved to `betadiv_output/`).
3. `Scripts/prepare_dataset_for_trait_analysis.R`
   - Helper script for preparing selected datasets for trait analyses.

HPC wrappers:
- `Scripts/run_gawdis.sh`
- `Scripts/run_betadiv.sh`

## Inputs

- `S1_Preprocessing/Processed/*_clean.rds`

## Outputs

- `betadiv_input/*.rds`
- `betadiv_output/*_beta_Output.rds`

## Reproducibility notes

- Run scripts in numeric order.
- Ensure package versions used for `gawdis`, `BAT`, and `hypervolume` are recorded in your environment.
