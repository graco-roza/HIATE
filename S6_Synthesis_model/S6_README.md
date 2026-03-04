# S6_Synthesis_model

This step compiles synthesis data and fits cross-dataset Bayesian models for direction, magnitude, shape, and convergence.

## Data files

- `data/synthesis_data.xlsx`: canonical synthesis table used by S6 modeling scripts.
- `data/BBGDM_SES.xlsx`: SES results used in convergence and extended summaries.

## Scripts

1. `scripts/1-compile_synthesis_data.R`
   - Compiles synthesis inputs and writes `data/synthesis_data.xlsx`.
2. `scripts/2-model_direction.R`
3. `scripts/3-model_magnitude.R`
4. `scripts/4-model_shape.R`
5. `scripts/5-model_convergence.R`
- `scripts/slopes_randomeffect.R`: additional slope/random-effect diagnostics.

## Functions

- `functions/helper_functions_S6.R`
- `functions/helper_functions_plot_extended.R`

## Outputs

- Model objects in `S7_Model_outputs_figures_and_tables/model/`
- Figures in `S7_Model_outputs_figures_and_tables/main_figures/` and `.../supplementary_figures/`

## Reproducibility notes

- Run scripts in numeric order.
- Ensure CmdStan + `cmdstanr` are configured before fitting models.
- Keep the same seed settings if exact posterior replication is required.
