#' ###########################################################################
#' SCRIPT NAME: Bayesian Modeling of Magnitude of Effects (3-model_magnitude.R)
#'
#' DESCRIPTION:
#'   This script performs Bayesian modeling to estimate the magnitude of species 
#'   and trait turnover across different ecosystems, land-use types, and other 
#'   predictors. It evaluates multiple model formulations, including interactions 
#'   and nested structures, and generates visualizations of the results.
#'
#' USAGE:
#'   - The script processes preprocessed data to run Bayesian models using the 
#'     `brms` package.
#'   - Several models with varying complexity (base, interaction, nested taxa) 
#'     are compared using Leave-One-Out Cross-Validation (LOO).
#'   - Generates figures for publication and saves outputs for further use.
#'
#' INPUTS:
#'   - Preprocessed synthesis data from `S6_Synthesis_model/data/synthesis_data.xlsx`
#'   - Helper functions for model processing and plotting
#'
#' OUTPUTS:
#'   - Bayesian models saved in `S7_Model_outputs_figures_and_tables/model/magnitude`
#'   - Figures saved in `S7_Model_outputs_figures_and_tables/main_figures` and
#'     `S7_Model_outputs_figures_and_tables/extended_data`
#'
#' AUTHOR: Caio Graco-Roza
#' LAST UPDATED: 2024-11-24
#'
#' NOTES:
#'   - Adjust the global settings for Bayesian models as needed (chains, iterations, etc.).
#'   - Ensure all required libraries are installed, and the input files are in the correct path.
#' ###########################################################################

# Load packages ####
library(tidyverse) # For data manipulation and visualization
library(ggdist)  # For visualizing distributions of data and model outputs
library(ggtext)  # For adding rich text elements to ggplot2 plots
library(brms) # For running Bayesian models
library(distributional) # For defining and visualizing probability distributions, including priors
library(marginaleffects) # For extracting marginal effects and predictions from Bayesian models
library(janitor) # For cleaning and organizing data tables, and standardizing column names
library(ggh4x) # For additional ggplot2 features like nested facets
library(cowplot) # For creating multi-panel plots

source("S6_Synthesis_model/functions/helper_functions_S6.R") # Load helper functions specific to synthesis models
source("S6_Synthesis_model/functions/helper_functions_plot_extended.R") # Load helper functions for extended plotting

# Set some global Stan options
base::options(mc.cores = 4,          # Set the number of parallel chains for Bayesian sampling
              brms.backend = "cmdstanr")   # Use cmdstanr backend for running brms models
showtext::showtext_auto() #show custom fonts on ggplot2

# All models will be run with same settings
CHAINS <- 4                          # Number of chains to run for the Bayesian model
ITER <- 5000                         # Total number of iterations per chain
WARMUP <- 2500                       # Number of iterations used for warmup (not included in posterior samples)
BAYES_SEED <- 1234                   # Seed for reproducibility in Bayesian modeling


model_data <- process_model_data("S6_Synthesis_model/data/synthesis_data.xlsx", "magnitude")
#Set priors ------------------------------

priors_magnitude <- c(
  prior(normal(0,1), class = b) #slopes
  ,prior(student_t(3, 0, .1), class = sd, lb=0) #standard deviation
)

#Set model formula --------------------------

formula_base <- bf(magnitude ~ species_number + spatial_extent + spatial_distance_min + absolute_latitude_mean + human_pressure_baseline + human_pressure_range + main_land_use_type + direction + ecosystem_type + (1 | taxa_coarse), decomp = "QR")
formula_interaction_ecosystem_land_use <- bf(magnitude ~ species_number + spatial_extent + spatial_distance_min + absolute_latitude_mean + human_pressure_baseline + human_pressure_range + ecosystem_type * main_land_use_type + direction + (1 | taxa_coarse), decomp = "QR")
formula_interaction_direction_land_use <- bf(magnitude ~ species_number + spatial_extent + spatial_distance_min + absolute_latitude_mean + human_pressure_baseline + human_pressure_range + direction * main_land_use_type + ecosystem_type + (1 | taxa_coarse), decomp = "QR")
formula_three_way_interaction <- bf(magnitude ~ species_number + spatial_extent + spatial_distance_min + absolute_latitude_mean + human_pressure_baseline + human_pressure_range + ecosystem_type * main_land_use_type * direction + (1 | taxa_coarse), decomp = "QR")
formula_base_nested <- bf(magnitude ~ species_number + spatial_extent + spatial_distance_min + absolute_latitude_mean + human_pressure_baseline + human_pressure_range + main_land_use_type + direction + ecosystem_type + (1 | taxa_coarse/taxa_fine), decomp = "QR")
formula_interaction_ecosystem_land_use_nested <- bf(magnitude ~ species_number + spatial_extent + spatial_distance_min + absolute_latitude_mean + human_pressure_baseline + human_pressure_range + ecosystem_type * main_land_use_type + direction + (1 | taxa_coarse/taxa_fine), decomp = "QR")
formula_interaction_direction_land_use_nested <- bf(magnitude ~ species_number + spatial_extent + spatial_distance_min + absolute_latitude_mean + human_pressure_baseline + human_pressure_range + direction * main_land_use_type + ecosystem_type + (1 | taxa_coarse/taxa_fine), decomp = "QR")
formula_three_way_interaction_nested <- bf(magnitude ~ species_number + spatial_extent + spatial_distance_min + absolute_latitude_mean + human_pressure_baseline + human_pressure_range + ecosystem_type * main_land_use_type * direction + (1 | taxa_coarse/taxa_fine), decomp = "QR")

#Run model ------------------------------------------

md_species_base <- brms::brm(
formula_base,
data = filter(model_data, facet == "Species"), #subset only the taxonomic parts
prior = priors_magnitude,
family=brms::hurdle_lognormal(),
chains = CHAINS,
iter = ITER,
warmup = WARMUP,
control = base::list(max_treedepth = 15, adapt_delta = .99),
save_pars = brms::save_pars(all = TRUE),
seed = BAYES_SEED,
file = "S7_Model_outputs_figures_and_tables/model/magnitude/magnitude_species_base", #here we save models automatically
file_refit = "on_change" #models are only refit if parameters change, otherwise load the saved ones
) 

# Interaction model: Ecosystem * Land Use
md_species_interaction_ecosystem_land_use <- update(md_species_base,
                                                    newdata = filter(model_data, facet == "Species"),
                                                    formula. = formula_interaction_ecosystem_land_use,
                                                    file = "S7_Model_outputs_figures_and_tables/model/magnitude/magnitude_species_interaction_ecosystem_land_use",
                                                    file_refit = "on_change"
)

# Interaction model: Direction * Land Use
md_species_interaction_direction_land_use <- update(md_species_base,
                                                    formula. = formula_interaction_direction_land_use,
                                                    newdata = filter(model_data, facet == "Species"),
                                                    file = "S7_Model_outputs_figures_and_tables/model/magnitude/magnitude_species_interaction_direction_land_use",
                                                    file_refit = "on_change"
)

# Three-way interaction model
md_species_three_way_interaction <- update(md_species_base,
                                           formula. = formula_three_way_interaction,
                                           newdata = filter(model_data, facet == "Species"),
                                           file = "S7_Model_outputs_figures_and_tables/model/magnitude/magnitude_species_three_way_interaction",
                                           file_refit = "on_change"
)

# Nested Base model
md_species_base_nested <- update(md_species_base,
                                 formula. = formula_base_nested,
                                 newdata = filter(model_data, facet == "Species"),
                                 file = "S7_Model_outputs_figures_and_tables/model/magnitude/magnitude_species_base_nested",
                                 file_refit = "on_change"
)

# Nested Interaction model: Ecosystem * Land Use
md_species_interaction_ecosystem_land_use_nested <- update(md_species_base,
                                                           formula. = formula_interaction_ecosystem_land_use_nested,
                                                           newdata = filter(model_data, facet == "Species"),
                                                           file = "S7_Model_outputs_figures_and_tables/model/magnitude/magnitude_species_interaction_ecosystem_land_use_nested",
                                                           file_refit = "on_change"
)

# Nested Interaction model: Direction * Land Use
md_species_interaction_direction_land_use_nested <- update(md_species_base,
                                                           formula. = formula_interaction_direction_land_use_nested,
                                                           newdata = filter(model_data, facet == "Species"),
                                                           file = "S7_Model_outputs_figures_and_tables/model/magnitude/magnitude_species_interaction_direction_land_use_nested",
                                                           file_refit = "on_change"
)

# Nested Three-way interaction model
md_species_three_way_interaction_nested <- update(md_species_base,
                                                  formula. = formula_three_way_interaction_nested,
                                                  newdata = filter(model_data, facet == "Species"),
                                                  file = "S7_Model_outputs_figures_and_tables/model/magnitude/magnitude_species_three_way_interaction_nested",
                                                  file_refit = "on_change"
)

#Perform Leave-one-out -------------------------------------

# md_species_base <- add_criterion(md_species_base, "loo", moment_match = FALSE)
# md_species_interaction_ecosystem_land_use <- add_criterion(md_species_interaction_ecosystem_land_use, "loo", moment_match = TRUE)
# md_species_interaction_direction_land_use <- add_criterion(md_species_interaction_direction_land_use, "loo", moment_match = TRUE)
# md_species_three_way_interaction <- add_criterion(md_species_three_way_interaction, "loo", moment_match = TRUE)
# md_species_base_nested <- add_criterion(md_species_base_nested, "loo", moment_match = TRUE)
# md_species_interaction_ecosystem_land_use_nested <- add_criterion(md_species_interaction_ecosystem_land_use_nested, "loo", moment_match = TRUE)
# md_species_interaction_direction_land_use_nested <- add_criterion(md_species_interaction_direction_land_use_nested, "loo", moment_match = TRUE)
# md_species_three_way_interaction_nested <- add_criterion(md_species_three_way_interaction_nested, "loo", moment_match = TRUE)

loo_species_comparison <- loo(
  md_species_base,
  md_species_interaction_ecosystem_land_use,
  md_species_interaction_direction_land_use,
  md_species_three_way_interaction,
  md_species_base_nested,
  md_species_interaction_ecosystem_land_use_nested,
  md_species_interaction_direction_land_use_nested,
  md_species_three_way_interaction_nested
)

# Run Models for Trait Facet ------------------------------------------

md_trait_base <- brms::brm(
  formula_base,
  data = filter(model_data, facet == "Traits"), # subset only the trait parts
  prior = priors_magnitude,
  family = brms::zero_one_inflated_beta(),
  chains = CHAINS,
  iter = ITER,
  warmup = WARMUP,
  control = list(max_treedepth = 15, adapt_delta = .99),
  save_pars = brms::save_pars(all = TRUE),
  seed = BAYES_SEED,
  file = "S7_Model_outputs_figures_and_tables/model/magnitude/magnitude_trait_base", 
  file_refit = "on_change"
)

md_trait_interaction_ecosystem_land_use <- update(md_trait_base, 
                                                  formula = formula_interaction_ecosystem_land_use, 
                                                  file = "S7_Model_outputs_figures_and_tables/model/magnitude/magnitude_trait_interaction_ecosystem_land_use"
)

md_trait_interaction_direction_land_use <- update(md_trait_base, 
                                                  formula = formula_interaction_direction_land_use, 
                                                  file = "S7_Model_outputs_figures_and_tables/model/magnitude/magnitude_trait_interaction_direction_land_use"
)

md_trait_three_way_interaction <- update(md_trait_base, 
                                         formula = formula_three_way_interaction, 
                                         file = "S7_Model_outputs_figures_and_tables/model/magnitude/magnitude_trait_three_way_interaction"
)

md_trait_base_nested <- update(md_trait_base, 
                               formula = formula_base_nested, 
                               newdata=filter(model_data, facet == "Traits"),
                               file = "S7_Model_outputs_figures_and_tables/model/magnitude/magnitude_trait_base_nested"
)

md_trait_interaction_ecosystem_land_use_nested <- update(md_trait_base, 
                                                         formula = formula_interaction_ecosystem_land_use_nested, 
                                                         newdata=filter(model_data, facet == "Traits"),
                                                         file = "S7_Model_outputs_figures_and_tables/model/magnitude/magnitude_trait_interaction_ecosystem_land_use_nested"
)

md_trait_interaction_direction_land_use_nested <- update(md_trait_base, 
                                                         formula = formula_interaction_direction_land_use_nested, 
                                                         newdata=filter(model_data, facet == "Traits"),
                                                         file = "S7_Model_outputs_figures_and_tables/model/magnitude/magnitude_trait_interaction_direction_land_use_nested"
)

md_trait_three_way_interaction_nested <- update(md_trait_base, 
                                                formula = formula_three_way_interaction_nested, 
                                                newdata=filter(model_data, facet == "Traits"),
                                                file = "S7_Model_outputs_figures_and_tables/model/magnitude/magnitude_trait_three_way_interaction_nested"
)

# Add LOO with moment matching for Trait Models -------------------------------------

loo_traits_comparison <- loo(
  md_trait_base,
  md_trait_interaction_ecosystem_land_use,
  md_trait_interaction_direction_land_use,
  md_trait_three_way_interaction,
  md_trait_base_nested,
  md_trait_interaction_ecosystem_land_use_nested,
  md_trait_interaction_direction_land_use_nested,
  md_trait_three_way_interaction_nested
)

loo_species_comparison
loo_traits_comparison


best_model_species <- md_species_base # md_species_base is chosen as the final model because it achieves similar predictive accuracy to more complex models but with fewer parameters
best_model_traits <- md_trait_base # md_traits_base is chosen as the final model because it achieves similar predictive accuracy but with fewer parameters

# Run the posterior predictive checks -----------------------------------------------------------------------------
cowplot::plot_grid(brms::pp_check(best_model_species, ndraws=100), brms::pp_check(best_model_traits, ndraws=100))
#'-----------------------------------------------------------------------------------------------------------------

# check R2 of each model ------------------------------------------------------------------------------------------
brms::bayes_R2(best_model_species) # Taxonomic  0.2646003
brms::bayes_R2(best_model_traits) # Functional 0.2970056
#'-----------------------------------------------------------------------------------------------------------------

#Make table ------------------------------------------------------------
# 
# hypothesis(best_model_species, "disturbanceurban > disturbanceforest")
# 
# emm <- emmeans(fit_magnitude_tax,specs=~ disturbance, level = 0.8)
# # Define the contrast for the hypothesis
# pairs(emm)
# # Test the contrast
# test(emm)

# confint(pairs(emmeans::emmeans(best_model_species, ~disturbance)), level=.8)
# contrast               estimate lower.HPD upper.HPD
# multiple - agriculture -0.50686    -0.848   -0.1797 ***
# multiple - forest       0.00573    -0.329    0.2991 
# multiple - urban       -0.78797    -1.538   -0.0355 ***
# agriculture - forest    0.51150     0.227    0.8066 ***
# agriculture - urban    -0.28780    -1.004    0.4405
# forest - urban         -0.79394    -1.574   -0.0419 ***


# confint(pairs(emmeans::emmeans(best_model_traits, ~disturbance)), level=.8)

# contrast               estimate lower.HPD upper.HPD
# multiple - agriculture   -0.318    -0.699    0.0328
# multiple - forest         0.193    -0.118    0.5357
# multiple - urban         -0.711    -1.477    0.0683
# agriculture - forest      0.514     0.216    0.8171 ***
# agriculture - urban      -0.395    -1.110    0.3985
# forest - urban           -0.907    -1.669   -0.1095 ***

#'---------------------------------------------------------------------------------------------------------------
# @Main text figures ===== #################################################################################
#'---------------------------------------------------------------------------------------------------------------

pred_direction <- extract_predictions(
  species_model  = best_model_species, 
  trait_model = best_model_traits,
  predictor = "direction",
  levels_list = c("Differentiation","Homogenisation") 
)

pred_main_lu <- extract_predictions(
  species_model = best_model_species, 
  trait_model = best_model_traits,
  predictor = "main_land_use_type",
  levels_list = rev(c("Agriculture", "Forest", "Urban", "Multiple"))
)

# Define colors outside the individual figures
#blue_green <- c(Freshwater = "#118ab2", Terrestrial = "#06d6a0")
land_use_colors <- c(Agriculture = "#a36627", Forest = "#848c04", Urban = "#1c1c0c", Multiple = "#dc7c5c")
direction_colors = c(Differentiation = "#ff7f32",Homogenisation = "#9f5cc0")

# For Species Replacement - Direction (figure_3a1)

figure_3a1 <- plot_magnitude(
  pred_data = filter(pred_direction, facet == "Species replacement"),  # Use the "Species replacement" facet
  colors = direction_colors, 
  x_limits = c(0, 1.2), 
  x_breaks = c(0, 1.2), 
  sample_sizes = c(Differentiation = 118, Homogenisation = 42), 
  predictor = "direction",  
  facet_var = "Species replacement",
  text_height = .6,
  justification = 0,
  height= 2
)

# For Trait Replacement - Direction (figure_3a2)
figure_3a2 <- plot_magnitude(
  pred_data = pred_direction |> filter(facet == "Trait replacement"),  # Use the "Trait replacement" facet
  colors = direction_colors, 
  x_limits = c(0, 0.5), 
  x_breaks = c(0, 0.5), 
  sample_sizes = c(Differentiation = 82, Homogenisation = 78), 
  predictor = "direction",  
  facet_var = "Trait replacement",
  text_height = .6,
 justification = 0,
 height= 2
)

# For Species Replacement - Disturbance (figure_3b1)
figure_3b1 <- plot_magnitude(
  pred_data = pred_main_lu |> filter(facet == "Species replacement"),  # Use the "Species replacement" facet
  colors = land_use_colors, 
  x_limits = c(0, 1.2), 
  x_breaks = c(0, 1.2), 
  sample_sizes = c(Agriculture = 38, Forest =  55, Urban =  15, Multiple = 52), 
  predictor = "main_land_use_type",  
  facet_var = "Species replacement",
  text_height = .55,
  justification = 0,
  height= 1.6
)
# For Trait Replacement - Disturbance (figure_3b2)
figure_3b2 <- plot_magnitude(
  pred_data = pred_main_lu |> filter(facet == "Trait replacement"),  # Use the "Trait replacement" facet
  colors = land_use_colors, 
  x_limits = c(0, 0.5), 
  x_breaks = c(0, 0.5), 
  sample_sizes =  c(Agriculture = 38, Forest =  55, Urban =  15, Multiple = 52), 
  predictor = "main_land_use_type",  
  facet_var = "Trait replacement",
  text_height = .55,
  justification = 0,
  height= 1.6
)

combined_figure_3 <- cowplot::plot_grid(
  figure_3a1 + labs(x="",tag="a"),figure_3b1 + labs(x="",tag="b"),  # First row
  figure_3a2, figure_3b2,  # Second row
  ncol = 2)

ggsave(plot=combined_figure_3
       ,filename=here::here("S7_Model_outputs_figures_and_tables", "main_figures", "Figure_3.pdf")
       ,width = 89, height=50, units="mm")

# @Extended Data ---------------------------------------------------------------------------------------------------
#Don't worry about the warnings, they are not relevant to us.

magnitude_figure_data <- generate_figure_data(best_model_species, best_model_traits, "magnitude", draw_min = -1, draw_max = 4)
extended_figure_3 <-plot_posterior_draws(magnitude_figure_data, "magnitude",trunc_upper = 3, trunc_lower = -1, tick_factor = 1, xlims = c(-1,4))


ggsave(filename = here::here("S7_Model_outputs_figures_and_tables", "extended_data", "Extended_Figure_3.pdf"),
       plot = extended_figure_3,
       device = cairo_pdf(),
       width =  80,  # Adjust width for shape figures
       height = 40, units = "mm")