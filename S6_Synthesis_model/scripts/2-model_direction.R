#' ###########################################################################
#' SCRIPT NAME: Bayesian Modeling of Response Direction (2-model_direction.R)
#'
#' DESCRIPTION:
#'   This script performs Bayesian modeling to analyze the direction of species 
#'   and trait turnover across different ecosystems and land-use types. 
#'   It also generates exploratory figures and saves model outputs for further
#'   analysis and visualization.
#'
#' USAGE:
#'   - The script reads preprocessed data and runs Bayesian models using the 
#'     `brms` package.
#'   - Several models with varying complexity (base, interaction, nested taxa) 
#'     are evaluated using Leave-One-Out Cross-Validation (LOO).
#'   - Generates the figures for inclusion in the publication.
#'
#' INPUTS:
#'   - Preprocessed synthesis data from `S6_Synthesis_model/data/synthesis_data.xlsx`
#'   - Helper functions for model processing and plotting
#'
#' OUTPUTS:
#'   - Bayesian models saved in `S7_Model_outputs_figures_and_tables/model/direction`
#'   - Figures saved in `S7_Model_outputs_figures_and_tables/main_figures` and
#'     `S7_Model_outputs_figures_and_tables/extended_data`
#'
#' AUTHOR: Caio Graco-Roza
#' LAST UPDATED: 2024-11-24
#'
#' NOTES:
#'   - Adjust the global settings for the Bayesian models as needed (chains, iterations, etc.).
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

# All models will be run with same settings
CHAINS <- 4                          # Number of chains to run for the Bayesian model
ITER <- 5000                         # Total number of iterations per chain
WARMUP <- 2500                       # Number of iterations used for warmup (not included in posterior samples)
BAYES_SEED <- 1234                   # Seed for reproducibility in Bayesian modeling


model_data <- process_model_data("S6_Synthesis_model/data/synthesis_data.xlsx", "direction")

# model_data |>
#   tabyl(direction,facet) |>
#   adorn_percentages("col") %>%
#   adorn_pct_formatting(digits = 2) %>%
#   adorn_ns() |>
#   kableExtra::kable(format = "rst")

# ===============  ============  ===========
# direction        Taxonomic     Functional 
# ===============  ============  ===========
# differentiation  73.75% (118)  55.62% (89)
# homogenisation   26.25%  (42)  44.38% (71)
# ===============  ============  ===========


# model_data |>
#   tabyl(direction,ecosystem_type, facet) |>
#   adorn_percentages("col") %>%
#   adorn_pct_formatting(digits = 2) %>%
#   adorn_ns() |>
#   kableExtra::kable(format = "rst")

#  Species
# ===============  ===========  ===========
# direction        terrestrial  freshwater 
# ===============  ===========  ===========
# differentiation  74.62% (97)  70.00% (21)
# homogenisation   25.38% (33)  30.00%  (9)
# ===============  ===========  ===========
# 
# Traits
# ===============  ===========  ===========
# direction        terrestrial  freshwater 
# ===============  ===========  ===========
# differentiation  53.08% (69)  66.67% (20)
# homogenisation   46.92% (61)  33.33% (10)
# ===============  ===========  ===========

# model_data |>
#   tabyl(direction,facet,main_land_use_type) |>
#   adorn_percentages("col") %>%
#   adorn_pct_formatting(digits = 2) %>%
#   adorn_ns()|>
#   bind_rows(.id="main_land_use_type")
  
#   
#                   Species	       Traits
# =============================================
# forest
# 
# differentiation	  69.09% (38)	   43.64% (24)
# homogenisation	  30.91% (17)	   56.36% (31)
# =============================================
# multiple
# 
# differentiation	  78.85% (41)	   67.31% (35)
# homogenisation	  21.15% (11)	   32.69% (17)
# ============================================
# agriculture     
# 
# differentiation	  78.95% (30)	   68.42% (26)
# homogenisation	  21.05% (8)	   31.58% (12)
# ============================================
# urban     
# 
# differentiation	  60.00% (9)	   26.67% (4)
# homogenisation	  40.00% (6)	   73.33% (11)
# ============================================

#Set priors --------------------------------------------------------------

priors_direction <- c(
  prior(normal(0,3), class = b) #slopes
 ,prior(student_t(3, 0, 1), class = sd, lb=0) #standard deviation
)

#Set model formula --------------------------

# Base model with main effects only
formula_base <- bf(
  direction_harrel_davis ~ species_number + spatial_extent + spatial_distance_min + 
    absolute_latitude_mean + human_pressure_baseline + human_pressure_range + 
    main_land_use_type + (1 | taxa_coarse) + buffer + ecosystem_type, 
  decomp = "QR"
)  # Base model

# Model with ecosystem and land use interaction
formula_interaction <- bf(
  direction_harrel_davis ~ species_number + spatial_extent + spatial_distance_min + 
    absolute_latitude_mean + human_pressure_baseline + human_pressure_range + 
    ecosystem_type * main_land_use_type + (1 | taxa_coarse) + buffer, 
  decomp = "QR"
)  # Adds ecosystem * land use interaction

# Model with slope for fine taxa nested within coarse taxa
formula_taxa_nested <- bf(
  direction_harrel_davis ~ species_number + spatial_extent + spatial_distance_min + 
    absolute_latitude_mean + human_pressure_baseline + human_pressure_range + 
    ecosystem_type + main_land_use_type + (1 | taxa_coarse/taxa_fine) + buffer, 
  decomp = "QR"
)  # Adds nested taxa slope

# Model with interaction and nested taxa slope
formula_interaction_nested <- bf(
  direction_harrel_davis ~ species_number + spatial_extent + spatial_distance_min + 
    absolute_latitude_mean + human_pressure_baseline + human_pressure_range + 
    ecosystem_type * main_land_use_type + (1 | taxa_coarse/taxa_fine) + buffer, 
  decomp = "QR"
)  # Adds interaction and nested taxa slope

#'----------------------------------------------------

#Run model ----

# Base Model with Main Effects Only
md_species_base <- brms::brm(
  formula_base,
  data = filter(model_data, facet == "Species"),
  prior = priors_direction,
  family = student(link = "identity"),
  chains = CHAINS,
  iter = ITER,
  warmup = WARMUP,
  control = list(max_treedepth = 12, adapt_delta = .99),
  save_pars = save_pars(all = TRUE),
  seed = BAYES_SEED,
  file = "S7_Model_outputs_figures_and_tables/model/direction/direction_species_base",
  file_refit = "on_change"
)

# Models with Increased Complexity, saved under distinct names
md_species_interaction <- update(
  md_species_base,
  newdata = filter(model_data, facet == "Species"),
  formula = formula_interaction,
  file = "S7_Model_outputs_figures_and_tables/model/direction/direction_species_interaction",
  file_refit = "on_change"
)

md_species_taxa_nested <- update(
  md_species_base,
  newdata = filter(model_data, facet == "Species"),
  formula = formula_taxa_nested,
  file = "S7_Model_outputs_figures_and_tables/model/direction/direction_species_taxa_nested",
  file_refit = "on_change"
)

md_species_interaction_nested <- update(
  md_species_base,
  newdata = filter(model_data, facet == "Species"),
  formula = formula_interaction_nested,
  file = "S7_Model_outputs_figures_and_tables/model/direction/direction_species_interaction_nested",
  file_refit = "on_change"
)

# # Add Model Selection Criteria
md_species_base <- add_criterion(md_species_base, "loo", moment_match = FALSE)
md_species_interaction <- add_criterion(md_species_interaction, "loo", moment_match = FALSE)
md_species_taxa_nested <- add_criterion(md_species_taxa_nested, "loo", moment_match = FALSE)
md_species_interaction_nested <- add_criterion(md_species_interaction_nested, "loo", moment_match = FALSE)

# Compare Models Using Leave-One-Out Cross-Validation
loo_compare_results <- loo(md_species_base, md_species_interaction, md_species_taxa_nested, md_species_interaction_nested)

# Base model for traits
md_traits_base <- brms::brm(
  formula_base,
  data = filter(model_data, facet == "Traits"),
  prior = priors_direction,
  chains = CHAINS,
  family = brms::student(),
  iter = ITER,
  warmup = WARMUP,
  control = list(max_treedepth = 12, adapt_delta = .99),
  save_pars = save_pars(all = TRUE),
  seed = BAYES_SEED,
  file = "S7_Model_outputs_figures_and_tables/model/direction/direction_traits_base",
  file_refit = "on_change"
)

# Model with interaction, saved as a unique name
md_traits_interaction <- update(
  md_traits_base,
  newdata = filter(model_data, facet == "Traits"),
  formula = formula_interaction,
  file = "S7_Model_outputs_figures_and_tables/model/direction/direction_traits_interaction",
  file_refit = "on_change"
)

# Model with nested taxa slope, saved as a unique name
md_traits_taxa_nested <- update(
  md_traits_base,
  newdata = filter(model_data, facet == "Traits"),
  formula = formula_taxa_nested,
  file = "S7_Model_outputs_figures_and_tables/model/direction/direction_traits_taxa_nested",
  file_refit = "on_change"
)

# Model with interaction and nested taxa slope, saved as a unique name
md_traits_interaction_nested <- update(
  md_traits_base,
  newdata = filter(model_data, facet == "Traits"),
  formula = formula_interaction_nested,
  file = "S7_Model_outputs_figures_and_tables/model/direction/direction_traits_interaction_nested",
  file_refit = "on_change"
)

# Add criteria for LOO
md_traits_base <- add_criterion(md_traits_base, "loo", moment_match = FALSE)
md_traits_interaction <- add_criterion(md_traits_interaction, "loo", moment_match = FALSE)
md_traits_taxa_nested <- add_criterion(md_traits_taxa_nested, "loo", moment_match = FALSE)
md_traits_interaction_nested <- add_criterion(md_traits_interaction_nested, "loo", moment_match = FALSE)

# LOO comparison
loo(md_traits_base, md_traits_interaction, md_traits_taxa_nested, md_traits_interaction_nested)

best_model_species <- md_species_base # md_species_base is chosen as the final model because it achieves similar predictive accuracy to more complex models but with fewer parameters
best_model_traits <- md_traits_base # md_traits_base is chosen as the final model because it achieves similar predictive accuracy but with fewer parameters

 # Run the posterior predictive checks -----------------------------------------------------------------------------
cowplot::plot_grid(brms::pp_check(best_model_species, ndraws=100), brms::pp_check(best_model_traits, ndraws=100))
#'-----------------------------------------------------------------------------------------------------------------

# check R2 of each model ------------------------------------------------------------------------------------------
brms::bayes_R2(best_model_species) # Taxonomic  0.11
brms::bayes_R2(best_model_traits) # Functional 0.07
#'-----------------------------------------------------------------------------------------------------------------

# Make summary table ----------------------------------------------------------------------------------------------

# confint(pairs(emmeans::emmeans(best_model_species, ~main_land_use_type)), level=.8)

#SPECIES TURNOVER
# contrast               estimate lower.HPD upper.HPD
# multiple - agriculture -0.00553   -0.0267    0.0155
# multiple - forest      -0.00706   -0.0288    0.0141
# multiple - urban       -0.08951   -0.1391   -0.0365 **
# agriculture - forest   -0.00098   -0.0213    0.0186
# agriculture - urban    -0.08463   -0.1332   -0.0348 **
# forest - urban         -0.08247   -0.1347   -0.0303 **

# confint(pairs(emmeans::emmeans(best_model_traits, ~main_land_use_type)), level=.8)

#TRAIT TURNOVER
# contrast               estimate lower.HPD upper.HPD
# multiple - agriculture  -0.0153   -0.0339   0.00261
# multiple - forest       -0.0266   -0.0434  -0.00878 **
# multiple - urban        -0.0494   -0.0892  -0.00918 **
# agriculture - forest    -0.0111   -0.0259   0.00260
# agriculture - urban     -0.0339   -0.0716   0.00525
# forest - urban          -0.0229   -0.0626   0.01604

#Exploratory figures -------------------------------------------------------------
#Species turnover

# Predictions ecosystem_type ----

pred_eco_type <- extract_predictions(
  species_model = best_model_species, 
  trait_model = best_model_traits,
  predictor = "ecosystem_type",
  levels_list = c("Freshwater","Terrestrial") 
)

pred_main_lu <- extract_predictions(
  species_model = best_model_species, 
  trait_model = best_model_traits,
  predictor = "main_land_use_type",
  levels_list = rev(c("Agriculture", "Forest", "Urban", "Multiple"))
)

# Define colors outside the individual figures
blue_green <- c(Freshwater = "#118ab2", Terrestrial = "#06d6a0")
land_use_colors <- c(Agriculture = "#a36627", Forest = "#848c04", Urban = "#1c1c0c", Multiple = "#dc7c5c")

# Example usage for two figures:
figure_2a <- plot_direction(pred_data = pred_eco_type, model_data = model_data, predictor_col = "ecosystem_type", color_map = blue_green, vjust_text = 0.65, height = 1, justification = 0.15,width=1.5)
figure_2b <- plot_direction(pred_data = pred_main_lu, model_data = model_data, predictor_col = "main_land_use_type", color_map = land_use_colors, vjust_text = 0.55, height=1, justification=0.7,width=5)


# Reduce spacing between Figure 2A and Figure 2B
Figure_2 <- plot_grid(
  figure_2a,           # Top plot
  figure_2b,           # Bottom plot
  ncol = 2,            # Arrange vertically
  rel_heights = c(0.5, 0.5),  # Adjust relative heights (equal height for both plots)
  align = "v",         # Align vertically
  labels = c("a", "b"),  # Labels for subfigures
  label_size = 5,      # Control size of labels
  label_fontface = "bold"  # Bold the labels
)

# Saving the plot for half-page width
ggsave(
  plot = Figure_2,
  filename = here::here("S7_Model_outputs_figures_and_tables", "main_figures", "Figure_2.pdf"),
  device = cairo_pdf(),
  width = 89,  # Half-page width in mm
  height = 40, units = "mm"
)

#Extended Data

direction_figure_data <- generate_figure_data(best_model_species, best_model_traits, "direction", draw_min = -3, draw_max =8 )

extended_figure_2<-plot_posterior_draws(direction_figure_data, "direction",trunc_upper = 7.5, trunc_lower = -2.5, tick_factor = 2.5)

ggsave(filename = here::here("S7_Model_outputs_figures_and_tables", "extended_data", "Extended_Figure_2.pdf"),
       plot = extended_figure_2,
       device = cairo_pdf(),
       width =  95,  # Adjust width for shape figures
       height = 40, units = "mm")
