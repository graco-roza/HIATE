#'###############################################################################
#' SCRIPT NAME: Modeling and Visualizing Relationship Shapes
#'
#' DESCRIPTION:
#'   This script models the shape of relationships between species and their 
#'   functional traits and environmental predictors. It includes testing various
#'   interaction effects and conducting K-fold cross-validation for model 
#'   comparison due to the multinomial nature of the response variable, which 
#'   precludes the use of standard LOO-CV (Leave-One-Out Cross-Validation).
#'   Additionally, the script generates main text and extended figures for 
#'   scientific reporting.
#'
#' USAGE:
#'   - This script loads data, defines model formulas, fits Bayesian models 
#'     using `brms`, evaluates models using K-fold cross-validation, and 
#'     produces visualizations of the results.
#'   - Save paths are defined for all models and outputs to facilitate re-use.
#'
#' INPUTS:
#'   - Synthesis data file: `S6_Synthesis_model/data/synthesis_data.xlsx`
#'   - Helper functions for modeling and plotting: 
#'       * `S6_Synthesis_model/functions/helper_functions_S6.R`
#'       * `S6_Synthesis_model/functions/helper_functions_plot_extended.R`
#'
#' OUTPUTS:
#'   - Fitted models saved in `S7_Model_outputs_figures_and_tables/model/shape/`
#'   - Main figures saved in `S7_Model_outputs_figures_and_tables/main_figures/`
#'   - Extended figures saved in `S7_Model_outputs_figures_and_tables/extended_data/`
#'   - K-fold cross-validation results for species and traits models saved in:
#'       * `kfold_results_species.RData`
#'       * `kfold_results_traits.RData`
#'
#' AUTHOR: Caio Graco-Roza
#' LAST UPDATED: 2024-11-24
#'
#' NOTES:
#'   - Multinomial models are fit with Bayesian methods using the `brms` package.
#'   - K-fold cross-validation is used instead of LOO due to the multinomial 
#'     response variable.
#'   - Visualization of the modeled relationship shapes includes interactions 
#'     between predictors such as land-use types, ecosystem types, and 
#'     homogenization/differentiation.
#'###############################################################################

# Load packages ####
library(rstan)
library(tidyverse) # For data manipulation and visualization
library(ggdist)  # For visualizing distributions of data and model outputs
library(ggtext)  # For adding rich text elements to ggplot2 plots
library(brms) # For running Bayesian models
library(distributional) # For defining and visualizing probability distributions, including priors
library(marginaleffects) # For extracting marginal effects and predictions from Bayesian models
library(future) # For running k-fold cross-validation in parallel using multiple cores
library(tidybayes) # For tidy workflows and working with Bayesian models
library(janitor) # For cleaning and organizing data tables, and standardizing column names
library(bayestestR) # For assessing Bayesian models and extracting posterior summaries
library(ggh4x) # For additional ggplot2 features like nested facets
library(ggalt) # For additional ggplot2 features such as enhanced line and ribbon plots
library(glue) # For constructing strings and paths dynamically

# Load helper scripts ####
# These scripts contain custom functions used in the analysis and plotting.
source("S6_Synthesis_model/functions/helper_functions_S6.R")
source("S6_Synthesis_model/functions/helper_functions_plot_extended.R")

# Set custom options ####
# Setting default directory for Stan files
options(cmdstanr_write_stan_file_dir = "S7_Model_outputs_figures_and_tables/model/shape")

# Global Stan options for Bayesian models ####
# Set the number of cores for parallel sampling and specify the CMDSTAN backend
base::options(mc.cores = 4,               # Use 5 cores for parallel processing
              brms.backend = "cmdstanr")  # Use cmdstanr for running Bayesian models

# Define global Bayesian model settings ####
CHAINS <- 4          # Number of chains for MCMC sampling
ITER <- 1500         # Total number of iterations per chain
WARMUP <- ITER/2       # Number of warmup iterations (not included in posterior)
BAYES_SEED <- 1234   # Seed for reproducibility in Bayesian sampling

# Load and process model data ####
# Load the synthesis data for modeling the shape of responses
model_data <- process_model_data("S6_Synthesis_model/data/synthesis_data.xlsx", "shape")

#Description of the shapes ----
# model_data |> 
#   select(dataset, direction, facet,Absent, Revlog, Exponential, Saturating, realm) |>
#   #mutate(N = Absent + Exponential + Saturating + Revlog)  |> 
#   pivot_longer(cols=!c(facet,dataset,direction,realm)) |> 
#   group_by(dataset,direction,facet) |> 
#   mutate(value = ifelse(dataset == "bexis_grassland_sweepnet_beetles" & facet == "Taxonomic" & name == "Saturating", 311,value)) |> 
#   slice_max(value) |> 
#   select(-value) |> 
#   ungroup() |> 
#   select(-direction) |> 
#   #pivot_wider(names_from=facet,values_from = name) |> 
#   tabyl(name,facet) |>
#   adorn_percentages("col") %>%
#   adorn_pct_formatting(digits = 2) %>%
#   adorn_ns() |>
#   kableExtra::kable(format = "rst")
# 
# # ===========  ===========  ===========
# # Shape        Taxonomic    Functional 
# # ===========  ===========  ===========
# # Absent       20.00% (32)  21.88% (35)
# # Exponential  24.38% (39)  34.38% (55)
# # Revlog       11.88% (19)  8.12% (13) 
# # Saturating   43.75% (70)  35.62% (57)
# # ===========  ===========  ===========
#   
# 
# model_data |> 
#   select(dataset, direction, facet,Absent, Revlog, Exponential, Saturating, realm) |>
#   pivot_longer(cols=!c(facet,dataset,direction,realm)) |> 
#   group_by(dataset,direction,facet) |> 
#   mutate(value = ifelse(dataset == "bexis_grassland_sweepnet_beetles" & facet == "Taxonomic" & name == "Saturating", 311,value)) |> 
#   slice_max(value) |> 
#   select(-value) |> 
#   tabyl(name,realm,facet) |>
#   adorn_percentages("col") %>%
#   adorn_pct_formatting(digits = 2) %>%
#   adorn_ns() |>
#   kableExtra::kable(format = "rst")
# 
# # ===========  ===========  ===========
# # Taxonomic    terrestrial  freshwater 
# # ===========  ===========  ===========
# # Absent       18.46% (24)  26.67%  (8)
# # Exponential  27.69% (36)  10.00%  (3)
# # Revlog       10.00% (13)  20.00%  (6)
# # Saturating   43.85% (57)  43.33% (13)
# # ===========  ===========  ===========
# # 
# # ===========  ===========  ===========
# # Functional   terrestrial  freshwater 
# # ===========  ===========  ===========
# # Absent       22.31% (29)  20.00%  (6)
# # Exponential  32.31% (42)  43.33% (13)
# # Revlog       9.23% (12)   3.33%  (1) 
# # Saturating   36.15% (47)  33.33% (10)
# # ===========  ===========  ===========
# 
# model_data |> 
#     select(dataset, direction, facet,Absent, Revlog, Exponential, Saturating, disturbance) |>
#     #mutate(N = Absent + Exponential + Saturating + Revlog)  |> 
#     pivot_longer(cols=!c(facet,dataset,direction,disturbance)) |> 
#     group_by(dataset,direction,facet) |> 
#     mutate(value = ifelse(dataset == "bexis_grassland_sweepnet_beetles" & facet == "Taxonomic" & name == "Saturating", 311,value)) |> 
#     slice_max(value) |> 
#     select(-value) |> 
#     tabyl(name,disturbance,facet) |>
#     adorn_percentages("col") %>%
#     adorn_pct_formatting(digits = 2) %>%
#     adorn_ns() |>
#     kableExtra::kable(format = "rst")
# 
# # ===========  ===========  ===========  ===========  ==========
# # Taxonomic    multiple     agriculture  forest       urban     
# # ===========  ===========  ===========  ===========  ==========
# # Absent       28.85% (15)  10.53%  (4)  23.64% (13)  0.00% (0) 
# # Exponential  21.15% (11)  23.68%  (9)  20.00% (11)  53.33% (8)
# # Revlog       11.54%  (6)  18.42%  (7)  5.45%  (3)   20.00% (3)
# # Saturating   38.46% (20)  47.37% (18)  50.91% (28)  26.67% (4)
# # ===========  ===========  ===========  ===========  ==========
# # 
# # ===========  ===========  ===========  ===========  ==========
# # Functional   multiple     agriculture  forest       urban     
# # ===========  ===========  ===========  ===========  ==========
# # Absent       32.69% (17)  15.79%  (6)  21.82% (12)  0.00% (0) 
# # Exponential  28.85% (15)  39.47% (15)  29.09% (16)  60.00% (9)
# # Revlog       3.85%  (2)   15.79%  (6)  3.64%  (2)   20.00% (3)
# # Saturating   34.62% (18)  28.95% (11)  45.45% (25)  20.00% (3)
# #===========  ===========  ===========  ===========  ==========
# #Description of the human pressure across main land-use type ----
# 
# "S6_Synthesis_model/data/synthesis_data.xlsx" |>
#   readxl::read_xlsx() |>
#   filter(!dataset %in% c("N67TTP","N78TTP","S47TTP")) |>
#   group_by(disturbance) |>
#   summarise(mean=mean(hfp.min),
#             sd=sd(hfp.min, na.rm=TRUE)) |>
#   kableExtra::kable(format = "rst")

# # ===========  =========  ========
# # disturbance       mean        sd
# # ===========  =========  ========
# # agriculture   5.788789  3.731489
# # forest        4.592429  3.521031
# # multiple      6.015048  5.462427
# # urban        31.598707  3.083456
# # ===========  =========  ========


# Base formula for shape, testing the effects of main predictors
formula_base <- bf(
  shape | trials(trials) ~ species_number + spatial_extent + spatial_distance_min + 
    absolute_latitude_mean + human_pressure_baseline + human_pressure_range + 
    main_land_use_type + direction + ecosystem_type + buffer + (1 | taxa_coarse),
  decomp = "QR"
)

priors_shape <- c(
  brms::prior(normal(0,1), class = Intercept,  dpar="murevlog")
 ,brms::prior(normal(0,1), class = b, dpar="murevlog")
 ,brms::prior(student_t(3, 0, 0.5), class=sd, dpar="murevlog")
 ,brms::prior(normal(0,1), class = Intercept,  dpar="musaturating")
 ,brms::prior(normal(0,1), class = b, dpar="musaturating")
 ,brms::prior(student_t(3, 0, 0.5), class=sd, dpar="musaturating")
 ,brms::prior(normal(0,1), class = Intercept,  dpar="muexponential")
 ,brms::prior(normal(0,1), class = b, dpar="muexponential")
 ,brms::prior(student_t(3, 0, 0.5), class=sd, dpar="muexponential")
)

# Base model for species
model_base_species <- brms::brm(
  formula = formula_base,
  data = model_data |> filter(facet == "Species"),
  family = brms::multinomial(link = "logit", refcat = "absent"),
  prior = priors_shape,
  chains = CHAINS,
  iter = ITER,
  warmup = WARMUP,
  control = base::list(max_treedepth = 15, adapt_delta = .99),
  save_pars = brms::save_pars(all = TRUE),
  seed = BAYES_SEED,
  file = "S7_Model_outputs_figures_and_tables/model/shape/model_base_species", # Automatically saves model
  file_refit = "on_change"
)

# Update models to test complexity
# model_direction_land_use_species <- update(
#   model_base_species,
#   newdata = model_data |> filter(facet == "Species"),
#   formula. = bf(update(shape_formula_base, direction * main_land_use_type)),
#   file = "S7_Model_outputs_figures_and_tables/model/shape/model_direction_land_use_species"
# )
# 
# model_direction_ecosystem_species <- update(
#   model_base_species,
#   newdata = model_data |> filter(facet == "Species"),
#   formula. = bf(update(shape_formula_base, direction * ecosystem_type)),
#   file = "S7_Model_outputs_figures_and_tables/model/shape/model_direction_ecosystem_species"
# )
# 
# model_ecosystem_land_use_species <- update(
#   model_base_species,
#   newdata = model_data |> filter(facet == "Species"),
#   formula. = bf(update(shape_formula_base, ecosystem_type * main_land_use_type)),
#   file = "S7_Model_outputs_figures_and_tables/model/shape/model_ecosystem_land_use_species"
# )
# 
# model_three_way_interaction_species <- update(
#   model_base_species,
#   newdata = model_data |> filter(facet == "Species"),
#   formula. = bf(update(shape_formula_base, direction * ecosystem_type * main_land_use_type)),
#   file = "S7_Model_outputs_figures_and_tables/model/shape/model_three_way_interaction_species"
# )
# 
# model_base_nested_species <- update(
#   model_base_species,
#   newdata = model_data |> filter(facet == "Species"),
#   formula. = bf(update(shape_formula_base, - (1|taxa_corase) + (1 | taxa_coarse / taxa_fine))), #remove the old structure and add the new random effect
#   file = "S7_Model_outputs_figures_and_tables/model/shape/model_base_nested_species"
# )
# 
# model_direction_land_use_nested_species <- update(
#   model_base_nested_species,
#   newdata = model_data |> filter(facet == "Species"),
#   formula. = bf(update(shape_formula_base, direction * main_land_use_type)),
#   file = "S7_Model_outputs_figures_and_tables/model/shape/model_direction_land_use_nested_species"
# )
# 
# model_direction_ecosystem_nested_species <- update(
#   model_base_nested_species,
#   newdata = model_data |> filter(facet == "Species"),
#   formula. = bf(update(shape_formula_base, direction * ecosystem_type)),
#   file = "S7_Model_outputs_figures_and_tables/model/shape/model_direction_ecosystem_nested_species"
# )
# 
# model_ecosystem_land_use_nested_species <- update(
#   model_base_nested_species,
#   newdata = model_data |> filter(facet == "Species"),
#   formula. = bf(update(shape_formula_base, ecosystem_type * main_land_use_type)),
#   file = "S7_Model_outputs_figures_and_tables/model/shape/model_ecosystem_land_use_nested_species"
# )
# 
# model_three_way_interaction_nested_species <- update(
#   model_base_nested_species,
#   newdata = model_data |> filter(facet == "Species"),
#   formula. = bf(update(shape_formula_base, direction * ecosystem_type * main_land_use_type)),
#   file = "S7_Model_outputs_figures_and_tables/model/shape/model_three_way_interaction_nested_species"
# )

# # K-fold Cross-Validation for Species Models
# plan(multisession)
# kf_base_species <- kfold(model_base_species, chains = 1, K = 3)
# plan(multisession)
# kf_direction_land_use_species <- kfold(model_direction_land_use_species, chains = 1, K = 3)
# plan(multisession)
# kf_direction_ecosystem_species <- kfold(model_direction_ecosystem_species, chains = 1, K = 3)
# plan(multisession)
# kf_ecosystem_land_use_species <- kfold(model_ecosystem_land_use_species, chains = 1, K = 3)
# plan(multisession)
# kf_three_way_interaction_species <- kfold(model_three_way_interaction_species, chains = 1, K = 3)
# plan(multisession)
# kf_base_nested_species <- kfold(model_base_nested_species, chains = 1, K = 3)
# plan(multisession)
# kf_direction_land_use_nested_species <- kfold(model_direction_land_use_nested_species, chains = 1, K = 3)
# plan(multisession)
# kf_direction_ecosystem_nested_species <- kfold(model_direction_ecosystem_nested_species, chains = 1, K = 3)
# plan(multisession)
# kf_ecosystem_land_use_nested_species <- kfold(model_ecosystem_land_use_nested_species, chains = 1, K = 3)
# plan(multisession)
# kf_three_way_interaction_nested_species <- kfold(model_three_way_interaction_nested_species, chains = 1, K = 3)
#
# save(
#   kf_base_species,
#   kf_direction_land_use_species,
#   kf_direction_ecosystem_species,
#   kf_ecosystem_land_use_species,
#   kf_three_way_interaction_species,
#   kf_base_nested_species,
#   kf_direction_land_use_nested_species,
#   kf_direction_ecosystem_nested_species,
#   kf_ecosystem_land_use_nested_species,
#   kf_three_way_interaction_nested_species,
#   file = "S7_Model_outputs_figures_and_tables/model/shape/kfold_results_species.RData"
# )
# 
# # Compare K-fold results for Species Models
# loo_compare(
#   kf_base_species,
#   kf_direction_land_use_species,
#   kf_direction_ecosystem_species,
#   kf_ecosystem_land_use_species,
#   kf_three_way_interaction_species,
#   kf_base_nested_species,
#   kf_direction_land_use_nested_species,
#   kf_direction_ecosystem_nested_species,
#   kf_ecosystem_land_use_nested_species,
#   kf_three_way_interaction_nested_species
# )


# Base model for traits
model_base_traits <- brms::brm(
  formula = formula_base,
  data = model_data |> filter(facet == "Traits"),
  family = brms::multinomial(link = "logit", refcat = "absent"),
  prior = priors_shape,
  chains = CHAINS,
  iter = ITER,
  warmup = WARMUP,
  control = base::list(max_treedepth = 15, adapt_delta = .99),
  save_pars = brms::save_pars(all = TRUE),
  seed = BAYES_SEED,
  file = "S7_Model_outputs_figures_and_tables/model/shape/model_base_traits", # Automatically saves model
  file_refit = "on_change"
)

# # Update models to test complexity
# model_direction_land_use_traits <- update(
#   model_base_traits,
#   newdata = model_data |> filter(facet == "Traits"),
#   formula. = bf(update(shape_formula_base, direction * main_land_use_type)),
#   file = "S7_Model_outputs_figures_and_tables/model/shape/model_direction_land_use_traits"
# )
# 
# model_direction_ecosystem_traits <- update(
#   model_base_traits,
#   newdata = model_data |> filter(facet == "Traits"),
#   formula. = bf(update(shape_formula_base, direction * ecosystem_type)),
#   file = "S7_Model_outputs_figures_and_tables/model/shape/model_direction_ecosystem_traits"
# )
# 
# model_ecosystem_land_use_traits <- update(
#   model_base_traits,
#   newdata = model_data |> filter(facet == "Traits"),
#   formula. = bf(update(shape_formula_base, ecosystem_type * main_land_use_type)),
#   file = "S7_Model_outputs_figures_and_tables/model/shape/model_ecosystem_land_use_traits"
# )
# 
# model_three_way_interaction_traits <- update(
#   model_base_traits,
#   newdata = model_data |> filter(facet == "Traits"),
#   formula. = bf(update(shape_formula_base, direction * ecosystem_type * main_land_use_type)),
#   file = "S7_Model_outputs_figures_and_tables/model/shape/model_three_way_interaction_traits"
# )
# 
# model_base_nested_traits <- update(
#   model_base_traits,
#   newdata = model_data |> filter(facet == "Traits"),
#   formula. = bf(update(shape_formula_base, - (1 | taxa_coarse) + (1 | taxa_coarse / taxa_fine))), # Remove old structure and add nested random effect
#   file = "S7_Model_outputs_figures_and_tables/model/shape/model_base_nested_traits"
# )
# 
# model_direction_land_use_nested_traits <- update(
#   model_base_nested_traits,
#   newdata = model_data |> filter(facet == "Traits"),
#   formula. = bf(update(shape_formula_base, direction * main_land_use_type)),
#   file = "S7_Model_outputs_figures_and_tables/model/shape/model_direction_land_use_nested_traits"
# )
# 
# model_direction_ecosystem_nested_traits <- update(
#   model_base_nested_traits,
#   newdata = model_data |> filter(facet == "Traits"),
#   formula. = bf(update(shape_formula_base, direction * ecosystem_type)),
#   file = "S7_Model_outputs_figures_and_tables/model/shape/model_direction_ecosystem_nested_traits"
# )
# 
# model_ecosystem_land_use_nested_traits <- update(
#   model_base_nested_traits,
#   newdata = model_data |> filter(facet == "Traits"),
#   formula. = bf(update(shape_formula_base, ecosystem_type * main_land_use_type)),
#   file = "S7_Model_outputs_figures_and_tables/model/shape/model_ecosystem_land_use_nested_traits"
# )
# 
# model_three_way_interaction_nested_traits <- update(
#   model_base_nested_traits,
#   newdata = model_data |> filter(facet == "Traits"),
#   formula. = bf(update(shape_formula_base, direction * ecosystem_type * main_land_use_type)),
#   file = "S7_Model_outputs_figures_and_tables/model/shape/model_three_way_interaction_nested_traits"
# )

# # K-fold Cross-Validation for Trait Models
# plan(multisession)
# kf_base_traits <- kfold(model_base_traits, chains = 1, K = 3)
# plan(multisession)
# kf_direction_land_use_traits <- kfold(model_direction_land_use_traits, chains = 1, K = 3)
# plan(multisession)
# kf_direction_ecosystem_traits <- kfold(model_direction_ecosystem_traits, chains = 1, K = 3)
# plan(multisession)
# kf_ecosystem_land_use_traits <- kfold(model_ecosystem_land_use_traits, chains = 1, K = 3)
# plan(multisession)
# kf_three_way_interaction_traits <- kfold(model_three_way_interaction_traits, chains = 1, K = 3)
# plan(multisession)
# kf_base_nested_traits <- kfold(model_base_nested_traits, chains = 1, K = 3)
# plan(multisession)
# kf_direction_land_use_nested_traits <- kfold(model_direction_land_use_nested_traits, chains = 1, K = 3)
# plan(multisession)
# kf_direction_ecosystem_nested_traits <- kfold(model_direction_ecosystem_nested_traits, chains = 1, K = 3)
# plan(multisession)
# kf_ecosystem_land_use_nested_traits <- kfold(model_ecosystem_land_use_nested_traits, chains = 1, K = 3)
# plan(multisession)
# kf_three_way_interaction_nested_traits <- kfold(model_three_way_interaction_nested_traits, chains = 1, K = 3)
# 
# 
# # Compare K-fold results for Trait Models
# loo_compare(
#   kf_base_traits,
#   kf_direction_land_use_traits,
#   kf_direction_ecosystem_traits,
#   kf_ecosystem_land_use_traits,
#   kf_three_way_interaction_traits,
#   kf_base_nested_traits,
#   kf_direction_land_use_nested_traits,
#   kf_direction_ecosystem_nested_traits,
#   kf_ecosystem_land_use_nested_traits,
#   kf_three_way_interaction_nested_traits
# )



best_model_species <- model_base_species
best_model_traits <- model_base_traits # this was not the case here, but we use this only for consistency and mention in the main paper that there interaction was found.

# Create facet labels for shapes ####
# This section prepares labels for each shape type, including species and trait counts.
facet_labels <- tribble(
  ~shape, ~species_count, ~trait_count,                # Column names
  "Exponential", 39, 55,                               # Exponential shape with counts
  "Saturating", 70, 57,                                # Saturating shape with counts
  "Reverse logistic", 19, 13,                          # Reverse logistic shape with counts
  "Absent", 32, 35                                     # Absent shape with counts
) |> 
  # Add formatted labels combining shape type and counts
  mutate(label = glue("{shape} <br>(Species = {species_count}; Trait = {trait_count})"))

# Convert labels to a named vector, where shape types are used as names for better access
facet_labels <- setNames(facet_labels$label, facet_labels$shape)

# Define shapes data for plotting ####
# This section defines the x and y coordinates for each shape and their corresponding labels.
shapes <- tribble(
  ~shape, ~x, ~y,                                      # Column names for shape, x-coordinates, and y-coordinates
  "Exponential", 0, 0, "Exponential", 1, .05,         # Exponential: low initial increase, accelerating
  "Exponential", 2, .2, "Exponential", 3, 1,
  "Saturating", 0, 0, "Saturating", 1, 1,             # Saturating: rapid initial increase, plateauing
  "Saturating", 2, 1, "Saturating", 3, 1,
  "Reverse logistic", 0, 0, "Reverse logistic", 1, .6,# Reverse logistic: rapid initial growth, deceleration
  "Reverse logistic", 2, .4, "Reverse logistic", 3, 1,
  "Absent", 0, 0, "Absent", 1, 0,                     # Absent: flat line at zero
  "Absent", 2, 0, "Absent", 3, 0
) |> 
  # Convert shape column to a factor for consistent plotting order
  mutate(shape = factor(shape, levels = c("Saturating", "Reverse logistic", "Exponential", "Absent")))

# Define text annotations for coefficients ####
# This section creates the coefficient labels for each shape at specific points along the x-axis.
text <- tribble(
  ~shape, ~x, ~y, ~label,                             # Column names for shape, x/y-coordinates, and label text
  "Exponential", .5, 1.2, "Coeff<sub>1</sub>",        # Exponential: Coefficient labels at key x positions
  "Exponential", 1.5, 1.2, "Coeff<sub>2</sub>",
  "Exponential", 2.5, 1.2, "Coeff<sub>3</sub>",
  "Saturating", .5, 1.2, "Coeff<sub>1</sub>",         # Saturating: Coefficient labels at key x positions
  "Saturating", 1.5, 1.2, "Coeff<sub>2</sub>",
  "Saturating", 2.5, 1.2, "Coeff<sub>3</sub>",
  "Reverse logistic", .5, 1.2, "Coeff<sub>1</sub>",   # Reverse logistic: Coefficient labels at key x positions
  "Reverse logistic", 1.5, 1.2, "Coeff<sub>2</sub>",
  "Reverse logistic", 2.5, 1.2, "Coeff<sub>3</sub>"
) |> 
  # Convert shape column to a factor to maintain plotting order
  mutate(shape = factor(shape, levels = c("Saturating", "Reverse logistic", "Exponential", "Absent")))

# Define color palette for shapes ####
# Each shape type is assigned a unique color for differentiation in the plot.
shape_colors <- c(
  "#b1c0bc",  # Light grey-green for "Saturating"
  "#5C4E63",  # Deep purple for "Reverse logistic"
  "#333f48",  # Dark grey for "Exponential"
  "#d3d3d3"   # Light grey for "Absent"
)


# Generate predicted probabilities for different predictors ####

# Extract predictions for the shape of relationships between species and trait replacement
# based on the "ecosystem_type" predictor.
pred_shape_ecosystem_type <- extract_shape_predictions(
  species_model = best_model_species,           # Bayesian model for species replacement
  trait_model = best_model_traits,             # Bayesian model for trait replacement
  predictor = "ecosystem_type",                # Predictor variable: "ecosystem_type"
  predictor_levels = c("Freshwater", "Terrestrial")  # Levels of the predictor to include in the predictions
)
# This generates the probabilities for each shape (e.g., Exponential, Saturating) 
# across different ecosystem types (Freshwater and Terrestrial).

# Extract predictions for the shape of relationships between species and trait replacement
# based on the "direction" predictor.
pred_shape_direction <- extract_shape_predictions(
  species_model = best_model_species,           # Bayesian model for species replacement
  trait_model = best_model_traits,             # Bayesian model for trait replacement
  predictor = "direction",                     # Predictor variable: "direction"
  predictor_levels = c("Differentiation", "Homogenisation")  # Levels of the predictor to include in the predictions
)
# This generates the probabilities for each shape across the two types of direction (Differentiation and Homogenisation).

# Extract predictions for the shape of relationships between species and trait replacement
# based on the "main_land_use_type" predictor.
pred_shape_land_use_type <- extract_shape_predictions(
  species_model = best_model_species,           # Bayesian model for species replacement
  trait_model = best_model_traits,             # Bayesian model for trait replacement
  predictor = "main_land_use_type",            # Predictor variable: "main_land_use_type"
  predictor_levels = c("Agriculture", "Forest", "Urban", "Multiple")  # Levels of the predictor to include in the predictions
)
# This generates the probabilities for each shape across the different land-use types 
# (Agriculture, Forest, Urban, and Multiple).


ecosystem_type_colors<- c("#118ab2", "#06d6a0")
direction_colors = c("#ff7f32", "#9f5cc0")  
landu_use_type_colors<-c("#a36627","#848c04","#1c1c0c","#dc7c5c")


# Create individual shape plots ####

# Generate Figure 4a: The shape scheme plot, which visually illustrates the different potential relationship shapes
# (e.g., Exponential, Saturating, Reverse Logistic, Absent) and their corresponding coefficient structures.
fig_4a <- plot_shape_scheme(
  shapes,          # Data representing the various shape curves
  text,            # Coefficient text labels for each shape
  shape_colors,    # Color palette for the shapes
  facet_labels     # Labels for each shape facet, including species and trait counts
) + 
  theme(plot.margin = margin(l = 2, t = 2, b = 2))  # Add margin to improve layout aesthetics

# Generate Figure 4b: The plot showing predicted probabilities for shape relationships based on "Ecosystem type."
fig_4b <- plot_shape(
  data = pred_shape_realm,          # Predictions for shape based on "Ecosystem type"
  predictor_label = "Ecosystem type",  # Label for the predictor
  fill_var = "ecosystem_type",      # Variable used to color the plot
  color_palette = realm_colors      # Color palette for ecosystem types (e.g., Freshwater, Terrestrial)
) + 
  theme(
    strip.text.y = element_blank(), # Remove facet labels on the y-axis
    axis.title.x = element_blank()  # Remove x-axis title to streamline the figure
  )

# Generate Figure 4c: The plot showing predicted probabilities for shape relationships based on "Direction."
fig_4c <- plot_shape(
  data = pred_shape_direction,    # Predictions for shape based on "Direction"
  predictor_label = "Direction",  # Label for the predictor
  fill_var = "direction",         # Variable used to color the plot
  color_palette = direction_colors  # Color palette for directions (e.g., Differentiation, Homogenisation)
) + 
  theme(
    axis.title.y = element_blank(),  # Remove y-axis title for a cleaner layout
    axis.text.y = element_blank(),   # Remove y-axis text to avoid redundancy
    axis.title.x = element_blank()   # Remove x-axis title
  )

# Generate Figure 4d: The plot showing predicted probabilities for shape relationships based on "Main land-use type."
fig_4d <- plot_shape(
  data = pred_shape_disturbance,          # Predictions for shape based on "Main land-use type"
  predictor_label = "Main land-use type", # Label for the predictor
  fill_var = "main_land_use_type",        # Variable used to color the plot
  color_palette = disturbance_colors      # Color palette for land-use types (e.g., Agriculture, Forest, Urban, Multiple)
)

# Combine Figures 4b and 4c into an "upper panel" ####
# The upper panel displays the results for "Ecosystem type" (4b) and "Direction" (4c) side-by-side.
upper_panel <- cowplot::plot_grid(
  fig_4b, fig_4c,              # The two plots to be combined
  ncol = 2,                    # Arrange in two columns
  labels = c("b", "c"),        # Add labels for subplots
  label_size = 5,              # Set size for subplot labels
  align = "h"                  # Horizontally align the plots
)

# Combine all panels into the final Figure 4 ####
# Combine Figure 4a as the first column and the upper panel with Figure 4d as the second column.
final_figure <- cowplot::plot_grid(
  fig_4a,                                      # The shape scheme plot (4a)
  cowplot::plot_grid(                         # Combine the upper panel and Figure 4d
    upper_panel,                              # The upper panel (4b and 4c)
    fig_4d,                                   # Figure 4d (Main land-use type)
    ncol = 1,                                # Arrange in one column
    rel_heights = c(0.5, 0.5),               # Adjust the relative heights of the subplots
    labels = c("", "d"),                     # Add a label for Figure 4d
    label_size = 5                           # Set size for subplot labels
  ),
  ncol = 2,                                  # Arrange in two columns
  rel_widths = c(0.3, 0.7),                  # Adjust the relative widths of the columns
  labels = "a",                              # Add a label for Figure 4a
  label_size = 5                             # Set size for the label
)

# Save the final figure as a PDF file ####
ggsave(
  plot = final_figure,                                      # The combined final figure
  filename = here::here("S7_Model_outputs_figures_and_tables", "main_figures", "Figure_4.pdf"),  # File path and name
  width = 89,                                               # Width of the figure in mm
  height = 80,                                              # Height of the figure in mm
  units = "mm"                                              # Units for dimensions (millimeters)
)




# Extended Data ---------------------------------------------------------------------------------------------------

# Generate the data required for plotting extended figure 3
# This data includes the posterior draws of the model parameters related to the shape of relationships 
# (e.g., Exponential, Saturating, Reverse Logistic, Absent) across biodiversity facets and predictors.
shape_figure_data <- generate_figure_data(
  best_model_species,        # Bayesian model for species replacement
  best_model_traits,        # Bayesian model for trait replacement
  "shape",    # The type of data being generated (in this case, relationship shapes)
  draw_min = -0.5,  # Minimum value for posterior draws
  draw_max = 2.3    # Maximum value for posterior draws
)

# Create Extended Figure 3
# This figure visualizes the posterior draws of the model parameters for different relationship shapes.
extended_figure_3 <- plot_figure(
  shape_figure_data,  # The data generated from the posterior draws
  "shape",            # Type of plot (relationship shapes)
  trunc_upper = 2.0,  # Upper truncation limit for the figure
  trunc_lower = -0.5, # Lower truncation limit for the figure
  tick_factor = 0.5   # Interval for tick marks on the axes
) +
  # Additional theming for figure aesthetics
  theme(
    strip.text.y = element_text(angle = 90, margin = margin(r = 2)), # Rotate facet labels for readability
    plot.margin = margin(3, 5, 3, 3)                                # Add margin around the plot
  )

# Save Extended Figure 3
# The figure is saved as a high-resolution PDF file for inclusion in the extended data of the manuscript.
ggsave(
  filename = here::here("S7_Model_outputs_figures_and_tables", "extended_data", "Extended_Figure_3.pdf"),  # File path and name
  plot = extended_figure_3,        # The figure to save
  width = 80,                      # Width of the figure in mm
  height = 80,                     # Height of the figure in mm
  units = "mm"                     # Units for figure dimensions (millimeters)
)

# pred_shape_interaction <-
#   epred_draws(
#     fsf3,  # Use the model that includes the interaction
#     newdata = insight::get_datagrid(
#       fsf3,  # Ensure you pass the correct model
#       at = c("disturbance", "realm", "trials=1"),  # Include both disturbance and direction
#       preserve_range = FALSE,
#       data = fsf3$data,
#       include_response = TRUE
#     ),
#     re_formula = NA
#   ) |> 
#   mutate(facet = "Trait replacement") |>  # Add any necessary transformations
#   dplyr::select(.epred, facet, realm, .category, disturbance, direction) |>  # Include direction and disturbance in the select
#   dplyr::mutate(
#     disturbance = factor(stringr::str_to_title(disturbance)),  # Convert disturbance to factor
#     direction = factor(stringr::str_to_title(direction))       # Convert direction to factor
#   )
# 
# 
# pred_shape_interaction |> 
#   mutate(predictor = "Main land-use type") |> 
#   ggplot(aes(x=.epred, y=.category, fill=disturbance, colour=disturbance)) +
#   stat_ccdfinterval(slab_alpha=.7, point_interval = "median_hdci", .width=c(.8,.9,.95), fatten_point=.6, interval_size_range=c(.3,.7))+
#   ggh4x::facet_nested(realm~predictor+disturbance)+
#   #scale_fill_manual(values=disturbance_colors)+
#   #scale_colour_manual(values=colorspace::darken(disturbance_colors,.3, space="HLS"))+
#   theme_pred()+
#   theme(panel.border = ggplot2::element_rect(colour="black", fill=NA),
#         legend.position="bottom",
#         panel.spacing.x = unit(.5, "lines"),
#         axis.title.y=element_blank(),
#         #axis.text.y = element_blank(),
#         axis.text.x = ggplot2::element_text(angle=90))+
#   labs(y="Relationship shape", x="Predicted probability", fill="",colour="")




