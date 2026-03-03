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
library(parameters) # For assessing Bayesian models and extracting posterior summaries
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
model_data <- process_model_data("S6_Synthesis_model/data/synthesis_data_complete.xlsx", "shape")


model_data <- model_data |>  filter(framework == "Podani", feature == "abun")



model_data |>
  group_by(facet,direction) |> 
  summarise(across(contains("_q_"),\(x) mean(x, na.rm = TRUE))) 

#   facet   direction       hfp_q_25 hfp_q_50 hfp_q_75 hfp_q_100 het_q_25 het_q_50 het_q_75 het_q_100
#   <fct>   <fct>              <dbl>    <dbl>    <dbl>     <dbl>    <dbl>    <dbl>    <dbl>     <dbl>
# 1 Species Differentiation     49.2     60.7     73.9       100     34.9     47.3     61.6       100
# 2 Species Homogenisation      46.2     55.9     68.0       100     47.7     62.5     75.4       100
# 3 Traits  Differentiation     51.5     62.8     73.2       100     24.5     34.4     52.2       100
# 4 Traits  Homogenisation      41.2     52.8     68.6       100     31.1     41.9     56.2       100
   
model_data |> 
  group_by(facet,direction) |> 
  summarise(across(matches("^(hfp|het)_b[1-4]$"),\(x) mean(x, na.rm = TRUE))) |> 
  group_by(facet,direction) 

#   facet   direction       hfp_b10 hfp_b20 hfp_b30 hfp_b40 het_b0.5 het_b1 het_b1.5 het_b2
#   <fct>   <fct>            <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>
# 1 Species Differentiation   40.6   72.2   87.7   95.6   41.2   68.4   92.6  100  
# 2 Species Homogenisation    52.4   68.7   79.1   87.6   48.4   78.3   97.9  100  
# 3 Traits  Differentiation   42.3   69.4   81.6   89.1   25.0   52.0   85.7   99.8
# 4 Traits  Homogenisation    39.4   69.0   83.5   94.4   35.0   63.9   92.3  100   


#Description of the shapes ----

model_data |>
  select(dataset, direction, facet, hfp_modal_shape,het_modal_shape) |>
  group_by(dataset,direction,facet) |>
  tabyl(hfp_modal_shape,direction,facet) |>
  adorn_percentages("col") %>%
  adorn_pct_formatting(digits = 2) %>%
  adorn_ns() |>
  kableExtra::kable(format = "rst")

#' ===============  ===============  ==============
#' hfp_modal_shape  Differentiation  Homogenisation
#' ===============  ===============  ==============
#' Absent           16.36% (18)      27.45% (14)   
#' Exponential      29.09% (32)      29.41% (15)   
#' Revlog           5.45%  (6)       5.88%  (3)    
#' Saturating       49.09% (54)      37.25% (19)   
#' ===============  ===============  ==============
#' 
#' ===============  ===============  ==============
#' hfp_modal_shape  Differentiation  Homogenisation
#' ===============  ===============  ==============
#' Absent           8.70%  (4)       33.04% (38)   
#' Exponential      28.26% (13)      30.43% (35)   
#' Revlog           2.17%  (1)       4.35%  (5)    
#' Saturating       60.87% (28)      32.17% (37)   
#' ===============  ===============  ==============


model_data |>
  select(dataset, direction, facet, hfp_modal_shape,het_modal_shape) |>
  group_by(dataset,direction,facet) |>
  tabyl(het_modal_shape,direction,facet) |>
  adorn_percentages("col") %>%
  adorn_pct_formatting(digits = 2) %>%
  adorn_ns() |>
  kableExtra::kable(format = "rst")

#' ===============  ===============  ==============
#' het_modal_shape  Differentiation  Homogenisation
#' ===============  ===============  ==============
#' Absent           9.09% (10)       21.57% (11)   
#' Exponential      41.82% (46)      25.49% (13)   
#' Revlog           7.27%  (8)       5.88%  (3)    
#' Saturating       41.82% (46)      47.06% (24)   
#' ===============  ===============  ==============
#' 
#' ===============  ===============  ==============
#' het_modal_shape  Differentiation  Homogenisation
#' ===============  ===============  ==============
#' Absent           15.22%  (7)      33.91% (39)   
#' Exponential      54.35% (25)      36.52% (42)   
#' Revlog           4.35%  (2)       5.22%  (6)    
#' Saturating       26.09% (12)      24.35% (28)   
#' ===============  ===============  ==============


model_data |>
  select(dataset, ecosystem_type, facet, hfp_modal_shape,het_modal_shape) |>
  group_by(dataset,ecosystem_type,facet) |>
  tabyl(hfp_modal_shape,ecosystem_type,facet) |>
  adorn_percentages("col") %>%
  adorn_pct_formatting(digits = 2) %>%
  adorn_ns() |>
  kableExtra::kable(format = "rst")

# Human footprint
# ===============  ===========  ===========
# species          Terrestrial  Freshwater 
# ===============  ===========  ===========
# Absent           19.85% (26)  20.00%  (6)
# Exponential      32.06% (42)  16.67%  (5)
# Revlog           4.58%  (6)   10.00%  (3)
# Saturating       43.51% (57)  53.33% (16)
# ===============  ===========  ===========
# 
# ===============  ===========  ===========
# traits           Terrestrial  Freshwater 
# ===============  ===========  ===========
# Absent           25.95% (34)  26.67%  (8)
# Exponential      33.59% (44)  13.33%  (4)
# Revlog           4.58%  (6)   0.00%  (0) 
# Saturating       35.88% (47)  60.00% (18)
# ===============  ===========  ===========


model_data |>
  select(dataset, ecosystem_type, facet, hfp_modal_shape,het_modal_shape) |>
  group_by(dataset,ecosystem_type,facet) |>
  tabyl(het_modal_shape,ecosystem_type,facet) |>
  adorn_percentages("col") %>%
  adorn_pct_formatting(digits = 2) %>%
  adorn_ns() |>
  kableExtra::kable(format = "rst")

#Habitat heterogeneity
# ===============  ===========  ===========
# species          Terrestrial  Freshwater 
# ===============  ===========  ===========
# Absent           14.50% (19)  6.67%  (2) 
# Exponential      37.40% (49)  33.33% (10)
# Revlog           7.63% (10)   3.33%  (1) 
# Saturating       40.46% (53)  56.67% (17)
# ===============  ===========  ===========
# 
# ===============  ===========  ===========
# traits           Terrestrial  Freshwater 
# ===============  ===========  ===========
# Absent           32.82% (43)  10.00%  (3)
# Exponential      41.98% (55)  40.00% (12)
# Revlog           5.34%  (7)   3.33%  (1) 
# Saturating       19.85% (26)  46.67% (14)
# ===============  ===========  ===========



#'Description of the human pressure across main land-use type ----
 
# "S6_Synthesis_model/data/synthesis_data.xlsx" |>
#   readxl::read_xlsx() |>
#   filter(!dataset %in% c("N67TTP","N78TTP","S47TTP")) |>
#   group_by(disturbance) |>
#   summarise(mean=mean(hfp.min),
#             sd=sd(hfp.min, na.rm=TRUE)) |>
#   kableExtra::kable(format = "rst")

#' ===========  =========  ========
#' disturbance       mean        sd
#' ===========  =========  ========
#' agriculture   5.788789  3.731489
#' forest        4.592429  3.521031
#' multiple      6.015048  5.462427
#' urban        31.598707  3.083456
#' ===========  =========  ========


# Base formula for shape, testing the effects of main predictors
formula_base_hfp <- bf(
  hfp_shape | trials(hfp_trials) ~ species_number + spatial_extent + spatial_distance_min + 
    absolute_latitude_mean + human_footprint_baseline + human_footprint_range + habitat_heterogeneity_baseline +habitat_heterogeneity_range +
    (main_land_use_type + ecosystem_type)*direction + hfp_buffer + het_buffer+ (1 |p| taxa_coarse),
  decomp = "QR"
) 

formula_base_het<- bf(
  het_shape | trials(het_trials) ~ species_number + spatial_extent + spatial_distance_min + 
    absolute_latitude_mean + human_footprint_baseline + human_footprint_range + habitat_heterogeneity_baseline +habitat_heterogeneity_range +
    (main_land_use_type + ecosystem_type)*direction  + hfp_buffer + het_buffer+ (1 |p| taxa_coarse),
  decomp = "QR"
)


colnames(model_data$hfp_shape)<-c("absent","exponential","saturating","revlog")
colnames(model_data$het_shape)<-c("absent","exponential","saturating","revlog")

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
model_base_hfp_species <- brms::brm(
  formula = formula_base_hfp,
  data = model_data |> filter(facet == "Species"),
  family = brms::multinomial(link = "logit", refcat ="absent"),
  prior = priors_shape,
  chains = CHAINS,
  iter = ITER,
  warmup = WARMUP,
  control = base::list(max_treedepth = 15, adapt_delta = .99),
  save_pars = brms::save_pars(all = TRUE),
  seed = BAYES_SEED,
  file = "S7_Model_outputs_figures_and_tables/model/Podani/Abundance/shape/model_base_hfp_species", # Automatically saves model
  file_refit = "on_change"
)

# Base model for species
model_base_het_species <- brms::brm(
  formula = formula_base_het,
  data = model_data |> filter(facet == "Species"),
  family = brms::multinomial(link = "logit", refcat ="absent"),
  prior = priors_shape,
  chains = CHAINS,
  iter = ITER,
  warmup = WARMUP,
  control = base::list(max_treedepth = 15, adapt_delta = .99),
  save_pars = brms::save_pars(all = TRUE),
  seed = BAYES_SEED,
  file = "S7_Model_outputs_figures_and_tables/model/Podani/Abundance/shape/model_base_het_species", # Automatically saves model
  file_refit = "on_change"
)

formula_base_hfp_traits <- bf(
  hfp_shape | trials(hfp_trials) ~ species_number + spatial_extent + spatial_distance_min + 
    absolute_latitude_mean + human_footprint_baseline + human_footprint_range + habitat_heterogeneity_baseline +habitat_heterogeneity_range +
    (main_land_use_type+ecosystem_type)*direction + hfp_buffer + het_buffer+ (1 |p| taxa_coarse),
  decomp = "QR"
) 

formula_base_het_traits<- bf(
  het_shape | trials(het_trials) ~ species_number + spatial_extent + spatial_distance_min + 
    absolute_latitude_mean + human_footprint_baseline + human_footprint_range + habitat_heterogeneity_baseline +habitat_heterogeneity_range +
    (main_land_use_type + ecosystem_type)*direction  + hfp_buffer + het_buffer+ (1 |p| taxa_coarse),
  decomp = "QR"
)

# Base model for species
model_base_hfp_traits <- brms::brm(
  formula = formula_base_hfp_traits,
  data = model_data |> filter(facet == "Traits"),
  family = brms::multinomial(link = "logit", refcat ="absent"),
  prior = priors_shape,
  chains = CHAINS,
  iter = ITER,
  warmup = WARMUP,
  control = base::list(max_treedepth = 15, adapt_delta = .99),
  save_pars = brms::save_pars(all = TRUE),
  seed = BAYES_SEED,
  file = "S7_Model_outputs_figures_and_tables/model/Podani/Abundance/shape/model_base_hfp_traits", # Automatically saves model
  file_refit = "on_change"
)

# Base model for species
model_base_het_traits <- brms::brm(
  formula = formula_base_het_traits,
  data = model_data |> filter(facet == "Traits"),
  family = brms::multinomial(link = "logit", refcat ="absent"),
  prior = priors_shape,
  chains = CHAINS,
  iter = ITER,
  warmup = WARMUP,
  control = base::list(max_treedepth = 15, adapt_delta = .99),
  save_pars = brms::save_pars(all = TRUE),
  seed = BAYES_SEED,
  file = "S7_Model_outputs_figures_and_tables/model/Podani/Abundance/shape/model_base_het_traits", # Automatically saves model
  file_refit = "on_change"
)


#Raw data observations 


model_data |>
  select(dataset, direction, facet, hfp_modal_shape,het_modal_shape) |>
  group_by(dataset,direction,facet) |>
  tabyl(hfp_modal_shape,direction,facet) |>
  adorn_percentages("col") %>%
  adorn_pct_formatting(digits = 2) %>%
  adorn_ns() |>
  kableExtra::kable(format = "rst")

# ===============  ===============  ==============
# hfp_modal_shape  Differentiation  Homogenisation
# ===============  ===============  ==============
# Absent           16.36% (18)      27.45% (14)   
# Exponential      29.09% (32)      29.41% (15)   
# Revlog           5.45%  (6)       5.88%  (3)    
# Saturating       49.09% (54)      37.25% (19)   
# ===============  ===============  ==============
# 
# ===============  ===============  ==============
# hfp_modal_shape  Differentiation  Homogenisation
# ===============  ===============  ==============
# Absent           8.70%  (4)       33.04% (38)   
# Exponential      28.26% (13)      30.43% (35)   
# Revlog           2.17%  (1)       4.35%  (5)    
# Saturating       60.87% (28)      32.17% (37)   
# ===============  ===============  ==============

model_data |>
  select(dataset, direction, facet, hfp_modal_shape,het_modal_shape) |>
  group_by(dataset,direction,facet) |>
  tabyl(het_modal_shape,direction,facet) |>
  adorn_percentages("col") %>%
  adorn_pct_formatting(digits = 2) %>%
  adorn_ns() |>
  kableExtra::kable(format = "rst")

# ===============  ===============  ==============
# het_modal_shape  Differentiation  Homogenisation
# ===============  ===============  ==============
# Absent           9.09% (10)       21.57% (11)   
# Exponential      41.82% (46)      25.49% (13)   
# Revlog           7.27%  (8)       5.88%  (3)    
# Saturating       41.82% (46)      47.06% (24)   
# ===============  ===============  ==============
# 
# ===============  ===============  ==============
# het_modal_shape  Differentiation  Homogenisation
# ===============  ===============  ==============
# Absent           15.22%  (7)      33.91% (39)   
# Exponential      54.35% (25)      36.52% (42)   
# Revlog           4.35%  (2)       5.22%  (6)    
# Saturating       26.09% (12)      24.35% (28)   
# ===============  ===============  ==============


model_data |>
  select(dataset, ecosystem_type, facet, hfp_modal_shape,het_modal_shape) |>
  group_by(dataset,ecosystem_type,facet) |>
  tabyl(hfp_modal_shape,ecosystem_type,facet) |>
  adorn_percentages("col") %>%
  adorn_pct_formatting(digits = 2) %>%
  adorn_ns() |>
  kableExtra::kable(format = "rst")

# Human footprint
# ===============  ===========  ===========
# species          Terrestrial  Freshwater 
# ===============  ===========  ===========
# Absent           19.85% (26)  20.00%  (6)
# Exponential      32.06% (42)  16.67%  (5)
# Revlog           4.58%  (6)   10.00%  (3)
# Saturating       43.51% (57)  53.33% (16)
# ===============  ===========  ===========
# 
# ===============  ===========  ===========
# traits           Terrestrial  Freshwater 
# ===============  ===========  ===========
# Absent           25.95% (34)  26.67%  (8)
# Exponential      33.59% (44)  13.33%  (4)
# Revlog           4.58%  (6)   0.00%  (0) 
# Saturating       35.88% (47)  60.00% (18)
# ===============  ===========  ===========


model_data |>
  select(dataset, ecosystem_type, facet, hfp_modal_shape,het_modal_shape) |>
  group_by(dataset,ecosystem_type,facet) |>
  tabyl(het_modal_shape,ecosystem_type,facet) |>
  adorn_percentages("col") %>%
  adorn_pct_formatting(digits = 2) %>%
  adorn_ns() |>
  kableExtra::kable(format = "rst")

#Habitat heterogeneity
# ===============  ===========  ===========
# species          Terrestrial  Freshwater 
# ===============  ===========  ===========
# Absent           14.50% (19)  6.67%  (2) 
# Exponential      37.40% (49)  33.33% (10)
# Revlog           7.63% (10)   3.33%  (1) 
# Saturating       40.46% (53)  56.67% (17)
# ===============  ===========  ===========
# 
# ===============  ===========  ===========
# traits           Terrestrial  Freshwater 
# ===============  ===========  ===========
# Absent           32.82% (43)  10.00%  (3)
# Exponential      41.98% (55)  40.00% (12)
# Revlog           5.34%  (7)   3.33%  (1) 
# Saturating       19.85% (26)  46.67% (14)
# ===============  ===========  ===========

#LAND USE  -------------------


model_data %>%
  split(.$facet) |> 
  map(~ .x |> 
  select(dataset, main_land_use_type, direction, hfp_modal_shape,het_modal_shape) |>
  group_by(dataset,main_land_use_type, direction) |>
  tabyl(hfp_modal_shape,main_land_use_type) |>
  adorn_percentages("col") %>%
  adorn_pct_formatting(digits = 2) %>%
  adorn_ns() |>
  kableExtra::kable(format = "rst"))

# Human footprint
# $Species
# 
# 
# ===============  ===========  ===========  ===========  ==========
# hfp_modal_shape  Multiple     Agriculture  Forest       Urban     
# ===============  ===========  ===========  ===========  ==========
# Absent           26.92% (14)  15.79%  (6)  19.64% (11)  6.67% (1) 
# Exponential      25.00% (13)  34.21% (13)  21.43% (12)  60.00% (9)
# Revlog           3.85%  (2)   5.26%  (2)   5.36%  (3)   13.33% (2)
# Saturating       44.23% (23)  44.74% (17)  53.57% (30)  20.00% (3)
# ===============  ===========  ===========  ===========  ==========
# 
# $Traits
# 
# 
# ===============  ===========  ===========  ===========  ==========
# hfp_modal_shape  Multiple     Agriculture  Forest       Urban     
# ===============  ===========  ===========  ===========  ==========
# Absent           23.08% (12)  15.79%  (6)  35.71% (20)  26.67% (4)
# Exponential      25.00% (13)  31.58% (12)  30.36% (17)  40.00% (6)
# Revlog           0.00%  (0)   7.89%  (3)   1.79%  (1)   13.33% (2)
# Saturating       51.92% (27)  44.74% (17)  32.14% (18)  20.00% (3)
# ===============  ===========  ===========  ===========  ==========

model_data %>%
  split(.$facet) |> 
  map(~ .x |> 
        select(dataset, main_land_use_type, direction, hfp_modal_shape,het_modal_shape) |>
        group_by(dataset,main_land_use_type, direction) |>
        tabyl(het_modal_shape,main_land_use_type) |>
        adorn_percentages("col") %>%
        adorn_pct_formatting(digits = 2) %>%
        adorn_ns() |>
        kableExtra::kable(format = "rst"))

# Habitat heterogeneity
# $Species
# 
# ===============  ===========  ===========  ===========  ==========
# het_modal_shape  Multiple     Agriculture  Forest       Urban     
# ===============  ===========  ===========  ===========  ==========
# Absent           7.69%  (4)   21.05%  (8)  10.71%  (6)  20.00% (3)
# Exponential      42.31% (22)  39.47% (15)  32.14% (18)  26.67% (4)
# Revlog           9.62%  (5)   5.26%  (2)   5.36%  (3)   6.67% (1) 
# Saturating       40.38% (21)  34.21% (13)  51.79% (29)  46.67% (7)
# ===============  ===========  ===========  ===========  ==========
# 
# $Traits
# 
# ===============  ===========  ===========  ===========  ==========
# het_modal_shape  Multiple     Agriculture  Forest       Urban     
# ===============  ===========  ===========  ===========  ==========
# Absent           19.23% (10)  26.32% (10)  30.36% (17)  60.00% (9)
# Exponential      42.31% (22)  47.37% (18)  41.07% (23)  26.67% (4)
# Revlog           1.92%  (1)   7.89%  (3)   7.14%  (4)   0.00% (0) 
# Saturating       36.54% (19)  18.42%  (7)  21.43% (12)  13.33% (2)
# ===============  ===========  ===========  ===========  ==========


#'########################################################
#Contrast----
#'########################################################
#Human footprint effect on species shape across land-use types for each direction

#Direction x Ecosystem type
grid_hfp_species <- datagrid(
  newdata = model_base_hfp_species$data,
  direction = unique,
  ecosystem_type = unique,
  main_land_use_type = unique,
  hfp_trials = 1
)

grid_hfp_traits <- datagrid(
  newdata = model_base_hfp_traits$data,
  direction = unique,
  ecosystem_type = unique,
  main_land_use_type = unique,
  hfp_trials = 1
)

grid_het_species <- datagrid(
  newdata = model_base_het_species$data,
  direction = unique,
  ecosystem_type = unique,
  main_land_use_type = unique,
  het_trials = 1
)

grid_het_traits <- datagrid(
  newdata = model_base_het_traits$data,
  direction = unique,
  ecosystem_type = unique,
  main_land_use_type = unique,
  het_trials = 1
)


species_hfp_ecosystem <-
  compare_direction(
    model_base_hfp_species,
    grid_hfp_species,
    ecosystem_type
  )


# |.category   |ecosystem_type |direction                        |     .epred|     .lower|     .upper|direction_effect |
# |:-----------|:--------------|:--------------------------------|----------:|----------:|----------:|:----------------|
# |absent      |Terrestrial    |Homogenisation - Differentiation |  0.0145930|  0.0093763|  0.0196569|hom              |
# |absent      |Freshwater     |Homogenisation - Differentiation | -0.0335501| -0.0441689| -0.0232587|diff             |
# |exponential |Terrestrial    |Homogenisation - Differentiation |  0.0017777| -0.0023964|  0.0066019|                 |
# |exponential |Freshwater     |Homogenisation - Differentiation |  0.0217147|  0.0138074|  0.0302083|hom              |
# |saturating  |Terrestrial    |Homogenisation - Differentiation | -0.0069087| -0.0105401| -0.0035134|diff             |
# |saturating  |Freshwater     |Homogenisation - Differentiation | -0.0560246| -0.0650703| -0.0474451|diff             |
# |revlog      |Terrestrial    |Homogenisation - Differentiation | -0.0093357| -0.0114582| -0.0072611|diff             |
# |revlog      |Freshwater     |Homogenisation - Differentiation |  0.0678206|  0.0575451|  0.0774371|hom              |
# 


  compare_direction(
    model_base_hfp_traits,
    grid_hfp_traits,
    ecosystem_type
  )
  
# |.category   |ecosystem_type |direction                        |     .epred|     .lower|     .upper|direction_effect |
# |:-----------|:--------------|:--------------------------------|----------:|----------:|----------:|:----------------|
# |absent      |Terrestrial    |Homogenisation - Differentiation |  0.1297643|  0.1237559|  0.1361611|hom              |
# |absent      |Freshwater     |Homogenisation - Differentiation |  0.1175499|  0.1068902|  0.1293772|hom              |
# |exponential |Terrestrial    |Homogenisation - Differentiation | -0.0540485| -0.0597937| -0.0487663|diff             |
# |exponential |Freshwater     |Homogenisation - Differentiation | -0.0267731| -0.0363384| -0.0173526|diff             |
# |saturating  |Terrestrial    |Homogenisation - Differentiation | -0.0506669| -0.0554339| -0.0462368|diff             |
# |saturating  |Freshwater     |Homogenisation - Differentiation |  0.0084066|  0.0009580|  0.0160050|hom              |
# |revlog      |Terrestrial    |Homogenisation - Differentiation | -0.0249504| -0.0277681| -0.0217400|diff             |
# |revlog      |Freshwater     |Homogenisation - Differentiation | -0.0986107| -0.1071452| -0.0906140|diff             |  


  compare_direction(
    model_base_hfp_species,
    grid_hfp_species,
    main_land_use_type
  )
  
# |.category   |main_land_use_type |direction                        |     .epred|     .lower|     .upper|direction_effect |
# |:-----------|:------------------|:--------------------------------|----------:|----------:|----------:|:----------------|
# |absent      |Multiple           |Homogenisation - Differentiation | -0.0820821| -0.0891421| -0.0741130|diff             |
# |absent      |Agriculture        |Homogenisation - Differentiation |  0.1726205|  0.1624728|  0.1832791|hom              |
# |absent      |Forest             |Homogenisation - Differentiation | -0.0385152| -0.0452161| -0.0321021|diff             |
# |absent      |Urban              |Homogenisation - Differentiation | -0.0897783| -0.1051568| -0.0740470|diff             |
# |exponential |Multiple           |Homogenisation - Differentiation |  0.2047712|  0.1969013|  0.2120669|hom              |
# |exponential |Agriculture        |Homogenisation - Differentiation | -0.1596950| -0.1668413| -0.1522925|diff             |
# |exponential |Forest             |Homogenisation - Differentiation | -0.0335945| -0.0393602| -0.0273748|diff             |
# |exponential |Urban              |Homogenisation - Differentiation |  0.0355400|  0.0242922|  0.0453964|hom              |
# |saturating  |Multiple           |Homogenisation - Differentiation | -0.0627596| -0.0676579| -0.0584054|diff             |
# |saturating  |Agriculture        |Homogenisation - Differentiation | -0.0691578| -0.0782208| -0.0593972|diff             |
# |saturating  |Forest             |Homogenisation - Differentiation |  0.0595440|  0.0508009|  0.0672164|hom              |
# |saturating  |Urban              |Homogenisation - Differentiation | -0.0538020| -0.0630679| -0.0445499|diff             |
# |revlog      |Multiple           |Homogenisation - Differentiation | -0.0596942| -0.0625606| -0.0570081|diff             |
# |revlog      |Agriculture        |Homogenisation - Differentiation |  0.0562150|  0.0464513|  0.0640954|hom              |
# |revlog      |Forest             |Homogenisation - Differentiation |  0.0125936|  0.0079453|  0.0171167|hom              |
# |revlog      |Urban              |Homogenisation - Differentiation |  0.1078730|  0.0932497|  0.1199378|hom              |


  compare_direction(
    model_base_hfp_traits,
    grid_hfp_traits,
    main_land_use_type
  )

# |.category   |main_land_use_type |direction                        |     .epred|     .lower|     .upper|direction_effect |
# |:-----------|:------------------|:--------------------------------|----------:|----------:|----------:|:----------------|
# |absent      |Multiple           |Homogenisation - Differentiation |  0.0919307|  0.0834153|  0.0998665|hom              |
# |absent      |Agriculture        |Homogenisation - Differentiation | -0.0498228| -0.0616165| -0.0403754|diff             |
# |absent      |Forest             |Homogenisation - Differentiation |  0.0353030|  0.0235604|  0.0476393|hom              |
# |absent      |Urban              |Homogenisation - Differentiation |  0.4169803|  0.4006188|  0.4311950|hom              |
# |exponential |Multiple           |Homogenisation - Differentiation | -0.0117352| -0.0187279| -0.0048990|diff             |
# |exponential |Agriculture        |Homogenisation - Differentiation |  0.0644279|  0.0576409|  0.0715204|hom              |
# |exponential |Forest             |Homogenisation - Differentiation |  0.1249567|  0.1173946|  0.1317219|hom              |
# |exponential |Urban              |Homogenisation - Differentiation | -0.3392340| -0.3578201| -0.3209937|diff             |
# |saturating  |Multiple           |Homogenisation - Differentiation | -0.0776354| -0.0847775| -0.0705897|diff             |
# |saturating  |Agriculture        |Homogenisation - Differentiation |  0.0329815|  0.0238065|  0.0404384|hom              |
# |saturating  |Forest             |Homogenisation - Differentiation | -0.1376579| -0.1480884| -0.1277050|diff             |
# |saturating  |Urban              |Homogenisation - Differentiation |  0.0980170|  0.0911167|  0.1052704|hom              |
# |revlog      |Multiple           |Homogenisation - Differentiation | -0.0024509| -0.0035650| -0.0012134|diff             |
# |revlog      |Agriculture        |Homogenisation - Differentiation | -0.0473418| -0.0522499| -0.0420682|diff             |
# |revlog      |Forest             |Homogenisation - Differentiation | -0.0222434| -0.0257988| -0.0190372|diff             |
# |revlog      |Urban              |Homogenisation - Differentiation | -0.1753327| -0.1933878| -0.1579573|diff             |
  

  compare_direction(
    model_base_het_species,
    grid_het_species,
    ecosystem_type
  )

# |.category   |ecosystem_type |direction                        |     .epred|     .lower|     .upper|direction_effect |
# |:-----------|:--------------|:--------------------------------|----------:|----------:|----------:|:----------------|
# |absent      |Terrestrial    |Homogenisation - Differentiation |  0.0652280|  0.0612218|  0.0691382|hom              |
# |absent      |Freshwater     |Homogenisation - Differentiation | -0.1051660| -0.1123138| -0.0971603|diff             |
# |exponential |Terrestrial    |Homogenisation - Differentiation | -0.0686782| -0.0719429| -0.0654108|diff             |
# |exponential |Freshwater     |Homogenisation - Differentiation | -0.0093360| -0.0134573| -0.0048322|diff             |
# |saturating  |Terrestrial    |Homogenisation - Differentiation |  0.0726678|  0.0674601|  0.0780856|hom              |
# |saturating  |Freshwater     |Homogenisation - Differentiation |  0.1696284|  0.1603364|  0.1789778|hom              |
# |revlog      |Terrestrial    |Homogenisation - Differentiation | -0.0691814| -0.0739944| -0.0639858|diff             |
# |revlog      |Freshwater     |Homogenisation - Differentiation | -0.0550260| -0.0610209| -0.0480886|diff             |

  compare_direction(
    model_base_het_traits,
    grid_het_traits,
    ecosystem_type
  )
  
# 
# |.category   |ecosystem_type |direction                        |     .epred|     .lower|     .upper|direction_effect |
# |:-----------|:--------------|:--------------------------------|----------:|----------:|----------:|:----------------|
# |absent      |Terrestrial    |Homogenisation - Differentiation |  0.1663330|  0.1626429|  0.1702721|hom              |
# |absent      |Freshwater     |Homogenisation - Differentiation |  0.1641000|  0.1586184|  0.1702847|hom              |
# |exponential |Terrestrial    |Homogenisation - Differentiation | -0.1479528| -0.1537858| -0.1422883|diff             |
# |exponential |Freshwater     |Homogenisation - Differentiation | -0.1264475| -0.1384821| -0.1155789|diff             |
# |saturating  |Terrestrial    |Homogenisation - Differentiation | -0.0287285| -0.0326276| -0.0244437|diff             |
# |saturating  |Freshwater     |Homogenisation - Differentiation |  0.0793943|  0.0731339|  0.0858702|hom              |
# |revlog      |Terrestrial    |Homogenisation - Differentiation |  0.0103399|  0.0066059|  0.0142431|hom              |
# |revlog      |Freshwater     |Homogenisation - Differentiation | -0.1173041| -0.1288549| -0.1057206|diff             |

  compare_direction(
    model_base_het_species,
    grid_het_species,
    main_land_use_type
  )
# 
# |.category   |main_land_use_type |direction                        |     .epred|     .lower|     .upper|direction_effect |
# |:-----------|:------------------|:--------------------------------|----------:|----------:|----------:|:----------------|
# |absent      |Multiple           |Homogenisation - Differentiation |  0.0571078|  0.0516947|  0.0623877|hom              |
# |absent      |Agriculture        |Homogenisation - Differentiation | -0.0157690| -0.0214336| -0.0098793|diff             |
# |absent      |Forest             |Homogenisation - Differentiation |  0.0235898|  0.0191326|  0.0275241|hom              |
# |absent      |Urban              |Homogenisation - Differentiation | -0.1450329| -0.1576416| -0.1316647|diff             |
# |exponential |Multiple           |Homogenisation - Differentiation | -0.0568103| -0.0620713| -0.0519335|diff             |
# |exponential |Agriculture        |Homogenisation - Differentiation |  0.0682738|  0.0605225|  0.0746600|hom              |
# |exponential |Forest             |Homogenisation - Differentiation | -0.1071743| -0.1112461| -0.1033917|diff             |
# |exponential |Urban              |Homogenisation - Differentiation | -0.0602871| -0.0639670| -0.0558663|diff             |
# |saturating  |Multiple           |Homogenisation - Differentiation | -0.0502681| -0.0573422| -0.0428202|diff             |
# |saturating  |Agriculture        |Homogenisation - Differentiation |  0.0586924|  0.0481402|  0.0689133|hom              |
# |saturating  |Forest             |Homogenisation - Differentiation |  0.1245412|  0.1164569|  0.1310962|hom              |
# |saturating  |Urban              |Homogenisation - Differentiation |  0.3516371|  0.3387868|  0.3651318|hom              |
# |revlog      |Multiple           |Homogenisation - Differentiation |  0.0499823|  0.0441378|  0.0555636|hom              |
# |revlog      |Agriculture        |Homogenisation - Differentiation | -0.1112534| -0.1176236| -0.1049267|diff             |
# |revlog      |Forest             |Homogenisation - Differentiation | -0.0409364| -0.0460819| -0.0349825|diff             |
# |revlog      |Urban              |Homogenisation - Differentiation | -0.1460342| -0.1576436| -0.1336546|diff             |  
#   
  

  compare_direction(
    model_base_het_traits,
    grid_het_traits,
    main_land_use_type
  )
# 
# |.category   |main_land_use_type |direction                        |     .epred|     .lower|     .upper|direction_effect |
# |:-----------|:------------------|:--------------------------------|----------:|----------:|----------:|:----------------|
# |absent      |Multiple           |Homogenisation - Differentiation | -0.0140324| -0.0189664| -0.0084836|diff             |
# |absent      |Agriculture        |Homogenisation - Differentiation |  0.0671294|  0.0633124|  0.0711681|hom              |
# |absent      |Forest             |Homogenisation - Differentiation |  0.1222769|  0.1181846|  0.1261528|hom              |
# |absent      |Urban              |Homogenisation - Differentiation |  0.4857321|  0.4734480|  0.4969369|hom              |
# |exponential |Multiple           |Homogenisation - Differentiation | -0.0679211| -0.0760746| -0.0604570|diff             |
# |exponential |Agriculture        |Homogenisation - Differentiation | -0.0429126| -0.0551163| -0.0332022|diff             |
# |exponential |Forest             |Homogenisation - Differentiation | -0.0833846| -0.0937016| -0.0708512|diff             |
# |exponential |Urban              |Homogenisation - Differentiation | -0.3545451| -0.3704187| -0.3403215|diff             |
# |saturating  |Multiple           |Homogenisation - Differentiation |  0.0649516|  0.0596988|  0.0704563|hom              |
# |saturating  |Agriculture        |Homogenisation - Differentiation |  0.0025416| -0.0032414|  0.0080451|                 |
# |saturating  |Forest             |Homogenisation - Differentiation |  0.0817073|  0.0750571|  0.0885240|hom              |
# |saturating  |Urban              |Homogenisation - Differentiation | -0.0476076| -0.0575962| -0.0368007|diff             |
# |revlog      |Multiple           |Homogenisation - Differentiation |  0.0171180|  0.0104011|  0.0226201|hom              |
# |revlog      |Agriculture        |Homogenisation - Differentiation | -0.0266595| -0.0375805| -0.0156651|diff             |
# |revlog      |Forest             |Homogenisation - Differentiation | -0.1203455| -0.1312816| -0.1080534|diff             |
# |revlog      |Urban              |Homogenisation - Differentiation | -0.0833028| -0.0937687| -0.0722282|diff             |



# Create individual shape plots ####

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
  "Absent", .5, .5, "Coeff<sub>1</sub>",         # Saturating: Coefficient labels at key x positions
  "Absent", 1.5, .5, "Coeff<sub>2</sub>",
  "Absent", 2.5, .5, "Coeff<sub>3</sub>"

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


# Generate Figure 4a: The shape scheme plot, which visually illustrates the different potential relationship shapes
# (e.g., Exponential, Saturating, Reverse Logistic, Absent) and their corresponding coefficient structures.
fig_4a <-
  ggplot(shapes) +
  aes(x = x, y = y, group = shape) +
  facet_wrap(~shape, ncol = 1) +
  geom_rect(aes(xmin = 0, xmax = 1, ymin = 0, ymax = 1.5), fill = shape_colors[1], colour = "white", alpha = .3) +
  geom_rect(aes(xmin = 1, xmax = 2, ymin = 0, ymax = 1.5), fill = shape_colors[2], colour = "white", alpha = .3) +
  geom_rect(aes(xmin = 2, xmax = 3, ymin = 0, ymax = 1.5), fill = shape_colors[3], colour = "white", alpha = .3) +
  ggtext::geom_richtext(
    data = text, 
    aes(x = x, y = y + .05, label = label),
    size = convert_size(10), colour = NA, fill = NA, text.colour = "black"
  ) +
  ggalt::geom_xspline(spline_shape = 1, size =0.5) +
    ggplot2::theme_void(base_family = "sans", base_size = 14)+ 
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      strip.background = ggplot2::element_rect(fill = NA, color = NA),
      strip.text = element_text(face="bold", margin=margin(b=3)),
      axis.title.x = ggtext::element_markdown(family = "sans", face = "bold", hjust = 0.5, size = 12, margin=margin(t=2)),
      axis.text.y = ggtext::element_markdown(family = "sans", size = 12, margin=margin(t=2),hjust=1),
      axis.title.y = ggtext::element_markdown(family = "sans", face = "bold", hjust = 0.5, size = 12, margin=margin(r=5), angle=90),
      axis.text.x = ggtext::element_markdown(family = "sans", size = 12, margin=margin(b=5)),
      panel.spacing.x = unit(0.2, "lines"),
      panel.spacing.y = unit(0.2, "lines"),
      axis.line.y = element_line(color = "grey20", linewidth=0.1),
      axis.ticks.y = element_line(color = "grey20", linewidth=0.1),
      axis.line.x = element_line(color = "grey20", linewidth=0.1),
      axis.ticks.x = element_line(color = "grey20", linewidth=0.1),
      axis.ticks.length = unit(0.3, "lines"),
      plot.margin = margin(3, 0, 3, 0, unit="mm"),
      legend.position="none",
      plot.tag = element_text(face="bold",size=14))+
  labs(x = "Human pressure gradient", y = "Species/Trait replacement", tag="a") +
  scale_y_continuous(breaks = seq(0, 1, 1), limits = c(0, 1.5)) +
  guides(x = ggplot2::guide_axis(cap = TRUE),
         y = ggplot2::guide_axis(cap = TRUE)) +
  scale_x_continuous(breaks = seq(0, 3, length.out = 5), labels = c("0", "", "", "", "1"), limits = c(0, 3)) +
  theme(plot.margin = margin(l = 2, t = 2, b = 2))  # Add margin to improve layout aesthetics


# Prepare the data
heatmap_raw_shapes <- model_data |> 
  select(dataset, facet, hfp_modal_shape, het_modal_shape) |> 
  pivot_longer(
    cols = c(hfp_modal_shape, het_modal_shape),
    names_to = "predictor", 
    values_to = "shape"
  ) |> 
  mutate(
    predictor = ifelse(predictor == "hfp_modal_shape", "Human Footprint", "Habitat Heterogeneity"),
    predictor = factor(predictor, levels = c("Human Footprint", "Habitat Heterogeneity"))
  ) |> 
  pivot_wider(
    names_from = facet, 
    values_from = shape
  ) |> 
  count(predictor, Traits, Species) |> 
  complete(predictor, Traits, Species, fill = list(n = 0)) |> 
  mutate(
    Traits = factor(Traits, levels = c("Absent", "Exponential", "Saturating","Revlog")),
    Species = factor(Species, levels = c("Absent", "Exponential", "Saturating","Revlog"))
  ) |> 
  mutate(across(Species:Traits, ~ gsub("Revlog","Reverse<br>logistic",.x)))
  

# Create square heatmap
fig_4b <- ggplot(heatmap_raw_shapes, aes(x = Traits, y = Species, fill = n)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = n), color = "black", size = 3.5) +
  facet_wrap(~ predictor, ncol=1) +
  scale_fill_gradient(
    low = "white", 
    high = "#1E88E5", 
    name = "Number of datasets"
  ) +
  labs(
    x = "Trait replacement",
    y = "Species replacement",
    tag="b"
  ) +
  ggplot2::theme_void(base_family = "sans", base_size = 12) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    plot.background = ggplot2::element_rect(fill = "white", color = NA),
    strip.background = ggplot2::element_rect(fill = NA, color = NA),
    strip.text = element_text(face="bold", margin=margin(b=3)),
    axis.title.x = ggtext::element_markdown(family = "sans", face = "bold", hjust = 0.5, size = 12, margin=margin(t=2)),
    axis.text.y = ggtext::element_markdown(family = "sans", size = 10, margin=margin(t=2),hjust=1),
    axis.title.y = ggtext::element_markdown(family = "sans", face = "bold", hjust = 0.5, size = 12, margin=margin(r=5), angle=90),
    axis.text.x = ggtext::element_markdown(family = "sans", size = 8, margin=margin(b=5), angle=30, hjust=1, vjust=.8),
    panel.spacing.x = unit(0.2, "lines"),
    panel.spacing.y = unit(0.2, "lines"),
    axis.line.x = element_line(color = "grey20", linewidth=0.1),
    axis.ticks.x = element_line(color = "grey20", linewidth=0.1),
    axis.ticks.length.x = unit(0.3, "lines"),
    plot.margin = margin(3, 0, 3, 0, unit="mm"),
    legend.position="none",
    plot.tag = element_text(face="bold",size=14)
  )


ecosystem_type_colors<- c("Freshwater" = "#118ab2","Terrestrial"= "#06d6a0")
direction_colors = c("Differentiation" = "#ff7f32","Homogenisation" = "#9f5cc0")  
land_use_colors<-c("Agriculture"="#a36627","Forest"="#848c04","Urban"="#1c1c0c","Multiple"="#dc7c5c")

#####
# HUMAN FOOTPRINT
####


## Ecosystem type × direction  ----

hfp_species_eco <-
  model_base_hfp_species |>
  epred_draws(
    newdata = grid_hfp_species,
    re_formula = NA,
  ) |>
  group_by(.draw, ecosystem_type, direction, .category) |>
  summarise(.epred = mean(.epred), .groups = "drop") |>
  mutate(
    facet = "Species",
    class = "Ecosystem type",
    Predictor = ecosystem_type,
    response = "Human Footprint",
  )

hfp_traits_eco <-
  model_base_hfp_traits |>
  epred_draws(
    newdata = grid_hfp_traits,
    re_formula = NA
  ) |>
  group_by(.draw, ecosystem_type, direction, .category) |>
  summarise(.epred = mean(.epred), .groups = "drop") |>
  mutate(
    facet = "Traits",
    class = "Ecosystem type",
    Predictor = ecosystem_type,
    response = "Human Footprint"
  )


## Land-use × direction ----

hfp_species_land <-
  model_base_hfp_species |>
  epred_draws(
    newdata = grid_hfp_species,
    re_formula = NA
  ) |>
  group_by(.draw, main_land_use_type, direction, .category) |>
  summarise(.epred = mean(.epred), .groups = "drop") |>
  mutate(
    facet = "Species",
    class = "Main land use type",
    Predictor = main_land_use_type,
    response = "Human Footprint"
  )

hfp_traits_land <-
  model_base_hfp_traits |>
  epred_draws(
    newdata = grid_hfp_traits,
    re_formula = NA
  ) |>
  group_by(.draw, main_land_use_type, direction, .category) |>
  summarise(.epred = mean(.epred), .groups = "drop") |>
  mutate(
    facet = "Traits",
    class = "Main land use type",
    Predictor = main_land_use_type,
    response = "Human Footprint"
  )

####
#HABITAT HETEROGENEITY 
####

## Ecosystem type × direction ----

het_species_eco <-
  model_base_het_species |>
  epred_draws(
    newdata = grid_het_species,
    re_formula = NA
  ) |>
  group_by(.draw, ecosystem_type, direction, .category) |>
  summarise(.epred = mean(.epred), .groups = "drop") |>
  mutate(
    facet = "Species",
    class = "Ecosystem type",
    Predictor = ecosystem_type,
    response = "Habitat Heterogeneity"
  )

het_traits_eco <-
  model_base_het_traits |>
  epred_draws(
    newdata = grid_het_traits,
    re_formula = NA
  ) |>
  group_by(.draw, ecosystem_type, direction, .category) |>
  summarise(.epred = mean(.epred), .groups = "drop") |>
  mutate(
    facet = "Traits",
    class = "Ecosystem type",
    Predictor = ecosystem_type,
    response = "Habitat Heterogeneity"
  )


## Land-use × direction ----

het_species_land <-
  model_base_het_species |>
  epred_draws(
    newdata = grid_het_species,
    re_formula = NA
  ) |>
  group_by(.draw, main_land_use_type, direction, .category) |>
  summarise(.epred = mean(.epred), .groups = "drop") |>
  mutate(
    facet = "Species",
    class = "Main land use type",
    Predictor = main_land_use_type,
    response = "Habitat Heterogeneity"
  )

het_traits_land <-
  model_base_het_traits |>
  epred_draws(
    newdata = grid_het_traits,
    re_formula = NA
  ) |>
  group_by(.draw, main_land_use_type, direction, .category) |>
  summarise(.epred = mean(.epred), .groups = "drop") |>
  mutate(
    facet = "Traits",
    class = "Main land use type",
    Predictor = main_land_use_type,
    response = "Habitat Heterogeneity"
  )

contrasts_shape <-
  bind_rows(
    hfp_species_eco,
    hfp_traits_eco,
    hfp_species_land,
    hfp_traits_land,
    het_species_eco,
    het_traits_eco,
    het_species_land,
    het_traits_land
  ) |>
  rename(
    draw = .draw,
    shape = .category,
    probability = .epred
  )


fig_4c <- contrasts_shape |>  
  select(draw, response, shape,facet,class,Predictor, direction, probability) |> 
  nest(data=c(draw, shape, direction,probability)) |> 
  mutate(results = purrr::map(data, regime_test,.80)) |>
  tidyr::unnest(results) |> 
  mutate(
    prob_lab = sprintf(" %.2f", regime_probability)
  ) |> 
 ggplot(aes(x = Predictor, y = shape)) +
  # centered label
  geom_label(
    aes(label = glue::glue("{prob_lab}"), colour = direction),
    hjust=.5,
    fill = "white",
  ) +
  ggh4x::facet_nested(
    response   ~ facet +class ,
    scales = "free_x",
    space  = "free_x",
    strip = strip_nested(
      text_y = elem_list_text(face = "bold", size = c(10, 9)),
      by_layer_y = TRUE,
      text_x = elem_list_text(face = "bold", size = c(10, 9)),
      by_layer_x = TRUE
    )
  ) +
  
  scale_colour_manual(
    values = direction_colors
  ) +
  
  #  scale_x_continuous(limits = c(0.84, 1.12), expand = c(0,0)) +
  theme_bw(base_family = "sans", base_size = 12) +
  theme(
    strip.text.x =  element_text(face="bold",margin=margin(t=5)),
    strip.text.y = element_text(size=10,face="bold",angle=270, margin=margin(l=5)),
    axis.title.x       = element_blank(),
    axis.text.y        = element_text(size = 10, hjust = 1),
    axis.text.x        = ggtext::element_markdown(size = 8, margin=margin(t=3),angle=30,vjust=1),
    axis.title.y       = element_blank(),
    axis.line.x        = element_line(color = "grey20", linewidth = 0.1),
    axis.ticks.x       = element_line(color = "grey20", linewidth = 0.1),
    axis.ticks.length  = unit(0.2, "lines"),
    plot.tag           = element_text(face = "bold", size = 14),
    panel.spacing.x = unit(2,"lines"),
    plot.margin = margin(0,0,0,0),
    legend.position="bottom",
    legend.title = element_text(face="bold"),
    strip.background = element_blank(),
    panel.border = element_blank(),
    panel.background = element_rect_round(
      fill = NA, 
      color = "black", 
      radius = unit(0.5, "cm") # Controla o arredondamento
    )
  ) +
  guides(
    x = guide_axis(cap = TRUE)
  ) + labs(
    colour="Dominant direction",
    tag="c")


Figure_4<- ((fig_4a/fig_4b) | (fig_4c)) + plot_layout(widths=c(.2,.8))

ggsave(filename = here::here("S7_Model_outputs_figures_and_tables", "main_figures", "Figure_4.pdf"),
       plot = Figure_4,
       device = cairo_pdf,
       width=11,height=9,units="in")


#supplemantary material


fig_new_supplemantary <- contrasts_shape |> 
  mutate(shape = str_to_sentence(shape)) |> 
  mutate(shape = gsub("Revlog","Reverse logistic",shape)) |> 
  mutate(
    Predictor = case_when(
      Predictor == "Freshwater"  ~ "Freshwater (n=30)",
      Predictor == "Terrestrial" ~ "Terrestrial (n=131)",
      Predictor == "Agriculture" ~ "Agriculture (n=38)",
      Predictor == "Forest"      ~ "Forest (n=56)",
      Predictor == "Urban"       ~ "Urban (n=15)",
      Predictor == "Multiple"    ~ "Multiple (n=52)",
      TRUE                        ~ Predictor
    ),
    Predictor = factor(
      Predictor,
      levels = c(
        "Freshwater (n=30)",
        "Terrestrial (n=131)",
        rev(c(
          "Agriculture (n=38)",
          "Forest (n=56)",
          "Urban (n=15)",
          "Multiple (n=52)"
        ))
      )
    )
  ) |> 
  mutate(response = factor(response,levels=c("Human Footprint","Habitat Heterogeneity"))) |> 
  ggplot(aes(y=shape,x=probability,colour=direction, fill=direction))+
  ggstats::geom_stripped_rows(
    aes(y = shape),
    colour=NA,
    odd  = "white",
    even = "grey90",
    xfrom = -Inf, xto = Inf
  ) +
  stat_ccdfinterval(position="dodge", .width=c(.8,.89,.95), slab_alpha=.5) +
  ggh4x::facet_nested(class+Predictor~facet+response, scales="free_y",space="free_y",
                      # use strip_themed() instead of strip_nested()
                      strip = strip_nested(
                        # supply one size for the "class" layer, one for the "Predictor" layer
                        text_y      = elem_list_text(size = c(10, 8)),
                        by_layer_y = TRUE,
                        text_x      = elem_list_text(size = c(10, 8)),
                        by_layer_x = TRUE,
                      ))+
  scale_colour_manual(values=direction_colors)+
  scale_fill_manual(values=direction_colors)+
  theme_void(base_family = "sans", base_size = 12) +
  theme(
    panel.grid.minor   = element_blank(),
    plot.background    = element_blank(),
    panel.background =  element_blank(),
    strip.text.x =  element_text(face="bold",margin=margin(t=5)),
    strip.text.y = element_text(size=10,face="bold",angle=270, margin=margin(l=5)),
    axis.title.x       = ggtext::element_markdown(face = "bold", size = 12, margin=margin(t=5)),
    axis.text.y        = element_text(size = 10, hjust = 1),
    axis.text.x        = ggtext::element_markdown(size = 8, margin=margin(t=3)),
    axis.title.y       = element_blank(),
    axis.line.x        = element_line(color = "grey20", linewidth = 0.1),
    axis.ticks.x       = element_line(color = "grey20", linewidth = 0.1),
    axis.ticks.length  = unit(0.2, "lines"),
    plot.tag           = element_text(face = "bold", size = 14),
    panel.spacing.x = unit(2,"lines"),
    plot.margin = margin(0,0,0,0),
    legend.position="bottom",
    legend.title = element_text(face="bold")
  ) +
  guides(
    x = guide_axis(cap = TRUE)
  ) + labs(x = "Posterior predictions (Probability of observing a shape)",
           colour="Relationship direction",
           fill="Relationship direction",
           tag="c")+
  scale_x_continuous(limits=c(0,NA),
                     labels = scales::label_number(accuracy = 0.01)
  )

posterior_plot_data_HFP <- generate_figure_data(model_base_hfp_species,model_base_hfp_traits, type="shape",draw_min=-5,draw_max=5) |> 
mutate(response = "Human Footprint")  

posterior_plot_data_HH <- generate_figure_data(model_base_het_species,model_base_het_traits, type="shape",draw_min=-5,draw_max=5) |> 
  mutate(response = "Habitat Heterogeneity") 

posterior_plot_data<- bind_rows(posterior_plot_data_HFP,posterior_plot_data_HH) |> 
  mutate(response = factor(response,levels=c("Human Footprint","Habitat Heterogeneity"))) |> 
  mutate(shape = str_to_sentence(shape))

# 1. Tag each row as either a draw or a text annotation
draws_df <- posterior_plot_data %>%
  mutate(layer = " ")

text_df <- posterior_plot_data %>%
  distinct(shape,parameter, facet, response, text, draw_min, draw_max) %>%
  # position text slightly beyond the max draw
  mutate(posterior = draw_max * 1.05,
         layer     = "  ")

plot_df <- bind_rows(draws_df, text_df) %>%
  mutate(layer = factor(layer, levels = c(" ", "  "))) 

 Extended_Figure_5 <- 
  ggplot(plot_df) +
  facet_nested(
    facet+response~shape+layer,
   scales      = "free_x",
   space       = "free_x",
  )+
  stat_pointinterval(
    data    = ~ .x |>  filter(layer == " "),
    aes(x = posterior, y = parameter, colour = significant),
    .width = c(.8, .89, .95)
  ) +
  # text labels
  geom_text(
    data    =~ .x |> filter(layer == "  "),
    aes(x = posterior, y = parameter, label = text),
    hjust   = 0,
    size    = convert_size(5)
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  # nested facets: rows = Posterior/Text, cols = your facet (Species vs Trait)
  # per‐strip x scales
  facetted_pos_scales(
    x = list(
      # for the Posterior strip
      layer == " " ~ scale_x_continuous(
        name   = "Posterior draws",
        expand = expansion(mult = 0.02)
      ),
      # for the Text strip
      layer == "  " ~ scale_x_continuous(
        limits = c(5.2,7),
        breaks=NULL,
        name   = NULL,
        expand = expansion(mult = c(0, 0.02))
      )
    )
  ) +
  scale_colour_manual(values = c("gray60", "black")) +
  theme_void(base_family = "sans", base_size = 12) +
  theme(
    panel.grid.minor   = element_blank(),
    plot.background    = element_blank(),
    panel.background =  element_blank(),
    strip.text.y = element_text(face="bold",angle=270,size=12),
    strip.text.x = element_text(face="bold",size=12),,
    axis.title.x       = ggtext::element_markdown(face = "bold", size = 12),
    axis.text.y        = ggtext::element_markdown(size = 10, hjust = 1),
    axis.text.x        = ggtext::element_markdown(size = 12),
    axis.title.y       = element_blank(),
    axis.line.x        = element_line(color = "grey20", linewidth = 0.1),
    axis.ticks.x       = element_line(color = "grey20", linewidth = 0.1),
    axis.ticks.length  = unit(0.1, "lines"),
    legend.position    = "none",
    plot.tag           = element_text(face = "bold", size = 14),
    panel.spacing.x = unit(0,"lines"),
    plot.margin = margin(2,2,2,2, unit="mm"),
    plot.title= element_text(face="bold",size=12)
  ) +
  guides(
    x = guide_axis(cap = TRUE)
  ) 

 ggsave(filename = here::here("S7_Model_outputs_figures_and_tables", "extended_data", "Extended_Figure_5.pdf"),
        plot = Extended_Figure_5,
        device = cairo_pdf,
        width=14,height=14,units="in")

 
 
 avg_predictions(best_model_species, variables="ecosystem_type", conf_level = .8)
 
 # ecosystem_type Estimate 10.0 % 90.0 %
 #     Terrestrial   -1.456 -1.710 -1.171
 #     Freshwater    -0.315 -0.889  0.304
 # 
 # Type:  response 
 
 
 avg_predictions(best_model_traits, variables="ecosystem_type", conf_level = .8)
 
 
 # ecosystem_type Estimate 10.0 % 90.0 %
 #   Terrestrial    0.622  0.457  0.787
 # Freshwater     0.245 -0.130  0.616
 # 
 # Type:  response 
 
 
 avg_predictions(model_base_hfp_species, variables = c("ecosystem_type","direction"),
                 newdata = datagrid(newdata=model_base_hfp_species$data,
                                    direction = unique,
                                    ecosystem_type = unique,
                                    main_land_use_type=unique,
                                    hfp_trials=1), re_formula=NA,conf_level=0.8) 
 
 # Group       direction ecosystem_type Estimate 10.0 % 90.0 %
 # absent      Differentiation    Terrestrial   0.3988 0.3401  0.459
 # absent      Homogenisation     Terrestrial   0.4173 0.3547  0.477
 # absent      Differentiation    Freshwater    0.3841 0.3270  0.440
 # absent      Homogenisation     Freshwater    0.3437 0.2863  0.398
 # exponential Differentiation    Terrestrial   0.2906 0.2178  0.356
 # exponential Homogenisation     Terrestrial   0.2905 0.2207  0.353
 # exponential Differentiation    Freshwater    0.1759 0.1228  0.222
 # exponential Homogenisation     Freshwater    0.1931 0.1438  0.245
 # revlog      Differentiation    Terrestrial   0.0923 0.0594  0.119
 # revlog      Homogenisation     Terrestrial   0.0805 0.0531  0.104
 # revlog      Differentiation    Freshwater    0.1531 0.1062  0.197
 # revlog      Homogenisation     Freshwater    0.2246 0.1625  0.277
 # saturating  Differentiation    Terrestrial   0.2098 0.1652  0.256
 # saturating  Homogenisation     Terrestrial   0.2039 0.1596  0.243
 # saturating  Differentiation    Freshwater    0.2767 0.2272  0.331
 # saturating  Homogenisation     Freshwater    0.2299 0.1866  0.276
 
 avg_predictions(model_base_hfp_traits, variables = c("ecosystem_type","direction"),
                 newdata = datagrid(newdata=model_base_hfp_traits$data,
                                    direction = unique,
                                    ecosystem_type = unique,
                                    main_land_use_type=unique,
                                    hfp_trials=1), re_formula=NA,conf_level=0.8) 
 
 #        Group       direction ecosystem_type Estimate 10.0 % 90.0 %
 # absent      Differentiation    Terrestrial   0.1781 0.1400 0.2147
 # absent      Homogenisation     Terrestrial   0.2704 0.2214 0.3243
 # absent      Differentiation    Freshwater    0.2157 0.1700 0.2630
 # absent      Homogenisation     Freshwater    0.3285 0.2746 0.3839
 # exponential Differentiation    Terrestrial   0.3241 0.2608 0.3969
 # exponential Homogenisation     Terrestrial   0.3140 0.2489 0.3849
 # exponential Differentiation    Freshwater    0.2055 0.1490 0.2580
 # exponential Homogenisation     Freshwater    0.2350 0.1791 0.2925
 # revlog      Differentiation    Terrestrial   0.1649 0.1071 0.2165
 # revlog      Homogenisation     Terrestrial   0.1288 0.0834 0.1784
 # revlog      Differentiation    Freshwater    0.2491 0.1785 0.3147
 # revlog      Homogenisation     Freshwater    0.0585 0.0352 0.0823
 # saturating  Differentiation    Terrestrial   0.3257 0.2736 0.3799
 # saturating  Homogenisation     Terrestrial   0.2768 0.2228 0.3326
 # saturating  Differentiation    Freshwater    0.3217 0.2681 0.3695
 # saturating  Homogenisation     Freshwater    0.3701 0.3031 0.4329
 
 avg_comparisons(model_base_hfp_species,  variables= c("direction"), by = c("ecosystem_type"),
                 newdata = datagrid(newdata=model_base_hfp_species$data,
                                    direction = unique,
                                    ecosystem_type = unique,
                                    main_land_use_type=unique,
                                    hfp_trials=1), re_formula=NA,conf_level=0.8) 
 
 #Species HFP
 #       Group      Term ecosystem_type  Estimate   10.0 %    90.0 %
 # absent      direction    Freshwater  -0.040387 -0.05482 -0.027027 * diff
 # absent      direction    Terrestrial  0.018210  0.01214  0.024648 * hom
 # exponential direction    Freshwater   0.016618  0.00845  0.025119 * hom
 # exponential direction    Terrestrial  0.000193 -0.00670  0.006323
 # saturating  direction    Freshwater  -0.046584 -0.05932 -0.034088 * diff
 # saturating  direction    Terrestrial -0.005567 -0.01128 -0.000458 * diff
 # revlog      direction    Freshwater   0.070476  0.05618  0.085277 * hom
 # revlog      direction    Terrestrial -0.011669 -0.01820 -0.007385 * diff
 
 avg_comparisons(model_base_hfp_traits, variables= c("direction"), by = c("ecosystem_type"),
                 newdata = datagrid(newdata=model_base_hfp_traits$data,
                                    direction = unique,
                                    ecosystem_type = unique,
                                    main_land_use_type=unique,
                                    hfp_trials=1), re_formula=NA,conf_level=0.8) 
 #Trait HFP
 #       Group      Term ecosystem_type Estimate   10.0 %   90.0 %
 # absent      direction    Freshwater   0.11242  0.09536  0.13146 * hom
 # absent      direction    Terrestrial  0.09294  0.07680  0.10929 * hom
 # exponential direction    Freshwater   0.02864  0.00386  0.05667 * hom
 # exponential direction    Terrestrial -0.00887 -0.02492  0.00657
 # saturating  direction    Freshwater   0.04722  0.02613  0.07507 * hom
 # saturating  direction    Terrestrial -0.04819 -0.05368 -0.04176 *diff
 # revlog      direction    Freshwater  -0.19043 -0.23985 -0.14804 * diff
 # revlog      direction    Terrestrial -0.03554 -0.04651 -0.02602 * diff
 
 avg_predictions(model_base_het_species, variables = c("ecosystem_type","direction"),
                 newdata = datagrid(newdata=model_base_het_species$data,
                                    direction = unique,
                                    ecosystem_type = unique,
                                    main_land_use_type=unique,
                                    het_trials=1), re_formula=NA,conf_level=0.8) 

 # Group       direction ecosystem_type Estimate 10.0 % 90.0 %
 # absent      Differentiation    Terrestrial   0.1386 0.1243  0.153
 # absent      Homogenisation     Terrestrial   0.2181 0.1993  0.236
 # absent      Differentiation    Freshwater    0.2444 0.2224  0.264
 # absent      Homogenisation     Freshwater    0.1381 0.1246  0.154
 # exponential Differentiation    Terrestrial   0.3814 0.2938  0.467
 # exponential Homogenisation     Terrestrial   0.2702 0.2010  0.333
 # exponential Differentiation    Freshwater    0.2232 0.1572  0.285
 # exponential Homogenisation     Freshwater    0.2092 0.1510  0.272
 # revlog      Differentiation    Terrestrial   0.1783 0.1311  0.223
 # revlog      Homogenisation     Terrestrial   0.1158 0.0842  0.145
 # revlog      Differentiation    Freshwater    0.1478 0.1062  0.185
 # revlog      Homogenisation     Freshwater    0.0914 0.0653  0.117
 # saturating  Differentiation    Terrestrial   0.2948 0.2209  0.365
 # saturating  Homogenisation     Terrestrial   0.3901 0.3268  0.456
 # saturating  Differentiation    Freshwater    0.3772 0.2993  0.443
 # saturating  Homogenisation     Freshwater    0.5566 0.4837  0.625
 
 
 avg_predictions(model_base_het_traits, variables = c("ecosystem_type","direction"),
                 newdata = datagrid(newdata=model_base_het_traits$data,
                                    direction = unique,
                                    ecosystem_type = unique,
                                    main_land_use_type=unique,
                                    het_trials=1), re_formula=NA,conf_level=0.8) 
 
 
 # Group       direction ecosystem_type Estimate 10.0 % 90.0 %
 # absent      Differentiation    Terrestrial   0.1586 0.1232  0.194
 # absent      Homogenisation     Terrestrial   0.3264 0.2790  0.374
 # absent      Differentiation    Freshwater    0.0963 0.0722  0.117
 # absent      Homogenisation     Freshwater    0.2581 0.2192  0.303
 # exponential Differentiation    Terrestrial   0.4638 0.3726  0.538
 # exponential Homogenisation     Terrestrial   0.3330 0.2688  0.400
 # exponential Differentiation    Freshwater    0.4208 0.3428  0.502
 # exponential Homogenisation     Freshwater    0.2980 0.2274  0.357
 # revlog      Differentiation    Terrestrial   0.1030 0.0696  0.131
 # revlog      Homogenisation     Terrestrial   0.1112 0.0800  0.142
 # revlog      Differentiation    Freshwater    0.2652 0.1960  0.326
 # revlog      Homogenisation     Freshwater    0.1374 0.0993  0.175
 # saturating  Differentiation    Terrestrial   0.2667 0.1890  0.345
 # saturating  Homogenisation     Terrestrial   0.2208 0.1579  0.286
 # saturating  Differentiation    Freshwater    0.2079 0.1441  0.275
 # saturating  Homogenisation     Freshwater    0.2978 0.2213  0.378
 
 
 
 #Direction x Ecosystem type
 avg_comparisons(model_base_het_species,  variables= c("direction"), by = c("ecosystem_type"),
                 newdata = datagrid(newdata=model_base_het_species$data,
                                    direction = unique,
                                    ecosystem_type = unique,
                                    main_land_use_type=unique,
                                    het_trials=1), re_formula=NA,conf_level=0.8) 
 
 #Species HH
 #       Group      Term ecosystem_type Estimate  10.0 %   90.0 %
 # absent      direction    Freshwater   -0.1064 -0.1180 -0.09395
 # absent      direction    Terrestrial   0.0803  0.0711  0.09040
 # exponential direction    Freshwater   -0.0136 -0.0230 -0.00489
 # exponential direction    Terrestrial  -0.1102 -0.1300 -0.08764
 # saturating  direction    Freshwater    0.1779  0.1657  0.18956
 # saturating  direction    Terrestrial   0.0939  0.0819  0.10284
 # revlog      direction    Freshwater   -0.0565 -0.0745 -0.04263
 # revlog      direction    Terrestrial  -0.0623 -0.0826 -0.04780
 
 avg_comparisons(model_base_het_traits, variables= c("direction"), by = c("ecosystem_type"),
                 newdata = datagrid(newdata=model_base_het_traits$data,
                                    direction = unique,
                                    ecosystem_type = unique,
                                    main_land_use_type=unique,
                                    het_trials=1), re_formula=NA,conf_level=0.8) 
 #Trait HFP
 #       Group      Term ecosystem_type Estimate   10.0 %  90.0 %
 # absent      direction    Freshwater   0.16138  0.13880  0.1845
 # absent      direction    Terrestrial  0.16764  0.15416  0.1789
 # exponential direction    Freshwater  -0.12193 -0.14373 -0.0956
 # exponential direction    Terrestrial -0.12959 -0.14769 -0.1075
 # saturating  direction    Freshwater   0.08835  0.06977  0.1072
 # saturating  direction    Terrestrial -0.04540 -0.06447 -0.0301
 # revlog      direction    Freshwater  -0.12788 -0.15949 -0.1008
 # revlog      direction    Terrestrial  0.00787  0.00303  0.0128
 
 
 
 #Main land use type 
 
 #model for direction alone, its a different one. (negative = diff, positive = hom)
 avg_predictions(best_model_species, variables = "main_land_use_type", conf_level = .8)
 
 # main_land_use_type Estimate 10.0 % 90.0 %
 #        Multiple       -1.52 -2.037 -1.022
 #        Agriculture    -1.30 -1.787 -0.765
 #        Forest         -1.89 -2.374 -1.448
 #        Urban           2.33  0.702  3.996
 
 avg_predictions(best_model_traits, variables = "main_land_use_type",  conf_level = .8)
 
 # main_land_use_type Estimate 10.0 % 90.0 %
 #        Multiple       0.522  0.244  0.801
 #        Agriculture    0.547  0.236  0.892
 #        Forest         0.829  0.566  1.083
 #        Urban         -0.365 -1.270  0.541
 
avg_predictions(model_base_hfp_species, variables = c("main_land_use_type","direction"),
                 newdata = datagrid(newdata=model_base_hfp_species$data,
                                    direction = unique,
                                    ecosystem_type = unique,
                                    main_land_use_type=unique,
                                    hfp_trials=1), re_formula=NA,conf_level=0.8) |> as_tibble() |> print(n=32)
 
# A tibble: 32 × 6
#    group       main_land_use_type direction     estimate conf.low conf.high
#    <chr>       <fct>              <fct>            <dbl>    <dbl>     <dbl>
#  1 absent      Multiple           Differentiat…   0.512    0.446     0.573 
#  2 absent      Multiple           Homogenisati…   0.435    0.356     0.499 
#  3 absent      Agriculture        Differentiat…   0.271    0.220     0.320 
#  4 absent      Agriculture        Homogenisati…   0.452    0.386     0.507 
#  5 absent      Forest             Differentiat…   0.320    0.271     0.378 
#  6 absent      Forest             Homogenisati…   0.275    0.224     0.322 
#  7 absent      Urban              Differentiat…   0.465    0.399     0.531 
#  8 absent      Urban              Homogenisati…   0.361    0.296     0.421 
#  9 exponential Multiple           Differentiat…   0.218    0.158     0.273 
# 10 exponential Multiple           Homogenisati…   0.406    0.323     0.489 
# 11 exponential Agriculture        Differentiat…   0.293    0.222     0.361 
# 12 exponential Agriculture        Homogenisati…   0.136    0.0959    0.176 
# 13 exponential Forest             Differentiat…   0.215    0.153     0.266 
# 14 exponential Forest             Homogenisati…   0.187    0.135     0.239 
# 15 exponential Urban              Differentiat…   0.207    0.150     0.263 
# 16 exponential Urban              Homogenisati…   0.237    0.176     0.303 
# 17 revlog      Multiple           Differentiat…   0.0913   0.0623    0.121 
# 18 revlog      Multiple           Homogenisati…   0.0276   0.0174    0.0371
# 19 revlog      Agriculture        Differentiat…   0.135    0.0931    0.175 
# 20 revlog      Agriculture        Homogenisati…   0.187    0.134     0.235 
# 21 revlog      Forest             Differentiat…   0.101    0.0661    0.130 
# 22 revlog      Forest             Homogenisati…   0.118    0.0827    0.154 
# 23 revlog      Urban              Differentiat…   0.162    0.106     0.206 
# 24 revlog      Urban              Homogenisati…   0.279    0.202     0.343 
# 25 saturating  Multiple           Differentiat…   0.171    0.134     0.210 
# 26 saturating  Multiple           Homogenisati…   0.127    0.0902    0.156 
# 27 saturating  Agriculture        Differentiat…   0.291    0.232     0.349 
# 28 saturating  Agriculture        Homogenisati…   0.216    0.169     0.260 
# 29 saturating  Forest             Differentiat…   0.356    0.293     0.418 
# 30 saturating  Forest             Homogenisati…   0.412    0.339     0.472 
# 31 saturating  Urban              Differentiat…   0.155    0.118     0.190 
# 32 saturating  Urban              Homogenisati…   0.112    0.0802    0.138 

 
 avg_predictions(model_base_hfp_traits, variables = c("main_land_use_type","direction"),
                 newdata = datagrid(newdata=model_base_hfp_traits$data,
                                    direction = unique,
                                    ecosystem_type = unique,
                                    main_land_use_type=unique,
                                    hfp_trials=1), re_formula=NA,conf_level=0.8) |>  as_tibble() |> print(n=32)
 
#     group       main_land_use_type direction       estimate conf.low conf.high
#     <chr>       <fct>              <fct>              <dbl>    <dbl>     <dbl>
#  1 absent      Multiple           Differentiation   0.233    0.191     0.282 
#  2 absent      Multiple           Homogenisation    0.307    0.253     0.359 
#  3 absent      Agriculture        Differentiation   0.234    0.184     0.288 
#  4 absent      Agriculture        Homogenisation    0.216    0.177     0.264 
#  5 absent      Forest             Differentiation   0.245    0.193     0.296 
#  6 absent      Forest             Homogenisation    0.287    0.237     0.340 
#  7 absent      Urban              Differentiation   0.0755   0.0514    0.101 
#  8 absent      Urban              Homogenisation    0.387    0.326     0.458 
#  9 exponential Multiple           Differentiation   0.286    0.222     0.354 
# 10 exponential Multiple           Homogenisation    0.292    0.228     0.357 
# 11 exponential Agriculture        Differentiation   0.195    0.139     0.242 
# 12 exponential Agriculture        Homogenisation    0.270    0.204     0.332 
# 13 exponential Forest             Differentiation   0.135    0.0948    0.172 
# 14 exponential Forest             Homogenisation    0.286    0.222     0.351 
# 15 exponential Urban              Differentiation   0.444    0.338     0.547 
# 16 exponential Urban              Homogenisation    0.250    0.195     0.313 
# 17 revlog      Multiple           Differentiation   0.0543   0.0326    0.0763
# 18 revlog      Multiple           Homogenisation    0.0508   0.0309    0.0722
# 19 revlog      Agriculture        Differentiation   0.210    0.137     0.278 
# 20 revlog      Agriculture        Homogenisation    0.0973   0.0595    0.133 
# 21 revlog      Forest             Differentiation   0.128    0.0810    0.176 
# 22 revlog      Forest             Homogenisation    0.0762   0.0441    0.104 
# 23 revlog      Urban              Differentiation   0.435    0.318     0.547 
# 24 revlog      Urban              Homogenisation    0.149    0.0890    0.202 
# 25 saturating  Multiple           Differentiation   0.421    0.344     0.487 
# 26 saturating  Multiple           Homogenisation    0.343    0.276     0.405 
# 27 saturating  Agriculture        Differentiation   0.351    0.291     0.416 
# 28 saturating  Agriculture        Homogenisation    0.408    0.336     0.472 
# 29 saturating  Forest             Differentiation   0.483    0.419     0.554 
# 30 saturating  Forest             Homogenisation    0.341    0.271     0.399 
# 31 saturating  Urban              Differentiation   0.0403   0.0278    0.0528
# 32 saturating  Urban              Homogenisation    0.201    0.156     0.244 
#  
 avg_comparisons(model_base_hfp_species,  variables= c("direction"), by = c("main_land_use_type"),
                 newdata = datagrid(newdata=model_base_hfp_species$data,
                                    direction = unique,
                                    ecosystem_type = unique,
                                    main_land_use_type=unique,
                                    hfp_trials=1), re_formula=NA,conf_level=0.8) 
 
 #Species HFP
#        Group      Term main_land_use_type Estimate   10.0 %  90.0 %
#  absent      direction        Multiple     -0.0786 -0.10598 -0.0525
#  absent      direction        Agriculture   0.1811  0.16140  0.2037
#  absent      direction        Forest       -0.0441 -0.05480 -0.0344
#  absent      direction        Urban        -0.1034 -0.12240 -0.0835
#  exponential direction        Multiple      0.1881  0.16400  0.2149
#  exponential direction        Agriculture  -0.1555 -0.18908 -0.1263
#  exponential direction        Forest       -0.0269 -0.03540 -0.0194
#  exponential direction        Urban         0.0295  0.01480  0.0431
#  saturating  direction        Multiple     -0.0435 -0.05460 -0.0318
#  saturating  direction        Agriculture  -0.0756 -0.09725 -0.0542
#  saturating  direction        Forest        0.0554  0.04517  0.0653
#  saturating  direction        Urban        -0.0409 -0.05248 -0.0278
#  revlog      direction        Multiple     -0.0641 -0.08189 -0.0429
#  revlog      direction        Agriculture   0.0502  0.03419  0.0655
#  revlog      direction        Forest        0.0162  0.00942  0.0227
#  revlog      direction        Urban         0.1146  0.09153  0.1389
# 
# Type:  response 
# Comparison: mean(Homogenisation) - mean(Differentiation)

avg_comparisons(model_base_hfp_traits, variables= c("direction"), by = c("main_land_use_type"),
                 newdata = datagrid(newdata=model_base_hfp_traits$data,
                                    direction = unique,
                                    ecosystem_type = unique,
                                    main_land_use_type=unique,
                                    hfp_trials=1), re_formula=NA,conf_level=0.8) 
 #Trait HFP
 #       Group      Term main_land_use_type Estimate   10.0 %    90.0 %
 # absent      direction        Multiple     0.07454  0.06387  8.67e-02
 # absent      direction        Agriculture -0.01839 -0.03411 -2.86e-03
 # absent      direction        Forest       0.04395  0.02691  6.33e-02
 # absent      direction        Urban        0.31225  0.27054  3.58e-01
 # exponential direction        Multiple     0.00600 -0.00327  1.59e-02
 # exponential direction        Agriculture  0.07406  0.05639  9.08e-02
 # exponential direction        Forest       0.15089  0.12467  1.79e-01
 # exponential direction        Urban       -0.19055 -0.25173 -1.17e-01
 # saturating  direction        Multiple    -0.07683 -0.08867 -6.62e-02
 # saturating  direction        Agriculture  0.05706  0.03344  7.64e-02
 # saturating  direction        Forest      -0.14209 -0.16458 -1.23e-01
 # saturating  direction        Urban        0.16053  0.12575  1.94e-01
 # revlog      direction        Multiple    -0.00354 -0.00729 -2.62e-05
 # revlog      direction        Agriculture -0.11245 -0.14578 -7.67e-02
 # revlog      direction        Forest      -0.05075 -0.07190 -3.27e-02
 # revlog      direction        Urban       -0.28408 -0.34437 -2.19e-01

# Type:  response 
# Comparison: mean(Homogenisation) - mean(Differentiation)
 
avg_predictions(model_base_het_species, variables = c("main_land_use_type","direction"),
                 newdata = datagrid(newdata=model_base_het_species$data,
                                    direction = unique,
                                    ecosystem_type = unique,
                                    main_land_use_type=unique,
                                    het_trials=1), re_formula=NA,conf_level=0.8) |>  as_tibble() |>  print(n=32)
 
# A tibble: 32 × 6
#    group       main_land_use_type direction     estimate conf.low conf.high
#    <chr>       <fct>              <fct>            <dbl>    <dbl>     <dbl>
#  1 absent      Multiple           Differentiat…   0.139    0.125     0.155 
#  2 absent      Multiple           Homogenisati…   0.194    0.172     0.214 
#  3 absent      Agriculture        Differentiat…   0.167    0.151     0.183 
#  4 absent      Agriculture        Homogenisati…   0.140    0.126     0.156 
#  5 absent      Forest             Differentiat…   0.118    0.105     0.128 
#  6 absent      Forest             Homogenisati…   0.158    0.142     0.174 
#  7 absent      Urban              Differentiat…   0.342    0.309     0.374 
#  8 absent      Urban              Homogenisati…   0.220    0.184     0.251 
#  9 exponential Multiple           Differentiat…   0.435    0.334     0.518 
# 10 exponential Multiple           Homogenisati…   0.352    0.264     0.433 
# 11 exponential Agriculture        Differentiat…   0.272    0.195     0.344 
# 12 exponential Agriculture        Homogenisati…   0.389    0.294     0.484 
# 13 exponential Forest             Differentiat…   0.350    0.261     0.435 
# 14 exponential Forest             Homogenisati…   0.178    0.120     0.237 
# 15 exponential Urban              Differentiat…   0.152    0.103     0.196 
# 16 exponential Urban              Homogenisati…   0.0398   0.0232    0.0548
# 17 revlog      Multiple           Differentiat…   0.128    0.0897    0.162 
# 18 revlog      Multiple           Homogenisati…   0.179    0.133     0.227 
# 19 revlog      Agriculture        Differentiat…   0.187    0.136     0.232 
# 20 revlog      Agriculture        Homogenisati…   0.0698   0.0498    0.0912
# 21 revlog      Forest             Differentiat…   0.163    0.120     0.206 
# 22 revlog      Forest             Homogenisati…   0.141    0.102     0.180 
# 23 revlog      Urban              Differentiat…   0.174    0.128     0.220 
# 24 revlog      Urban              Homogenisati…   0.0244   0.0154    0.0329
# 25 saturating  Multiple           Differentiat…   0.292    0.221     0.367 
# 26 saturating  Multiple           Homogenisati…   0.267    0.205     0.333 
# 27 saturating  Agriculture        Differentiat…   0.366    0.290     0.440 
# 28 saturating  Agriculture        Homogenisati…   0.396    0.310     0.482 
# 29 saturating  Forest             Differentiat…   0.363    0.281     0.441 
# 30 saturating  Forest             Homogenisati…   0.517    0.445     0.596 
# 31 saturating  Urban              Differentiat…   0.323    0.262     0.386 
# 32 saturating  Urban              Homogenisati…   0.714    0.660     0.760 
#  
 
avg_predictions(model_base_het_traits, variables = c("main_land_use_type","direction"),
                 newdata = datagrid(newdata=model_base_het_traits$data,
                                    direction = unique,
                                    ecosystem_type = unique,
                                    main_land_use_type=unique,
                                    het_trials=1), re_formula=NA,conf_level=0.8) |> as_tibble() |>  print(n=32)
 
 
# A tibble: 32 × 6
#    group       main_land_use_type direction     estimate conf.low conf.high
#    <chr>       <fct>              <fct>            <dbl>    <dbl>     <dbl>
#  1 absent      Multiple           Differentiat…   0.199   0.156      0.241 
#  2 absent      Multiple           Homogenisati…   0.177   0.138      0.216 
#  3 absent      Agriculture        Differentiat…   0.0835  0.0616     0.102 
#  4 absent      Agriculture        Homogenisati…   0.158   0.120      0.189 
#  5 absent      Forest             Differentiat…   0.0418  0.0309     0.0514
#  6 absent      Forest             Homogenisati…   0.172   0.133      0.208 
#  7 absent      Urban              Differentiat…   0.186   0.141      0.229 
#  8 absent      Urban              Homogenisati…   0.660   0.598      0.732 
#  9 exponential Multiple           Differentiat…   0.407   0.320      0.476 
# 10 exponential Multiple           Homogenisati…   0.340   0.259      0.407 
# 11 exponential Agriculture        Differentiat…   0.474   0.389      0.552 
# 12 exponential Agriculture        Homogenisati…   0.429   0.347      0.504 
# 13 exponential Forest             Differentiat…   0.456   0.370      0.534 
# 14 exponential Forest             Homogenisati…   0.365   0.284      0.437 
# 15 exponential Urban              Differentiat…   0.430   0.343      0.510 
# 16 exponential Urban              Homogenisati…   0.127   0.0896     0.164 
# 17 revlog      Multiple           Differentiat…   0.138   0.0975     0.174 
# 18 revlog      Multiple           Homogenisati…   0.149   0.103      0.188 
# 19 revlog      Agriculture        Differentiat…   0.222   0.163      0.278 
# 20 revlog      Agriculture        Homogenisati…   0.191   0.136      0.237 
# 21 revlog      Forest             Differentiat…   0.279   0.211      0.346 
# 22 revlog      Forest             Homogenisati…   0.145   0.102      0.184 
# 23 revlog      Urban              Differentiat…   0.0967  0.0659     0.125 
# 24 revlog      Urban              Homogenisati…   0.0126  0.00831    0.0163
# 25 saturating  Multiple           Differentiat…   0.248   0.180      0.325 
# 26 saturating  Multiple           Homogenisati…   0.324   0.247      0.413 
# 27 saturating  Agriculture        Differentiat…   0.211   0.145      0.278 
# 28 saturating  Agriculture        Homogenisati…   0.212   0.144      0.274 
# 29 saturating  Forest             Differentiat…   0.213   0.149      0.285 
# 30 saturating  Forest             Homogenisati…   0.308   0.227      0.390 
# 31 saturating  Urban              Differentiat…   0.277   0.189      0.350 
# 32 saturating  Urban              Homogenisati…   0.192   0.133      0.253
  
 

avg_comparisons(model_base_het_species,  variables= c("direction"), by = c("main_land_use_type"),
                 newdata = datagrid(newdata=model_base_het_species$data,
                                    direction = unique,
                                    ecosystem_type = unique,
                                    main_land_use_type=unique,
                                    het_trials=1), re_formula=NA,conf_level=0.8) 
 
 #Species HH
 #       Group      Term main_land_use_type Estimate   10.0 %   90.0 %
 # absent      direction        Multiple      0.0546  0.04544  0.06249
 # absent      direction        Agriculture  -0.0270 -0.03785 -0.01643
 # absent      direction        Forest        0.0413  0.03124  0.05121
 # absent      direction        Urban        -0.1218 -0.15249 -0.09666
 # exponential direction        Multiple     -0.0811 -0.09207 -0.06922
 # exponential direction        Agriculture   0.1163  0.09256  0.14361
 # exponential direction        Forest       -0.1720 -0.20193 -0.14375
 # exponential direction        Urban        -0.1116 -0.14227 -0.07951
 # saturating  direction        Multiple     -0.0240 -0.03570 -0.01048
 # saturating  direction        Agriculture   0.0284  0.00391  0.04959
 # saturating  direction        Forest        0.1516  0.13826  0.16459
 # saturating  direction        Urban         0.3863  0.36592  0.40652
 # revlog      direction        Multiple      0.0504  0.03788  0.06221
 # revlog      direction        Agriculture  -0.1171 -0.14226 -0.08649
 # revlog      direction        Forest       -0.0221 -0.03604 -0.00874
 # revlog      direction        Urban        -0.1500 -0.18316 -0.10668

avg_comparisons(model_base_het_traits, variables= c("direction"), by = c("main_land_use_type"),
                 newdata = datagrid(newdata=model_base_het_traits$data,
                                    direction = unique,
                                    ecosystem_type = unique,
                                    main_land_use_type=unique,
                                    het_trials=1), re_formula=NA,conf_level=0.8) 
 #Trait HFP
#        Group      Term main_land_use_type  Estimate   10.0 %   90.0 %
#  absent      direction        Multiple    -0.020951 -0.02881 -0.01152
#  absent      direction        Agriculture  0.074590  0.05722  0.08816
#  absent      direction        Forest       0.130404  0.10272  0.15947
#  absent      direction        Urban        0.473465  0.44364  0.50649
#  exponential direction        Multiple    -0.065827 -0.07624 -0.05529
#  exponential direction        Agriculture -0.042621 -0.05617 -0.02994
#  exponential direction        Forest      -0.091995 -0.11550 -0.06952
#  exponential direction        Urban       -0.301024 -0.35326 -0.25082
#  saturating  direction        Multiple     0.075802  0.06195  0.09005
#  saturating  direction        Agriculture -0.000541 -0.00877  0.00793
#  saturating  direction        Forest       0.093735  0.07356  0.11281
#  saturating  direction        Urban       -0.083087 -0.11739 -0.04294
#  revlog      direction        Multiple     0.011107  0.00126  0.02049
#  revlog      direction        Agriculture -0.030908 -0.04517 -0.01789
#  revlog      direction        Forest      -0.133768 -0.16148 -0.10276
#  revlog      direction        Urban       -0.084176 -0.10863 -0.05579
# 
# Type:  response 
# Comparison: mean(Homogenisation) - mean(Differentiation)

