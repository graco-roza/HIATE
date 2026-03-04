#'###############################################################################
#' SCRIPT NAME: Modeling and Visualizing Relationship Shapes
#'
#' DESCRIPTION:
#'   This script models the shape of relationships between species and their 
#'   functional traits and environmental predictors. It includes testing various
#'   interaction effects and conducting K-fold cross-validation for model 
#'   comparison due to the multinomial nature of the response variable, which 
#'   precludes the use of standard LOO-CV (Leave-One-Out Cross-Validation).
#'   Additionally, the script generates main text and supplementary figures for 
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
#'   - Supplementary figures saved in `S7_Model_outputs_figures_and_tables/supplementary_figures/`
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
model_data <- process_model_data("S6_Synthesis_model/data/synthesis_data.xlsx", "shape")


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

compare_shape(
  model_base_hfp_species,
  grid_hfp_species
)

# |.category                |     .epred|     .lower|     .upper|
# |:------------------------|----------:|----------:|----------:|
# |exponential - absent     | -0.1498571| -0.2530423| -0.0502986|
# |revlog - absent          | -0.2462761| -0.3295768| -0.1729792|
# |revlog - exponential     | -0.0996516| -0.1842164| -0.0185743|
# |revlog - saturating      | -0.0919647| -0.1557749| -0.0256309|
# |saturating - absent      | -0.1551497| -0.2306059| -0.0681635|
# |saturating - exponential | -0.0081839| -0.1001608|  0.0797620|

compare_shape(
  model_base_hfp_traits,
  grid_hfp_traits
)

# |.category                |     .epred|     .lower|     .upper|
# |:------------------------|----------:|----------:|----------:|
# |exponential - absent     |  0.0194905| -0.0555846|  0.0982362|
# |revlog - absent          | -0.0978465| -0.1900776| -0.0207790|
# |revlog - exponential     | -0.1192972| -0.2161064| -0.0231896|
# |revlog - saturating      | -0.1740355| -0.2501104| -0.0993366|
# |saturating - absent      |  0.0734508|  0.0004250|  0.1666688|
# |saturating - exponential |  0.0535876| -0.0503984|  0.1613838|

compare_shape(
  model_base_het_species,
  grid_het_species
)

# |.category                |     .epred|     .lower|     .upper|
# |:------------------------|----------:|----------:|----------:|
# |exponential - absent     |  0.0854804|  0.0170618|  0.1534742|
# |revlog - absent          | -0.0517063| -0.0928446| -0.0103914|
# |revlog - exponential     | -0.1353334| -0.2218299| -0.0470283|
# |revlog - saturating      | -0.2689518| -0.3569478| -0.1835213|
# |saturating - absent      |  0.2191866|  0.1395471|  0.2835925|
# |saturating - exponential |  0.1335188|  0.0001094|  0.2641538|

compare_shape(
  model_base_het_traits,
  grid_het_traits
)

# |.category                |     .epred|     .lower|     .upper|
# |:------------------------|----------:|----------:|----------:|
# |exponential - absent     |  0.1648017|  0.0733641|  0.2545385|
# |revlog - absent          | -0.0555173| -0.1004087| -0.0141831|
# |revlog - exponential     | -0.2214001| -0.3132475| -0.1220499|
# |revlog - saturating      | -0.0929100| -0.1839139|  0.0032229|
# |saturating - absent      |  0.0365118| -0.0474154|  0.1311994|
# |saturating - exponential | -0.1307026| -0.2573557|  0.0019898|


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
  labs(x = "Human pressure gradient", y = "Community turnover", tag="a") +
  scale_y_continuous(breaks = seq(0, 1, 1), limits = c(0, 1.5)) +
  guides(x = ggplot2::guide_axis(cap = TRUE),
         y = ggplot2::guide_axis(cap = TRUE)) +
  scale_x_continuous(breaks = seq(0, 3, length.out = 5), labels = c("0", "", "", "", "1"), limits = c(0, 3)) +
  theme(plot.margin = margin(l = 2, t = 2, b = 2))  # Add margin to improve layout aesthetics


compare_shape <- function(model, grid, width = 0.80){
  
  epred_draws(model, newdata = grid, re_formula = NA) %>%
    group_by(.draw, .category) %>%
    summarise(.epred = mean(.epred), .groups = "drop") %>%
    compare_levels(.epred, by = .category) %>%
    # posterior summary
    median_hdci(.epred, .width = .8) |> 
    select(-.point, -.interval,-.width) |> 
    kableExtra::kable()
}


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
    x = "Functional turnover",
    y = "Taxonomic turnover",
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


direction_colors = c("Differentiation" = "#ff7f32","Homogenisation" = "#9f5cc0")  

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
    facet = "Taxonomic turnover",
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
    facet = "Functional turnover",
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
    facet = "Taxonomic turnover",
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
    facet = "Functional turnover",
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
    facet = "Taxonomic turnover",
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
    facet = "Functional turnover",
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
    facet = "Taxonomic turnover",
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
    facet = "Functional turnover",
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
  ) |> 
  select(draw,facet,response,class,Predictor,direction,shape,probability) |> 
  mutate(facet = factor(facet, levels=c("Taxonomic turnover", "Functional turnover")))

contrasts_shape |> 
  filter(class == "Ecosystem type") |> 
  group_by(draw,facet,response,shape) |> 
  summarise(probability=mean(probability)) |> 
  ungroup() |> 
  group_split(facet,response) |> 
  map(~ .x |> group_by(facet,response,shape) |>   median_hdci(probability,width=.8))

# Summary:
# Species × Habitat heterogeneity: mostly saturating.
# Species × Human footprint: mostly absent.
# Traits × Habitat heterogeneity: mostly exponential.
# Traits × Human footprint: mostly saturating.

fig_4c<-contrasts_shape |> 
  filter(class == "Ecosystem type") |> 
   mutate(shape = str_to_sentence(shape)) |> 
   mutate(shape = gsub("Revlog","Reverse logistic",shape)) |> 
  mutate(response = factor(response,levels=c("Human Footprint","Habitat Heterogeneity"))) |> 
  group_by(draw,facet,response, shape) |> 
  summarise(probability=mean(probability)) |> 
  ungroup() |> 
ggplot(aes(y=shape,x=probability)) + 
   ggstats::geom_stripped_rows(
     aes(y = shape),
     colour=NA,
     odd  = "white",
     even = "grey97",
     xfrom = -Inf, xto = Inf
   ) +
  stat_slabinterval(.width=c(.8,.89,.95), fill = "black", slab_alpha=.5,
                    point_interval="median_hdci", 
                    normalize="panels",
                    point_size = 1.2,      # smaller median point
                    interval_size_range = c(.5,1))+    # thinner interval lines
  facet_grid(response~facet)+
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
  ) + labs(x = "Probability of turnover shape",
           tag="c")+
  scale_x_continuous(limits=c(0,NA),
                     labels = scales::label_number(accuracy = 0.01)
  )


annotations_p2 <- data.frame(
  label = c(
    "<span>&larr; Dissimilarity</span>",
    "<span>Similarity &rarr;</span>"
  ),
  x      = rep(c(-0.01,  0.01),2),   # tweak these to your actual x‐limits
  y      = c( 4.7,     4.7,4.7,4.7),   # y position above your points
  angle  = c( 0,     0,0,0   ),
  hjust  = c( 1,     0,1,0   ),
  facet  = c("Taxonomic turnover", "Taxonomic turnover", "Functional turnover","Functional turnover")  # or repeat for "Trait replacement"
) |> 
  mutate(facet = factor(facet, levels=c("Taxonomic turnover", "Functional turnover")))
  


fig_4d<- contrasts_shape |> 
   mutate(shape = str_to_sentence(shape)) |> 
   mutate(shape = gsub("Revlog","Reverse logistic",shape)) |> 
   filter(class == "Ecosystem type") |> 
   mutate(response = factor(response,levels=c("Human Footprint","Habitat Heterogeneity"))) |> 
   group_by(draw,facet,response, direction, shape) |> 
   summarise(probability=mean(probability)) |> 
   ungroup() |> 
   pivot_wider(values_from = "probability",names_from = "direction") |> 
   mutate(p_direction = Homogenisation-Differentiation) |> 
   ggplot(aes(y=shape,x=p_direction)) + 
   ggstats::geom_stripped_rows(
     aes(y = shape),
     colour=NA,
     odd  = "white",
     even = "grey97",
     xfrom = -Inf, xto = Inf
   ) +
   stat_slab(
     aes(fill = after_stat(ifelse(x >= 0,
                                  "Homogenisation",
                                  "Differentiation"))),
     alpha = 1,
     colour = NA,
     normalize="panels",
     height=0.6
   ) +
   stat_pointinterval(
     point_interval = "median_hdci",
     .width = c(.8, .89, .95),
     point_size = 1.2,
     interval_size_range = c(.5,1),
     colour = "black"
   )+
  geom_richtext(
    data = annotations_p2,
    aes(x = x, y = y, label = label, angle = angle, hjust = hjust), 
    size=convert_size(8),
    label.color = NA,
    fill=NA,
    inherit.aes=FALSE
  ) +
  geom_vline(xintercept=0, linetype="11", linewidth=.5, colour = "gray40")+
  scale_fill_manual(
    name = "Turnover direction",
    values = direction_colors,
    breaks = c("Homogenisation", "Differentiation"),
    labels = c("Similarity increase", "Dissimilarity increase")
  )+
     facet_grid(response~facet)+
   theme_void(base_family = "sans", base_size = 12) +
   theme(
     panel.grid.minor   = element_blank(),
     plot.background    = element_blank(),
     panel.background =  element_blank(),
     strip.text.x =  element_text(face="bold",margin=margin(b=10)),
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
     panel.spacing.y = unit(1, "lines"),
     plot.margin = margin(0,0,0,0),
     legend.position = c(0.5,.06),
     legend.direction = "vertical",
     legend.background = element_blank(),
     legend.box.background = element_blank(),
     legend.title = element_blank(),
     legend.text = element_text(size = 8),
     legend.title.position = "left",
     legend.key.size = unit(0.6, "lines"),
     legend.spacing.x = unit(0.4, "lines"),
     legend.spacing.y = unit(0.2, "lines")
   ) +
   guides(
     x = guide_axis(cap = TRUE),
     fill = guide_legend(
       title.position = "top",
       keywidth = unit(0.6, "lines"),
       keyheight = unit(0.6, "lines")
     )
   ) + labs(x = "Association of turnover shape with direction",
            tag="d")+
   scale_x_continuous(limits=c(NA,NA),
                      labels = scales::label_number(accuracy = 0.01)
   )+
  coord_cartesian(clip = "off")


Figure_4<- ((fig_4a/fig_4b) | (fig_4c/fig_4d)) + plot_layout(widths=c(.25,.75))

ggsave(filename = here::here("S7_Model_outputs_figures_and_tables", "main_figures", "Figure_4.pdf"),
       plot = Figure_4,
       device = cairo_pdf,
       width=11,height=9,units="in")


# fig_4c<-contrasts_shape |> 
#   mutate(shape = str_to_sentence(shape)) |> 
#   mutate(shape = gsub("Revlog","Reverse logistic",shape)) |> 
#   mutate(
#     Predictor = case_when(
#       Predictor == "Freshwater"  ~ "Freshwater (n=30)",
#       Predictor == "Terrestrial" ~ "Terrestrial (n=131)",
#       Predictor == "Agriculture" ~ "Agriculture (n=38)",
#       Predictor == "Forest"      ~ "Forest (n=56)",
#       Predictor == "Urban"       ~ "Urban (n=15)",
#       Predictor == "Multiple"    ~ "Multiple (n=52)",
#       TRUE                        ~ Predictor
#     ),
#     Predictor = factor(
#       Predictor,
#       levels = c(
#         "Freshwater (n=30)",
#         "Terrestrial (n=131)",
#         rev(c(
#           "Agriculture (n=38)",
#           "Forest (n=56)",
#           "Urban (n=15)",
#           "Multiple (n=52)"
#         ))
#       )
#     )
#   ) |> 
#   mutate(response = factor(response,levels=c("Human Footprint","Habitat Heterogeneity"))) |> 
#   ggplot(aes(y=shape,x=probability,colour=direction, fill=direction))+
#   ggstats::geom_stripped_rows(
#     aes(y = shape),
#     colour=NA,
#     odd  = "white",
#     even = "grey90",
#     xfrom = -Inf, xto = Inf
#   ) +
#   stat_ccdfinterval(position="dodge", .width=c(.8,.89,.95), slab_alpha=.5) +
#   ggh4x::facet_nested(class+Predictor~facet+response, scales="free_y",space="free_y",
#                       # use strip_themed() instead of strip_nested()
#                       strip = strip_nested(
#                         # supply one size for the "class" layer, one for the "Predictor" layer
#                         text_y      = elem_list_text(size = c(10, 8)),
#                         by_layer_y = TRUE,
#                         text_x      = elem_list_text(size = c(10, 8)),
#                         by_layer_x = TRUE,
#                       ))+
#   scale_colour_manual(values=direction_colors)+
#   scale_fill_manual(values=direction_colors)+
#   theme_void(base_family = "sans", base_size = 12) +
#   theme(
#     panel.grid.minor   = element_blank(),
#     plot.background    = element_blank(),
#     panel.background =  element_blank(),
#     strip.text.x =  element_text(face="bold",margin=margin(t=5)),
#     strip.text.y = element_text(size=10,face="bold",angle=270, margin=margin(l=5)),
#     axis.title.x       = ggtext::element_markdown(face = "bold", size = 12, margin=margin(t=5)),
#     axis.text.y        = element_text(size = 10, hjust = 1),
#     axis.text.x        = ggtext::element_markdown(size = 8, margin=margin(t=3)),
#     axis.title.y       = element_blank(),
#     axis.line.x        = element_line(color = "grey20", linewidth = 0.1),
#     axis.ticks.x       = element_line(color = "grey20", linewidth = 0.1),
#     axis.ticks.length  = unit(0.2, "lines"),
#     plot.tag           = element_text(face = "bold", size = 14),
#     panel.spacing.x = unit(2,"lines"),
#     plot.margin = margin(0,0,0,0),
#     legend.position="bottom",
#     legend.title = element_text(face="bold")
#   ) +
#   guides(
#     x = guide_axis(cap = TRUE)
#   ) + labs(x = "Posterior predictions (Probability of observing a shape)",
#            colour="Relationship direction",
#            fill="Relationship direction",
#            tag="c")+
#   scale_x_continuous(limits=c(0,NA),
#                      labels = scales::label_number(accuracy = 0.01)
#   )
# 
# 
# 
# 
# ggsave(filename = here::here("S7_Model_outputs_figures_and_tables", "main_figures", "Figure_4.pdf"),
#        plot = Figure_4,
#        device = cairo_pdf,
#        width=11,height=9,units="in")


#Interactions figure 

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

 supplementary_figure_shape_posteriors <- 
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

 ggsave(filename = here::here("S7_Model_outputs_figures_and_tables", "supplementary_figures", "Supplementary_Figure_shape_posteriors.pdf"),
        plot = supplementary_figure_shape_posteriors,
        device = cairo_pdf,
        width=14,height=14,units="in")

 
 
 
 
 
 

 