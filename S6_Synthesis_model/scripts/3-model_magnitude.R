#' ###########################################################################
#' SCRIPT NAME: Bayesian Modeling of Magnitude of Effects (3-model_magnitude.R)
#'
#' DESCRIPTION:
#'   This script performs Bayesian modeling to estimate the magnitude of species 
#'   and Functional turnover across different ecosystems, land-use types, and other 
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


model_data <- process_model_data("S6_Synthesis_model/data/synthesis_data_complete.xlsx", "magnitude")
#Set priors ------------------------------

model_data<- model_data |> filter(framework == "Podani",
                                  feature == "abun")

# Statistical Test for Magnitude Difference ---------------------------------

# Conduct a Wilcoxon signed-rank test to check if the magnitudes differ significantly 
# between Species and Traits across datasets. This test accounts for paired observations.
# 

#check if Human footprint magnitude is larger in species or traits

model_data |> 
  select(dataset,facet,hfp_magnitude) |> 
  group_by(facet) |> 
  summarise(mean_hfp_magnitude = mean(hfp_magnitude), sd_hfp_magnitude=sd(hfp_magnitude))

# # A tibble: 2 × 3
#   facet   mean_hfp_magnitude sd_hfp_magnitude
#   <fct>                <dbl>            <dbl>
# 1 Species              0.200            0.300
# 2 Traits               0.173            0.267

model_data |>
  select(dataset, facet, hfp_magnitude) |>
  pivot_wider(names_from = facet, values_from = hfp_magnitude) |>
  mutate(higher = Species < Traits) |>
  janitor::tabyl(higher) |>
  adorn_pct_formatting(digits = 2) |>
  kableExtra::kable(format = "rst")


#How many times traits are more impacted than species? (human footprint)
# ======  ===  =======
# higher    n  percent
# ======  ===  =======
# FALSE    96  59.63% 
# TRUE     65  40.37% 
# ======  ===  =======

model_data |> 
  select(dataset,facet,hfp_magnitude) |> 
  group_by(dataset) |> 
  filter(n_distinct(facet) == 2)  |>  
  ungroup() |> 
  rstatix::wilcox_test(formula=hfp_magnitude~facet,
                       paired = TRUE,
                       alternative = "two.sided")

#Are differences significant in human footprint effect on species and traits?
# # A tibble: 1 × 7
#   .y.           group1  group2    n1    n2 statistic     p
# * <chr>         <chr>   <chr>  <int> <int>     <dbl> <dbl>
# 1 hfp_magnitude Species Traits   161   161      6395  0.17

#check if Habitat heterogeneity magnitude is larger in species or traits

model_data |> 
  select(dataset,facet,het_magnitude) |> 
  group_by(facet) |> 
  summarise(mean_het_magnitude = mean(het_magnitude), sd_het_magnitude=sd(het_magnitude))

# # A tibble: 2 × 3
#   facet   mean_het_magnitude sd_het_magnitude
#   <fct>                <dbl>            <dbl>
# 1 Species              0.167            0.197
# 2 Traits               0.220            0.371

model_data |>
  select(dataset, facet, het_magnitude) |>
  pivot_wider(names_from = facet, values_from = het_magnitude) |>
  mutate(higher = Species < Traits) |>
  janitor::tabyl(higher) |>
  adorn_pct_formatting(digits = 2) |>
  kableExtra::kable(format = "rst")

#How many times traits are more impacted than species? (habitat heterogeneity)
# ======  ===  =======
# higher    n  percent
# ======  ===  =======
# FALSE    93  57.76% 
# TRUE     68  42.24% 
# ======  ===  =======

model_data |> 
  select(dataset,facet,het_magnitude) |> 
  group_by(dataset) |> 
  filter(n_distinct(facet) == 2)  |>  
  ungroup() |> 
  rstatix::wilcox_test(formula=het_magnitude~facet,
                       paired = TRUE,
                       alternative = "two.sided")

#Are differences significant in habitat heterogeneity effect on species and traits?
# # A tibble: 1 × 7
#   .y.           group1  group2    n1    n2 statistic     p
# * <chr>         <chr>   <chr>  <int> <int>     <dbl> <dbl>
# 1 het_magnitude Species Traits   161   161      6151  0.93

model_data |> 
  group_split(facet) |> 
  map(~ .x |>  select(hfp_magnitude, het_magnitude) |>  cor(method="s"))

#Check if human footprint effect somehow correlates with habitat heterogeneity effect
# species
# hfp_magnitude het_magnitude
# hfp_magnitude    1.00000000    0.09657465
# het_magnitude    0.09657465    1.00000000
# 
# traits
# hfp_magnitude het_magnitude
# hfp_magnitude     1.0000000     0.1030768
# het_magnitude     0.1030768     1.0000000

priors_magnitude <- c(
  # slope priors
  prior(normal(0, 1), class = "b", resp = "hfpmagnitude"),
  prior(normal(0, 1), class = "b", resp = "hetmagnitude"),
  
  # group-level (sd) priors
  prior(student_t(3, 0, 0.1), class = "sd", resp = "hfpmagnitude"),
  prior(student_t(3, 0, 0.1), class = "sd", resp = "hetmagnitude")
  
)

#Set model formula --------------------------



formula_base <- bf(mvbind(hfp_magnitude,het_magnitude) ~ species_number + spatial_extent + spatial_distance_min + absolute_latitude_mean + habitat_heterogeneity_baseline + habitat_heterogeneity_range + human_footprint_baseline + human_footprint_range + (main_land_use_type + ecosystem_type)*direction + het_buffer + hfp_buffer+  (1 | taxa_coarse), decomp = "QR")
formula_base_trait <- bf(mvbind(hfp_magnitude,het_magnitude) ~ species_number + spatial_extent + spatial_distance_min + absolute_latitude_mean + habitat_heterogeneity_baseline + habitat_heterogeneity_range + human_footprint_baseline + human_footprint_range + (main_land_use_type + ecosystem_type)*direction+ het_buffer + hfp_buffer+   (1 | taxa_coarse), sigma ~ hypervolume_quality, decomp = "QR")

#Run model ------------------------------------------

md_species_base <- brms::brm(
formula_base + set_rescor(FALSE),
data = filter(model_data, facet == "Species"), #subset only the taxonomic parts
prior = priors_magnitude,
family=brms::hurdle_lognormal(),
chains = CHAINS,
iter = ITER,
warmup = WARMUP,
control = base::list(max_treedepth = 15, adapt_delta = .99),
save_pars = brms::save_pars(all = TRUE),
seed = BAYES_SEED,
file = "S7_Model_outputs_figures_and_tables/model/Podani/Abundance/magnitude/magnitude_species_base", #here we save models automatically
file_refit = "on_change" #models are only refit if parameters change, otherwise load the saved ones
) 

# Run Models for Trait Facet ------------------------------------------

md_trait_base <- brms::brm(
  formula_base_trait + set_rescor(FALSE),
  data = filter(model_data, facet == "Traits"), # subset only the trait parts
  prior = priors_magnitude,
  family = brms::hurdle_lognormal(),
  chains = CHAINS,
  iter = ITER,
  warmup = WARMUP,
  control = list(max_treedepth = 15, adapt_delta = .99),
  save_pars = brms::save_pars(all = TRUE),
  seed = BAYES_SEED,
  file = "S7_Model_outputs_figures_and_tables/model/Podani/Abundance/magnitude/magnitude_trait_base", 
  file_refit = "on_change"
)

best_model_species <- md_species_base # md_species_base is chosen as the final model because it achieves similar predictive accuracy to more complex models but with fewer parameters
best_model_traits <- md_trait_base # md_traits_base is chosen as the final model because it achieves similar predictive accuracy but with fewer parameters

# Run the posterior predictive checks -----------------------------------------------------------------------------
cowplot::plot_grid(brms::pp_check(best_model_species, ndraws=100, resp="hfpmagnitude")+coord_cartesian(xlim=c(0,3)),
                   brms::pp_check(best_model_traits, ndraws=100,resp="hfpmagnitude")+coord_cartesian(xlim=c(0,3)),
                   brms::pp_check(best_model_species, ndraws=100, resp="hetmagnitude")+coord_cartesian(xlim=c(0,3)),
                   brms::pp_check(best_model_traits, ndraws=100,resp="hetmagnitude")+coord_cartesian(xlim=c(0,3))
                   )
#'-----------------------------------------------------------------------------------------------------------------

# check R2 of each model ------------------------------------------------------------------------------------------
brms::bayes_R2(best_model_species)
#                 Estimate  Est.Error      Q2.5     Q97.5
# R2hfpmagnitude 0.3254577 0.08753501 0.1643733 0.4935231
# R2hetmagnitude 0.4298391 0.06844250 0.2799503 0.5415490

brms::bayes_R2(best_model_traits) 
#                 Estimate  Est.Error      Q2.5     Q97.5
# R2hfpmagnitude 0.3998435 0.07739485 0.2346578 0.5293771
# R2hetmagnitude 0.3917507 0.05826916 0.2665438 0.4880509
#'-----------------------------------------------------------------------------------------------------------------


bayestestR::ci(best_model_species, ci=.8, method="HDI")

# Parameter                                             |     Response |        80% HDI
# -------------------------------------------------------------------------------------
# (Intercept)                                           | hfpmagnitude | [-2.90, -1.47]
# species_number                                        | hfpmagnitude | [-0.14,  0.12]
# spatial_extent                                        | hfpmagnitude | [-0.35,  0.06]
# spatial_distance_min                                  | hfpmagnitude | [-0.08,  0.15]
# absolute_latitude_mean                                | hfpmagnitude | [-0.12,  0.17]
# habitat_heterogeneity_baseline                        | hfpmagnitude | [ 0.33,  1.77]
# habitat_heterogeneity_range                           | hfpmagnitude | [-0.49,  0.44]
# human_footprint_baseline                              | hfpmagnitude | [-0.06,  0.00]
# human_footprint_range                                 | hfpmagnitude | [-0.01,  0.02]
# main_land_use_typeAgriculture                         | hfpmagnitude | [ 0.57,  1.32]
# main_land_use_typeForest                              | hfpmagnitude | [ 0.16,  0.87]
# main_land_use_typeUrban                               | hfpmagnitude | [-0.12,  1.78]
# ecosystem_typeFreshwater                              | hfpmagnitude | [ 0.04,  0.78]
# directionHomogenisation                               | hfpmagnitude | [ 0.94,  2.05]
# het_buffer1000                                        | hfpmagnitude | [-0.99, -0.36]
# het_buffer1500                                        | hfpmagnitude | [-0.98, -0.28]
# het_buffer2000                                        | hfpmagnitude | [-1.02, -0.36]
# hfp_buffer1500                                        | hfpmagnitude | [-0.28,  0.25]
# hfp_buffer2000                                        | hfpmagnitude | [-0.03,  0.54]
# main_land_use_typeAgriculture:directionHomogenisation | hfpmagnitude | [-2.55, -0.87]
# main_land_use_typeForest:directionHomogenisation      | hfpmagnitude | [-1.78, -0.53]
# main_land_use_typeUrban:directionHomogenisation       | hfpmagnitude | [-1.07,  0.65]
# ecosystem_typeFreshwater:directionHomogenisation      | hfpmagnitude | [-0.92,  0.47]
# 
# # Fixed effects () (hetmagnitude)
# 
# Parameter                                             |     Response |        80% HDI
# -------------------------------------------------------------------------------------
# (Intercept)                                           | hetmagnitude | [-3.71, -2.06]
# species_number                                        | hetmagnitude | [-0.06,  0.24]
# spatial_extent                                        | hetmagnitude | [-0.24,  0.22]
# spatial_distance_min                                  | hetmagnitude | [ 0.05,  0.31]
# absolute_latitude_mean                                | hetmagnitude | [-0.13,  0.22]
# habitat_heterogeneity_baseline                        | hetmagnitude | [ 0.28,  1.92]
# habitat_heterogeneity_range                           | hetmagnitude | [-0.13,  0.94]
# human_footprint_baseline                              | hetmagnitude | [-0.04,  0.03]
# human_footprint_range                                 | hetmagnitude | [-0.03,  0.00]
# main_land_use_typeAgriculture                         | hetmagnitude | [-0.48,  0.39]
# main_land_use_typeForest                              | hetmagnitude | [-0.44,  0.38]
# main_land_use_typeUrban                               | hetmagnitude | [-2.29,  0.13]
# ecosystem_typeFreshwater                              | hetmagnitude | [-1.04, -0.22]
# directionHomogenisation                               | hetmagnitude | [-0.12,  1.14]
# het_buffer1000                                        | hetmagnitude | [ 0.22,  0.94]
# het_buffer1500                                        | hetmagnitude | [ 0.44,  1.32]
# het_buffer2000                                        | hetmagnitude | [ 0.00,  0.81]
# hfp_buffer1500                                        | hetmagnitude | [-0.18,  0.40]
# hfp_buffer2000                                        | hetmagnitude | [-0.45,  0.17]
# main_land_use_typeAgriculture:directionHomogenisation | hetmagnitude | [-0.92,  0.88]
# main_land_use_typeForest:directionHomogenisation      | hetmagnitude | [-0.69,  0.81]
# main_land_use_typeUrban:directionHomogenisation       | hetmagnitude | [-0.66,  1.44]
# ecosystem_typeFreshwater:directionHomogenisation      | hetmagnitude | [-0.40,  1.19]

bayestestR::ci(best_model_traits, ci=.8, method="HDI")

# Highest Density Interval () (hfpmagnitude)
# 
# Parameter                                             |     Response |        80% HDI
# -------------------------------------------------------------------------------------
# (Intercept)                                           | hfpmagnitude | [-5.85, -3.66]
# species_number                                        | hfpmagnitude | [-0.50, -0.17]
# spatial_extent                                        | hfpmagnitude | [-0.15,  0.35]
# spatial_distance_min                                  | hfpmagnitude | [-0.01,  0.26]
# absolute_latitude_mean                                | hfpmagnitude | [-0.26,  0.13]
# habitat_heterogeneity_baseline                        | hfpmagnitude | [-0.91,  0.69]
# habitat_heterogeneity_range                           | hfpmagnitude | [-0.53,  0.68]
# human_footprint_baseline                              | hfpmagnitude | [-0.03,  0.03]
# human_footprint_range                                 | hfpmagnitude | [ 0.00,  0.04]
# main_land_use_typeAgriculture                         | hfpmagnitude | [ 0.09,  1.62]
# main_land_use_typeForest                              | hfpmagnitude | [-0.33,  1.15]
# main_land_use_typeUrban                               | hfpmagnitude | [ 0.27,  2.50]
# ecosystem_typeFreshwater                              | hfpmagnitude | [-0.19,  1.01]
# directionHomogenisation                               | hfpmagnitude | [ 1.69,  2.91]
# het_buffer1000                                        | hfpmagnitude | [-0.66,  0.19]
# het_buffer1500                                        | hfpmagnitude | [-0.70,  0.20]
# het_buffer2000                                        | hfpmagnitude | [-0.22,  0.61]
# hfp_buffer1500                                        | hfpmagnitude | [-0.06,  0.66]
# hfp_buffer2000                                        | hfpmagnitude | [-0.02,  0.65]
# main_land_use_typeAgriculture:directionHomogenisation | hfpmagnitude | [-1.12,  0.58]
# main_land_use_typeForest:directionHomogenisation      | hfpmagnitude | [-1.17,  0.42]
# main_land_use_typeUrban:directionHomogenisation       | hfpmagnitude | [-2.48, -0.54]
# ecosystem_typeFreshwater:directionHomogenisation      | hfpmagnitude | [-1.54,  0.01]
# 
# # Fixed effects () (hetmagnitude)
# 
# Parameter                                             |     Response |        80% HDI
# -------------------------------------------------------------------------------------
#   (Intercept)                                           | hetmagnitude | [-6.64, -4.52]
# species_number                                        | hetmagnitude | [-0.40, -0.08]
# spatial_extent                                        | hetmagnitude | [-0.35,  0.10]
# spatial_distance_min                                  | hetmagnitude | [-0.08,  0.20]
# absolute_latitude_mean                                | hetmagnitude | [-0.46, -0.09]
# habitat_heterogeneity_baseline                        | hetmagnitude | [ 0.52,  1.96]
# habitat_heterogeneity_range                           | hetmagnitude | [ 0.29,  1.39]
# human_footprint_baseline                              | hetmagnitude | [-0.04,  0.03]
# human_footprint_range                                 | hetmagnitude | [-0.02,  0.02]
# main_land_use_typeAgriculture                         | hetmagnitude | [ 0.67,  2.08]
# main_land_use_typeForest                              | hetmagnitude | [ 0.16,  1.61]
# main_land_use_typeUrban                               | hetmagnitude | [-0.91,  1.37]
# ecosystem_typeFreshwater                              | hetmagnitude | [ 0.37,  1.48]
# directionHomogenisation                               | hetmagnitude | [ 1.64,  2.83]
# het_buffer1000                                        | hetmagnitude | [ 0.02,  0.80]
# het_buffer1500                                        | hetmagnitude | [-0.02,  0.78]
# het_buffer2000                                        | hetmagnitude | [ 0.49,  1.28]
# hfp_buffer1500                                        | hetmagnitude | [-0.48,  0.24]
# hfp_buffer2000                                        | hetmagnitude | [-0.70, -0.08]
# main_land_use_typeAgriculture:directionHomogenisation | hetmagnitude | [-1.36,  0.24]
# main_land_use_typeForest:directionHomogenisation      | hetmagnitude | [-1.28,  0.30]
# main_land_use_typeUrban:directionHomogenisation       | hetmagnitude | [-2.19,  0.39]
# ecosystem_typeFreshwater:directionHomogenisation      | hetmagnitude | [-0.63,  0.76]



options("marginaleffects_posterior_interval" = "hdi",
        "marginaleffects_posterior_center" = "median")

#Contrast

avg_predictions(best_model_species,  variables= c("direction"),
                newdata = datagrid(
                                   direction = unique,
                                   ecosystem_type = unique,
                                   main_land_use_type=unique),
                re_formula=NA,
                conf_level=0.8) 


 #        Group       direction Estimate 10.0 % 90.0 %
 # hetmagnitude Differentiation    0.150 0.0963  0.196
 # hetmagnitude Homogenisation     0.334 0.2034  0.462
 # hfpmagnitude Differentiation    0.174 0.1035  0.247
 # hfpmagnitude Homogenisation     0.356 0.1795  0.532



avg_comparisons(best_model_species,  variables= c("direction"),
                newdata = datagrid(
                  direction = unique,
                  ecosystem_type = unique,
                  main_land_use_type=unique),
                re_formula=NA,
                conf_level=0.8) 

 #        Group Estimate 10.0 % 90.0 %
 # hetmagnitude    0.181 0.0690  0.290
 # hfpmagnitude    0.178 0.0245  0.336

# Comparison: mean(Homogenisation) - mean(Differentiation)

avg_predictions(best_model_traits, variables= c("direction"), 
                newdata = datagrid(newdata=best_model_traits$data,
                                   direction = unique,
                                   ecosystem_type = unique,
                                   main_land_use_type=unique),
                re_formula=NA,
                conf_level=0.8)

 #        Group       direction Estimate 10.0 % 90.0 %
 # hetmagnitude Differentiation    0.145 0.0750  0.217
 # hetmagnitude Homogenisation     0.852 0.4863  1.219
 # hfpmagnitude Differentiation    0.077 0.0326  0.119
 # hfpmagnitude Homogenisation     0.274 0.1614  0.372

avg_comparisons(best_model_traits, variables= c("direction"), 
                newdata = datagrid(newdata=best_model_traits$data,
                                   direction = unique,
                                   ecosystem_type = unique,
                                   main_land_use_type=unique), 
                re_formula=NA,
                conf_level=0.8) 

#Traits
 #        Group Estimate 10.0 % 90.0 %
 # hetmagnitude    0.700 0.3633  1.043
 # hfpmagnitude    0.189 0.0933  0.287

avg_comparisons(best_model_species,  variables= c("ecosystem_type"),
                newdata = datagrid(newdata=best_model_species$data,
                                   direction = unique,
                                   ecosystem_type = unique,
                                   main_land_use_type=unique), 
                re_formula=NA,
                conf_level=0.8) 

 #        Group Estimate  10.0 % 90.0 %
 # hetmagnitude  -0.0739 -0.2049  0.044
 # hfpmagnitude   0.0833 -0.0521  0.221

avg_comparisons(best_model_traits, variables= c("ecosystem_type"), 
                newdata = datagrid(newdata=best_model_traits$data,
                                   direction = unique,
                                   ecosystem_type = unique,
                                   main_land_use_type=unique),
                re_formula=NA,
                conf_level=0.8) 

 #        Group Estimate 10.0 % 90.0 %
 # hetmagnitude   0.4409  0.130 0.7337
 # hfpmagnitude  -0.0196 -0.101 0.0547

avg_comparisons(best_model_species,  variables= list(main_land_use_type = "pairwise"),
                newdata = datagrid(newdata=best_model_species$data,
                                   direction = unique,
                                   ecosystem_type = unique,
                                   main_land_use_type=unique),
                re_formula=NA,
                conf_level=0.8) 

 #        Group               Contrast Estimate  10.0 % 90.0 %
 # hetmagnitude Agriculture - Multiple -0.02051 -0.1744 0.1365
 # hetmagnitude Forest - Agriculture    0.01138 -0.1334 0.1512
 # hetmagnitude Forest - Multiple      -0.00937 -0.1437 0.1160
 # hetmagnitude Urban - Agriculture    -0.10512 -0.3220 0.1111
 # hetmagnitude Urban - Forest         -0.11844 -0.3190 0.1071
 # hetmagnitude Urban - Multiple       -0.12757 -0.3393 0.0739
 # hfpmagnitude Agriculture - Multiple -0.02434 -0.1325 0.0781
 # hfpmagnitude Forest - Agriculture   -0.02455 -0.0956 0.0457
 # hfpmagnitude Forest - Multiple      -0.04898 -0.1402 0.0369
 # hfpmagnitude Urban - Agriculture     0.23699 -0.0971 0.6366
 # hfpmagnitude Urban - Forest          0.26304 -0.0859 0.6571
 # hfpmagnitude Urban - Multiple        0.21277 -0.1420 0.6275

avg_comparisons(best_model_traits, variables= list(main_land_use_type = "pairwise"), 
                newdata = datagrid(newdata=best_model_traits$data,
                                   direction = unique,
                                   ecosystem_type = unique,
                                   main_land_use_type=unique),
                re_formula=NA,
                conf_level=0.8) 

 #        Group               Contrast Estimate   10.0 %   90.0 %
 # hetmagnitude Agriculture - Multiple   0.4502  0.10484  0.79886
 # hetmagnitude Forest - Agriculture    -0.2560 -0.56182  0.05626
 # hetmagnitude Forest - Multiple        0.1924 -0.02904  0.37240
 # hetmagnitude Urban - Agriculture     -0.5352 -1.08492 -0.05938
 # hetmagnitude Urban - Forest          -0.2970 -0.70103  0.09999
 # hetmagnitude Urban - Multiple        -0.1143 -0.40860  0.22498
 # hfpmagnitude Agriculture - Multiple   0.1047  0.00498  0.20491
 # hfpmagnitude Forest - Agriculture    -0.0871 -0.17514  0.00228
 # hfpmagnitude Forest - Multiple        0.0163 -0.03979  0.07114
 # hfpmagnitude Urban - Agriculture     -0.0534 -0.24556  0.12247
 # hfpmagnitude Urban - Forest           0.0328 -0.12058  0.21039
 # hfpmagnitude Urban - Multiple         0.0498 -0.09666  0.21822

#Habitat heterogeneity effect on species shape across land-use types for each direction

#Direction x Ecosystem type
avg_comparisons(best_model_species,  variables= c("direction"), by = c("ecosystem_type"),
                newdata = datagrid(newdata=best_model_species$data,
                                   direction = unique,
                                   ecosystem_type = unique,
                                   main_land_use_type=unique), 
                re_formula=NA,
                conf_level=0.8) 
#Species
 #        Group ecosystem_type Estimate   10.0 % 90.0 %
 # hetmagnitude    Terrestrial    0.166  0.04172  0.284
 # hfpmagnitude    Terrestrial    0.168  0.05693  0.282
 # hetmagnitude    Freshwater     0.182  0.00218  0.360
 # hfpmagnitude    Freshwater     0.178 -0.07281  0.445

avg_comparisons(best_model_traits, variables= c("direction"), by = c("ecosystem_type"),
                newdata = datagrid(newdata=best_model_traits$data,
                                   direction = unique,
                                   ecosystem_type = unique,
                                   main_land_use_type=unique),
                re_formula=NA,
                conf_level=0.8) 
#Trait 
 #        Group ecosystem_type Estimate  10.0 % 90.0 %
 # hetmagnitude    Terrestrial    0.381  0.2308  0.529
 # hfpmagnitude    Terrestrial    0.247  0.1445  0.339
 # hetmagnitude    Freshwater     1.004  0.3727  1.579
 # hfpmagnitude    Freshwater     0.128 -0.0144  0.272

#Comparison: mean(Homogenisation) - mean(Differentiation)

#Direction x Main land use type
avg_comparisons(best_model_species,  variables= c("direction"), by = c("main_land_use_type"),
                newdata = datagrid(newdata=best_model_species$data,
                                   direction = unique,
                                   ecosystem_type = unique,
                                   main_land_use_type=unique), 
                re_formula=NA,
                conf_level=.8) 

#Species 
 #        Group main_land_use_type Estimate   10.0 % 90.0 %
 # hetmagnitude        Multiple      0.1793  0.01170 0.3358
 # hfpmagnitude        Multiple      0.2656  0.10986 0.4072
 # hetmagnitude        Agriculture   0.1483 -0.07092 0.3574
 # hfpmagnitude        Agriculture  -0.0665 -0.20247 0.0592
 # hetmagnitude        Forest        0.1723 -0.00219 0.3527
 # hfpmagnitude        Forest        0.0411 -0.05321 0.1360
 # hetmagnitude        Urban         0.1214 -0.01852 0.3278
 # hfpmagnitude        Urban         0.4183 -0.00681 0.9539

avg_comparisons(best_model_traits, variables= c("direction"), by = c("main_land_use_type"),
                newdata = datagrid(newdata=best_model_traits$data,
                                   direction = unique,
                                   ecosystem_type = unique,
                                   main_land_use_type=unique),
                re_formula=NA,
                conf_level=0.8) 

#Traits
 #        Group main_land_use_type Estimate  10.0 % 90.0 %
 # hetmagnitude        Multiple      0.5613  0.3156  0.816
 # hfpmagnitude        Multiple      0.1868  0.0977  0.267
 # hetmagnitude        Agriculture   1.0525  0.3648  1.742
 # hfpmagnitude        Agriculture   0.3022  0.1228  0.480
 # hetmagnitude        Forest        0.7635  0.3305  1.142
 # hfpmagnitude        Forest        0.1818  0.0874  0.279
 # hetmagnitude        Urban         0.2330 -0.0824  0.700
 # hfpmagnitude        Urban         0.0602 -0.1018  0.252

#'---------------------------------------------------------------------------------------------------------------
# @Main text figures ===== #################################################################################
#'---------------------------------------------------------------------------------------------------------------

direction_colors = c("Increasing dissimilarity" = "#ff7f32","Increasing similarity" = "#9f5cc0")  

contrasts_magnitude<-predictions(best_model_species,  by = c("ecosystem_type","direction"),
                                 newdata = datagrid(newdata=best_model_species$data,
                                                    direction = unique,
                                                    ecosystem_type = unique,
                                                    main_land_use_type=unique), re_formula=NA) |> 
  posterior_draws() |> 
  mutate(facet="Taxonomic turnover") |> 
  bind_rows(
    predictions(best_model_traits,  by = c("ecosystem_type","direction"),
                newdata = datagrid(newdata=best_model_traits$data,
                                   direction = unique,
                                   ecosystem_type = unique,
                                   main_land_use_type=unique), re_formula=NA) |> 
      posterior_draws() |> mutate(facet="Functional turnover")) |> 
  mutate(class = "Ecosystem type") |> 
  rename(Predictor = "ecosystem_type") |> 
  select(draw,group,facet,class,Predictor,direction) |> 
  bind_rows(
    predictions(best_model_species,  by = c("main_land_use_type","direction"),
                newdata = datagrid(newdata=best_model_species$data,
                                   direction = unique,
                                   ecosystem_type = unique,
                                   main_land_use_type=unique), re_formula=NA) |> 
      posterior_draws() |> 
      mutate(facet="Taxonomic turnover") |> 
      bind_rows(predictions(best_model_traits,  by = c("main_land_use_type","direction"),
                            newdata = datagrid(newdata=best_model_traits$data,
                                               direction = unique,
                                               ecosystem_type = unique,
                                               main_land_use_type=unique), re_formula=NA) |> 
                  posterior_draws() |> mutate(facet="Functional turnover")) |> 
      mutate(class = "Main land use type") |>
      rename(Predictor = "main_land_use_type") |> 
      select(draw,group,facet,class,Predictor,direction)
  )

Fig_3<-
  contrasts_magnitude |> 
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
  mutate(
    direction = recode(
      direction,
      "Differentiation" = "Increasing dissimilarity",
      "Homogenisation"  = "Increasing similarity"
    ),
    direction = factor(
      direction,
      levels = c("Increasing dissimilarity", "Increasing similarity")
    )
  ) |> 
  mutate(
    group = case_when(
      str_detect(group, "^hfpmagnitude") ~ "Human Footprint",
      str_detect(group, "^hetmagnitude") ~ "Habitat Heterogeneity"
    )) |> 
  mutate(group = factor(group,levels=c("Human Footprint","Habitat Heterogeneity"))) |> 
  mutate(facet = factor(facet,levels=c("Taxonomic turnover","Functional turnover"))) |> 
  ggplot(aes(y=Predictor,x=draw,colour=direction))+
    ggstats::geom_stripped_rows(
      aes(y = Predictor),
      colour=NA,
      odd  = "white",
      even = "grey90",
      xfrom = -Inf, xto = Inf
    ) +
  stat_pointinterval(position="dodge", .width=c(.8,.89,.95))+
  ggh4x::facet_nested(class~facet+group, scales="free_y",space="free_y")+
  scale_colour_manual(values=direction_colors)+
  theme_void(base_family = "sans", base_size = 14) +
  theme(
    panel.grid.minor   = element_blank(),
    plot.background    = element_blank(),
    panel.background =  element_blank(),
    strip.text.x =  element_text(face="bold",margin=margin(t=5)),
    strip.text.y = element_text(face="bold",angle=270, margin=margin(l=5)),
    axis.title.x       = ggtext::element_markdown(face = "bold", size = 12, margin=margin(t=5)),
    axis.text.y        = element_text(size = 12, hjust = 1),
    axis.text.x        = ggtext::element_markdown(size = 10, margin=margin(t=3)),
    axis.title.y       = element_blank(),
    axis.line.x        = element_line(color = "grey20", linewidth = 0.1),
    axis.ticks.x       = element_line(color = "grey20", linewidth = 0.1),
    axis.ticks.length  = unit(0.2, "lines"),
    legend.position    = "bottom",
    legend.title = element_text(face = "bold", size = 11),
    legend.text  = element_text(size = 10),
    legend.key.width  = unit(1.2, "lines"),
    legend.key.height = unit(0.8, "lines"),
    legend.spacing.x  = unit(0.5, "lines"),
    legend.spacing.y  = unit(0.2, "lines"),
    legend.margin     = margin(t = 5, b = 3),
    legend.box.margin = margin(t = -5),
    plot.tag           = element_text(face = "bold", size = 14),
    panel.spacing.x = unit(2,"lines"),
    plot.margin = margin(0,0,0,0)
  ) +
  guides(
    x = guide_axis(cap = TRUE)
  ) + labs(x = "<b>Cumulative fitted turnover</b>",
           colour="Turnover direction",
           )+
  scale_x_continuous(limits=c(0,NA),
                     labels = scales::label_number(accuracy = 0.01)
  ) +
    coord_cartesian(xlim=c(0,3), clip="off")



ggsave(Fig_3, 
       filename="S7_Model_outputs_figures_and_tables/main_figures/Figure_3.pdf", 
       device=cairo_pdf, 
       width=10,
       height=5,
       units="in")

library(bayesplot)
library(dplyr)
library(stringr)
library(ggplot2)
library(forcats)


# 1. extract and filter
posterior_plot_data <- generate_figure_data(best_model_species,best_model_traits, type="magnitude",draw_min=-5,draw_max=5) |> 
  mutate(
    response = case_when(
      str_detect(response, "^hfpmagnitude") ~ "Human Footprint",
      str_detect(response, "^hetmagnitude") ~ "Habitat Heterogeneity"
    )) |> 
  mutate(response = factor(response,levels=c("Human Footprint","Habitat Heterogeneity")))

Extended_data_Fig4 <- ggplot(posterior_plot_data, aes(x = posterior, y = parameter, colour = significant)) +
stat_pointinterval(.width=c(.8,.89,.95)) +
  geom_text(
    data = ~ .x |>  distinct(across(any_of(c("response","shape", "parameter", "facet"))), .keep_all = TRUE), # Add text annotations from get_support
    aes(y = parameter, x = 4, label = text),
    hjust = 0, size = convert_size(8)
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_grid(response ~ facet, space = "free_y") +
  labs(x = "Posterior draws", y = NULL) +
  theme_void(base_family = "sans", base_size = 12) +
  theme(
    panel.grid.minor   = element_blank(),
    plot.background    = element_blank(),
    panel.background =  element_blank(),
    strip.text.y = element_text(face="bold",angle=90,size=12),
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
  ) +
  coord_cartesian(xlim=c(-5,10), clip="off")+
  scale_colour_manual(values=rev(c("black","gray60")))

ggsave(Extended_data_Fig4, filename="S7_Model_outputs_figures_and_tables/extended_data/Extended_Figure_4.pdf", device=cairo_pdf, width=10.5,height=7.5,units="in")

