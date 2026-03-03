#' -----------------------------------------------------------------------------------------------------------------
#' Script: Get Best Predictors
#' Author: Caio Graco-Roza
#' Last Updated: 2024-11-20
#' -----------------------------------------------------------------------------------------------------------------
#' Description:
#' This script is designed to identify the best environmental predictors for taxonomic and functional beta diversity.
#' - Reads pre-computed beta diversity metrics and predictor variables.
#' - Explores all combinations of predictors to find the best based on coefficient size and R² values.
#' - Outputs the selected predictors and their relationship directions.
#' -----------------------------------------------------------------------------------------------------------------
#' Input:
#' - Beta diversity metrics from `S2_get_beta_diversity/betadiv_output/`.
#' - Predictor variables from `S1_Preprocessing/Processed/`.
#' 
#' Output:
#' - Best predictors based on coefficients and R² saved in `S3_get_best_predictors/best_predictors/`.
#' - Relationship directions for taxonomic and functional beta diversity in `S3_get_best_predictors/relationship_direction/`.
#' 
#' Notes:
#' - This script is designed to run on an HPC cluster with an array job configuration.
#' - To run locally, set `ii <- 1` to process the first dataset manually.
#' -----------------------------------------------------------------------------------------------------------------

# Packages .............................................................................................................
if(!require("pacman")) {install.packages("pacman")}
pacman::p_load("tidyverse", "gdm", "surveillance", "plotrix", "glue", "magrittr")

# Only for HPC cluster (Doesn't affect anything in local machines) ---------------------------------------------------
# Grab the array ID value from the environment variable passed from sbatch
#Dont change this because it is mostly needed for the HPC
beta_types <- c("Podani_pa","Podani_abun", "Baselga_abun", "Baselga_pa")
#combined_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))  # Combined ID from 1 to 660
ii <- ifelse(combined_id > 165, ((combined_id-1) %% 165) + 1,combined_id)        # ii cycles from 1 to 165
beta_index <- 1#((combined_id - 1) %/% 165) + 1  # beta_index cycles 1,2,3,4
beta_type <- beta_types[beta_index]
# -----------------------------------------------------------------------------------------------------------------
# @ Run Analysis ######
# -----------------------------------------------------------------------------------------------------------------

# Source custom functions for predictor selection
source("S3_get_best_predictors/Functions/functions_get_best_predictors.R")

# Get dataset names from beta diversity output
files <- tools::file_path_sans_ext(list.files(glue::glue("S2_get_beta_diversity/betadiv_output/{beta_type}")))
focal_dataset <- "Predicts_vk1_2011_edenius_1"#gsub("_beta_Output", "", files[ii]) # Choose one dataset
write(focal_dataset, stderr())

# Read beta diversity metrics and predictor variables
dissim_data <- glue::glue("S2_get_beta_diversity/betadiv_output/{beta_type}/{focal_dataset}_beta_Output.rds") %>% read_rds()
predictors <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds() %>% pluck("predictors")

# Define predictor groups
hfp_var <- predictors %>% dplyr::select(contains("hfp")) %>% names()
modis_var <- predictors %>% dplyr::select(contains("modis")) %>% names()
het_var <- predictors %>% dplyr::select(contains("het")) %>% names()
comb <- t(expand.grid(hfp_var, modis_var, het_var)) %>% data.frame
metrics <- dissim_data %>% map(~ .x %>% pluck(2) %>% round(3))

# Variable Selection for Taxonomic Beta Diversity ----------------------------------------
print("Taxonomic best variables...")
var_select_tax <- lapply(1:ncol(comb), function(j) {
  jj <- comb %>% pull(j)
  res <- get_bbgdm_predictors(
    metrics$Taxonomic,
    pred = predictors %>% dplyr::select(site, x, y, all_of(jj), Temp, Prec)
  )
  return(res)
})

best_coeff_tax_group <- var_select_tax %>%
  imap_dfr(~ bind_rows(
    as_tibble(.x$coefficient$differentiation) %>% 
      mutate(direction = "differentiation",
             r2 = .x$r2$differentiation,
             group = .y),
    as_tibble(.x$coefficient$homogenization) %>% 
      mutate(direction = "homogenization",
             r2 = .x$r2$homogenization,
             group = .y)
  ) %>%
    rename(predictors = var,
           coefficients = effect_size) %>%
    select(group, direction, predictors, coefficients, r2)
  ) |> 
  group_by(group, direction) %>%
  summarise(
    r2 = unique(r2),
    sum_buffer_coef = sum(coefficients[str_detect(predictors, "_")]),
    .groups = "drop"
  ) |>
  slice_max(r2)  |> 
  slice_max(sum_buffer_coef) |> 
  slice_min(group) |> 
  pull(group)

# Prepare predictors (best ones from variable selection)
pred_tax_best <- predictors %>% dplyr::select(site, x, y,
                                               all_of(comb %>% data.frame() %>% pull(best_coeff_tax_group)),
                                               Temp, Prec
)

# Variable Selection for Functional Beta Diversity --------------------------------------
print("Functional best variables...")
var_select_fun <- lapply(1:ncol(comb), function(j) {
  jj <- comb %>% data.frame() %>% pull(j)
  res <- get_bbgdm_predictors(
    metrics$Functional,
    pred = predictors %>% dplyr::select(site, x, y, all_of(jj), Temp, Prec)
  )
  return(res)
})

best_coeff_fun_group <- var_select_fun %>%
  imap_dfr(~ bind_rows(
    as_tibble(.x$coefficient$differentiation) %>% 
      mutate(direction = "differentiation",
             r2 = .x$r2$differentiation,
             group = .y),
    as_tibble(.x$coefficient$homogenization) %>% 
      mutate(direction = "homogenization",
             r2 = .x$r2$homogenization,
             group = .y)
  ) 
  ) |> 
  rename(predictors = var,
         coefficients = effect_size) %>%
  select(group, direction, predictors, coefficients, r2) |> 
  group_by(group, direction) %>%
  summarise(
    r2 = unique(r2),
    sum_buffer_coef = sum(coefficients[str_detect(predictors, "_")]),
    .groups = "drop"
  ) |>
  slice_max(r2)  |> 
  slice_max(sum_buffer_coef) |> 
  slice_min(group) |> 
  pull(group)

# Prepare predictors (best ones from variable selection)
  pred_fun_best <- predictors %>% dplyr::select(site, x, y,
                                                all_of(comb %>% data.frame() %>% pull(best_coeff_fun_group)),
                                                Temp, Prec
  )
  

pred_best <- list(Taxonomic = pred_tax_best, Functional = pred_fun_best)


saveRDS(pred_best, glue::glue("S3_get_best_predictors/best_predictors/{beta_type}/{focal_dataset}_pred.rds"))