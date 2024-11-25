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
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
# Coerce the value to an integer
ii <- as.numeric(slurm_arrayid) # e.g., ii <- 1 for local testing

# -----------------------------------------------------------------------------------------------------------------
# @ Run Analysis ######
# -----------------------------------------------------------------------------------------------------------------

# Source custom functions for predictor selection
source("S3_get_best_predictors/Functions/functions_get_best_predictors.R")

# Get dataset names from beta diversity output
files <- tools::file_path_sans_ext(list.files("S2_get_beta_diversity/betadiv_output"))
focal_dataset <- gsub("_beta_Output", "", files[ii]) # Choose one dataset
write(focal_dataset, stderr())

# Read beta diversity metrics and predictor variables
dissim_data <- glue::glue("S2_get_beta_diversity/betadiv_output/{focal_dataset}_beta_Output.rds") %>% read_rds()
predictors <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds() %>% pluck("predictors")

# Define predictor groups
hfp_var <- predictors %>% dplyr::select(contains("hfp")) %>% names()
modis_var <- predictors %>% dplyr::select(contains("modis")) %>% names()
het_var <- predictors %>% dplyr::select(contains("het")) %>% names()
comb <- t(expand.grid(hfp_var, modis_var, het_var)) %>% data.frame
metrics <- dissim_data %>% map(~ .x %>% pluck("Brepl") %>% round(3))

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

best_coeff_tax <- lapply(var_select_tax, function(x) x$coefficient %>% bind_rows(.id = "type")) %>%
  bind_rows(.id = "group") %>%
  as_tibble() %>%
  filter(!var %in% c("Temp", "Prec", "Geographic")) %>%
  dplyr::mutate(across(c(group, var, type), as.factor)) %>%
  group_by(group, type, var) %>%
  dplyr::summarise(total_effect = sum(effect_size)) %>%
  ungroup() %>%
  slice_max(total_effect)

best_r2_tax <- sapply(var_select_tax, function(x) unlist(x$r2)) %>%
  data.frame() %>%
  set_names(paste0(1:48)) %>%
  rownames_to_column("type") %>%
  pivot_longer(cols = !type, names_to = "group") %>%
  dplyr::mutate(across(c(group, type), as.factor)) %>%
  slice_max(value)

# Prepare predictors (best ones from variable selection)
pred_tax_coeff <- predictors %>% dplyr::select(site, x, y,
                                               all_of(comb %>% data.frame() %>% pull(as.numeric(best_coeff_tax$group)[1])),
                                               Temp, Prec
)

pred_tax_r2 <- predictors %>% dplyr::select(site, x, y,
                                            all_of(comb %>% data.frame() %>% pull(as.numeric(best_r2_tax$group)[1])),
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

best_coeff_fun <- lapply(var_select_fun, function(x) x$coefficient %>% bind_rows(.id = "type")) %>%
  bind_rows(.id = "group") %>%
  as_tibble() %>%
  dplyr::mutate(across(c(group, var, type), as.factor)) %>%
  group_by(group, type) %>%
  dplyr::summarise(total_effect = sum(effect_size)) %>%
  ungroup() %>%
  slice_max(total_effect)

best_r2_fun <- sapply(var_select_fun, function(x) unlist(x$r2)) %>%
  data.frame() %>%
  set_names(paste0(1:48)) %>%
  rownames_to_column("type") %>%
  pivot_longer(cols = !type, names_to = "group") %>%
  dplyr::mutate(across(c(group, type), as.factor)) %>%
  slice_max(value)

# Prepare predictors (best ones from variable selection)
pred_fun_coeff <- predictors %>% dplyr::select(site, x, y,
                                               all_of(comb %>% data.frame() %>% pull(as.numeric(best_coeff_fun$group)[1])),
                                               Temp, Prec
)

pred_fun_r2 <- predictors %>% dplyr::select(site, x, y,
                                            all_of(comb %>% data.frame() %>% pull(as.numeric(best_r2_fun$group)[1])),
                                            Temp, Prec
)

# Save Results ------------------------------------------------------------------------
pred_coeff <- list(Taxonomic = pred_tax_coeff, Functional = pred_fun_coeff)
pred_r2 <- list(Taxonomic = pred_tax_r2, Functional = pred_fun_r2)

saveRDS(pred_coeff, glue::glue("S3_get_best_predictors/best_predictors/coeff/{focal_dataset}_pred.rds"))
saveRDS(pred_r2, glue::glue("S3_get_best_predictors/best_predictors/r2/{focal_dataset}_pred.rds"))

model_direction_coeff <- list(Taxonomic = data.frame(dataset = focal_dataset, type = best_coeff_tax$type[1]),
                              Functional = data.frame(dataset = focal_dataset, type = best_coeff_fun$type[1])
)
model_direction_r2 <- list(Taxonomic = data.frame(dataset = focal_dataset, type = best_r2_tax$type[1]),
                           Functional = data.frame(dataset = focal_dataset, type = best_r2_fun$type[1])
)

saveRDS(model_direction_coeff, glue::glue("S3_get_best_predictors/relationship_direction/coeff/{focal_dataset}_dir.rds"))
saveRDS(model_direction_r2, glue::glue("S3_get_best_predictors/relationship_direction/r2/{focal_dataset}_dir.rds"))
