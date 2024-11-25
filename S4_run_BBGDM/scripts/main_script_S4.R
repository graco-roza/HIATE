#' -----------------------------------------------------------------------------------------------------------------
#' Script: Run BBGDM Analysis
#' Author: Caio Graco-Roza
#' Last Updated: 2024-11-24
#' -----------------------------------------------------------------------------------------------------------------
#' Description:
#' This script performs Bayesian Bootstrap Generalized Dissimilarity Model (BBGDM) analysis for a selected dataset.
#' It integrates beta diversity dissimilarity metrics and environmental predictors to generate differentiation and 
#' homogenization models. The script is designed to run on HPC clusters and can also be executed locally.
#' -----------------------------------------------------------------------------------------------------------------
#' Input:
#' - Beta diversity dissimilarity metrics from `S2_get_beta_diversity/betadiv_output/`.
#' - Environmental predictor coefficients from `S3_get_best_predictors/best_predictors/coeff/`.
#' 
#' Output:
#' - BBGDM model results saved in `S4_run_BBGDM/bbgdm_output/` as RDS files.
#' 
#' Notes:
#' - The script uses SLURM array job variables (`SLURM_ARRAY_TASK_ID`) for parallel execution on HPC clusters.
#' - Local execution is also supported by manually setting the variable `ii`.
#' -----------------------------------------------------------------------------------------------------------------

# Packages .............................................................................................................
if(!require("pacman")) {install.packages("pacman")}
pacman::p_load("tidyverse","gdm", "surveillance","plotrix", "glue", "magrittr","bbgdm")


#Only for HPC cluster (Doesn't affect anything in local machines) ---------------------------------------------------
# grab the array id value from the environment variable passed from sbatch
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
# coerce the value to an integer
ii <- as.numeric(slurm_arrayid) # ii <- 1
#' =================================================================================================================

#' -----------------------------------------------------------------------------------------------------------------
# @ run analysis ######
#' -----------------------------------------------------------------------------------------------------------------
files <- tools::file_path_sans_ext(list.files("S3_get_best_predictors/best_predictors/r2")) #get a vector with dataset names
focal_dataset <- gsub("_pred","",files[ii]) #chose one dataset 
write(focal_dataset,stderr())

source("S4_run_BBGDM/functions/helper_functions_S4.R")

#read dissimilarity data
dissim_data <- glue::glue("S2_get_beta_diversity/betadiv_output/{focal_dataset}_beta_Output.rds") %>%  read_rds() #read the dataset beta diversity file
#read predictor data
predictors <- glue::glue("S3_get_best_predictors/best_predictors/coeff/{focal_dataset}_pred.rds") %>%  read_rds() #read the dataset pre processed file

#extract replacement portion of beta diversity
metrics <- dissim_data %>% map(~ .x %>%  pluck("Brepl") %>% round(3))

tax_bbgdm <- run_bbgdm(metrics$Taxonomic,predictors$Taxonomic)
fun_bbgdm <- run_bbgdm(metrics$Functional,predictors$Functional)

results_bbgdm <-list(Taxonomic = tax_bbgdm, Functional = fun_bbgdm)

saveRDS(results_bbgdm, glue::glue("S4_run_BBGDM/bbgdm_output/{focal_dataset}_bbgdm.rds"))

