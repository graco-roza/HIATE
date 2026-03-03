# ==============================================================================
# Script: 1-run_null_analysis.R
# Author: Caio Graco-Roza
# Date: 2024-11-24
# Description: 
#   This script estimates taxonomic and functional beta diversity by applying 
#   abundance weights to species. It leverages null model analyses using 
#   Bayesian Bootstrap Generalized Dissimilarity Modeling (BBGDM).
#
# Workflow:
# 1. Load and preprocess community, trait, and predictor data for a selected dataset.
# 2. Execute null model analyses to calculate beta diversity components:
#    - Total beta diversity
#    - Replacement beta diversity
#    - Richness beta diversity
# 3. Save the results for further statistical evaluation and model comparison.
#
# Usage:
# - Designed to run in high-performance computing (HPC) environments using SLURM arrays.
# - For local execution, manually set `ii <- 1` to process a single dataset.
#
# Key Steps:
# - Load community matrix and synthetic trait data from S2.
# - Incorporate best predictors from S3 for BBGDM.
# - Run `run_null()` function (defined in `functions_run_null_bbgdm.R`) 
#   iteratively for multiple randomizations.
# - Save outputs in `S5_Null_models/null_output/coeff/`.
#
# Dependencies:
# - Libraries: BAT, hypervolume, tidyverse, magrittr, glue, readxl, future, doSNOW.
# - Functions: `functions_run_null_bbgdm.R` (source script).
#
# Notes:
# - Ensure community data matches available trait data to avoid missing species.
# - Adjust randomization iterations (`replicate(1000, ...)`) as needed based on 
#   computational resources.
# - This script can be executed locally or via SLURM for batch processing.
# ==============================================================================

# Packages .............................................................................................................
if(!require("pacman")) {install.packages("pacman")}

pacman::p_load(
  BAT #make taxonomic beta diversity
  ,hypervolume #make functional beta diversity
  ,tidyverse #manage data
  ,magrittr #use set names
  ,glue #create strings
  ,readxl #read rds files
  ,future #set up cores available
  ,doSNOW # make parallel procedure
  ,doRNG
)

source("S5_run_BBGDM/functions_run_null_bbgdm.R")


#' This script is designed to:
#' - estimate taxonomic and functional beta diversity by applying abundance weights to species.

#Only for HPC cluster (Doesn't affect anything in local machines) ---------------------------------------------------
# grab the array id value from the environment variable passed from sbatch
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
# coerce the value to an integer
ii <- as.numeric(slurm_arrayid) # ii <- 1
#' =================================================================================================================

#' -----------------------------------------------------------------------------------------------------------------
# @ run analysis ######
#' -----------------------------------------------------------------------------------------------------------------
#' Reference: 

files <- tools::file_path_sans_ext(list.files("S2_get_beta_diversity/betadiv_input")) #get a vector with dataset names

focal_dataset <- gsub("_beta_Input","",files[ii]) #chose one dataset 
write(focal_dataset, stderr())

species_diff <- glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds") %>%  read_rds() #read the dataset pre processed file

comm <- species_diff %>%  pluck("comm") #extract the community data
traits <- species_diff %>%  pluck("trait_syndrome") #extract the synthetic traits (PCoA Axes)
predictors <- lapply(c("Podani_abun","Podani_pa","Baselga_abun","Baselga_pa"),function(x) readRDS(glue::glue("S3_get_best_predictors/best_predictors/{x}/{focal_dataset}_pred.rds"))$Functional) #read the dataset pre processed file
names(predictors)<-c("Podani_abun","Podani_pa","Baselga_abun","Baselga_pa")

#make sure the community only has species with trait values. 
comm <- comm %>%  select(any_of(rownames(traits)))
res <- transpose(replicate(1000, run_null(comm = comm, traits = traits[sample(1:nrow(traits)),], pred = predictors), simplify=FALSE))

Podani_abun_null <- res$Podani_abun
Podani_pa_null <- res$Podani_pa
Baselga_abun_null <- res$Baselga_abun
Baselga_pa_null <- res$Baselga_pa

saveRDS(Podani_abun_null, glue::glue("S5_Null_models/null_output/Podani_abun/{focal_dataset}_null_bbgdm.rds"))
saveRDS(Podani_pa_null, glue::glue("S5_Null_models/null_output/Podani_pa/{focal_dataset}_null_bbgdm.rds"))
saveRDS(Baselga_abun_null, glue::glue("S5_Null_models/null_output/Baselga_abun/{focal_dataset}_null_bbgdm.rds"))
saveRDS(Baselga_pa_null, glue::glue("S5_Null_models/null_output/Baselga_pa/{focal_dataset}_null_bbgdm.rds"))