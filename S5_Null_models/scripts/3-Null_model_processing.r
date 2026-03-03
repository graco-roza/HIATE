# ==============================================================================
# Script: 3-Null_model_processing.R (HPC version)
# Author: Caio Graco-Roza
# Description:
#   This script processes ONE dataset per HPC array job.
#   Each array job gets its dataset index from SLURM_ARRAY_TASK_ID.
# ==============================================================================

library(tidyverse)
source("S5_Null_models/functions/functions_run_null_bbgdm.R")
source("S5_Null_models/functions/functions_prepare_null_results.R")
# ---- Paths ----
null_model_folder <- "S5_Null_models/null_output/Podani_abun"
observed_folder   <- "S4_run_BBGDM/bbgdm_output/Podani_abun"
out_folder        <- "S5_Null_models/post_processed/Podani_abun"
metadata<-readxl::read_excel("S6_Synthesis_model/data/synthesis_data_complete.xlsx") |> 
  filter(facet=="Functional",beta_type == "Podani",metric_type == "abun", ) |> 
  mutate(across(hfp_q_25:het_n_Uncertain, ~ifelse(is.na(.x),0,.x))) 

# Make sure output folder exists
if (!dir.exists(out_folder)) dir.create(out_folder, recursive = TRUE)

# ---- Get dataset index from SLURM ----
task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))  # default to 1 if not set
null_model_files <- list.files(null_model_folder, pattern = "_null_bbgdm.rds", full.names = FALSE)

# Safety check: bail out if index exceeds available files
if (task_id > length(null_model_files)) {
  stop(glue::glue("Task ID {task_id} exceeds number of files {length(null_model_files)}"))
}

# ---- Pick dataset ----
filename <- gsub("_null_bbgdm","",tools::file_path_sans_ext(null_model_files[task_id]))
null_model_file <- file.path(null_model_folder, paste0(filename,"_null_bbgdm.rds"))
observed_file   <- file.path(observed_folder, paste0(filename,"_bbgdm.rds"))

cat(glue::glue("Processing dataset {task_id}: {filename}\n"))

tryCatch({
  # Load models
  null_model <- readRDS(null_model_file)
  observed_model <- readRDS(observed_file)
  
  # Ensure null_model list has consistent names
  names(null_model) <- paste0("estimated_", seq_along(null_model))
  
  # Run SES analyses
  SES_direction <- estimate_SES_direction(null_model = null_model, observed_model = observed_model, param = FALSE) %>%
    mutate(dataset = filename)
  SES_magnitude <- estimate_SES_magnitude(null_model = null_model, observed_model = observed_model, param = FALSE) %>%
    mutate(dataset = filename)
  SES_shape     <- estimate_SES_shape(null_model = null_model, observed_model = observed_model, param = FALSE) %>%
    mutate(dataset = filename)
  
  # Combine results
  out <- list(
    SES_direction = SES_direction,
    SES_magnitude = SES_magnitude,
    SES_shape     = SES_shape
  )
  out <- Reduce(full_join,out)
  # Save individual dataset output
  saveRDS(out, file = file.path(out_folder, paste0(filename, "_SES_results.rds")))
  
}, error = function(e) {
  message(glue::glue("Dataset {task_id} - {filename} FAILED: {conditionMessage(e)}"))
})
