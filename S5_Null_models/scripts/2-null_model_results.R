# ==============================================================================
# Script: 2-Null_model_results.R
# Author: Caio Graco-Roza
# Date: 2024-11-24
# Description:
#   This script compiles and analyzes the results of null model simulations and 
#   compares them to observed beta diversity models. It calculates standardized 
#   effect sizes (SES) for direction, magnitude, and shape metrics across datasets.
#
# Workflow:
# 1. Load null model and observed model results from their respective directories.
# 2. Perform SES calculations using `estimate_SES_direction`, `estimate_SES_magnitude`, 
#    and `estimate_SES_shape` functions.
# 3. Combine SES results across datasets into a unified format.
# 4. Save final SES results in an Excel file for further analysis.
#
# Key Steps:
# - Extract filenames from null and observed model folders.
# - Process each pair of null and observed model files iteratively.
# - Use SES estimation functions to evaluate model performance.
# - Combine results into a single structured output.
#
# Usage:
# - Designed for high-performance computing environments but can also run locally.
# - Ensure that the file paths in `null_model_folder` and `observed_folder` are accurate.
# - For local execution, simply adjust the file paths and run the script.
#
# Dependencies:
# - Libraries: parallel, doSNOW, tidyverse, xlsx.
# - Functions: `estimate_SES_direction`, `estimate_SES_magnitude`, and `estimate_SES_shape`.
# - Input Files: Null model results (`*_null_bbgdm.rds`) and observed model results (`*_bbgdm.rds`).
#
# Notes:
# - Ensure the `functions_run_null_bbgdm.R` script is correctly sourced.
# - Output is saved as an Excel file in the specified location.
# - If running in parallel, uncomment and configure the cluster setup section.
# ==============================================================================
library(parallel)
library(doSNOW)

# Specify the folders and file extensions
null_model_folder <- "S5_Null_models/null_output"
observed_folder <- "S4_run_BBGDM/bbgdm_output"

# Get a list of all null model files
null_model_files <- list.files(null_model_folder, pattern = "_null_bbgdm.rds", full.names = FALSE)

# numCores<-1
# cl <- makeCluster(numCores)
# registerDoSNOW(cl)
# 
# # progress bar ------------------------------------------------------------
# library(progress)
iterations <- length(null_model_files)         # used for the foreach loop  
# 
# pb <- progress_bar$new(
#   format = ":letter [:bar] :elapsed | eta: :eta",
#   total = iterations,    # 100 
#   width = 100)
# 
# progress_letter <- seq_along(null_model_files)  # token reported in progress bar
# 
# # allowing progress bar to be used in foreach -----------------------------
# progress <- function(n){
#   pb$tick(tokens = list(letter = progress_letter[n]))
# } 
# 
# opts <- list(progress = progress)
# 
# # foreach loop ------------------------------------------------------------
# library(foreach)

# output<- foreach(i = 1:iterations,
#                  .combine = rbind,
#                 # .options.snow = opts,
#                  .multicombine = FALSE,
#                  .errorhandling = 'pass') %do% {
output <- list() 
for (i in 1:iterations){
  print(i)
  i<-1
                   require(tidyverse)
                   source("S5_Null_models/functions/functions_run_null_bbgdm.R")

                   # Extract the filename without the extension
                   filename <- gsub("_null_bbgdm","",tools::file_path_sans_ext(null_model_files[i]))
                   # Generate the corresponding observed file name
                   null_model_file<-file.path(null_model_folder, paste0(filename,"_null_bbgdm.rds"))
                   observed_file <- file.path(observed_folder, paste0(filename,"_bbgdm.rds"))
                   
                   # Process the null model file and observed file accordingly
                   # Replace the placeholder code below with your desired operations
                   # For example, you can load the files using `readRDS()` and perform analyses
                   
                   # Load the null model file
                   null_model <- readRDS(null_model_file)
                   
                   # Load the observed file
                   observed_model <- readRDS(observed_file)
                   
                   # Perform desired operations with the files
                   
                   names(null_model) <- paste0("estimated_",rep(1:999))
                   
                   # # Simulated Scenario to Compare Approaches
                   # 
                   # # Set the random seed for reproducibility
                   # set.seed(123)
                   # 
                   # # Simulated data
                   # n_total <- 1000
                   # n1 <- 200
                   # n2 <- 800
                   # 
                   # # Generate random values for differentiation and homogenisation
                   # differentiation <- c(runif(n1, min = 0.5, max = 0.8), rep(NA, n_total - n1))
                   # homogenisation <- c(runif(n2, min = 0.1, max = 0.4), rep(NA, n_total - n2))
                   # 
                   # # Calculate the mean values using the weighted approach
                   # mean_weighted_diff <- (n1 / (n1 + n2)) * mean(differentiation, na.rm = TRUE) - (n2 / (n1 + n2)) * mean(homogenisation, na.rm = TRUE)
                   # 
                   # # Calculate the standard errors
                   # se_weighted_diff <- sqrt((n1^2 / ((n1 + n2)^2 * (n1 + n2 - 1))) * var(differentiation, na.rm = TRUE) +
                   #                            (n2^2 / ((n1 + n2)^2 * (n1 + n2 - 1))) * var(homogenisation, na.rm = TRUE))
                   # 
                   # # Calculate the confidence interval using a desired confidence level (e.g., 95%)
                   # confidence_level <- 0.95
                   # z_value <- qnorm((1 + confidence_level) / 2)
                   # ci_lower <- mean_weighted_diff - z_value * se_weighted_diff
                   # ci_upper <- mean_weighted_diff + z_value * se_weighted_diff
                   # 
                   # # Print the results with error estimation
                   # cat("Mean Difference (Weighted Approach):", mean_weighted_diff, "\n")
                   # cat("Standard Error (Weighted Approach):", se_weighted_diff, "\n")
                   # cat("Confidence Interval (Weighted Approach):", ci_lower, "-", ci_upper, "\n")
                   # 
                   # # Explanation of the Example:
                   # # This updated code includes the calculation of the standard error and confidence interval for the weighted difference value.
                   # # The standard error quantifies the variability or uncertainty around the mean difference estimate.
                   # # The confidence interval provides a range of plausible values for the true mean difference, based on the desired confidence level.
                   # # By incorporating the error estimation, you can better understand the precision of the difference estimate and the level of uncertainty associated with it.
                   
                   SES_direction <- estimate_SES_direction(null_model = null_model,observed_model = observed_model,param=TRUE) |> mutate(dataset=filename)
                   SES_magnitude <- estimate_SES_magnitude(null_model = null_model,observed_model = observed_model, param=TRUE) |> mutate(dataset=filename)
                   SES_shape <- estimate_SES_shape(null_model = null_model,observed_model = observed_model,param=TRUE) |> mutate(dataset=filename)
                   
                   out<- Reduce(full_join,list(SES_direction,SES_magnitude,SES_shape)) |>  relocate(dataset,direction)
                   output[[i]]<- SES_direction
                   
                 }

#stopCluster(cl) 

#Save final results.
output |> 
  bind_rows() |> 
  xlsx::write.xlsx("S6_Synthesis_model/BBGDM_SES.xlsx")



