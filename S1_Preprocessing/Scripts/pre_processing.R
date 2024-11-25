###############################################################################
# SCRIPT NAME: Pre-processing Community and Trait Data
#
# DESCRIPTION:
#   This script handles the data cleaning and transformation of community and
#   trait datasets. It is designed to prepare the data for downstream analyses,
#   including the estimation of species dissimilarities based on selected 
#   functional traits. The pre-processing steps include:
#     - Cleaning and filtering community and trait datasets
#     - Estimating predictors (e.g., human footprint, climate variables)
#     - Saving cleaned datasets for further analysis
#
# USAGE:
#   This script is executed as part of a batch process using SLURM on an HPC.
#   It can also be run manually for specific datasets.
#
# INPUTS:
#   - Raw data files: located in `S1_Preprocessing/Raw_data`
#   - Metadata: `S1_Preprocessing/Miscellaneous/dataset_info_all.xlsx`
#   - Auxiliary functions: `S1_Preprocessing/Functions/auxiliary_functions.R`
#
# OUTPUTS:
#   - Cleaned datasets: saved in `S1_Preprocessing/pre_processed`
#
# AUTHOR: Caio Graco-Roza
# DATE CREATED: 2021-05-14
# LAST UPDATED: 2024-11-20
#
# NOTES:
#   - Ensure all required packages are installed before execution.
#   - This script leverages SLURM to process datasets in parallel.
#   - Key cleaning steps:
#     1. Retain sites with coordinates.
#     2. Retain species with trait data.
#     3. Remove sites with fewer than 2 species (needed for hypervolume estimation).
#     4. Remove species with zero abundance across all sites.
#     5. Ensure consistency between community, trait, and coordinate datasets.
###############################################################################

# Packages .............................................................................................................
if(!require("pacman")) {install.packages("pacman")}
if(!require("MODISTools")) devtools::install_github("khufkens/MODISTools")

pacman::p_load(
  BAT,
  hypervolume,
  tidyverse,
  glue,
  gawdis,
  readxl,
  vegan,
  terra,
  geodata,
  MODISTools)

#' This script is designed to:
#' - data cleaning and data transformation of the community and trait datasets
#' - estimation of species dissimilarities based on the selected functional traits

#' This is the cleaning procedure applied to all datasets as performed in the function clean_data()
#' - keep only sites for which there are coordinates
#' - keep only species for which we have traits
#' - remove sites with < 2 species (needed for hypervolume estimation)
#' - remove species with zero abundance in all sites (avoid computational time for unused species)
#' - In the trait table, keep only species which are present in the community data
#' - In the coordinates, keep only the ones for the sites with communities

# grab the array id value from the environment variable passed from sbatch
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')

# coerce the value to an integer
ii <- as.numeric(slurm_arrayid)

#' -----------------------------------------------------------------------------------------------------------------
# @ run analysis ######
#' -----------------------------------------------------------------------------------------------------------------
 

files<-list.files("S1_Preprocessing/Raw_data")

focal_dataset <- tools::file_path_sans_ext(files[ii])


  # read dataset info
  data_info <- readxl::read_xlsx("S1_Preprocessing/Miscellaneous/dataset_info_all.xlsx") 
  source("S1_Preprocessing/Functions/auxiliary_functions.R")
  
  tryCatch({
    write("clean data", stdout())
    write(focal_dataset, stdout())
    
    # clean data
    data_clean <- clean_data(focal_dataset)
    
    # get predictors
    data_predictors <- data_clean %>% 
      pluck("coord") %>%
      drop_na()  %>% 
      mutate(site = as.character(site)) %>%
      get_predictors(focal_dataset)

    data_clean$predictors <- data_predictors 
    
    # save pre-processed data
    saveRDS(object = data_clean, file = glue::glue("S1_Preprocessing/pre_processed/{focal_dataset}_clean.rds"))
    
    # remove objects from memory
    rm(list = setdiff(ls(), "files"))
  }, error = function(e) {
    write(glue::glue("Error occurred in processing {focal_dataset}: {e}\n"), stdout())
  })


# -------------------------------------------------------------------------------------------------------------------




