###############################################################################
# SCRIPT NAME: pre_processing.R
#
# DESCRIPTION:
#   Cleans and harmonizes community/trait datasets and extracts site-level
#   predictors (human footprint, MODIS, climate) for downstream analyses.
#
# USAGE:
#   - SLURM mode: set `SLURM_ARRAY_TASK_ID` to run one dataset per job array.
#   - Local mode: run without `SLURM_ARRAY_TASK_ID` to process all datasets.
#
# INPUTS:
#   - `S1_Preprocessing/Raw_data/*.xlsx`
#   - `S1_Preprocessing/Miscellaneous/dataset_info_all.xlsx`
#   - `S1_Preprocessing/Functions/auxiliary_functions.R`
#
# OUTPUTS:
#   - `S1_Preprocessing/Processed/<dataset>_clean.rds`
#
# AUTHOR: Caio Graco-Roza
# DATE CREATED: 2021-05-14
# LAST UPDATED: 2026-03-03
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
slurm_arrayid <- Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "")

<<<<<<< Updated upstream
files <- list.files("S1_Preprocessing/Raw_data", pattern = "\\.xlsx$", ignore.case = TRUE)
dataset_ids <- tools::file_path_sans_ext(files)
=======
# coerce the value to an integer
ii <- 1 #as.numeric(slurm_arrayid) #if using HPC
>>>>>>> Stashed changes

if (slurm_arrayid != "") {
  ii <- as.numeric(slurm_arrayid)
  if (is.na(ii) || ii < 1 || ii > length(dataset_ids)) {
    stop(glue::glue("Invalid SLURM_ARRAY_TASK_ID: {slurm_arrayid}. Valid range: 1-{length(dataset_ids)}"))
  }
  focal_datasets <- dataset_ids[ii]
} else {
  focal_datasets <- dataset_ids
}

# @ run analysis ######
<<<<<<< Updated upstream
data_info <- readxl::read_xlsx("S1_Preprocessing/Miscellaneous/dataset_info_all.xlsx")
source("S1_Preprocessing/Functions/auxiliary_functions.R")

for (focal_dataset in focal_datasets) {
=======
#' -----------------------------------------------------------------------------------------------------------------

for (ii in length(list.files("S1_Preprocessing/Raw_data"))){
  
files<-list.files("S1_Preprocessing/Raw_data")

focal_dataset <- tools::file_path_sans_ext(files[ii])


  # read dataset info
  data_info <- readxl::read_xlsx("S1_Preprocessing/Miscellaneous/dataset_info_all.xlsx") 
  source("S1_Preprocessing/Functions/auxiliary_functions.R")
  
>>>>>>> Stashed changes
  tryCatch({
    write("clean data", stdout())
    write(focal_dataset, stdout())

    # clean data
    data_clean <- clean_data(focal_dataset)

    # get predictors
    data_predictors <- data_clean %>%
      pluck("coord") %>%
      drop_na() %>%
      mutate(site = as.character(site)) %>%
      get_predictors(focal_dataset)

    data_clean$predictors <- data_predictors

    # save pre-processed data
    saveRDS(
      object = data_clean,
      file = glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds")
    )
  }, error = function(e) {
    write(glue::glue("Error occurred in processing {focal_dataset}: {e}\n"), stdout())
  })
}
# -------------------------------------------------------------------------------------------------------------------
