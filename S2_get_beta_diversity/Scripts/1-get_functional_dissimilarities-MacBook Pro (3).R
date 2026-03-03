#'-----------------------------------------------------------------------------
# SCRIPT NAME: Functional Dissimilarities Estimation
#
# DESCRIPTION:
#   This script calculates functional dissimilarities between species using
#   the `gawdis` function, based on selected functional traits.
#
# USAGE:
#   - This script is designed to process individual datasets included in the
#     analysis. Each dataset corresponds to a specific set of functional traits
#     and species.
#   - To run the analysis for a specific dataset, uncomment the relevant lines
#     and execute the script.
#
# INPUTS:
#   - Functional trait data for species
#   - Dataset files corresponding to the selected traits
#
# OUTPUTS:
#   - Functional dissimilarity matrices for each processed dataset
#
# AUTHOR: Caio Graco-Roza
# LAST UPDATED: 2024-11-20

# NOTES:
#   - This script is modular, allowing users to process datasets one at a time
#     by uncommenting specific sections.
#   - Ensure all required libraries are installed, and the input files are
#     correctly formatted.
#'-----------------------------------------------------------------------------

# Packages .............................................................................................................
if (!require("pacman")) {
  install.packages("pacman")
}
pacman::p_load(tidyverse, magrittr, glue, gawdis, readxl, vegan, FD)

calculate_trait_space_quality <- function(dist_matrix, pcoa_vectors, axes) {
  # Original and reduced-space distances
  original_dist <- as.matrix(dist_matrix)
  pcoa_dist <- as.matrix(dist(pcoa_vectors))
  # Calculate Mean Squared Deviation
  msd <- mean((original_dist - pcoa_dist)^2, na.rm = TRUE)  # using lower triangle optional
  # Calculate Mean Absolute Deviation
  mad <- mean(abs(original_dist - pcoa_dist), na.rm = TRUE)
  return(list(MSD = msd, MAD = mad))
}

# Function to extract code blocks from an R script
extract_dataset_chunks <- function(file_path) {
  # Read all lines
  lines <- readLines(file_path)
  
  # Find all header line triplets
  header_starts <- c()
  dataset_names <- c()
  
  for (i in 1:(length(lines) - 2)) {
    if (
      grepl("^#'#{78,}$", lines[i]) &&
      grepl("^#  Dataset [0-9]+: .+ ####$", lines[i + 1]) &&
      grepl("^#'#{78,}$", lines[i + 2])
    ) {
      header_starts <- c(header_starts, i)
      dataset_names <- c(dataset_names, sub("^#  Dataset [0-9]+: (.+) ####$", "\\1", lines[i + 1]))
    }
  }
  
  # Extract code chunks between header blocks
  chunks <- lapply(seq_along(header_starts), function(i) {
    start <- header_starts[i] + 3
    end <- if (i < length(header_starts)) header_starts[i + 1] - 1 else length(lines)
    lines[start:end]
  })
  
  # Return named list
  names(chunks) <- dataset_names
  return(chunks)
}

blocks<-extract_dataset_chunks("S2_get_beta_diversity/Scripts/prepare_dataset_for_trait_analysis.R")

  tryCatch(
    {
      eval(parse(text = blocks[[ii]]))
      cat(glue("✅ Success: {names(blocks)[ii]}\n\n"))
      TRUE
    },
    error = function(e) {
      cat(glue("❌ Error in block {ii} ({names(blocks)[ii]}): {e$message}\n\n"))
      failed <<- c(failed, name)
      FALSE
    }
  )

