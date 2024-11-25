###############################################################################
#' SCRIPT NAME: 1-compile_synthesis_data.R
#'
#' DESCRIPTION:
#'   This script compiles, processes, and synthesizes data from previous steps 
#'   of the analysis workflow to generate a comprehensive dataset. The output 
#'   integrates key metrics such as species richness, spatial extent, human 
#'   footprint metrics, effect sizes, and shapes of relationships, along with 
#'   metadata for each dataset.
#'
#' USAGE:
#'   - This script processes all datasets included in the analysis workflow.
#'   - It integrates outputs from Steps 1-5, consolidating information into 
#'     a single dataset for further analysis and visualization.
#'
#' INPUTS:
#'   - BBGDM results: Located in `S4_run_BBGDM/bbgdm_output`.
#'   - Processed datasets: Located in `S1_Preprocessing/Processed`.
#'   - Metadata file: `S1_Preprocessing/Miscellaneous/dataset_info_all.xlsx`.
#'
#' OUTPUTS:
#'   - Synthesized dataset saved as `synthesis_data.xlsx` in 
#'     `S6_Synthesis_model/data/`.
#'
#' AUTHOR: Caio Graco-Roza
#' LAST UPDATED: 2024-11-24
#'
#' NOTES:
#'   - Ensure all required libraries are installed before running the script.
#'   - Excludes datasets that do not cross land-use gradients (`N67TTP` and 
#'     `N78TTP`).
#'   - Designed to be modular and scalable for additional datasets or metrics.
###############################################################################

# Synthesis analysis -----------------------------------------------------------------------------
library(tidyverse) #data wrangling
library(magrittr) #manipulate data
library(readxl) # read excel files
library(here) #make file paths 
library(GeoRange) # spatial extent
library(fossil) # min and max spatial distance
library(qdapRegex)
library(glue)

source("S4_run_BBGDM/functions/helper_functions_S4.R")

files <- tools::file_path_sans_ext(list.files("S4_run_BBGDM/bbgdm_output")) # List file names in the directory without extensions

dataset_features <- "S1_Preprocessing/Processed" |> # Define the directory for processed datasets
  list.files() |> # List all files in the directory
  purrr::map( # Iterate through each file
    ~ here::here("S1_Preprocessing/Processed", .x) |> # Generate the full file path
      readr::read_rds() |> # Read the RDS file
      pluck("predictors") # Extract the "predictors" element
  ) |>
  set_names(tools::file_path_sans_ext(list.files("S1_Preprocessing/Processed"))) |> # Set names to match file base names
  bind_rows(.id = "dataset") |> # Combine all rows into a single data frame, with a column for dataset name
  mutate(dataset = str_remove(dataset, "_clean")) |> # Remove "_clean" from dataset names
  left_join( # Join with additional metadata
    "S1_Preprocessing/Miscellaneous/dataset_info_all.xlsx" |> # Define the metadata file
      readxl::read_xlsx(), # Read the Excel file
    by = c("dataset" = "dataset_name") # Match datasets by name
  ) 

bbgdms_output <- purrr::map(files, ~ readr::read_rds(glue::glue("S4_run_BBGDM/bbgdm_output/{.x}.rds"))) |> # Read BBGDM results for each file
  set_names(gsub("_bbgdm", "", files)) |> # Set names by removing "_bbgdm" from file names
  discard_at(c("N67TTP", "N78TTP")) # Remove datasets that do not cross the land-use gradient

#'---------------------------------------------------------------------------------------------------------------
#  synthesis BBGDM ###########################################################################################
#'---------------------------------------------------------------------------------------------------------------

# Extract the response features for all datasets -----------------------------------------------------------------------
mean_r2 <- bbgdms_output |> 
  enframe(name = "dataset", value = "dataset_value") |> # Convert list of datasets into a tibble with "dataset" as a column
  mutate(dataset_value = purrr::map(dataset_value, ~ enframe(.x, name = "facet", value = "facet_value"))) |> # Extract each dataset's facets as tibbles
  unnest(dataset_value, keep_empty = TRUE) |> # Expand nested tibbles for each dataset's facets
  mutate(facet_value = purrr::map(facet_value, ~ enframe(.x, name = "direction", value = "direction_value"))) |> # Extract each facet's directions as tibbles
  unnest(facet_value, keep_empty = TRUE) |> # Expand nested tibbles for each facet's directions
  dplyr::mutate(r2 = purrr::map(direction_value, ~ .x |>  pluck("r2")),.keep="unused") |> # Extract R² values for each direction
  unnest(r2) |> # Expand nested R² values into rows
  mutate(across(c(r2), ~ ifelse(.x <= 0, NA, .x/100))) |> # Set non-positive R² values to NA and scale valid values by dividing by 100
  group_by(dataset,facet) |> # Group data by dataset and facet
  summarise(direction_r2 = gammaEffectSize(y=r2[direction == "homogenization"],x=r2[direction == "differentiation"],prob=.5)) # Compute gamma effect size between homogenization and differentiation R² values

mean_effect_size <- bbgdms_output |> 
  enframe(name = "dataset", value = "dataset_value") |> # Convert list of datasets into a tibble with "dataset" as a column
  mutate(dataset_value = purrr::map(dataset_value, ~ enframe(.x, name = "facet", value = "facet_value"))) |> # Extract each dataset's facets as tibbles
  unnest(dataset_value, keep_empty = TRUE) |> # Expand nested tibbles for each dataset's facets
  mutate(facet_value = purrr::map(facet_value, ~ enframe(.x, name = "direction", value = "direction_value"))) |> # Extract each facet's directions as tibbles
  unnest(facet_value, keep_empty = TRUE) |> # Expand nested tibbles for each facet's directions
  mutate(coef = purrr::map(direction_value, ~ .x |>  pluck("coefs") |> bind_rows() |> rownames_to_column("predictor")),.keep="unused") |> # Extract and bind coefficients for each direction into a single data frame
  unnest(coef) |> # Expand coefficients into rows
  mutate(predictor = str_remove(predictor,"\\.{3}.*")) |> # Remove trailing "...*" from predictor names
  filter(grepl("^hfp",predictor)) |> # Keep only predictors starting with "hfp"
  mutate(effect_size = (coefficient.1+coefficient.2+coefficient.3), .keep="unused") |> # Compute combined effect size from three coefficients
  group_by(dataset,facet,direction,predictor) |> # Group data by dataset, facet, direction, and predictor
  summarise(magnitude = median(effect_size)*(n()/1000)) |> # Calculate median magnitude, scaling by dataset size
  separate(predictor, into=c("predictor","buffer")) # Split predictor names into "predictor" and "buffer" components

shape <- bbgdms_output |> 
  enframe(name = "dataset", value = "dataset_value") |> # Convert list of datasets into a tibble with "dataset" as a column
  mutate(dataset_value = purrr::map(dataset_value, ~ enframe(.x, name = "facet", value = "facet_value"))) |> # Extract each dataset's facets as tibbles
  unnest(dataset_value, keep_empty = TRUE) |> # Expand nested tibbles for each dataset's facets
  mutate(facet_value = purrr::map(facet_value, ~ enframe(.x, name = "direction", value = "direction_value"))) |> # Extract each facet's directions as tibbles
  unnest(facet_value, keep_empty = TRUE) |> # Expand nested tibbles for each facet's directions
  mutate(coef = purrr::map(direction_value, ~ .x |>  pluck("coefs") |> bind_rows() |> rownames_to_column("predictor")),.keep="unused") |> # Extract and bind coefficients for each direction into a single data frame
  unnest(coef) |> # Expand coefficients into rows
  mutate(predictor = str_remove(predictor,"\\.{3}.*")) |> # Remove trailing "...*" from predictor names
  filter(grepl("^hfp",predictor)) |> # Keep only predictors starting with "hfp"
  dplyr::mutate(shape = dplyr::case_when( # Classify the shape of predictor effects
    coefficient.1 + coefficient.2 + coefficient.3 == 0 ~ "Absent",
    coefficient.1 <= coefficient.2 & coefficient.2 < coefficient.3 ~ "Exponential",
    coefficient.1 > coefficient.2 & coefficient.2 > coefficient.3 ~ "Saturating",
    coefficient.1 < coefficient.2 & coefficient.2 > coefficient.3 ~ "Saturating",
    coefficient.1 < coefficient.2 & coefficient.2 < coefficient.3 ~ "Saturating",
    coefficient.1 > coefficient.2 & coefficient.2 == coefficient.3 ~ "Saturating",
    coefficient.1 > coefficient.2 & coefficient.2 < coefficient.3 & coefficient.3 > 0 ~ "Revlog"))  |> 
  
  group_by(dataset,facet,direction) |> # Group data by dataset, facet, and direction
  mutate(shape = factor(shape, levels = c("Absent", "Saturating", "Exponential", "Revlog"))) |> # Convert shape categories to ordered factors
  count(shape) |>  # Count occurrences of each shape category within groups
  complete(shape) |> # Ensure all shape categories are present, even if counts are zero
  mutate(n = ifelse(is.na(n), 0, n)) |> # Replace NA counts with zeros
  pivot_wider(names_from = shape, values_from = n, values_fill = 0) # Reshape the data to have one column per shape category

species_number <- purrr::map(gsub("_bbgdm", "", files), # Remove the "_bbgdm" suffix from dataset filenames.
                             ~ readr::read_rds(glue::glue( # Load the processed data file for each dataset.
                               "S1_Preprocessing/Processed/{.x}_clean.rds"
                             )) |>
                               pluck("comm") |> # Extract the "comm" (community matrix) component.
                               ncol()) |> # Count the number of species (columns in the community matrix).
  set_names(gsub("_bbgdm", "", files)) |> # Assign dataset names to the resulting list.
  bind_rows() |>  # Combine the results into a single data frame.
  pivot_longer(everything(), names_to="dataset",values_to = "species.number") |> # Reshape data into long format with dataset names and species counts.
  filter(!dataset %in% c("N67TTP", "N78TTP")) # Exclude datasets that do not cross the land use gradient.

# Get extent as the area of a minimum convex hull enclosing all coordinates (square kilometers).
spatial.area <- purrr::map(gsub("_bbgdm", "", files), # Remove the "_bbgdm" suffix from dataset filenames.
                           ~ readr::read_rds(glue::glue( # Load the processed data file for each dataset.
                             "S1_Preprocessing/Processed/{.x}_clean.rds"
                           )) |>
                             pluck("coord") |> # Extract the "coord" (coordinates) component.
                             drop_na(x, y)) |> # Remove rows with missing x or y coordinates.
  purrr::map(~GeoRange::CHullAreaEarth(.x$x, .x$y)) |> # Calculate the area of the convex hull enclosing the coordinates.
  set_names(gsub("_bbgdm", "", files)) |> # Assign dataset names to the results.
  bind_rows() |>  # Combine the results into a single data frame.
  pivot_longer(everything(), names_to="dataset",values_to = "spatial.extent") |> # Reshape data into long format with dataset names and spatial extents.
  filter(!dataset %in% c("N67TTP", "N78TTP")) # Exclude datasets that do not cross the land use gradient.

spatial.summary <- purrr::map(gsub("_bbgdm", "", files), # Remove the "_bbgdm" suffix from dataset filenames.
                              ~ readr::read_rds(glue::glue( # Load the processed data file for each dataset.
                                "S1_Preprocessing/Processed/{.x}_clean.rds"
                              )) |>
                                pluck("coord") |> # Extract the "coord" (coordinates) component.
                                drop_na(x, y)) |> # Remove rows with missing x or y coordinates.
  purrr::map(~ summary(fossil::earth.dist(data.frame(.x$x, .x$y)))[c(1, 4, 6)]) |> # Calculate min, mean, and max distances between points.
  set_names(gsub("_bbgdm", "", files)) |> # Assign dataset names to the results.
  bind_rows(.id="dataset") |> # Combine the results into a single data frame with dataset names.
  set_names(c("dataset", "spatial.min", "spatial.mean", "spatial.max")) |> # Rename columns for clarity.
  dplyr::mutate(across(contains("spatial"), as.numeric)) |> # Ensure spatial values are numeric.
  filter(!dataset %in% c("N67TTP", "N78TTP")) # Exclude datasets that do not cross the land use gradient.

# Get mean latitude value of the dataset.
latitude.mean <- purrr::map(gsub("_bbgdm", "", files), # Remove the "_bbgdm" suffix from dataset filenames.
                            ~ readr::read_rds(glue::glue( # Load the processed data file for each dataset.
                              "S1_Preprocessing/Processed/{.x}_clean.rds"
                            )) |>
                              pluck("coord") |> # Extract the "coord" (coordinates) component.
                              select(x, y)) |> # Select x and y columns for the coordinates.
  set_names(gsub("_bbgdm", "", files)) |> # Assign dataset names to the results.
  bind_rows(.id="dataset") |> # Combine the results into a single data frame with dataset names.
  group_by(dataset) |> # Group by dataset for summary calculation.
  summarise(latitude.mean = mean(y, na.rm=TRUE)) |> # Calculate the mean latitude for each dataset, ignoring missing values.
  filter(!dataset %in% c("N67TTP", "N78TTP")) # Exclude datasets that do not cross the land use gradient.


hfp_min_max <- 
  purrr::map(gsub("_bbgdm", "", files), # Remove the "_bbgdm" suffix from dataset filenames.
             ~ readr::read_rds(glue::glue( # Load the processed data file for each dataset.
               "S1_Preprocessing/Processed/{.x}_clean.rds"
             )) |>
               pluck("predictors") |> # Extract the "predictors" component of the dataset.
               summarise(across(contains("hfp"), # Select all columns containing "hfp" (Human Footprint data).
                                list(min=~min(.x, na.rm=TRUE), # Calculate the minimum value for each "hfp" column.
                                     mean=~mean(.x, na.rm=TRUE), # Calculate the mean value for each "hfp" column.
                                     max=~max(.x, na.rm=TRUE)), # Calculate the maximum value for each "hfp" column.
                                .names="{.fn}.{.col}"))) |> # Assign new column names using the summary function names.
  set_names(gsub("_bbgdm", "", files)) |> # Assign dataset names to the results.
  bind_rows(.id="dataset") |> # Combine the results into a single data frame with dataset names.
  pivot_longer( # Reshape data into a long format with separate rows for each Human Footprint statistic.
    cols = !dataset, # All columns except "dataset".
    names_to = "predictor", # Column containing the predictor name.
    values_to = "value" # Column containing the computed value.
  ) |> 
  separate(predictor, into = c("fun", "predictor"), sep="\\.") |> # Split the "predictor" column into function (min/mean/max) and actual predictor name.
  pivot_wider(names_from="fun", values_from = "value", names_prefix = "hfp.") |> # Reshape data back into wide format with separate columns for min, mean, and max.
  separate(predictor, into=c("predictor", "buffer")) |> # Further split the "predictor" column into the actual predictor and its buffer.
  mutate(hfp.range = hfp.max - hfp.min) |> # Calculate the range of Human Footprint values (max - min).
  select(-predictor) # Drop the "predictor" column, leaving only relevant data.

synthesis_data <- 
  mean_r2 |> 
  mutate(direction = ifelse(sign(direction_r2) > 0 , "homogenization", "differentiation")) |> # Assign "homogenization" or "differentiation" based on the sign of `direction_r2`.
  relocate(dataset, facet, direction) |> # Reorganize columns to bring `dataset`, `facet`, and `direction` to the front.
  ungroup() |> # Remove grouping to allow subsequent operations to act on the full dataset.
  left_join(mean_effect_size, by = c("dataset", "facet", "direction")) |> # Merge effect size data on matching `dataset`, `facet`, and `direction`.
  left_join(shape, by = c("dataset", "facet", "direction")) |> # Merge shape data on matching `dataset`, `facet`, and `direction`.
  left_join(hfp_min_max, by = c("dataset", "buffer")) |> # Merge Human Footprint summary stats on matching `dataset` and `buffer`.
  left_join(species_number, by = c("dataset")) |> # Add species number data for each dataset.
  left_join(spatial.area, by = c("dataset")) |> # Add spatial extent (area) data for each dataset.
  left_join(spatial.summary, by = c("dataset")) |> # Add spatial summary stats (min, mean, max distances) for each dataset.
  left_join(latitude.mean, by = c("dataset")) |> # Add mean latitude data for each dataset.
  left_join(dataset_features |> distinct(dataset, disturbance, taxa, group, system), by = c("dataset")) |> # Add metadata (disturbance, taxa, biotic group, system) for each dataset.
  mutate(direction = gsub("z", "s", direction)) |> # Replace "z" with "s" in the `direction` column (e.g., fix typos or conventions).
  rename(realm = system, biotic.group = group) # Rename columns for consistency: `system` to `realm`, `group` to `biotic.group`.

writexl::write_xlsx(data.frame(synthesis_data), here::here("S6_Synthesis_model","data","synthesis_data.xlsx"))
