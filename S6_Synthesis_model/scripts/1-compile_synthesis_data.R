###############################################################################
#' SCRIPT NAME: 1-compile_synthesis_data.R
#'
#' DESCRIPTION:
#'   Compiles and synthesizes multiple data sources into one overarching dataset,
#'   including BBGDM outputs, effect sizes, shape classifications, and metadata.
#'   Designed to feed directly into downstream analyses and figures.
#'
#' USAGE:
#'   Run this script after all preprocessing (S1) and BBGDM model fitting (S4).
#'   It reads in raw model outputs, joins with metadata, computes summary
#'   statistics, and writes a combined `synthesis_data` object for further use.
#'
#' INPUTS:
#'   - BBGDM RDS output files in `S4_run_BBGDM/bbgdm_output/`
#'   - Preprocessed datasets in `S1_Preprocessing/Processed`
#'   - Metadata spreadsheet `dataset_info_all.xlsx`
#'
#' OUTPUTS:
#'   - A single `synthesis_data` tibble containing:
#'     * R² metrics, effect sizes, and shape classifications
#'     * Species counts, spatial extents, latitudinal zones
#'     * Human Footprint & Environmental Human Pressure (HET) summarizations
#'
#' AUTHOR: Caio Graco-Roza
#' LAST UPDATED: 2024-11-24
###############################################################################

# Load required packages for data wrangling and spatial calculations
library(tidyverse)     # core tidy data tools
library(magrittr)      # piping and data manipulation helpers
library(readxl)        # read Excel metadata
library(here)          # project‐relative file paths
library(GeoRange)      # convex hull area for spatial extent
library(geosphere)     # geographic distance calculations
library(gdm)           # BBGDM model utilities
library(glue)

# Source custom helper functions (e.g. classify_shapes, ispline_curve_row)
source("S4_run_BBGDM/functions/helper_functions_S4.R")

# 1. Discover all BBGDM output files and strip extensions for dataset names
files_all <- list.files(
  "S4_run_BBGDM/bbgdm_output", 
  recursive = TRUE,
  full.names = TRUE
)
files <- tools::file_path_sans_ext(files_all)

# 2. Read predictor metadata for each processed dataset
dataset_features <- list.files("S1_Preprocessing/Processed") %>%
  purrr::map(~ here("S1_Preprocessing/Processed", .x) %>%
               read_rds() %>%
               pluck("predictors")) %>%
  set_names(tools::file_path_sans_ext(list.files("S1_Preprocessing/Processed"))) %>%
  bind_rows(.id = "dataset") %>%
  mutate(dataset = str_remove(dataset, "_clean")) %>%
  left_join(
    read_xlsx("S1_Preprocessing/Miscellaneous/dataset_info_all.xlsx"),
    by = c("dataset" = "dataset_name")
  )

# 3. Extract dataset identifiers from file paths
basic_structure <- tibble(fullpath = files) %>%
  mutate(
    dataset     = basename(fullpath),                # raw dataset name
    beta_folder = basename(dirname(fullpath)),       # e.g. "Baselga_abun"
    beta_type   = str_extract(beta_folder, "^(Baselga|Podani)"),
    metric_type = str_extract(beta_folder, "(abun|pa)$")
  ) %>%
  select(dataset, beta_type, metric_type)

# --- Helper function: read BBGDM outputs and unwrap list structure ------------
read_bbgdm <- function(files, .pattern) {
  # .pattern = "r2" or "coefs" or etc.
  purrr::map(files, ~ read_rds(glue::glue("{.x}.rds"))) %>%
    set_names(gsub("_bbgdm$", "", files)) %>%
    enframe(name = "dataset", value = "data_list") %>%
    mutate(
      data_list = purrr::map(data_list, ~ enframe(.x, name = "facet", value = "facet_list")),
      facet_list = purrr::map(data_list, ~ unnest(.x, cols = "facet_list")),
      dir_list   = purrr::map(facet_list, ~ unnest(enframe(.x$facet_list, name="direction", value="direction_list")))
    ) %>%
    unnest(dir_list)
}

# 4A. Compute R² metrics per dataset/facet/direction ---------------------------
mean_r2 <- purrr::map(files, ~ readr::read_rds(glue::glue("{.x}.rds"))) |> # Read BBGDM results for each file
  set_names(gsub("_bbgdm", "", files)) |> 
  enframe(name = "dataset", value = "dataset_value") |> # Convert list of datasets into a tibble with "dataset" as a column
  mutate(dataset_value = purrr::map(dataset_value, ~ enframe(.x, name = "facet", value = "facet_value"))) |> # Extract each dataset's facets as tibbles
  unnest(dataset_value, keep_empty = TRUE) |> # Expand nested tibbles for each dataset's facets
  mutate(facet_value = purrr::map(facet_value, ~ enframe(.x, name = "direction", value = "direction_value"))) |> # Extract each facet's directions as tibbles
  unnest(facet_value, keep_empty = TRUE) |> # Expand nested tibbles for each facet's directions
  dplyr::mutate(r2 = purrr::map(direction_value, ~ .x |>  pluck("r2")),.keep="unused") |> # Extract R² values for each direction
  unnest(r2) |> # Expand nested R² values into rows
  mutate(across(c(r2), ~ ifelse(.x <= 0, NA, .x/100))) |> # Set non-positive R² values to NA and scale valid values by dividing by 100
  group_by(dataset,facet, direction) |> # Group data by dataset and facet
  mutate(iteration = row_number()) |> 
  pivot_wider(names_from = direction, values_from = r2) |>
  summarise(
    harrel_davis =  harrel_davis(homogenization, differentiation),
    median_overlap     = median_overlap(homogenization, differentiation)$hybrid,
    median_diff =  median_overlap(homogenization, differentiation)$r2_diff
  )

# 4B. Compute HFP effect sizes (median of sum of three coefs) --------------------
hfp_effect_size <- purrr::map(files, ~ readr::read_rds(glue::glue("{.x}.rds"))) |> # Read BBGDM results for each file
  set_names(gsub("_bbgdm", "", files)) |> 
  enframe(name = "dataset", value = "dataset_value") |> # Convert list of datasets into a tibble with "dataset" as a column
  mutate(dataset_value = purrr::map(dataset_value, ~ enframe(.x, name = "facet", value = "facet_value"))) |> # Extract each dataset's facets as tibbles
  unnest(dataset_value, keep_empty = TRUE) |> # Expand nested tibbles for each dataset's facets
  mutate(facet_value = purrr::map(facet_value, ~ enframe(.x, name = "direction", value = "direction_value"))) |> # Extract each facet's directions as tibbles
  unnest(facet_value, keep_empty = TRUE) |> # Expand nested tibbles for each facet's directions
  mutate(coef = purrr::map(direction_value, ~ .x |>  pluck("coefs") |> bind_rows() |> rownames_to_column("predictor")),.keep="unused") |> # Extract and bind coefficients for each direction into a single data frame
  unnest(coef) |> # Expand coefficients into rows
  mutate(predictor = str_remove(predictor,"\\.{3}.*")) |> # Remove trailing "...*" from predictor names
  mutate(magnitude = (coefficient.1+coefficient.2+coefficient.3), .keep="unused") |> # Compute combined effect size from three coefficients
  group_by(dataset,facet,direction,predictor) |> # Group data by dataset, facet, direction, and predictor
  summarise(magnitude = median(magnitude)*(n()/1000)) |> # Calculate median magnitude, scaling by dataset size
  group_by(dataset,facet,direction) |> # Group data by dataset, facet, direction, and predictor
  mutate(rel_magnitude = magnitude/sum(magnitude)) |> # Calculate median magnitude, scaling by dataset size
  filter(grepl("^hfp",predictor)) |> # Keep only predictors starting with "hfp"
  separate(predictor, into=c("predictor","hfp_buffer")) |> select(-predictor) |>  # Split predictor names into "predictor" and "buffer" components
  rename(hfp_magnitude = magnitude, hfp_magnitude_rel = rel_magnitude)


# 4C. Compute HET effect sizes (same logic) -------------------------------------
het_effect_size <- purrr::map(files, ~ readr::read_rds(glue::glue("{.x}.rds"))) |> # Read BBGDM results for each file
  set_names(gsub("_bbgdm", "", files)) |> 
  enframe(name = "dataset", value = "dataset_value") |> # Convert list of datasets into a tibble with "dataset" as a column
  mutate(dataset_value = purrr::map(dataset_value, ~ enframe(.x, name = "facet", value = "facet_value"))) |> # Extract each dataset's facets as tibbles
  unnest(dataset_value, keep_empty = TRUE) |> # Expand nested tibbles for each dataset's facets
  mutate(facet_value = purrr::map(facet_value, ~ enframe(.x, name = "direction", value = "direction_value"))) |> # Extract each facet's directions as tibbles
  unnest(facet_value, keep_empty = TRUE) |> # Expand nested tibbles for each facet's directions
  mutate(coef = purrr::map(direction_value, ~ .x |>  pluck("coefs") |> bind_rows() |> rownames_to_column("predictor")),.keep="unused") |> # Extract and bind coefficients for each direction into a single data frame
  unnest(coef) |> # Expand coefficients into rows
  mutate(predictor = str_remove(predictor,"\\.{3}.*")) |> # Remove trailing "...*" from predictor names
  mutate(magnitude = (coefficient.1+coefficient.2+coefficient.3), .keep="unused") |> # Compute combined effect size from three coefficients
  group_by(dataset,facet,direction,predictor) |> # Group data by dataset, facet, direction, and predictor
  summarise(magnitude = median(magnitude)*(n()/1000)) |> # Calculate median magnitude, scaling by dataset size
  group_by(dataset,facet,direction) |> # Group data by dataset, facet, direction, and predictor
  mutate(rel_magnitude = magnitude/sum(magnitude)) |> # Calculate median magnitude, scaling by dataset size
  filter(grepl("^het",predictor)) |> # Keep only predictors starting with "hfp"
  separate(predictor, into=c("predictor","het_buffer")) |> select(-predictor) |>  # Split predictor names into "predictor" and "buffer" components
  rename(het_magnitude = magnitude, het_magnitude_rel = rel_magnitude)

# First, load every file into a tibble and extract metadata
dissimilarity_values<-  tibble(file = list.files(
  "~/Library/CloudStorage/OneDrive-UniversityofHelsinki/Ongoing manuscripts/HIATES/Analysis/HIATE/S2_get_beta_diversity/betadiv_output",
  recursive = TRUE, full.names = TRUE
)) |>
  mutate(
    # Extract dataset name from the file basename
    dataset = basename(file) |> str_remove("_beta_Output\\.rds"),
    # Get the directory part to extract more info (e.g., "Baselga_abun", "Podani_pa", etc.)
    dir_info = dirname(file) |> basename(),
    # Extract framework: assume it is either "Baselga" or "Podani"
    framework = str_extract(dir_info, "Baselga|Podani"),
    # Extract feature: assume it is either "abun" or "pa"
    metric_type = str_extract(dir_info, "abun|pa"),
    # Load the RDS object; each object is assumed to be a list with facets ("taxonomic", "functional")
    obj = map(file, readRDS)
  ) |> 
  mutate(
    facet_data = map(obj, ~ enframe(.x, name = "facet", value = "beta"))
  ) |>
  unnest(facet_data) |>
  select(dataset, metric_type, framework, facet, beta) |> 
  pivot_wider(
    names_from = framework,
    values_from = beta
  ) |>
  arrange(dataset, facet, metric_type) |> 
  mutate(
    Baselga = map(Baselga, ~ .x[[2]]),
    Podani  = map(Podani,  ~ .x[[2]])
  )  |> 
  # rowwise() |> 
  # mutate(cor = cor(Baselga,Podani)) |> 
  ungroup() |> 
  # Convert the dist objects to numeric vectors
  mutate(
    Baselga_values = map(Baselga, ~ as.numeric(.x)),
    Podani_values = map(Podani, ~ as.numeric(.x)),
  ) |>
  # Pivot the two columns into a long format
  pivot_longer(cols = c("Baselga_values", "Podani_values"),
               names_to = "beta_type", values_to = "values") |>
  # Clean up beta_type names (remove _values suffix)
  mutate(beta_type = sub("_values", "", beta_type)) |>
  # Calculate median and quantiles for each row
  mutate(
    beta_median = map_dbl(values, ~ Hmisc::smedian.hilow(.x, na.rm = TRUE)[1]),
    beta_range = map_dbl(values, ~ abs(max(.x,na.rm=TRUE) - min(.x,na.rm=TRUE)))
  ) |> 
  select(dataset,beta_type,metric_type,facet,beta_median,beta_range) 


species_number <- purrr::map(gsub("_bbgdm", "", unique(basic_structure$dataset)), # Remove the "_bbgdm" suffix from dataset filenames.
                             ~ readr::read_rds(glue::glue( # Load the processed data file for each dataset.
                               "S1_Preprocessing/Processed/{.x}_clean.rds"
                             )) |>
                               pluck("comm") |> # Extract the "comm" (community matrix) component.
                               ncol()) |> # Count the number of species (columns in the community matrix).
  set_names(gsub("_bbgdm", "", unique(basic_structure$dataset))) |> # Assign dataset names to the resulting list.
  bind_rows() |>  # Combine the results into a single data frame.
  pivot_longer(everything(), names_to="dataset",values_to = "species.number") 

hyper_quality <- purrr::map(gsub("_bbgdm", "", unique(basic_structure$dataset)), # Remove the "_bbgdm" suffix from dataset filenames.
                            ~ readr::read_rds(glue::glue( # Load the processed data file for each dataset.
                              "S2_get_beta_diversity/betadiv_input/{.x}_beta_Input.rds"
                            )) |>
                              pluck("quality") |>  pluck("MAD")) |> 
  set_names(gsub("_bbgdm", "", unique(basic_structure$dataset))) |>  # Assign dataset names to the results.
  bind_rows() |> 
  pivot_longer(everything(), names_to="dataset",values_to = "hyper_quality") 

# Get extent as the area of a minimum convex hull enclosing all coordinates (square kilometers).
spatial.area <- purrr::map(gsub("_bbgdm", "", unique(basic_structure$dataset)), # Remove the "_bbgdm" suffix from dataset filenames.
                           ~ readr::read_rds(glue::glue( # Load the processed data file for each dataset.
                             "S1_Preprocessing/Processed/{.x}_clean.rds"
                           )) |>
                             pluck("coord") |> # Extract the "coord" (coordinates) component.
                             drop_na(x, y)) |> # Remove rows with missing x or y coordinates.
  purrr::map(~GeoRange::CHullAreaEarth(.x$x, .x$y)) |> # Calculate the area of the convex hull enclosing the coordinates.
  set_names(gsub("_bbgdm", "", unique(basic_structure$dataset))) |> # Assign dataset names to the results.
  bind_rows() |>  # Combine the results into a single data frame.
  pivot_longer(everything(), names_to="dataset",values_to = "spatial.extent") 


spatial.summary <-
  purrr::map(
    gsub("_bbgdm", "", unique(basic_structure$dataset)),
    ~ {
      coords <- read_rds(glue::glue("S1_Preprocessing/Processed/{.x}_clean.rds")) |>
        pluck("coord") |>
        drop_na(x, y)
      
      if (nrow(coords) < 2) {
        tibble(dataset = .x, min = NA, median = NA, max = NA)
      } else {
        dists <- distm(coords[, c("x", "y")], fun = distHaversine)
        dists_vec <- dists[lower.tri(dists)]
        
        tibble(
          dataset=.x,
          spatial_min = min(dists_vec),
          spatial_median = median(dists_vec),
          spatial_max = max(dists_vec)
        )
      }
    }
  ) |>
  bind_rows() # Calculate min, mean, and max distances between points.

# Get mean latitude value of the dataset.
latitude.mean <- purrr::map(gsub("_bbgdm", "", unique(basic_structure$dataset)), # Remove the "_bbgdm" suffix from dataset filenames.
                            ~ readr::read_rds(glue::glue( # Load the processed data file for each dataset.
                              "S1_Preprocessing/Processed/{.x}_clean.rds"
                            )) |>
                              pluck("coord") |> # Extract the "coord" (coordinates) component.
                              dplyr::select(x, y)) |> # Select x and y columns for the coordinates.
  set_names(gsub("_bbgdm", "", unique(basic_structure$dataset))) |> # Assign dataset names to the results.
  bind_rows(.id="dataset") |> # Combine the results into a single data frame with dataset names.
  group_by(dataset) |> # Group by dataset for summary calculation.
  summarise(latitude.mean = mean(y, na.rm=TRUE)) |> 
  mutate(latitudinal.zone = if_else(
    abs(latitude.mean) <= 23.437,   # inside the tropics
    "Tropical",
    "Temperate"
  ),
  latitudinal.zone = factor(latitudinal.zone, levels = c("Tropical","Temperate"))
  )


# 6. Summarize HFP/HET predictor ranges (min/median/max) -----------------------
hfp_min_max <- 
  purrr::map(gsub("_bbgdm", "", unique(basic_structure$dataset)), # Remove the "_bbgdm" suffix from dataset filenames.
             ~ readr::read_rds(glue::glue( # Load the processed data file for each dataset.
               "S1_Preprocessing/Processed/{.x}_clean.rds"
             )) |>
               pluck("predictors") |> # Extract the "predictors" component of the dataset.
               summarise(across(contains("hfp"), # Select all columns containing "hfp" (Human Footprint data).
                                list(min=~min(.x, na.rm=TRUE), # Calculate the minimum value for each "hfp" column.
                                     median=~median(.x, na.rm=TRUE), # Calculate the mean value for each "hfp" column.
                                     max=~max(.x, na.rm=TRUE)), # Calculate the maximum value for each "hfp" column.
                                .names="{.fn}.{.col}"))) |> # Assign new column names using the summary function names.
  set_names(gsub("_bbgdm", "", unique(basic_structure$dataset))) |> # Assign dataset names to the results.
  bind_rows(.id="dataset") |> # Combine the results into a single data frame with dataset names.
  pivot_longer( # Reshape data into a long format with separate rows for each Human Footprint statistic.
    cols = !dataset, # All columns except "dataset".
    names_to = "predictor", # Column containing the predictor name.
    values_to = "value" # Column containing the computed value.
  ) |> 
  separate(predictor, into = c("fun", "predictor"), sep="\\.") |> # Split the "predictor" column into function (min/mean/max) and actual predictor name.
  pivot_wider(names_from="fun", values_from = "value", names_prefix = "hfp_") |> # Reshape data back into wide format with separate columns for min, mean, and max.
  separate(predictor, into=c("predictor", "hfp_buffer")) |> # Further split the "predictor" column into the actual predictor and its buffer.
  mutate(hfp_range = hfp_max - hfp_min) |> # Calculate the range of Human Footprint values (max - min).
  select(-predictor) # Drop the "predictor" column, leaving only relevant data.





het_min_max <- 
  purrr::map(gsub("_bbgdm", "", unique(basic_structure$dataset)), # Remove the "_bbgdm" suffix from dataset filenames.
             ~ readr::read_rds(glue::glue( # Load the processed data file for each dataset.
               "S1_Preprocessing/Processed/{.x}_clean.rds"
             )) |>
               pluck("predictors") |> # Extract the "predictors" component of the dataset.
               summarise(across(contains("het"), # Select all columns containing "hfp" (Human Footprint data).
                                list(min=~min(.x, na.rm=TRUE), # Calculate the minimum value for each "hfp" column.
                                     median=~median(.x, na.rm=TRUE), # Calculate the mean value for each "hfp" column.
                                     max=~max(.x, na.rm=TRUE)), # Calculate the maximum value for each "hfp" column.
                                .names="{.fn}.{.col}"))) |> # Assign new column names using the summary function names.
  set_names(gsub("_bbgdm", "", unique(basic_structure$dataset))) |> # Assign dataset names to the results.
  bind_rows(.id="dataset") |> # Combine the results into a single data frame with dataset names.
  pivot_longer( # Reshape data into a long format with separate rows for each Human Footprint statistic.
    cols = !dataset, # All columns except "dataset".
    names_to = "predictor", # Column containing the predictor name.
    values_to = "value" # Column containing the computed value.
  ) |> 
  separate(predictor, into = c("fun", "predictor"), sep="\\.") |> # Split the "predictor" column into function (min/mean/max) and actual predictor name.
  pivot_wider(names_from="fun", values_from = "value", names_prefix = "het_") |> # Reshape data back into wide format with separate columns for min, mean, and max.
  separate(predictor, into=c("predictor", "het_buffer")) |> # Further split the "predictor" column into the actual predictor and its buffer.
  mutate(het_range = het_max - het_min) |> # Calculate the range of Human Footprint values (max - min).
  select(-predictor) # Drop the "predictor" column, leaving only relevant data.


# 7. Final synthesis: join everything into one tibble ---------------------------
combined_data <-
  mean_r2 |> 
  mutate(direction = ifelse(sign(median_diff) > 0 , "homogenization", "differentiation")) |> # Assign "homogenization" or "differentiation" based on the sign of `direction_r2`.
  ungroup() |> # Remove grouping to allow subsequent operations to act on the full dataset.
  left_join(hfp_effect_size, by = c("dataset", "facet", "direction")) |> # Merge effect size data on matching `dataset`, `facet`, and `direction`.
  left_join(het_effect_size, by = c("dataset", "facet", "direction")) |> # Merge effect size data on matching `dataset`, `facet`, and `direction`.
  mutate(
    beta_folder = basename(dirname(dataset)),       # e.g., "Baselga_abun"
    beta_type = str_extract(beta_folder, "^(Baselga|Podani)"),
    metric_type = str_extract(beta_folder, "(abun|pa)$")) |> 
  mutate(dataset = basename(dataset)) |> 
  left_join(hfp_min_max, by = c("dataset", "hfp_buffer")) |> # Merge Human Footprint summary stats on matching `dataset` and `buffer`.
  left_join(het_min_max, by = c("dataset", "het_buffer")) |> # Merge Human Footprint summary stats on matching `dataset` and `buffer`.
  left_join(species_number, by = c("dataset")) |> # Add species number data for each dataset.
  left_join(hyper_quality, by = c("dataset")) |> # Add the quality of hypervolume data for each dataset.
  left_join(spatial.area, by = c("dataset")) |> # Add spatial extent (area) data for each dataset.
  left_join(spatial.summary, by = c("dataset")) |> # Add spatial summary stats (min, mean, max distances) for each dataset.
  left_join(latitude.mean, by = c("dataset")) |> # Add mean latitude data for each dataset.
  left_join(dataset_features |> distinct(dataset, disturbance, taxa, group, system), by = c("dataset")) |> # Add metadata (disturbance, taxa, biotic group, system) for each dataset.
  mutate(direction = gsub("z", "s", direction)) |> # Replace "z" with "s" in the `direction` column (e.g., fix typos or conventions).
  rename(realm = system, biotic.group = group)  |> 
  relocate(dataset,facet,beta_type, metric_type) |> 
  left_join(dissimilarity_values) |> 
  mutate(beta_median = ifelse(direction == "homogenisation", 1-beta_median,beta_median))

glimpse(combined_data) 

combined_data |> 
  filter(hfp_magnitude> 0) |> 
  drop_na() |> 
  ggplot(aes(x=hfp_magnitude,y=beta_median))+
  geom_point(aes(color=direction))+
  geom_smooth(method="lm")+ 
  facet_grid(beta_type~metric_type)+
  coord_cartesian(xlim=c(0,1), ylim=c(0,1))


# -----------------------------------------------------------------------------
# 8. Merge BBGDM model outputs with synthesized metadata and classify shapes -----
# -----------------------------------------------------------------------------

models_collected<-purrr::map(files, ~ readr::read_rds(glue::glue("{.x}.rds"))) |> # Read BBGDM results for each file
  set_names(gsub("_bbgdm", "", files)) |>
  enframe(name = "dataset", value = "dataset_value") |> # Convert list of datasets into a tibble with "dataset" as a column
  mutate(dataset_value = purrr::map(dataset_value, ~ enframe(.x, name = "facet", value = "facet_value"))) |> # Extract each dataset's facets as tibbles
  unnest(dataset_value, keep_empty = TRUE) |> # Expand nested tibbles for each dataset's facets
  mutate(facet_value = purrr::map(facet_value, ~ enframe(.x, name = "direction", value = "direction_value"))) |> # Extract each facet's directions as tibbles
  unnest(facet_value, keep_empty = TRUE) |> # Expand nested tibbles for each facet's directions
  mutate(coef = purrr::map(direction_value, ~ .x |>  pluck("coefs") |> bind_rows() |> rownames_to_column("predictor")),.keep="unused") |> # Extract and bind coefficients for each direction into a single data frame
  unnest(coef) |> # Expand coefficients into rows
  mutate(predictor = str_remove(predictor,"\\.{3}.*")) |> 
  filter(grepl("hfp|het", predictor)) |> 
  mutate(
    beta_folder = basename(dirname(dataset)),       # e.g., "Baselga_abun"
    dataset = gsub("_bbgdm","",basename(dataset)),                    # e.g., "barbaro_birds_1"
    beta_type = str_extract(beta_folder, "^(Baselga|Podani)"),
    metric_type = str_extract(beta_folder, "(abun|pa)$"),
    direction = gsub("z","s",direction)
  ) 



library("furrr")
plan(multisession, workers = 8)



synthesis_data<-models_collected |> 
  right_join(
    combined_data |> 
      select(dataset,facet, beta_type,metric_type, direction,harrel_davis,
             hfp_min,hfp_median,hfp_max,
             het_min,het_median,het_max)) |> 
  separate(predictor, into=c("predictor","buffer")) |> 
  drop_na() |> 
  select(-buffer) |>
  dplyr::group_split(dataset,facet,beta_type,metric_type,direction) |> 
  furrr::future_map( 
    .f = function(df) {
      # 1) row‐wise: classify & build list‐columns
      df <- df %>%
        rowwise() %>%
        mutate(
          shape = classify_predictor_row(
            predictor, coefficient.1, coefficient.2, coefficient.3,
            hfp_min, hfp_median, hfp_max,
            het_min, het_median, het_max,
            rel_tol = 0.05
          ),
          props = list(pct_by_props_predictor(
            predictor, coefficient.1, coefficient.2, coefficient.3,
            hfp_min, hfp_median, hfp_max,
            het_min, het_median, het_max
          )),
          bands = list(pct_by_bands_predictor(
            predictor, coefficient.1, coefficient.2, coefficient.3,
            hfp_min, hfp_median, hfp_max,
            het_min, het_median, het_max
          ))
        ) %>%
        ungroup()
      
      # 2) expand the list‐columns
      df <- df %>%
        tidyr::unnest_wider(props) %>%
        tidyr::unnest_wider(bands)
      
      # 3) compute raw counts & modal shape per predictor×group
      counts_modal <- df %>%
        group_by(dataset, facet, direction, beta_type, metric_type, predictor) %>%
        summarise(
          n_Absent      = sum(shape == "Absent"),
          n_Linear      = sum(shape == "Linear"),
          n_Exponential = sum(shape == "Exponential"),
          n_Saturating  = sum(shape == "Saturating"),
          n_Revlog      = sum(shape == "Revlog"),
          n_Uncertain   = sum(shape == "Uncertain"),
          .groups = "drop"
        ) %>%
        rowwise() %>%
        mutate(
          modal_shape = case_match(
            which.max(c_across(n_Absent:n_Uncertain)),
            1 ~ "Absent", 2 ~ "Linear", 3 ~ "Exponential",
            4 ~ "Saturating", 5 ~ "Revlog", 6 ~ "Uncertain"
          )
        ) %>%
        ungroup()
      
      # 4) filter the full df to only those rows matching the modal shape
      df_modal <- dplyr::inner_join(
        df,
        counts_modal,
        by = c("dataset","facet","direction","beta_type","metric_type","predictor")
      ) %>%
        filter(shape == modal_shape)
      
      # 5) compute medians *only* on the modal‐shape draws
      medians <- df_modal %>%
        group_by(dataset, facet, direction, beta_type, metric_type, predictor) %>%
        summarise(
          q_25  = median(q_25,  na.rm = TRUE),
          q_50  = median(q_50,  na.rm = TRUE),
          q_75  = median(q_75,  na.rm = TRUE),
          q_100 = median(q_100, na.rm = TRUE),
          b1    = median(b1,    na.rm = TRUE),
          b2    = median(b2,    na.rm = TRUE),
          b3    = median(b3,    na.rm = TRUE),
          b4    = median(b4,    na.rm = TRUE),
          .groups = "drop"
        )
      
      # 6) join counts+modal with medians, then pivot-wider
      summary <- dplyr::left_join(counts_modal, medians,
                                  by = c("dataset","facet","direction",
                                         "beta_type","metric_type","predictor"))
      
      summary %>%
        tidyr::pivot_wider(
          names_from  = predictor,
          values_from = c(q_25, q_50, q_75, q_100,
                          b1, b2, b3, b4,
                          n_Absent, n_Linear, n_Exponential,
                          n_Saturating, n_Revlog, n_Uncertain,
                          modal_shape),
          names_glue  = "{predictor}_{.value}"
        )
    },
    .options = furrr_options(packages = c("dplyr","tidyr","gdm")),.progress = TRUE
  ) %>%
  bind_rows() %>%
  relocate(dataset, facet, direction, beta_type, metric_type) %>%
  relocate(matches("^hfp_"), .after = metric_type) %>%
  relocate(matches("^het_"), .after = last_col())

final_data<- synthesis_data |> 
  left_join(combined_data) 

writexl::write_xlsx(data.frame(final_data), here::here("S6_Synthesis_model","data","synthesis_data_relative_magnitude.xlsx"))
