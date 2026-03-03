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

files_all <- list.files("S4_run_BBGDM/bbgdm_output", recursive = TRUE, full.names = TRUE)

files <- tools::file_path_sans_ext(files_all)

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

basic_structure <- tibble(fullpath = files) %>%                   # Use your subset (or all) of files
  mutate(
    dataset = basename(fullpath),                    # e.g., "barbaro_birds_1"
    beta_folder = basename(dirname(fullpath)),       # e.g., "Baselga_abun"
    beta_type = str_extract(beta_folder, "^(Baselga|Podani)"),
    metric_type = str_extract(beta_folder, "(abun|pa)$")
  ) |> 
  select(dataset, beta_type, metric_type)
#'---------------------------------------------------------------------------------------------------------------
#  synthesis BBGDM ###########################################################################################
#'---------------------------------------------------------------------------------------------------------------

# Extract the response features for all datasets -----------------------------------------------------------------------
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
  group_by(dataset,facet) |> # Group data by dataset and facet
  summarise(direction_r2 = gammaEffectSize(y=r2[direction == "homogenization"],x=r2[direction == "differentiation"],prob=.5)) # Compute gamma effect size between homogenization and differentiation R² values

mean_effect_size <- purrr::map(files, ~ readr::read_rds(glue::glue("{.x}.rds"))) |> # Read BBGDM results for each file
  set_names(gsub("_bbgdm", "", files)) |> 
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

shape <- purrr::map(files, ~ readr::read_rds(glue::glue("{.x}.rds"))) |> # Read BBGDM results for each file
  set_names(gsub("_bbgdm", "", files)) |> 
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

species_number <- purrr::map(gsub("_bbgdm", "", unique(basic_structure$dataset)), # Remove the "_bbgdm" suffix from dataset filenames.
                             ~ readr::read_rds(glue::glue( # Load the processed data file for each dataset.
                               "S1_Preprocessing/Processed/{.x}_clean.rds"
                             )) |>
                               pluck("comm") |> # Extract the "comm" (community matrix) component.
                               ncol()) |> # Count the number of species (columns in the community matrix).
  set_names(gsub("_bbgdm", "", unique(basic_structure$dataset))) |> # Assign dataset names to the resulting list.
  bind_rows() |>  # Combine the results into a single data frame.
  pivot_longer(everything(), names_to="dataset",values_to = "species.number") 

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

spatial.summary <- purrr::map(gsub("_bbgdm", "", unique(basic_structure$dataset)), # Remove the "_bbgdm" suffix from dataset filenames.
                              ~ readr::read_rds(glue::glue( # Load the processed data file for each dataset.
                                "S1_Preprocessing/Processed/{.x}_clean.rds"
                              )) |>
                                pluck("coord") |> # Extract the "coord" (coordinates) component.
                                drop_na(x, y)) |> # Remove rows with missing x or y coordinates.
  purrr::map(~ summary(fossil::earth.dist(data.frame(.x$x, .x$y)))[c(1, 4, 6)]) |> # Calculate min, mean, and max distances between points.
  set_names(gsub("_bbgdm", "", unique(basic_structure$dataset))) |> # Assign dataset names to the results.
  bind_rows(.id="dataset") |> # Combine the results into a single data frame with dataset names.
  set_names(c("dataset", "spatial.min", "spatial.mean", "spatial.max")) |> # Rename columns for clarity.
  dplyr::mutate(across(contains("spatial"), as.numeric))
# Get mean latitude value of the dataset.

latitude.mean <- purrr::map(gsub("_bbgdm", "", unique(basic_structure$dataset)), # Remove the "_bbgdm" suffix from dataset filenames.
                            ~ readr::read_rds(glue::glue( # Load the processed data file for each dataset.
                              "S1_Preprocessing/Processed/{.x}_clean.rds"
                            )) |>
                              pluck("coord") |> # Extract the "coord" (coordinates) component.
                              select(x, y)) |> # Select x and y columns for the coordinates.
  set_names(gsub("_bbgdm", "", unique(basic_structure$dataset))) |> # Assign dataset names to the results.
  bind_rows(.id="dataset") |> # Combine the results into a single data frame with dataset names.
  group_by(dataset) |> # Group by dataset for summary calculation.
  summarise(latitude.mean = mean(y, na.rm=TRUE)) 

hfp_min_max <- 
  purrr::map(gsub("_bbgdm", "", unique(basic_structure$dataset)), # Remove the "_bbgdm" suffix from dataset filenames.
             ~ readr::read_rds(glue::glue( # Load the processed data file for each dataset.
               "S1_Preprocessing/Processed/{.x}_clean.rds"
             )) |>
               pluck("predictors") |> # Extract the "predictors" component of the dataset.
               summarise(across(contains("hfp"), # Select all columns containing "hfp" (Human Footprint data).
                                list(min=~min(.x, na.rm=TRUE), # Calculate the minimum value for each "hfp" column.
                                     mean=~mean(.x, na.rm=TRUE), # Calculate the mean value for each "hfp" column.
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
  pivot_wider(names_from="fun", values_from = "value", names_prefix = "hfp.") |> # Reshape data back into wide format with separate columns for min, mean, and max.
  separate(predictor, into=c("predictor", "buffer")) |> # Further split the "predictor" column into the actual predictor and its buffer.
  mutate(hfp.range = hfp.max - hfp.min) |> # Calculate the range of Human Footprint values (max - min).
  select(-predictor) # Drop the "predictor" column, leaving only relevant data.

synthesis_data <- 
  mean_r2 |> 
  mutate(direction = ifelse(sign(direction_r2) > 0 , "homogenization", "differentiation")) |> # Assign "homogenization" or "differentiation" based on the sign of `direction_r2`.
 # relocate(dataset, facet, direction) |> # Reorganize columns to bring `dataset`, `facet`, and `direction` to the front.
  ungroup() |> # Remove grouping to allow subsequent operations to act on the full dataset.
  left_join(mean_effect_size, by = c("dataset", "facet", "direction")) |> # Merge effect size data on matching `dataset`, `facet`, and `direction`.
  left_join(shape, by = c("dataset", "facet", "direction")) |> # Merge shape data on matching `dataset`, `facet`, and `direction`.
  mutate(
    beta_folder = basename(dirname(dataset)),       # e.g., "Baselga_abun"
    beta_type = str_extract(beta_folder, "^(Baselga|Podani)"),
    metric_type = str_extract(beta_folder, "(abun|pa)$")) |> 
  mutate(dataset = basename(dataset)) |> 
  left_join(hfp_min_max, by = c("dataset", "buffer")) |> # Merge Human Footprint summary stats on matching `dataset` and `buffer`.
  left_join(species_number, by = c("dataset")) |> # Add species number data for each dataset.
  left_join(spatial.area, by = c("dataset")) |> # Add spatial extent (area) data for each dataset.
  left_join(spatial.summary, by = c("dataset")) |> # Add spatial summary stats (min, mean, max distances) for each dataset.
  left_join(latitude.mean, by = c("dataset")) |> # Add mean latitude data for each dataset.
  left_join(dataset_features |> distinct(dataset, disturbance, taxa, group, system), by = c("dataset")) |> # Add metadata (disturbance, taxa, biotic group, system) for each dataset.
  mutate(direction = gsub("z", "s", direction)) |> # Replace "z" with "s" in the `direction` column (e.g., fix typos or conventions).
  rename(realm = system, biotic.group = group) 
  

# Rename columns for consistency: `system` to `realm`, `group` to `biotic.group`.

#writexl::write_xlsx(data.frame(synthesis_data), here::here("S6_Synthesis_model","data","synthesis_data.xlsx"))
# library(ggtext)
# library(ggh4x)
# annotations_df <- data.frame(
#   label = c(
#     "<span>&larr; Differentiation</span>",
#     "<span>Homogenisation &rarr;</span>",
#     "<span>&larr; Differentiation</span>",
#     "<span>Homogenisation &rarr;</span>"
#   ),
#   x = c(-0.03, 0.03, -7, -7),
#   y = c(-10, -10, -0.2, 0.2),
#   angle = c(0, 0, 90, 90),
#   hjust = c(1, 0, 1, 0))
# 
# synthesis_data |> 
#   select(dataset,facet,beta_type,metric_type,direction_r2) |> 
#   pivot_wider(values_from="direction_r2", names_from="facet") |> 
#   ggplot(aes(y=Taxonomic,x=Functional)) +
#   geom_point(size=2,alpha=0.7)+
#   facet_wrap(beta_type~metric_type) +
#   geom_vline(xintercept = 0)+
#   geom_hline(yintercept = 0)+
#   #coord_cartesian(xlim = c(-7, NA), ylim = c(-10, 7)) +
#   geom_richtext(
#     data = annotations_df,
#     aes(x = x, y = y, label = label, angle = angle, hjust = hjust), 
#     size=convert_size(7),
#     label.color = NA
#   ) +
#   geom_hline(yintercept = 0,
#              linetype = 3,
#              linewidth = .1) +
#   geom_vline(xintercept = 0,
#              linetype = 3,
#              linewidth = .1) +
#   ggplot2::theme_void(base_family = "sans", base_size = 8) +
#   ggplot2::theme(
#     panel.grid.minor = ggplot2::element_blank(),
#     plot.background = ggplot2::element_rect(fill = "white", color = NA),
#     strip.background = ggplot2::element_rect(fill = NA, color = NA),
#     strip.text = ggplot2::element_text(
#       face = "bold",
#       margin = margin(b = 2)
#     ),
#     axis.title.x = ggtext::element_markdown(
#       family = "sans",
#       face = "bold",
#       hjust = 0.5,
#       margin = margin(t = 2)
#     ),
#     axis.title.y = ggtext::element_markdown(
#       family = "sans",
#       face = "bold",
#       angle = 90,
#       margin = margin(r = 4)
#     ),
#     axis.text.x = element_markdown(
#       family = "sans",
#       color = "grey30",
#       margin = margin(t = 2)
#     ),
#     axis.text.y = element_markdown(
#       family = "sans",
#       color = "grey30",
#     ),
#     panel.spacing.x = unit(0.2, "lines"),
#     panel.spacing.y = unit(0.2, "lines"),
#     axis.line.x = element_line(color = "grey70", linewidth = 0.1),
#     axis.ticks.x = element_line(color = "grey70", linewidth = 0.1),
#     axis.ticks.length.x = unit(0.1, "lines"),
#     axis.line.y = element_line(color = "grey70", linewidth = 0.1),
#     axis.ticks.y = element_line(color = "grey70", linewidth = 0.1),
#     axis.ticks.length.y = unit(0.1, "lines"),
#     plot.margin =  margin(3, 3, 3, 3),
#     legend.position = "none",
#   ) +
#   labs(x = "Trait replacement [Harrel-Davis]", y = "Species replacement [Harrel-Davis]") +
#   scale_x_continuous(breaks=c(-6.5, -4.5, -2.5, 0, 2.5))+
#   guides(
#     x = guide_axis_truncated(trunc_lower = -6.5, trunc_upper = 2.5),
#     y = guide_axis_truncated(trunc_lower = -10, trunc_upper = 5)
#   )
# 
# foo<-synthesis_data |> 
#   filter(beta_type == "Podani", metric_type == "abun")
# 
# foo2<-synthesis_old |> 
#   filter(dataset %in% unique(foo$dataset))
# 
# synthesis_data |> 
# select(dataset,facet,beta_type,metric_type, magnitude) |> 
#   pivot_wider(names_from=beta_type,values_from=magnitude) |> 
#   ggplot(aes(y=log10(Baselga+0.01),x=log10(Podani+0.01))) +
#   geom_point() +
#   geom_smooth(method="lm")+
#   facet_grid(facet~metric_type)+
#   labs(title="Magnitude of effect")
# 
# ggplot(aes(x=direction,y=magnitude))+
#   geom_boxplot() +
#   facet_grid(facet~beta_type*metric_type)+
#   theme_bw()
# 
# 


# homog <- mean_r2 |> 
#   filter(id == "N25FBD", 
#          beta_type == "Podani", 
#          metric_type == "abun", 
#          facet == "Functional") |> 
#   pull(homogenization)
# 
# diffe <- mean_r2 |> 
#   filter(id == "N25FBD", 
#          beta_type == "Podani", 
#          metric_type == "pa", 
#          facet == "Functional") |> 
#   pull(differentiation)