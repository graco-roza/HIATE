###############################################################################
# SCRIPT NAME: Helper Functions for Synthesis Analysis
#
# DESCRIPTION:
#   This script contains helper functions to process, analyze, and visualize
#   synthesis data for the study. These functions include data wrangling, 
#   statistical computations, and plotting routines, tailored to support the 
#   synthesis of beta-diversity metrics across datasets.
#
# USAGE:
#   - This script is meant to be sourced by other scripts in the synthesis 
#     analysis workflow.
#   - Functions defined here can be called directly in analysis scripts for 
#     specific tasks such as extracting predictions, calculating effect sizes,
#     or plotting results.
#
# INPUTS:
#   - Processed synthesis datasets in .xlsx format.
#   - Beta-diversity models and predictions for taxonomic and functional facets.
#
# OUTPUTS:
#   - Data transformations, statistical summaries, and visualizations, including:
#     - Direction, magnitude, and shape of beta-diversity responses.
#     - Plots of relationships between human pressures and biodiversity metrics.
#
# AUTHOR: Caio Graco-Roza
# LAST UPDATED: 2024-11-24
#
# NOTES:
#   - Ensure all required libraries are installed and the input files are 
#     correctly formatted before sourcing this script.
#   - Functions are modular and can be reused or extended for related analyses.
###############################################################################

require(Hmisc) # For advanced statistical functions, including bootstrapped confidence intervals.
require(dplyr) # For data manipulation, including filtering, renaming, joining, and summarizing datasets.
require(readxl) # For reading Excel files, particularly synthesis data in .xlsx format.
require(tidyr) # For reshaping data, such as pivoting between long and wide formats.
require(glue) # For creating dynamic strings, especially for file paths and labeling.

#' Process Model Data
#'
#' Processes synthesis data for analysis, standardizing columns and filtering datasets.
#'
#' @param data_path Path to the synthesis data file (Excel format).
#' @param model_type The type of model to process. Options: "direction", "magnitude", "shape", "convergence".
#' @return A processed tibble ready for modeling or analysis.
#' @export
process_model_data <- function(data_path, model_type) {
  # data_path="S7_Synthesis_model/data/synthesis_data.xlsx"
  # model_type = "convergence"
  # Define base column renaming for consistency with text labels
  base_col_rename <- c(
    "species_number" = "species.number",
    "spatial_distance_min" = "spatial.min",
    "spatial_extent" = "spatial.extent",
    "absolute_latitude_mean" = "latitude.mean",
    "direction_harrel_davis" = "direction_r2",
    "taxa_coarse" = "biotic.group",
    "taxa_fine" = "taxa",
    "ecosystem_type" = "realm",
    "main_land_use_type" = "disturbance",
    "buffer" = "buffer",
    "human_pressure_baseline" = "hfp.min",
    "human_pressure_max" = "hfp.max",
    "human_pressure_range" = "hfp.range"
  )
  
  # Define model-specific columns to select based on model type
  specific_cols <- switch(
    model_type,
    "direction" = c("direction_harrel_davis"),
    "magnitude" = c("magnitude"),
    "shape" = c("Absent", "Exponential", "Saturating", "Revlog"),
    "convergence" = c("magnitude", "direction_harrel_davis", "Absent", "Exponential", "Saturating", "Revlog")
  )
  
  # Read and process data
  processed_data <- data_path |>
    readxl::read_xlsx() |>
    dplyr::filter(!dataset %in% c("N67TTP", "N78TTP", "S47TTP")) |> 
    dplyr::rename(all_of(base_col_rename)) |> 
    dplyr::select(
      dataset, facet, ecosystem_type, taxa_coarse, taxa_fine, buffer, species_number,
      spatial_extent, spatial_distance_min, absolute_latitude_mean, direction,
      human_pressure_baseline, human_pressure_max, human_pressure_range, main_land_use_type, all_of(specific_cols)
    ) |>
    
    # Apply general transformations and formatting
    dplyr::mutate(
      facet = dplyr::recode(facet, "Taxonomic" = "Species", "Functional" = "Traits"),
      facet = factor(facet, levels = c("Species", "Traits")),
      ecosystem_type = gsub("aquatic", "freshwater", ecosystem_type),
      ecosystem_type = factor(str_to_title(ecosystem_type), levels = c("Terrestrial", "Freshwater")),
      taxa_coarse = factor(taxa_coarse, levels = c("invertebrate", "vertebrate", "plant", "microorganism")),
      main_land_use_type = factor(str_to_title(main_land_use_type), levels = c("Multiple", "Agriculture", "Forest", "Urban")),
      species_number = log10(species_number + 1),
      spatial_extent = log10(spatial_extent + 1),
      spatial_distance_min = log10(spatial_distance_min + 1),
      human_pressure_range = human_pressure_max - human_pressure_baseline,
      absolute_latitude_mean = abs(absolute_latitude_mean),
      direction = gsub("z", "s", direction),
      direction = factor(str_to_title(direction), levels = c("Differentiation", "Homogenisation")),
      buffer = factor(buffer, levels = c("1000", "1500", "2000"))
    ) |>
    
    # Standardize numeric columns
    dplyr::mutate(dplyr::across(
      c("species_number", "spatial_distance_min", "spatial_extent", 
        "absolute_latitude_mean", "human_pressure_baseline", "human_pressure_range"), 
      ~ (.x - mean(.x, na.rm = TRUE)) / sd(.x, na.rm = TRUE)
    )) |> 
    select(-human_pressure_max) |> 
    janitor::clean_names()
  
  # Calculate the most frequent shape for each dataset and facet
  if (model_type == "convergence") {
    shape_data <- processed_data |>
      dplyr::select(dataset, facet, absent, exponential, saturating, revlog) |>
      tidyr::pivot_longer(cols = c(absent, exponential, saturating, revlog), names_to = "shape", values_to = "count") |>
      dplyr::group_by(dataset, facet) |>
      dplyr::summarize(shape = shape[which.max(count)], .groups = "drop") |> 
      dplyr::mutate(shape = factor(str_to_title(shape), levels = c("Absent", "Exponential", "Saturating", "Revlog")))
    
    # Merge the shape column back into the main data
    processed_data <- processed_data |>
      dplyr::left_join(shape_data, by = c("dataset", "facet"))
  }

  if (model_type == "shape") {
    processed_data <- processed_data |> 
      mutate(shape = cbind(absent, exponential, saturating, revlog),
             trials = absent + exponential + saturating + revlog) 
  }
  
  return(processed_data)
}


#' Count Datasets Matching a Direction and Facet
#'
#' Counts the number of distinct datasets for a specific facet and direction.
#'
#' @param df A data frame containing the synthesis data.
#' @param facet The facet to filter (e.g., "Species", "Traits").
#' @param direction The direction to filter (e.g., "Homogenisation").
#' @return An integer count of datasets matching the criteria.
#' @export
count_directions <- function(df, facet, direction) {
  df %>% 
    filter(facet == !!facet, direction == !!direction) %>% 
    distinct(dataset) %>% 
    nrow()
}

#' Calculate Effect Size
#'
#' Calculates the bootstrapped mean and confidence intervals for effect size.
#'
#' @param bbgdm_summary A data frame summarizing BBGDM results.
#' @param facet The facet to filter (e.g., "Species", "Traits").
#' @param direction The direction to filter (e.g., "Homogenisation").
#' @return A list containing the mean and confidence intervals.
#' @export
calculate_effect_size <- function(bbgdm_summary, facet, direction) {
  result <- bbgdm_summary %>%
    group_by(dataset,facet) %>% 
    mutate(effect_size = magnitude/sum(magnitude)) %>% 
    filter(predictor == "hfp") %>% 
    select(facet,direction,effect_size) %>% 
    filter(facet == !!facet, direction == !!direction) %>%
    pull(effect_size) %>% 
    Hmisc::smean.cl.boot(conf.int=.95, B=1000, na.rm=TRUE, reps=FALSE)
  
  return(result)
}

#' Convert Text Size
#'
#' Converts ggplot2 text size to physical dimensions for consistent scaling.
#'
#' @param size Numeric size in ggplot2 units.
#' @return The converted size.
#' @export
convert_size <- function(size) {
  size <- size * 0.3528 
  return(size)
}

#' Pooled Variance Calculation
#'
#' Computes pooled variance across two groups using a custom function.
#'
#' @param x Numeric vector of group 1.
#' @param y Numeric vector of group 2.
#' @param FUN Function to compute variance (e.g., `var` or `hdmad`).
#' @return Pooled variance.
#' @export
pooled <- function(x, y, FUN) {
  nx <- length(x)
  ny <- length(y)
  sqrt(((nx - 1) * FUN(x) ^ 2 + (ny - 1) * FUN(y) ^ 2) / (nx + ny - 2))
}

#' Harrell-Davis Median
#'
#' Computes the Harrell-Davis median of a numeric vector.
#'
#' @param x Numeric vector.
#' @return The median value.
#' @export
hdmedian <- function(x) as.numeric(hdquantile(x, 0.5))

#' Harrell-Davis MAD
#'
#' Computes the median absolute deviation using the Harrell-Davis method.
#'
#' @param x Numeric vector.
#' @return The MAD value.
#' @export
hdmad <- function(x) 1.4826 * hdmedian(abs(x - hdmedian(x)))

#' Pooled Harrell-Davis MAD
#'
#' Computes the pooled MAD for two numeric vectors using the Harrell-Davis method.
#'
#' @param x Numeric vector of group 1.
#' @param y Numeric vector of group 2.
#' @return The pooled MAD value.
#' @export
phdmad <- function(x, y) pooled(x, y, hdmad)

#' Gamma Effect Size
#'
#' Computes the gamma effect size using Harrell-Davis quantiles.
#'
#' @param x Numeric vector for group 1.
#' @param y Numeric vector for group 2.
#' @param prob Quantile probability for comparison.
#' @return The gamma effect size.
#' @export
gammaEffectSize <- function(x, y, prob){
  if(length(na.exclude(y)) < 300) y<-c(na.exclude(y),rep(0,length(na.exclude(x)) - length(na.exclude(y))))
  if(length(na.exclude(x)) < 300) x<-c(na.exclude(x),rep(0,length(na.exclude(y)) - length(na.exclude(x))))
  res<-as.numeric((hdquantile(na.exclude(y), prob) - hdquantile(na.exclude(x), prob)) / phdmad(na.exclude(x), na.exclude(y)))
  return(res)
}

#' Extract Predictions
#'
#' Extracts posterior predictions for species and trait models.
#'
#' @param species_model Bayesian model for species replacement.
#' @param trait_model Bayesian model for trait replacement.
#' @param predictor The predictor variable for marginal effects.
#' @param ndraws Number of posterior draws.
#' @param levels_list Optional levels for factor predictors.
#' @return A tibble with posterior predictions for each facet.
#' @export
extract_predictions <- function(species_model, trait_model, predictor, ndraws = 1000, levels_list = NULL) {
  
  # Generate predictions for taxonomic facet (Species replacement)
  species_pred <- predictions(species_model, by = predictor, ndraws = ndraws) |>
    marginaleffects::posteriordraws() |>
    mutate(facet = "Species replacement")  # Directly assign facet label
  
  # Generate predictions for functional facet (Trait replacement)
  trait_pred <- predictions(trait_model, by = predictor, ndraws = ndraws) |>
    marginaleffects::posteriordraws() |>
    mutate(facet = "Trait replacement")  # Directly assign facet label
  
  # Combine both predictions
  pred_combined <- bind_rows(species_pred, trait_pred) |>
    dplyr::select(draw, facet, !!predictor) |>
    dplyr::mutate(facet = factor(facet, levels = c("Species replacement", "Trait replacement")))  # Use final facet levels
  
  # If levels are provided, ensure factors are set correctly
  if (!is.null(levels_list)) {
    pred_combined <- pred_combined |> dplyr::mutate(across({{predictor}}, ~ factor(stringr::str_to_sentence(.x), levels = levels_list)))
  }
  
  return(pred_combined)
}

#' Plot Direction of Relationships
#'
#' Creates a plot showing the direction of relationships for each predictor.
#'
#' @param pred_data Posterior prediction data.
#' @param model_data Synthesis model data.
#' @param predictor_col Column name for predictors.
#' @param color_map Color mapping for predictors.
#' @param x_limits Limits for the x-axis.
#' @param vjust_text Vertical adjustment for text annotations.
#' @param height Height of the slabs.
#' @param justification Horizontal adjustment for slabs.
#' @param width Width of the point intervals.
#' @return A ggplot object.
#' @export
plot_direction <- function(pred_data, model_data, predictor_col, color_map, x_limits = c(-3, 3), vjust_text = 0.75, height=1, justification = 0.15, width=1) {
  
  
  theme_direction <- function() {
    ggplot2::theme_void(base_family = "sans", base_size = 3) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        plot.background = ggplot2::element_rect(fill = "white", color = NA),
        strip.background = ggplot2::element_rect(fill = NA, color = NA),
        strip.text = ggplot2::element_text(size = 4, face = "bold", margin=margin(b=2)),
        axis.title.x = ggtext::element_markdown(family = "sans", face = "bold", hjust = 0.5, size = 3, margin=margin(t=2)),
        axis.title.y = element_blank(),
        axis.text.x = element_markdown(family = "sans", color = "grey30", size = 3, margin = margin(t = 2)),
        axis.text.y = element_blank(),
        panel.spacing.x = unit(0.2, "lines"),
        panel.spacing.y = unit(0.2, "lines"),
        axis.line.x = element_line(color = "grey70", linewidth=0.1),
        axis.ticks.x = element_line(color = "grey70", linewidth=0.1),
        axis.ticks.length.x = unit(0.1, "lines"),
        plot.margin = margin(2, 0, 2, 0),
        legend.position="none",
        plot.tag = element_text(face="bold",size=5)
      )
  }
  
  # Helper function to add annotations (arrows and text)
  add_direction_annotations <- function() {
    list(
      annotate("text", family = "sans", label = "Differentiation", size = convert_size(2.5), x = -0.2, y = 1.5, hjust = 1, col = "gray30"),
      annotate("text", family = "sans", label = "Homogenisation", size = convert_size(2.5), x = 0.2, y = 1.5, hjust = 0, col = "gray30"),
      annotate("segment", linewidth = 0.1, x = -1.4, xend = -2, y = 1.5, yend = 1.5, arrow = arrow(type = "closed", length = unit(0.02, "inches")), col = "gray30"),
      annotate("segment", linewidth = 0.1, x = 1.6, xend = 2.2, y = 1.5, yend = 1.5, arrow = arrow(type = "closed", length = unit(0.02, "inches")), col = "gray30")
    )
  }
  
  # Helper function to extract dataset counts for each category, split by facet
  get_data_counts <- function(model_data, predictor_col, color_map) {
    # Count the data and organize it by `predictor_col` and `facet`
    counts <- model_data %>%
      group_by(facet, !!sym(predictor_col)) %>%
      count(!!sym(predictor_col)) %>%
      ungroup() 
    
    # Use the provided `color_map` if it’s supplied; otherwise, default to factor levels in `predictor_col`
    if (is.null(color_map)) {
      unique_levels <- levels(pred_data[[predictor_col]])  # Unique levels as factors
      default_colors <- RColorBrewer::brewer.pal(length(unique_levels), "Set1")  # Create colors if none provided
      color_map <- setNames(default_colors, unique_levels)  # Map colors to levels
    } else {
      # Ensure color mapping aligns with levels order in predictor column
      color_map <- color_map[levels(pred_data[[predictor_col]])]
    }
    
    # Map sample sizes using the counts
    sample_sizes <- setNames(counts$n, counts[[predictor_col]])
    
    list(color_map = color_map, sample_sizes = sample_sizes)
  }
  
  # Get counts data
  counts_data <- get_data_counts(model_data, predictor_col, color_map)
  
  # Helper function to prepare text for rich text annotations
  prepare_text_data <- function(draft_plot, model_data, color_map, sample_sizes) {

    
    text_prep <- ggplot2::layer_data(draft_plot, 3) |>
      dplyr::select(-dist) |>
      dplyr::filter(.width == 0.95) |>
      dplyr::mutate(
        # Assign facet variable dynamically based on PANEL
        facet = ifelse(PANEL == 1, "Species replacement", "Trait replacement"),
        type = purrr::map_chr(colour, ~ {
          category <- names(color_map)[color_map == .x]  # Match color to category
          sample_size <- sample_sizes[category]  # Retrieve sample size for category
          glue::glue("{category} (n = {sample_size})")  # Generate label
        })
      ) |>
      dplyr::mutate(label = glue::glue("<b style='color:{colour};'>{type}</b>"))
    
    return(text_prep)
  }
  
  # Generate draft plot
  draft_plot <- 
    ggplot(pred_data, aes(x = draw)) +
    geom_vline(xintercept = 0, linewidth = 0.1, linetype = "11", colour="gray60") +
    facet_wrap(~facet, ncol = 1, scales = "free") +
    labs(x = "Direction of Relationship (Harrell-Davis)", colour = "") +
    scale_x_continuous(limits = x_limits, breaks = pretty(seq(-2, 2, 0.2))) +
    scale_fill_manual(values = color_map) +
    scale_colour_manual(values = color_map) +
    ggdist::stat_slab(aes(fill = !!sym(predictor_col)), justification = {justification}, colour = "white", alpha = 0.5, height = {height}, size = 0.1, show.legend = FALSE) +
    ggdist::stat_pointinterval(aes(colour = !!sym(predictor_col)), point_interval = "median_hdci", .width = c(0.8, 0.9, 0.95), interval_size_range = c(0.1, 0.4), fatten_point = 0.1, position = position_dodge(width = {width}, preserve = "single")) +
    add_direction_annotations() +
    theme_direction()
  
  # Prepare text data for rich text annotations
  text_data <- prepare_text_data(
    draft_plot,
    model_data = model_data,
    color_map = counts_data$color_map,
    sample_sizes = counts_data$sample_sizes
  )
  
  # Add text annotations to draft plot
  final_plot <- draft_plot +
    geom_richtext(
      data = text_data,
      aes(y = ymax, x = xmax, label = label),
      hjust = 0.05,
      vjust = vjust_text,
      colour = NA,
      fill = NA,
      size = convert_size(2.5)
    )+
    guides(x = guide_axis_truncated(trunc_lower = -2,trunc_upper = 2)) 
  
  return(final_plot)
}

#' Plot Magnitude of Effects
#'
#' Creates a plot showing the magnitude of effects for predictors.
#'
#' @param pred_data Posterior prediction data.
#' @param colors Color mapping for predictors.
#' @param x_limits Limits for the x-axis.
#' @param x_breaks Break points for the x-axis.
#' @param predictor Column name for predictors.
#' @param facet_var Facet variable for the plot.
#' @param sample_sizes Sample sizes for each predictor level.
#' @param text_height Vertical adjustment for text annotations.
#' @param justification Horizontal adjustment for slabs.
#' @param height Height of the slabs.
#' @return A ggplot object.
#' @export
plot_magnitude <- function(pred_data, colors, x_limits, x_breaks, predictor, facet_var, sample_sizes, text_height = 1, justification, height){
  
  # Define a reusable base theme with common settings
  my_base_theme <- function() {
    ggplot2::theme_void(base_family = "sans", base_size = 3) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        plot.background = ggplot2::element_rect(fill = "white", color = NA),
        strip.background = ggplot2::element_rect(fill = NA, color = NA),
        strip.text = ggplot2::element_text(size = 4, face = "bold", margin=margin(b=2)),
        axis.title.x = ggtext::element_markdown(family = "sans", face = "bold", hjust = 0.5, size = 3, margin=margin(t=2)),
        axis.text.x = element_markdown(family = "sans", color = "grey30", size = 3, margin = margin(t = 2)),
        axis.text.y = element_blank(),
        panel.spacing.x = unit(0.2, "lines"),
        panel.spacing.y = unit(0.2, "lines"),
        axis.line.x = element_line(color = "grey70", linewidth=0.1),
        axis.ticks.x = element_line(color = "grey70", linewidth=0.1),
        axis.ticks.length.x = unit(0.1, "lines"),
        plot.margin = margin(2, 0, 2, 0),
        legend.position="none",
        plot.tag = element_text(face="bold",size=5)
      )
  }
  
  prepare_text_data <- function(draft_plot, facet_var, color_map, sample_sizes) {
    ggplot2::layer_data(draft_plot, 2) |>
      dplyr::select(-dist) |>
      dplyr::filter(.width == 0.95) |>
      dplyr::mutate(
        facet = facet_var,
        type = purrr::map_chr(colour, ~ {
          category <- names(color_map)[color_map == .x]  # Get the category from the color
          sample_size <- sample_sizes[category]  # Match sample size based on category
          glue::glue("{category} (n = {sample_size})")  # Generate the text label
        })
      ) |>
      dplyr::mutate(label = glue::glue("<b style='color:{colour};'>{type}</b>"))
  }
  
  # text_data_3a1 <- prepare_text_data(
  #   draft_plot = figure_3a1_draft,  # Pass the draft plot to this function
  #   facet_var = "Species replacement",  # This is the facet value for this plot
  #   color_map = direction_colors,  # Color mapping for the plot
  #   sample_sizes = c(118, 42),  # Sample sizes for Differentiation and Homogenisation
  #   text_mapping = c("Differentiation", "Homogenisation")  # Labels for the two groups
  # )
  
  
  generate_draft_plot <- function(pred_data, colors, x_limits, x_breaks, predictor, justification, height) {
    
    ggplot2::ggplot(pred_data, ggplot2::aes(x = draw)) +
      ggplot2::facet_wrap(~facet, ncol = 1, scales = "free_y") +
      my_base_theme() +
      ggplot2::theme(
        strip.background = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        legend.position = "none",
        legend.direction = "horizontal",
        panel.grid.major.x = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(size = 5, face = "bold", hjust = .5),
        plot.title.position = "plot"
      ) +
      ggplot2::labs(x = "Magnitude of effect (weighted effect size)", colour = "") +
      ggplot2::scale_fill_manual(values = colors) +
      ggplot2::scale_colour_manual(values = colors) +
      ggplot2::scale_x_continuous(breaks = base::seq(x_breaks[1], x_breaks[2], abs(x_breaks[1]-x_breaks[2])*.2), limits = c(x_limits[1], x_limits[2])) +
      ggplot2::scale_y_discrete(breaks = NULL) +
      ggdist::stat_slab(aes(fill = !!sym(predictor), colour = !!sym(predictor)),
                        justification = !!justification, colour = "white", alpha = .5, height = {height}, size = .1, show.legend = FALSE) +
      ggdist::stat_pointinterval(ggplot2::aes(colour = !!sym(predictor)),
                                 point_interval = "median_hdci", .width = c(.8, .9, .95),
                                 interval_size_range = c(.1, .5), fatten_point = .1,
                                 position = position_dodge(width = 4, preserve = "single")) +
      guides(x = guide_axis_truncated(trunc_lower = x_breaks[1], trunc_upper = x_breaks[2])) 
  }
  
  
  # Draft plot for Species Replacement - Direction (facet_var is 'direction')
  # figure_3a1_draft <- generate_draft_plot(
  #   pred_data = pred_mag_dir |> filter(facet == "Species replacement"),  # Use the "Species replacement" facet
  #   colors = direction_colors,
  #   x_limits = c(0, 1.3),
  #   x_breaks = c(0, 1.5),
  #   predictor = "direction"  # Use 'direction' for this plot
  # )
  
  # # Draft plot for Species Replacement - Direction (facet_var is 'direction')
  draft_plot <- generate_draft_plot(
    pred_data = pred_data,  # Use the "Species replacement" facet
    colors = colors,
    x_limits =x_limits,
    x_breaks = x_breaks,
    predictor = predictor,  # Use 'direction' for this plot
    justification= justification,
    height = height
  )

  # Prepare the text to overlay on the plot
  text_data <- prepare_text_data(draft_plot, facet_var, colors, sample_sizes)
  # Add the text to the plot
  final_plot <- draft_plot +
  ggnewscale::new_scale_colour()  +
    geom_text(
      data = text_data, 
      aes(y = ymax, x = xmax, label = type, colour = colour), 
      hjust = -.05, vjust = text_height,  
      size.unit = "pt",
      size= 3,
      fontface="bold"
    ) + scale_colour_identity()
  
 return(final_plot)
}

#' Extract Shape Predictions
#'
#' Extracts predictions for shape models.
#'
#' @param species_model Bayesian model for species replacement.
#' @param trait_model Bayesian model for trait replacement.
#' @param predictor The predictor variable.
#' @param predictor_levels Levels for predictor factors.
#' @param facet_labels Labels for facets.
#' @return A tibble with shape predictions.
#' @export
extract_shape_predictions <- function(species_model, trait_model, predictor, predictor_levels, facet_labels = c("Species replacement", "Trait replacement")) {
  
  # Generate predictions for the "Species replacement" model
  pred_species <- epred_draws(
    species_model, 
    newdata = insight::get_datagrid(species_model, at = c(predictor, "trials=1"), preserve_range = FALSE, data = species_model$data, include_response = TRUE),
    re_formula = NA
  ) |>
    ungroup() |> 
    dplyr::mutate(facet = facet_labels[1])
  
  # Generate predictions for the "Trait replacement" model
  pred_trait <- epred_draws(
    trait_model, 
    newdata = insight::get_datagrid(trait_model, at = c(predictor, "trials=1"), preserve_range = FALSE, data = trait_model$data, include_response = TRUE),
    re_formula = NA
  ) |>
    ungroup() |> 
    dplyr::mutate(facet = facet_labels[2])
  
  # Combine predictions and process them
  predictions <- dplyr::bind_rows(pred_species, pred_trait) |>
    dplyr::select(.epred, facet, !!sym(predictor), .category) |>
    dplyr::mutate(
      facet = factor(facet, levels = facet_labels),
      !!sym(predictor) := factor(
        stringr::str_to_title(!!sym(predictor)), 
        levels = stringr::str_to_title(predictor_levels)
      )
    )
  
  return(predictions)
}

#' Plot Shape Scheme
#'
#' Creates a plot illustrating the shape scheme (e.g., Absent, Exponential, Saturating, Revlog) for species/trait replacement.
#'
#' @param shapes_data A data frame containing information about shapes, including coordinates for drawing rectangles.
#' @param text_data A data frame containing rich text annotations to overlay on the plot.
#' @param shape_colors A vector of colors corresponding to the shapes.
#' @param facet_labels A named vector or list for relabeling the facets (shapes).
#' @return A ggplot object representing the shape scheme with annotated text and graphical elements.
#' @export
plot_shape_scheme <- function(shapes_data, text_data, shape_colors, facet_labels) {
  ggplot(shapes_data) +
    aes(x = x, y = y, group = shape) +
    facet_wrap(~shape, ncol = 1, labeller = labeller(shape = facet_labels)) +
    geom_rect(aes(xmin = 0, xmax = 1, ymin = 0, ymax = 1.5), fill = shape_colors[1], colour = "white", alpha = .3) +
    geom_rect(aes(xmin = 1, xmax = 2, ymin = 0, ymax = 1.5), fill = shape_colors[2], colour = "white", alpha = .3) +
    geom_rect(aes(xmin = 2, xmax = 3, ymin = 0, ymax = 1.5), fill = shape_colors[3], colour = "white", alpha = .3) +
    ggtext::geom_richtext(
      data = text_data, 
      aes(x = x, y = y + .05, label = label),
      size = convert_size(3), colour = NA, fill = NA, text.colour = "black"
    ) +
    ggalt::geom_xspline(spline_shape = 1, size =0.5) +
    theme_shape_plot() +  # Use a base shape theme if you've defined one
    theme(
      axis.line.y = element_line(color = "grey70", linewidth = 0.1),
      axis.ticks.y = element_line(color = "grey70", linewidth = 0.1),
      axis.ticks.length.y = unit(0.1, "lines"),
      panel.spacing.y = unit(0, "lines"),
      strip.text = element_textbox_simple(halign = 0.5, size = convert_size(8), face = "bold", lineheight = 1.2)
    ) +
    labs(x = "Human pressure gradient", y = "Species/Trait replacement") +
    scale_y_continuous(breaks = seq(0, 1, .25), limits = c(0, 1.5)) +
    guides(x = guide_axis_truncated(trunc_lower = 0, trunc_upper = 3),
           y = guide_axis_truncated(trunc_lower = 0, trunc_upper = 1)) +
    scale_x_continuous(breaks = seq(0, 3, length.out = 5), labels = c(0, .25, .5, .75, 1), limits = c(0, 3)) 
}

#' Create Shape Plot
#'
#' Generates a shape plot showing predicted probabilities.
#'
#' @param data Prediction data.
#' @param predictor_label Label for the predictor axis.
#' @param fill_var Column used for fill color.
#' @param color_palette Color palette for the plot.
#' @return A ggplot object.
#' @export
plot_shape <- function(data, predictor_label, fill_var, color_palette) {

 data |> 
         mutate(predictor := {{predictor_label}},
                group := !!sym(fill_var),
                .category = gsub("Revlog","Reverse<br>logistic", .category)) |> 
  ggplot(aes(x = .epred, y = .category, fill = group, colour = group)) +
  ggh4x::facet_nested(facet ~ predictor + group) +
  stat_ccdfinterval(slab_alpha = .7, point_interval = "median_hdci", .width = c(.8, .9, .95),
                    fatten_point = .3, interval_size_range = c(.1, .4)) +
  scale_fill_manual(values = color_palette) +
  scale_colour_manual(values = colorspace::darken(color_palette, .3, space = "HLS")) +
  theme_shape_plot() +
  labs(y = "Relationship shape", x = "Predicted probability") +
  guides(x = guide_axis_truncated(trunc_lower = 0, trunc_upper = .75),
         y = guide_axis_truncated(trunc_lower = function(x) {x - 0.45},
                                  trunc_upper = function(x) {x + 0.45})) +
  scale_x_continuous(breaks = seq(0, .75, length.out = 4), limits = c(0, .75))

}


# Custom Themes ---------------------------------------------------------------------------------------------------
theme_shape_plot <- function() {
  ggplot2::theme_void(base_family = "sans", base_size = 3) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      strip.background = ggplot2::element_rect(fill = NA, color = NA),
      strip.text = ggplot2::element_text(size = 2.8, face = "bold", margin=margin(b=2)),
      strip.text.x = element_text(margin = margin(t=2)),
      strip.text.y = element_text(angle = 90, margin=margin(r=3)),
      axis.title.x = ggtext::element_markdown(family = "sans", face = "bold", hjust = 0.5, size = 3, margin=margin(t=2)),
      axis.title.y = ggtext::element_markdown(family = "sans", face = "bold", angle=90, size = 3, margin=margin(r=4)),
      axis.text.x = ggtext::element_markdown(family = "sans", color = "grey30", size = 2.5, angle = 90, margin = margin(t = 2)),
      axis.text.y = ggtext::element_markdown(family = "sans", color = "grey30", hjust=1, size = 2.5, margin = margin(r = 3)),
      panel.spacing.x = unit(0.5, "lines"),
      panel.spacing.y = unit(0.2, "lines"),
      axis.line.x = element_line(color = "grey70", linewidth = 0.1),
      axis.ticks.x = element_line(color = "grey70", linewidth = 0.1),
      axis.line.y = element_line(color = "grey70", linewidth=0.1),
      axis.ticks.y = element_line(color = "grey70", linewidth = 0.1),
      axis.ticks.length.x = unit(0.1, "lines"),
      axis.ticks.length.y = unit(0.1, "lines"),
      plot.margin = margin(2, 0, 2, 0),
      legend.position = "none"
    )
}

theme_convergence_plot <- function() {
  ggplot2::theme_void(base_family = "sans", base_size = 3) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      strip.background = ggplot2::element_rect(fill = NA, color = NA),
      strip.text = ggplot2::element_text(size = 2.8, face = "bold", margin=margin(b=2)),
      strip.text.x = element_text(margin = margin(t=2)),
      strip.text.y = element_text(angle = 90, margin=margin(r=3)),
      axis.title.x = ggtext::element_markdown(family = "sans", face = "bold", hjust = 0.5, size = 3, margin=margin(t=2)),
      axis.title.y = ggtext::element_markdown(family = "sans", face = "bold", angle=90, size = 3, margin=margin(r=4)),
      axis.text.x = ggtext::element_markdown(family = "sans", color = "grey30", size = 2.5, angle = 90, margin = margin(t = 2)),
      axis.text.y = ggtext::element_markdown(family = "sans", color = "grey30", hjust=1, size = 2.5, margin = margin(r = 3)),
      panel.spacing.x = unit(0.5, "lines"),
      panel.spacing.y = unit(0.2, "lines"),
      axis.line.x = element_line(color = "grey70", linewidth = 0.1),
      axis.ticks.x = element_line(color = "grey70", linewidth = 0.1),
      axis.line.y = element_line(color = "grey70", linewidth=0.1),
      axis.ticks.y = element_line(color = "grey70", linewidth = 0.1),
      axis.ticks.length.x = unit(0.1, "lines"),
      axis.ticks.length.y = unit(0.1, "lines"),
      plot.margin = margin(2, 0, 2, 0),
      legend.position = "none"
    )
}




