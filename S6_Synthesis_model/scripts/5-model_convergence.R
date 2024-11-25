#' ---------------------------------------------------------------------------------------------------------------
#' Script Name: Testing Dependencies Between Species and Traits Replacement ####
#' Last Updated: 24 November 2024
#' Author: Caio Graco-Roza
#' ---------------------------------------------------------------------------------------------------------------
#' 
#' Purpose:
#' This script performs statistical modeling and visualization to assess dependencies between species and traits
#' replacement across multiple datasets. It relies on null model results to evaluate convergence, divergence, and
#' magnitude relationships. The script generates visualizations (e.g., Figures 5a, 5b, 5c) and extended data figures.
#' 
#' Key Steps:
#' 1. **Data Preparation:** Load and process synthesis and SES data for dependency analysis.
#' 2. **Statistical Testing:**
#'    - Assess differences in magnitude and shapes between species and traits.
#'    - Use Bayesian models to analyze directional SES, magnitude SES, and shape divergence.
#' 3. **Visualization:** Generate comprehensive plots to display model results, including directional trends, 
#'    magnitude comparisons, and shape divergence.
#' 4. **Model Comparison:** Use Leave-One-Out Cross-Validation (LOO) to compare model fits and select the best model.
#' 
#' Dependencies:
#' - Null model results from previous scripts (S5).
#' - Custom helper functions from `helper_functions_S6.R` and `helper_functions_plot_extended.R`.
#' - Libraries for Bayesian modeling (brms), data visualization (ggplot2, ggdist), and data manipulation (tidyverse).
#' 
#' Outputs:
#' - Main Figures: Figures 5a, 5b, and 5c.
#' - Extended Figures: Detailed SES results for directional, magnitude, and shape analyses.
#' - Saved Bayesian model objects for reproducibility.
#'
#' ---------------------------------------------------------------------------------------------------------------

#Load required libraries
library(tidyverse)         # For data manipulation and visualization
library(ggdist)            # For visualizing distributions of data and model outputs
library(brms)              # For running Bayesian models
library(marginaleffects)   # For extracting marginal effects and predictions from Bayesian models
library(ggtext)            # For adding rich text elements to ggplot2 plots
library(ggblend)           # For blending graphics layers in ggplot2 plots
library(ggforestplot)      # For creating forest plots; install via devtools::install_github("NightingaleHealth/ggforestplot")
library(janitor)           # For cleaning and organizing data tables, and standardizing column names
library(ggalluvial)        # For visualizing alluvial diagrams (e.g., flows in Sankey diagrams)
library(ggh4x)             # For additional ggplot2 features like nested facets
library(patchwork)         # For combining multiple ggplots into complex layouts
library(cowplot)           # For creating multi-panel plots and additional ggplot2 functionalities
library("rstan")

# Set global Stan options ####
# Define options for running Bayesian models with brms and cmdstanr backend
base::options(
  mc.cores = 4,                # Use 4 cores for parallel processing
  brms.backend = "cmdstanr"    # Use cmdstanr as the backend for brms
  # ggblend.check_blend = FALSE # Optionally, disable ggblend checks (uncomment if needed)
)

# Load helper functions ####
# These scripts contain custom functions used for modeling and plotting
source("S6_Synthesis_model/functions/helper_functions_S6.R")  # Functions specific to synthesis models
source("S6_Synthesis_model/functions/helper_functions_plot_extended.R")  # Extended plotting functions

# Enable custom fonts ####
# Automatically handle custom fonts in ggplot2 plots
showtext::showtext_auto()

# Define global Bayesian model settings ####
CHAINS <- 4          # Number of MCMC chains
ITER <- 5000         # Total number of iterations per chain
WARMUP <- 2500       # Number of warmup iterations (not included in posterior sampling)
BAYES_SEED <- 1234   # Seed for reproducibility in Bayesian sampling

# Load and preprocess data ####
# Process data for testing dependency between species and traits
# This relies on null model results from S5
model_data <- process_model_data("S6_Synthesis_model/data/synthesis_data.xlsx", "convergence")


# Check for convergence in direction 
#   model_data |>
#   select(dataset,facet,direction) |>
#   pivot_wider(names_from = facet, values_from = direction) |>
#   janitor::tabyl(Species, Traits) |>
#   adorn_percentages("row") %>%
#   adorn_pct_formatting(digits = 2) %>%
#   adorn_ns() |>
#   kableExtra::kable(format = "rst")

# ===============  ===============  ==============
# species/traits   differentiation  homogenisation
# ===============  ===============  ==============
# differentiation  60.17% (71)      39.83% (47)   
# homogenisation   42.86% (18)      57.14% (24)   
# ===============  ===============  ==============

#make figure 5a 
annotations_df <- data.frame(
  label = c(
    "<span>&larr; Differentiation</span>",
    "<span>Homogenisation &rarr;</span>",
    "<span>&larr; Differentiation</span>",
    "<span>Homogenisation &rarr;</span>"
  ),
  x = c(-0.03, 0.03, -7, -7),
  y = c(-10, -10, -0.2, 0.2),
  angle = c(0, 0, 90, 90),
  hjust = c(1, 0, 1, 0))

Figure_5a <-
  model_data |>
  select(dataset, facet, direction_harrel_davis) |>
  pivot_wider(names_from = facet, values_from = direction_harrel_davis) |>
  filter(Traits > -7) |>
  ggplot(aes(y = Species, x = Traits)) +
  geom_point(size = 2, alpha = .7) +
  coord_cartesian(xlim = c(-7, NA), ylim = c(-10, 7)) +
  geom_richtext(
    data = annotations_df,
    aes(x = x, y = y, label = label, angle = angle, hjust = hjust), 
    size=convert_size(5),
    label.color = NA
  ) +
  geom_hline(yintercept = 0,
             linetype = 3,
             linewidth = .1) +
  geom_vline(xintercept = 0,
             linetype = 3,
             linewidth = .1) +
  ggplot2::theme_void(base_family = "sans", base_size = 6) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    plot.background = ggplot2::element_rect(fill = "white", color = NA),
    strip.background = ggplot2::element_rect(fill = NA, color = NA),
    strip.text = ggplot2::element_text(
      face = "bold",
      margin = margin(b = 2)
    ),
    axis.title.x = ggtext::element_markdown(
      family = "sans",
      face = "bold",
      hjust = 0.5,
      margin = margin(t = 2)
    ),
    axis.title.y = ggtext::element_markdown(
      family = "sans",
      face = "bold",
      angle = 90,
      margin = margin(r = 4)
    ),
    axis.text.x = element_markdown(
      family = "sans",
      color = "grey30",
      margin = margin(t = 2)
    ),
    axis.text.y = element_markdown(
      family = "sans",
      color = "grey30",
    ),
    panel.spacing.x = unit(0.2, "lines"),
    panel.spacing.y = unit(0.2, "lines"),
    axis.line.x = element_line(color = "grey70", linewidth = 0.1),
    axis.ticks.x = element_line(color = "grey70", linewidth = 0.1),
    axis.ticks.length.x = unit(0.1, "lines"),
    axis.line.y = element_line(color = "grey70", linewidth = 0.1),
    axis.ticks.y = element_line(color = "grey70", linewidth = 0.1),
    axis.ticks.length.y = unit(0.1, "lines"),
    plot.margin =  margin(3, 3, 3, 3),
    legend.position = "none",
  ) +
  labs(x = "Trait replacement [Harrel-Davis]", y = "Species replacement [Harrel-Davis]") +
  scale_x_continuous(breaks=c(-6.5, -4.5, -2.5, 0, 2.5))+
  guides(
    x = guide_axis_truncated(trunc_lower = -6.5, trunc_upper = 2.5),
    y = guide_axis_truncated(trunc_lower = -10, trunc_upper = 5)
  )

# Create annotations for directional labels ####
# The annotations_df data frame defines directional labels and their positions for the plot
annotations_df <- data.frame(
  label = c(
    "<span>&larr; Differentiation</span>",  # Left arrow label for differentiation (horizontal)
    "<span>Homogenisation &rarr;</span>",  # Right arrow label for homogenization (horizontal)
    "<span>&larr; Differentiation</span>",  # Upward arrow label for differentiation (vertical)
    "<span>Homogenisation &rarr;</span>"   # Downward arrow label for homogenization (vertical)
  ),
  x = c(-0.03, 0.03, -7, -7),              # X coordinates for placement
  y = c(-10, -10, -0.2, 0.2),              # Y coordinates for placement
  angle = c(0, 0, 90, 90),                 # Angles of the labels (horizontal or vertical)
  hjust = c(1, 0, 1, 0)                    # Horizontal justification for alignment
)

# Create Figure 5a ####
# The plot visualizes the relationship between species and trait replacement
Figure_5a <-
  model_data |> 
  # Select necessary columns for the plot
  select(dataset, facet, direction_harrel_davis) |>
  # Pivot data to a wider format, with separate columns for species and traits
  pivot_wider(names_from = facet, values_from = direction_harrel_davis) |>
  # Filter out extreme values for traits to avoid distortions
  filter(Traits > -7) |>
  # Begin ggplot object with trait replacement on x-axis and species replacement on y-axis
  ggplot(aes(y = Species, x = Traits)) +
  # Add points to represent datasets, with moderate size and transparency
  geom_point(size = 2, alpha = .7) +
  # Set the coordinate limits for the plot, focusing on the relevant range
  coord_cartesian(xlim = c(-7, NA), ylim = c(-10, 7)) +
  # Add annotations for directional labels using geom_richtext
  geom_richtext(
    data = annotations_df,
    aes(x = x, y = y, label = label, angle = angle, hjust = hjust), 
    size = convert_size(5),        # Convert font size to match plot aesthetics
    label.color = NA               # Remove the label box border
  ) +
  # Add a horizontal line at y = 0 to represent the baseline for species replacement
  geom_hline(yintercept = 0, linetype = 3, linewidth = .1) +
  # Add a vertical line at x = 0 to represent the baseline for trait replacement
  geom_vline(xintercept = 0, linetype = 3, linewidth = .1) +
  # Apply a minimal theme with custom settings for plot appearance
  ggplot2::theme_void(base_family = "sans", base_size = 6) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),              # Remove minor grid lines
    plot.background = ggplot2::element_rect(fill = "white", color = NA),  # White background
    strip.background = ggplot2::element_rect(fill = NA, color = NA),      # Remove strip background
    strip.text = ggplot2::element_text(face = "bold", margin = margin(b = 2)),  # Bold strip text
    axis.title.x = ggtext::element_markdown(                 # Custom x-axis title
      family = "sans",
      face = "bold",
      hjust = 0.5,
      margin = margin(t = 2)
    ),
    axis.title.y = ggtext::element_markdown(                 # Custom y-axis title
      family = "sans",
      face = "bold",
      angle = 90,
      margin = margin(r = 4)
    ),
    axis.text.x = element_markdown(                          # Custom x-axis text
      family = "sans",
      color = "grey30",
      margin = margin(t = 2)
    ),
    axis.text.y = element_markdown(                          # Custom y-axis text
      family = "sans",
      color = "grey30",
    ),
    panel.spacing.x = unit(0.2, "lines"),                    # Spacing between panels (if faceted)
    panel.spacing.y = unit(0.2, "lines"),                    # Spacing between panels (if faceted)
    axis.line.x = element_line(color = "grey70", linewidth = 0.1),  # X-axis line
    axis.ticks.x = element_line(color = "grey70", linewidth = 0.1), # X-axis ticks
    axis.ticks.length.x = unit(0.1, "lines"),                # X-axis tick length
    axis.line.y = element_line(color = "grey70", linewidth = 0.1),  # Y-axis line
    axis.ticks.y = element_line(color = "grey70", linewidth = 0.1), # Y-axis ticks
    axis.ticks.length.y = unit(0.1, "lines"),                # Y-axis tick length
    plot.margin =  margin(3, 3, 3, 3),                       # Margins around the plot
    legend.position = "none"                                 # Remove legend
  ) +
  # Add axis labels with titles for replacement metrics
  labs(x = "Trait replacement [Harrel-Davis]", y = "Species replacement [Harrel-Davis]") +
  # Define x-axis breaks and truncation limits for improved readability
  scale_x_continuous(breaks = c(-6.5, -4.5, -2.5, 0, 2.5)) +
  guides(
    x = guide_axis_truncated(trunc_lower = -6.5, trunc_upper = 2.5),  # Truncated x-axis
    y = guide_axis_truncated(trunc_lower = -10, trunc_upper = 5)     # Truncated y-axis
  )

# Data Preparation for Magnitude Comparison ---------------------------------

# Create a dataset for comparing magnitude between biodiversity facets (Species and Traits).
# The data is reshaped to a wider format where Species and Traits are columns.
magnitude_comparison <- model_data |>
  select(dataset, facet, magnitude) |>
  pivot_wider(names_from = facet, values_from = magnitude)

# Statistical Test for Magnitude Difference ---------------------------------

# Conduct a Wilcoxon signed-rank test to check if the magnitudes differ significantly 
# between Species and Traits across datasets. This test accounts for paired observations.
# 
# Uncomment to run the test:
# rstatix::wilcox_test(
#   data = model_data,
#   magnitude ~ facet,
#   paired = TRUE,
#   alternative = "two.sided"
# ) |>
#   kableExtra::kable(format = "rst")

# Example output of the Wilcoxon test:
# =========  =======  ======  ===  ===  =========  =====
# .y.        group1   group2   n1   n2  statistic      p
# =========  =======  ======  ===  ===  =========  =====
# magnitude  Species  Traits  160  160       5875  0.587
# =========  =======  ======  ===  ===  =========  =====

# Additional Analysis of Magnitude Comparison ---------------------------------

# Analyze the proportion of datasets where the magnitude for Traits is greater than Species.
# Uncomment to compute the proportions:
# model_data |> 
#   select(dataset, facet, magnitude) |> 
#   pivot_wider(names_from = facet, values_from = magnitude) |> 
#   mutate(higher = Traits > Species) |> 
#   janitor::tabyl(higher) |> 
#   adorn_pct_formatting(digits = 2) |>
#   kableExtra::kable(format = "rst")

# Example output:
# ======  ===  =======
# higher    n  percent
# ======  ===  =======
# FALSE    92  57.50% 
# TRUE     68  42.50% 
# ======  ===  =======

# Visualization: Create Figure 5b --------------------------------------------

# Figure 5b visualizes the magnitude comparison between Species and Traits replacements.
# A segment connects the corresponding magnitudes for Species and Traits to highlight their relationship.

Figure_5b <- 
  model_data |> 
  # Convert the facet variable (Species, Traits) to numeric for positioning.
  mutate(facet = as.numeric(factor(facet, levels = c("Species", "Traits")))) |> 
  ggplot(aes(x = facet, y = magnitude)) +
  
  # Add connecting segments for pairs of Species and Traits magnitudes where Species magnitude < 0.7.
  geom_segment(
    data = magnitude_comparison |> filter(Species < 0.7), 
    aes(x = 1.97, y = Traits, xend = 1.03, yend = Species),
    linewidth = 0.1, colour = "gray60", alpha = 0.7
  ) +
  
  # Add density slabs for the Trait magnitudes, positioned to the left and right of the y-axis.
  ggdist::stat_slab(
    data = magnitude_comparison |> filter(Species < 0.7), 
    aes(y = Traits, x = 0.9), 
    width = 0.3, side = "left"
  ) +
  ggdist::stat_slab(
    data = magnitude_comparison |> filter(Species < 0.7), 
    aes(y = Traits, x = 2.1), 
    width = 0.3, side = "right"
  ) +
  
  # Add scatter points for individual magnitude values.
  geom_point(size = 0.5, colour = "gray30") +
  
  # Styling and layout adjustments for the plot.
  ggplot2::theme_void(base_family = "sans", base_size = 6) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    plot.background = ggplot2::element_rect(fill = "white", color = NA),
    strip.background = ggplot2::element_rect(fill = NA, color = NA),
    strip.text = ggplot2::element_text(face = "bold", margin = margin(b = 2)),
    strip.text.x = element_text(margin = margin(t = 2)),
    strip.text.y = element_text(angle = 90, margin = margin(r = 3)),
    axis.title.x = ggtext::element_markdown(
      family = "sans", face = "bold", hjust = 0.5, margin = margin(b = 2)
    ),
    axis.title.y = ggtext::element_markdown(
      family = "sans", face = "bold", angle = 90, margin = margin(r = 4)
    ),
    axis.text.x = ggtext::element_markdown(
      family = "sans", color = "grey30", margin = margin(b = 2)
    ),
    axis.text.y = ggtext::element_markdown(
      family = "sans", color = "grey30", hjust = 1, margin = margin(r = 3)
    ),
    panel.spacing.x = unit(0.5, "lines"),
    panel.spacing.y = unit(0.2, "lines"),
    axis.line.x = element_line(color = "grey70", linewidth = 0.1),
    axis.ticks.x = element_line(color = "grey70", linewidth = 0.1),
    axis.line.y = element_line(color = "grey70", linewidth = 0.1),
    axis.ticks.y = element_line(color = "grey70", linewidth = 0.1),
    axis.ticks.length.x = unit(0.1, "lines"),
    axis.ticks.length.y = unit(0.1, "lines"),
    plot.margin = margin(3, 3, 3, 7),
    legend.position = "none"
  ) +
  
  # Adjust coordinate limits and axis scales.
  coord_cartesian(ylim = c(0, 0.7), xlim = c(0.5, 2.5)) +
  labs(
    y = "Magnitude of effect", 
    x = NA  # No x-axis label
  ) +
  scale_x_continuous(
    breaks = seq(1, 2), 
    labels = c("Species<br>replacement", "Trait<br>replacement"), 
    position = "top"
  ) +
  scale_y_continuous(breaks = seq(0, 0.7, by = 0.1)) +
  
  # Add truncated axes to improve visual focus.
  guides(
    y = guide_axis_truncated(trunc_lower = 0, trunc_upper = 0.7),
    x = guide_axis_truncated(
      trunc_lower = function(x) { x - 0.2 },
      trunc_upper = function(x) { x + 0.2 }
    )
  )

# Statistical Test for Shape Convergence -------------------------------------
# Analyze if shapes of replacement converge between Species and Traits.
# Uncomment to check the convergence statistics:
# model_data |> 
#   select(dataset, facet, shape) |> 
#   group_by(dataset, facet) |>
#   pivot_wider(names_from = facet, values_from = shape) |> 
#   mutate(convergence = Species == Traits) |> 
#   janitor::tabyl(Species, Traits) |> 
#   adorn_percentages("row") %>%
#   adorn_pct_formatting(digits = 2) %>%
#   adorn_ns() |>
#   kableExtra::kable(format = "rst")

# Example output:
# ===========  ===========  ===========  ===========  ==========
# Species      Absent       Exponential  Saturating   Revlog    
# ===========  ===========  ===========  ===========  ==========
# Absent       46.88% (15)  18.75%  (6)  31.25% (10)  3.12% (1) 
# Exponential  17.50%  (7)  42.50% (17)  35.00% (14)  5.00% (2) 
# Saturating   17.39% (12)  34.78% (24)  37.68% (26)  10.14% (7)
# Revlog       5.26%  (1)   42.11%  (8)  36.84%  (7)  15.79% (3)
# ===========  ===========  ===========  ===========  ==========

# Visualization: Create Figure 5c -------------------------------------------

# Create an alluvial diagram to represent convergence and divergence in replacement shapes
# between Species and Traits.
Figure_5c <- 
  model_data |> 
  # Process data to map replacement shapes between Species and Traits.
  select(dataset, facet, shape) |> 
  mutate(shape = gsub("Revlog", "Reverse logistic", shape)) |>  # Standardize naming
  pivot_wider(names_from = facet, values_from = shape) |>      # Reshape for comparison
  group_by(Species) |> 
  count(Traits) |>                                             # Count occurrences of shape pairs
  mutate(Convergence = factor(Species == Traits)) |>           # Create a convergence flag
  mutate(Species = paste(Species, " ")) |>                     # Add padding to Species labels
  ggplot(
    aes(
      y = n,                                  # Number of datasets as y-axis
      axis1 = Species,                       # First axis: Species shapes
      axis2 = Traits                         # Second axis: Traits shapes
    )) +
  
  # Add alluvial flows (connections between Species and Traits shapes)
  geom_alluvium(aes(fill = Convergence), curve_type = "cubic", show.legend = FALSE) +
  
  # Add strata for Species and Traits shapes with counts as labels
  geom_stratum(size = 0.1) +
  geom_text(
    stat = "stratum",
    aes(label = gsub("   ", " ", paste0(stratum, " (n = ", after_stat(count), ")"))),
    color = "black", size = convert_size(5)
  ) +
  
  # Style adjustments
  ggplot2::theme_void(base_family = "sans", base_size = 6) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    plot.background = ggplot2::element_rect(fill = "white", color = NA),
    axis.title.y = ggtext::element_markdown(
      family = "sans", face = "bold", angle = 90, margin = margin(r = 4)
    ),
    axis.text.y = ggtext::element_markdown(
      family = "sans", color = "grey30", hjust = 1, margin = margin(r = 3)
    ),
    axis.text.x = ggtext::element_markdown(
      family = "sans", color = "grey30", margin = margin(b = 2)
    ),
    axis.line.x = element_line(color = "grey70", linewidth = 0.1),
    axis.ticks.x = element_line(color = "grey70", linewidth = 0.1),
    axis.line.y = element_line(color = "grey70", linewidth = 0.1),
    axis.ticks.y = element_line(color = "grey70", linewidth = 0.1),
    axis.ticks.length.x = unit(0.1, "lines"),
    axis.ticks.length.y = unit(0.1, "lines"),
    plot.margin = margin(7, 3, 3, 3),
    legend.position = "none"
  ) +
  
  # Configure x-axis as discrete with replacement facets
  scale_x_discrete(
    limits = c("Species", "Traits"),
    labels = c("Species replacement", "Trait replacement"),
    expand = c(0.15, 0.05),
    position = "top"
  ) +
  
  # Define fill colors for convergence and divergence
  scale_fill_manual(values = c("red", "black")) +
  
  # Add y-axis label and configure breaks
  labs(y = "Number of datasets") +
  scale_y_continuous(breaks = seq(0, 160, by = 20)) +
  
  # Truncate axes to improve focus
  guides(
    y = guide_axis_truncated(trunc_lower = 0, trunc_upper = 160),
    x = guide_axis_truncated(
      trunc_lower = function(x) { x - 0.1 },
      trunc_upper = function(x) { x + 0.1 }
    )
  )

# Combine Figures 5a, 5b, and 5c into a Composite Figure ---------------------

# Combine Figure 5a and Figure 5b in a single row, followed by Figure 5c below.
Figure_5 <- plot_grid(
  plot_grid(Figure_5a, Figure_5b, ncol = 2, rel_widths = c(0.65, 0.35),  
            labels = c("a", "b"), label_fontface = "bold", label_size = 8),
  Figure_5c,
  ncol = 1,
  rel_heights = c(1, 1),
  labels = c("", "c"),
  label_fontface = "bold",
  label_size = 8
)

# Save the Composite Figure -------------------------------------------------

# Save the combined Figure 5 as a PDF file with appropriate dimensions.
ggsave(
  plot = Figure_5,
  filename = here::here("S8_Model_outputs_figures_and_tables", "main_figures", "Figure_5.pdf"),
  width = 130, height = 120, units = "mm"
)

#' ---------------------------------------------------------------------------------------------------------------
# Modelling SES for Direction ####
#' ---------------------------------------------------------------------------------------------------------------

# Statistical Analysis of Direction SES -------------------------------------------------------

# Uncomment this block to check how datasets exhibit convergent or divergent directions
# based on SES and significance of p-values. This analysis summarizes the counts and
# proportions of datasets for each direction.

# model_data |> 
#   select(dataset, facet, direction) |> 
#   pivot_wider(names_from = facet, values_from = direction) |> 
#   rename(direction = Traits) |> 
#   left_join(
#     ses_data |>  
#       select(dataset, direction, direction_SES, ends_with("pvalue")) |>  
#       mutate(direction = str_to_title(gsub("z", "s", direction)))
#   ) |> 
#   mutate(
#     convergent = ifelse(Species == direction, "Convergent", "Divergent"),
#     signal_dir = case_when(
#       sign(direction_SES) == 1 & direction_pvalue < 0.05 ~ "Positive",
#       sign(direction_SES) == -1 & direction_pvalue < 0.05 ~ "Negative",
#       TRUE ~ "Random"
#     )
#   ) |> 
#   janitor::tabyl(signal_dir, Species, convergent) |> 
#   janitor::adorn_percentages("col") %>%
#   janitor::adorn_pct_formatting(digits = 2) %>%
#   janitor::adorn_ns() |> 
#   kableExtra::kable(format = "rst")

# Example Output:
# ==========  ===============  ==============
# Convergent  Differentiation  Homogenisation
# ==========  ===============  ==============
# Negative    32.39% (23)      16.67%  (4)   
# Positive    23.94% (17)      45.83% (11)   
# Random      43.66% (31)      37.50%  (9)   
# ==========  ===============  ==============

# ==========  ===============  ==============
# Divergent   Differentiation  Homogenisation
# ==========  ===============  ==============
# Negative    19.15%  (9)      22.22% (4)    
# Positive    36.17% (17)      27.78% (5)    
# Random      44.68% (21)      50.00% (9)    
# ==========  ===============  ==============

# Data Preparation for Direction SES ----------------------------------------------------------

# Load SES data for direction analysis and process it to align with model input.
ses_data <- 
  "S6_Synthesis_model/data/BBGDM_SES.xlsx" |>  # Path to SES data file
  readxl::read_xlsx() |>                       # Read the data
  filter(!dataset %in% c("N67TTP", "N78TTP", "S47TTP")) |>  # Exclude specific datasets
  select(-1) |>                                # Remove the first column
  mutate(direction = str_to_title(direction)) # Standardize direction names

# Merge SES data with model data for traits to analyze directional response
data_direction <-  
  model_data |> 
  filter(facet == "Traits") |>  # Focus on Traits facet only
  left_join(
    ses_data |>  
      select(dataset, direction, direction_SES, direction_pvalue) |>  
      mutate(direction = str_to_title(gsub("z", "s", direction)))  # Adjust direction column format
  ) |> 
  drop_na() |>  # Remove missing values
  # Categorize directional SES into "Smaller", "Higher", or "Random" categories
  mutate(
    direction_ses = factor(
      case_when(
        direction_pvalue <= 0.05 & direction_SES < 0 ~ "Smaller",
        direction_pvalue <= 0.05 & direction_SES > 0 ~ "Higher",
        TRUE ~ "Random"
      ),
      levels = c("Random", "Smaller", "Higher")
    ),
    .keep = "unused"
  ) |> 
  # Define the response variable based on direction SES and Harrel-Davis estimators
  mutate(
    response = case_when(
      direction_harrel_davis < 0 & direction_ses == "Smaller" ~ "homogenisation",
      direction_harrel_davis < 0 & direction_ses == "Higher" ~ "differentiation",
      direction_harrel_davis > 0 & direction_ses == "Smaller" ~ "homogenisation",
      direction_harrel_davis > 0 & direction_ses == "Higher" ~ "differentiation",
      TRUE ~ "Random"
    ),
    response = factor(response)
  )

# Model Formula and Priors for Direction Analysis ----------------------------------------------

# Define the base formula for the Bayesian model
formula_base <- bf(
  response ~ species_number + spatial_extent + spatial_distance_min + 
    absolute_latitude_mean + human_pressure_baseline + human_pressure_range + 
    main_land_use_type + (1 | taxa_coarse) + buffer + ecosystem_type, 
  decomp = "QR"
)

# Specify priors for the model
priors_direction <- c(
  # Priors for standard deviations
  prior(student_t(3, 0, 1), class = "sd", lb = 0, dpar = "mudifferentiation"),
  prior(student_t(3, 0, 1), class = "sd", lb = 0, dpar = "muhomogenisation"),
  prior(student_t(3, 0, 1), class = "sd", group = "taxa_coarse", lb = 0, dpar = "mudifferentiation"),
  prior(student_t(3, 0, 1), class = "sd", group = "taxa_coarse", lb = 0, dpar = "muhomogenisation"),
  
  # Priors for slopes
  prior(normal(0, 3), class = "b", dpar = "mudifferentiation"),
  prior(normal(0, 3), class = "b", dpar = "muhomogenisation")
)

# Fit the Base Model ---------------------------------------------------------------------------

# Run the Bayesian categorical model for direction
md_direction_base <- brms::brm(
  formula_base,
  data = data_direction,
  prior = priors_direction,
  chains = CHAINS,
  family = brms::categorical(refcat = "Random"),  # Reference category is "Random"
  iter = ITER,
  warmup = WARMUP,
  control = base::list(max_treedepth = 12, adapt_delta = .99),
  save_pars = save_pars(all = TRUE),
  seed = BAYES_SEED,
  file = "S7_Model_outputs_figures_and_tables/model/convergence/direction_base",
  file_refit = "on_change"
)

# Additional Models with Interaction Terms ----------------------------------------------------

# Add ecosystem type and land use interaction
md_direction_interaction_ecosystem_land_use <- update(
  md_direction_base,
  newdata = data_direction,
  formula. = . ~ . + ecosystem_type * main_land_use_type,
  file = "S7_Model_outputs_figures_and_tables/model/convergence/direction_interaction_ecosystem_land_use"
)

# Base model with nested taxa random effects
md_direction_base_nested <- update(
  md_direction_base,
  newdata = data_direction,
  formula. = . ~ . - (1 | taxa_coarse) + (1 | taxa_coarse / taxa_fine),
  file = "S7_Model_outputs_figures_and_tables/model/convergence/direction_base_nested"
)

# Interaction model with nested taxa random effects
md_direction_interaction_ecosystem_land_use_nested <- update(
  md_direction_base,
  newdata = data_direction,
  formula. = . ~ . + ecosystem_type * main_land_use_type - (1 | taxa_coarse) + (1 | taxa_coarse / taxa_fine),
  file = "S7_Model_outputs_figures_and_tables/model/convergence/direction_interaction_ecosystem_land_use_nested"
)

# Model Comparison Using Leave-One-Out Cross-Validation ----------------------------------------

# Compare models using LOO cross-validation with moment matching
loo_values <- loo(
  md_direction_base, 
  md_direction_interaction_ecosystem_land_use, 
  md_direction_base_nested, 
  md_direction_interaction_ecosystem_land_use_nested,
  moment_match = TRUE
)

# Visualization of SES Results -----------------------------------------------------------------

# Select the best model based on LOO and generate data for visualization
best_model <- md_direction_base  # Select the base model as the best for simplicity

# Generate posterior draws for visualization
direction_figure_data <- generate_ses_figure_data(
  fit_model = best_model,
  model_type = "direction",
  draw_min = -5,
  draw_max = 10
)

# Plot posterior distributions and save the figure
Extended_Figure_6 <- plot_posterior_draws(
  data = direction_figure_data,
  facet_var = "response",
  xlim = c(-5, 12),
  trunc_lower = -5,
  trunc_upper = 5,
  tick_factor = 2.5
)

# Save the extended figure to the specified directory
ggsave(
  filename = here::here("S7_Model_outputs_figures_and_tables", "extended_data", "Extended_Figure_6.pdf"),
  plot = Extended_Figure_6,
  device = cairo_pdf,
  width = 89,
  height = 40,
  units = "mm"
)

#' ---------------------------------------------------------------------------------------------------------------
# Modelling Magnitude Divergence ####
#' ---------------------------------------------------------------------------------------------------------------

# Data Preparation for Magnitude Divergence ----------------------------------------------------

# Prepare data by filtering for traits facet and merging SES data
data_magnitude <- model_data |> 
  filter(facet == "Traits") |>  # Focus on "Traits" facet only
  left_join(
    ses_data |>  
      select(dataset, direction, magnitude_SES, magnitude_pvalue) |>  # Select relevant SES columns
      mutate(direction = gsub("z", "s", direction))  # Standardize direction names
  ) |> 
  drop_na() |>  # Remove rows with missing values
  # Categorize magnitude SES into "Weaker", "Stronger", or "Random" categories
  mutate(
    response = factor(
      case_when(
        magnitude_pvalue <= 0.05 & magnitude_SES < 0 ~ "Weaker",
        magnitude_pvalue <= 0.05 & magnitude_SES > 0 ~ "Stronger",
        TRUE ~ "Random"
      ),
      levels = c("Random", "Weaker", "Stronger")
    )
  )

# Define Priors for Magnitude Divergence -------------------------------------------------------

# Specify priors for Bayesian categorical models
priors_magnitude <- c(
  # Priors for slopes
  prior(normal(0, 3), class = b, dpar = "muStronger"),
  prior(normal(0, 3), class = b, dpar = "muWeaker"),
  
  # Priors for standard deviations
  prior(student_t(3, 0, 0.5), class = sd, lb = 0, dpar = "muStronger"),
  prior(student_t(3, 0, 0.5), class = sd, lb = 0, dpar = "muWeaker")
)

# Fit Base Model for Magnitude Divergence ------------------------------------------------------

# Base model with main predictors and direction
md_magnitude_base <- brms::brm(
  update(formula_base, . ~ . + direction),  # Update formula to include direction
  data = data_magnitude,                   # Use prepared magnitude data
  prior = priors_magnitude,                # Apply specified priors
  chains = CHAINS,                         # Number of chains
  family = brms::categorical(refcat = "Random"),  # Set reference category to "Random"
  iter = ITER,                             # Total iterations
  warmup = WARMUP,                         # Warmup iterations
  control = base::list(max_treedepth = 12, adapt_delta = .99),  # Control parameters
  save_pars = save_pars(all = TRUE),       # Save parameters for reproducibility
  seed = BAYES_SEED,                       # Set seed for reproducibility
  file = "S7_Model_outputs_figures_and_tables/model/divergence_magnitude",  # File path for saved model
  file_refit = "on_change"                 # Refit only if parameters change
)

# Models with Interaction Terms ---------------------------------------------------------------

# Add interactions between ecosystem type and land use
md_magnitude_interaction_ecosystem_land_use <- update(
  md_magnitude_base,
  newdata = data_magnitude,
  formula. = . ~ . + ecosystem_type * main_land_use_type,
  file = "S7_Model_outputs_figures_and_tables/model/convergence/magnitude_interaction_ecosystem_land_use",
  file_refit = "on_change"
)

# Add interactions between ecosystem type and direction
md_magnitude_interaction_ecosystem_direction <- update(
  md_magnitude_base,
  newdata = data_magnitude,
  formula. = . ~ . + ecosystem_type * direction,
  file = "S7_Model_outputs_figures_and_tables/model/convergence/magnitude_interaction_ecosystem_direction",
  file_refit = "on_change"
)

# Add interactions between land use and direction
md_magnitude_interaction_land_use_direction <- update(
  md_magnitude_base,
  newdata = data_magnitude,
  formula. = . ~ . + main_land_use_type * direction,
  file = "S7_Model_outputs_figures_and_tables/model/convergence/magnitude_interaction_land_use_direction",
  file_refit = "on_change"
)

# Add three-way interaction
md_magnitude_three_way_interaction <- update(
  md_magnitude_base,
  newdata = data_magnitude,
  formula. = . ~ . + ecosystem_type * main_land_use_type * direction,
  file = "S7_Model_outputs_figures_and_tables/model/convergence/direction_three_way_interaction",
  file_refit = "on_change"
)

# Models with Nested Random Effects -----------------------------------------------------------

# Base model with nested random effects
md_magnitude_base_nested <- update(
  md_magnitude_base,
  newdata = data_magnitude,
  formula. = . ~ . - (1 | taxa_coarse) + (1 | taxa_coarse / taxa_fine),  # Remove old random effects and add nested
  file = "S7_Model_outputs_figures_and_tables/model/convergence/magnitude_base_nested",
  file_refit = "on_change"
)

# Nested model with ecosystem type and land use interaction
md_magnitude_interaction_ecosystem_land_use_nested <- update(
  md_magnitude_base_nested,
  newdata = data_magnitude,
  formula. = . ~ . + ecosystem_type * main_land_use_type,
  file = "S7_Model_outputs_figures_and_tables/model/convergence/magnitude_interaction_ecosystem_land_use_nested",
  file_refit = "on_change"
)

# Nested model with ecosystem type and direction interaction
md_magnitude_interaction_ecosystem_direction_nested <- update(
  md_magnitude_base_nested,
  newdata = data_magnitude,
  formula. = . ~ . + ecosystem_type * direction,
  control = base::list(max_treedepth = 15, adapt_delta = .99),  # Adjust control settings
  file = "S7_Model_outputs_figures_and_tables/model/convergence/magnitude_interaction_ecosystem_direction_nested",
  file_refit = "on_change"
)

# Nested model with land use and direction interaction
md_magnitude_interaction_land_use_direction_nested <- update(
  md_magnitude_base_nested,
  newdata = data_magnitude,
  formula. = . ~ . + main_land_use_type * direction,
  control = base::list(max_treedepth = 15, adapt_delta = .99),  # Adjust control settings
  file = "S7_Model_outputs_figures_and_tables/model/convergence/magnitude_interaction_land_use_direction_nested",
  file_refit = "on_change"
)

# Nested model with three-way interaction
md_magnitude_three_way_interaction_nested <- update(
  md_magnitude_base_nested,
  newdata = data_magnitude,
  formula. = . ~ . + ecosystem_type * main_land_use_type * direction,
  control = base::list(max_treedepth = 15, adapt_delta = .99),  # Adjust control settings
  file = "S7_Model_outputs_figures_and_tables/model/convergence/magnitude_three_way_interaction_nested",
  file_refit = "on_change"
)

# Model Comparison Using LOO Cross-Validation -------------------------------------------------

# Compare all models using LOO with moment matching
loo_magnitude <- loo(
  md_magnitude_base,
  md_magnitude_interaction_ecosystem_land_use,
  md_magnitude_interaction_ecosystem_direction,
  md_magnitude_interaction_land_use_direction,
  md_magnitude_three_way_interaction,
  md_magnitude_base_nested,
  md_magnitude_interaction_ecosystem_land_use_nested,
  md_magnitude_interaction_ecosystem_direction_nested,
  md_magnitude_interaction_land_use_direction_nested,
  md_magnitude_three_way_interaction_nested,
  moment_match = TRUE
)

# Visualization of Magnitude Divergence Results -----------------------------------------------

# Select the best model (base model in this case)
best_model_magnitude <- md_magnitude_base

# Generate posterior draws for visualization
magnitude_figure_data <- generate_ses_figure_data(
  fit_model = best_model_magnitude,
  model_type = "magnitude",
  draw_min = -5,
  draw_max = 10
)

# Plot posterior draws for magnitude divergence
Extended_Figure_7 <- plot_posterior_draws(
  data = magnitude_figure_data,
  facet_var = "response",
  xlim = c(-5, 12),
  trunc_lower = -5,
  trunc_upper = 5,
  tick_factor = 2.5
)

# Save the figure to the specified directory
ggsave(
  filename = here::here("S7_Model_outputs_figures_and_tables", "extended_data", "Extended_Figure_7.pdf"),
  plot = Extended_Figure_7,
  device = cairo_pdf,
  width = 89,
  height = 40,
  units = "mm"
)

#' ---------------------------------------------------------------------------------------------------------------
# Modelling Shape Divergence ####
#' ---------------------------------------------------------------------------------------------------------------

# Data Preparation for Shape Divergence -------------------------------------------------------

# Extract data specific to shape divergence analysis
data_shape <- extract_model_ses_data(
  ses_data = ses_data,         # SES data containing null model results
  model_data = model_data,     # Original model data
  model_type = "shape"         # Specify analysis type as "shape"
)

# Define the Base Formula ---------------------------------------------------------------------

# Base formula specifying predictors and response categories
base_formula_string <- "mvbind(Absent, Revlog, Saturating, Exponential) ~ species_number + spatial_extent + spatial_distance_min + absolute_latitude_mean + human_pressure_baseline + human_pressure_range + ecosystem_type + main_land_use_type + (1 | taxa_coarse) + buffer + direction"

# Convert the base formula string to a brms-friendly formula
formula_shape_base <- bf(as.formula(base_formula_string), decomp = "QR")

# Define Interaction and Nested Terms --------------------------------------------------------

# Specify interaction terms for testing combined effects
interaction_ecosystem_land_use <- "ecosystem_type * main_land_use_type"
interaction_ecosystem_direction <- "ecosystem_type * direction"
interaction_land_use_direction <- "main_land_use_type * direction"
interaction_all_three <- "ecosystem_type * main_land_use_type * direction"

# Specify nested term for hierarchical models
nested_term <- "(1 | taxa_coarse/taxa_fine)"

# Construct Interaction and Nested Formulas --------------------------------------------------

# Formulas for models with interaction terms
formula_interaction_ecosystem_land_use <- bf(as.formula(glue::glue("{base_formula_string} + {interaction_ecosystem_land_use}")), decomp = "QR")
formula_interaction_ecosystem_direction <- bf(as.formula(glue::glue("{base_formula_string} + {interaction_ecosystem_direction}")), decomp = "QR")
formula_interaction_land_use_direction <- bf(as.formula(glue::glue("{base_formula_string} + {interaction_land_use_direction}")), decomp = "QR")
formula_three_way_interaction <- bf(as.formula(glue::glue("{base_formula_string} + {interaction_all_three}")), decomp = "QR")

# Formulas for nested models
formula_base_nested <- bf(as.formula(glue::glue("{base_formula_string} - (1 | taxa_coarse) + {nested_term}")), decomp = "QR")
formula_interaction_ecosystem_land_use_nested <- bf(as.formula(glue::glue("{base_formula_string} + {interaction_ecosystem_land_use} - (1 | taxa_coarse) + {nested_term}")), decomp = "QR")
formula_interaction_ecosystem_direction_nested <- bf(as.formula(glue::glue("{base_formula_string} + {interaction_ecosystem_direction} - (1 | taxa_coarse) + {nested_term}")), decomp = "QR")
formula_interaction_land_use_direction_nested <- bf(as.formula(glue::glue("{base_formula_string} + {interaction_land_use_direction} - (1 | taxa_coarse) + {nested_term}")), decomp = "QR")
formula_three_way_interaction_nested <- bf(as.formula(glue::glue("{base_formula_string} + {interaction_all_three} - (1 | taxa_coarse) + {nested_term}")), decomp = "QR")

# Fit Base Model ------------------------------------------------------------------------------

# Fit the base model with main predictors and random effects
md_shape_base <- brms::brm(
  formula = formula_shape_base,
  data = data_shape,
  chains = CHAINS,                           # Number of MCMC chains
  family = brms::categorical(refcat = "Random"),  # Reference category is "Random"
  iter = ITER,                               # Total number of iterations
  warmup = WARMUP,                           # Number of warmup iterations
  control = base::list(max_treedepth = 12, adapt_delta = .99),  # Control parameters
  save_pars = save_pars(all = TRUE),         # Save parameters for reproducibility
  seed = BAYES_SEED,                         # Set seed for reproducibility
  file = "S7_Model_outputs_figures_and_tables/model/convergence/shape_base",  # File path for saving the model
  file_refit = "on_change"                   # Refit only if parameters change
)

# Fit Models with Interactions and Nested Structures ------------------------------------------

# Using the helper function `run_and_save_model` to simplify model fitting and saving
md_shape_interaction_ecosystem_land_use <- run_and_save_model(data = data_shape, formula = formula_interaction_ecosystem_land_use)  # Fails
md_shape_interaction_land_use_direction <- run_and_save_model(data = data_shape, formula = formula_interaction_land_use_direction)  # Successful
md_shape_three_way_interaction <- run_and_save_model(data = data_shape, formula = formula_three_way_interaction)  # Fails

# Nested models
md_shape_base_nested <- run_and_save_model(data = data_shape, formula = formula_base_nested)  # Successful
md_shape_interaction_ecosystem_land_use_nested <- run_and_save_model(data = data_shape, formula = formula_interaction_ecosystem_land_use_nested)  # Fails
md_shape_interaction_land_use_direction_nested <- run_and_save_model(data = data_shape, formula = formula_interaction_land_use_direction_nested)  # Successful
md_shape_three_way_interaction_nested <- run_and_save_model(data = data_shape, formula = formula_three_way_interaction_nested)  # Fails

# Model Comparison Using LOO Cross-Validation -------------------------------------------------

# Compare models using LOO for model selection
loo(
  md_shape_base,
  md_shape_base_nested,
  md_shape_interaction_land_use_direction_nested,
  md_shape_interaction_land_use_direction
)

# Select the best model (base model in this case)
best_shape_model <- md_shape_base

# Visualization of Shape Divergence Results ---------------------------------------------------

# Generate posterior draws for visualization
shape_figure_data <- generate_ses_figure_data(
  fit_model = best_shape_model,
  model_type = "shape",
  draw_min = -8,
  draw_max = 10
) |>  
  # Format shape levels for better readability
  mutate(
    shape = factor(
      gsub("Revlog", "Reverse logistic", shape), 
      levels = c("Absent", "Exponential", "Reverse logistic", "Saturating")
    )
  )

# Plot posterior draws for shape divergence
Extended_Figure_8 <- plot_posterior_draws(
  data = shape_figure_data,
  facet_var = c("response", "shape"),
  xlim = c(-8, 12),
  trunc_lower = -7.5,
  trunc_upper = 5,
  tick_factor = 2.5,
  add_annotations = NULL  # No annotations for this plot
)

# Save the plot to the specified directory
ggsave(
  filename = here::here("S7_Model_outputs_figures_and_tables", "extended_data", "Extended_Figure_8.pdf"),
  plot = Extended_Figure_8 + theme(plot.margin = margin(3, 3, 3, 3)),
  device = cairo_pdf(),
  width = 89,
  height = 120,
  units = "mm"
)

