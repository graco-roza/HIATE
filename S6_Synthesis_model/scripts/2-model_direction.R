#' ###########################################################################
#' SCRIPT NAME: Bayesian Modeling of Response Direction (2-model_direction.R)
#'
#' DESCRIPTION:
#'   This script performs Bayesian modeling to analyze the direction of species 
#'   and Functional turnover across different ecosystems and land-use types. 
#'   It also generates exploratory figures and saves model outputs for further
#'   analysis and visualization.
#'
#' USAGE:
#'   - The script reads preprocessed data and runs Bayesian models using the 
#'     `brms` package.
#'   - Several models with varying complexity (base, interaction, nested taxa) 
#'     are evaluated using Leave-One-Out Cross-Validation (LOO).
#'   - Generates the figures for inclusion in the publication.
#'
#' INPUTS:
#'   - Preprocessed synthesis data from `S6_Synthesis_model/data/synthesis_data.xlsx`
#'   - Helper functions for model processing and plotting
#'
#' OUTPUTS:
#'   - Bayesian models saved in `S7_Model_outputs_figures_and_tables/model/direction`
#'   - Figures saved in `S7_Model_outputs_figures_and_tables/main_figures` and
#'     `S7_Model_outputs_figures_and_tables/extended_data`
#'
#' AUTHOR: Caio Graco-Roza
#' LAST UPDATED: 2024-11-24
#'
#' NOTES:
#'   - Adjust the global settings for the Bayesian models as needed (chains, iterations, etc.).
#'   - Ensure all required libraries are installed, and the input files are in the correct path.
#' ###########################################################################

# Load packages ####
library(tidyverse) # For data manipulation and visualization
library(ggdist)  # For visualizing distributions of data and model outputs
library(ggtext)  # For adding rich text elements to ggplot2 plots
library(brms) # For running Bayesian models
library(distributional) # For defining and visualizing probability distributions, including priors
library(marginaleffects) # For extracting marginal effects and predictions from Bayesian models
library(janitor) # For cleaning and organizing data tables, and standardizing column names
library(ggh4x) # For additional ggplot2 features like nested facets
library(cowplot) # For creating multi-panel plots
library(patchwork)
library(rstan)
showtext::showtext_auto()

source("S6_Synthesis_model/functions/helper_functions_S6.R") # Load helper functions specific to synthesis models
source("S6_Synthesis_model/functions/helper_functions_plot_extended.R") # Load helper functions for extended plotting

# Set some global Stan options
base::options(mc.cores = 4,          # Set the number of parallel chains for Bayesian sampling
              brms.backend = "cmdstanr",
              "marginaleffects_posterior_interval" = "hdi",
              "marginaleffects_posterior_center" = "median")   # Use cmdstanr backend for running brms models


# All models will be run with same settings
CHAINS <- 4                          # Number of chains to run for the Bayesian model
ITER <- 5000                         # Total number of iterations per chain
WARMUP <- 2500                       # Number of iterations used for warmup (not included in posterior samples)
BAYES_SEED <- 1234                   # Seed for reproducibility in Bayesian modeling


model_data <- process_model_data("S6_Synthesis_model/data/synthesis_data_complete.xlsx", "direction")

model_data<- model_data |> filter(framework == "Podani",
                     feature == "abun") 

model_data |> 
select(dataset,direction,facet) |> 
  pivot_wider(names_from = facet, values_from = direction) |> 
  count(Traits, Species)

# # A tibble: 4 × 3
#   Traits          Species             n
#   <fct>           <fct>           <int>
# 1 Differentiation Differentiation    32
# 2 Differentiation Homogenisation     14
# 3 Homogenisation  Differentiation    78
# 4 Homogenisation  Homogenisation     37

model_data |>
  tabyl(direction,facet) |>
  adorn_percentages("col") %>%
  adorn_pct_formatting(digits = 2) %>%
  adorn_ns() |>
  kableExtra::kable(format = "rst")

# ===============  ============  ============
# direction        Species       Traits      
# ===============  ============  ============
# Differentiation  68.32% (110)  28.57%  (46)
# Homogenisation   31.68%  (51)  71.43% (115)
# ===============  ============  ============

model_data |>
  tabyl(direction,ecosystem_type, facet) |>
  adorn_percentages("col") %>%
  adorn_pct_formatting(digits = 2) %>%
  adorn_ns() |>
  kableExtra::kable(format = "rst")

#Species
# ===============  ===========  ===========
# direction        Terrestrial  Freshwater 
# ===============  ===========  ===========
# Differentiation  67.94% (89)  70.00% (21)
# Homogenisation   32.06% (42)  30.00%  (9)
# ===============  ===========  ===========
# 
# Traits
# ===============  ===========  ===========
# direction        Terrestrial  Freshwater 
# ===============  ===========  ===========
# Differentiation  25.19% (33)  43.33% (13)
# Homogenisation   74.81% (98)  56.67% (17)
# ===============  ===========  ===========

model_data |>
  tabyl(direction,facet,main_land_use_type) |>
  adorn_percentages("col") %>%
  adorn_pct_formatting(digits = 2) %>%
  adorn_ns()
  
# $Forest
#        direction     Species      Traits
#  Differentiation 67.86% (38) 12.50%  (7)
#   Homogenisation 32.14% (18) 87.50% (49)
# 
# $Multiple
#        direction     Species      Traits
#  Differentiation 71.15% (37) 34.62% (18)
#   Homogenisation 28.85% (15) 65.38% (34)
# 
# $Agriculture
#        direction     Species      Traits
#  Differentiation 73.68% (28) 42.11% (16)
#   Homogenisation 26.32% (10) 57.89% (22)
# 
# $Urban
#        direction    Species      Traits
#  Differentiation 46.67% (7) 33.33%  (5)
#   Homogenisation 53.33% (8) 66.67% (10)


#Set priors --------------------------------------------------------------
priors_direction <- c(
  prior(normal(0,3), class = b) #slopes
 ,prior(student_t(3, 0, 1), class = sd, lb=0) #standard deviation
)

# #Set model formula --------------------------
# 
# Base model with main effects only
formula_base <- bf(
  direction_harrel_davis ~ species_number + spatial_extent + spatial_distance_min +
    absolute_latitude_mean  + human_footprint_baseline + human_footprint_range + habitat_heterogeneity_baseline + habitat_heterogeneity_range + 
  +  main_land_use_type+ecosystem_type + (1 | taxa_coarse)  + hfp_buffer  + het_buffer,
  decomp = "QR"
)  # Base model
  
  formula_with_sigma <- bf(
    direction_harrel_davis ~ species_number + spatial_extent + spatial_distance_min +
      absolute_latitude_mean  + human_footprint_baseline + human_footprint_range + habitat_heterogeneity_baseline + habitat_heterogeneity_range + 
      +  main_land_use_type+ecosystem_type + (1 | taxa_coarse)  + hfp_buffer  + het_buffer,
    sigma ~ hypervolume_quality,
    decomp = "QR"
  )

#Run model ----
# Base Model with Main Effects Only

md_species_base <- brms::brm(
  formula_base,
  data = filter(model_data, facet == "Species"),
  prior = priors_direction,
  family = student(),
  chains = CHAINS,
  iter = ITER,
  warmup = WARMUP,
  control = list(max_treedepth = 12, adapt_delta = .99),
  save_pars = save_pars(all = TRUE),
  seed = BAYES_SEED,
  file = "S7_Model_outputs_figures_and_tables/model/Podani/Abundance/direction/direction_species_base",
  file_refit = "on_change"
)

# Base model for traits

md_traits_base <- brms::brm(
 formula_with_sigma,    # σ increases as quality (MAD) increases),
  data = filter(model_data, facet == "Traits"),
  prior = priors_direction,
  chains = CHAINS,
  family = "student",
  iter = ITER,
  warmup = WARMUP,
  control = list(max_treedepth = 12, adapt_delta = .99),
  save_pars = save_pars(all = TRUE),
  seed = BAYES_SEED,
  file = "S7_Model_outputs_figures_and_tables/model/Podani/Abundance/direction/direction_traits_base",
  file_refit = "on_change"
)

best_model_species <- md_species_base # md_species_base is chosen as the final model because it achieves similar predictive accuracy to more complex models but with fewer parameters
best_model_traits <- md_traits_base # md_traits_base is chosen as the final model because it achieves similar predictive accuracy but with fewer parameters

 # Run the posterior predictive checks -----------------------------------------------------------------------------
cowplot::plot_grid(brms::pp_check(best_model_species, ndraws=100)+coord_cartesian(xlim=c(-20,20), clip="on"), brms::pp_check(best_model_traits, ndraws=100)+coord_cartesian(xlim=c(-20,20), clip="on"))
#'-----------------------------------------------------------------------------------------------------------------

# check R2 of each model ------------------------------------------------------------------------------------------
brms::bayes_R2(best_model_species) # Taxonomic  0.2081661
brms::bayes_R2(best_model_traits) # Functional 0.1590586
#'-----------------------------------------------------------------------------------------------------------------


# Make summary table ----------------------------------------------------------------------------------------------


avg_predictions(best_model_species, 
                newdata = datagrid(ecosystem_type=unique),
                variables="ecosystem_type",
                re_formula=NA,
                conf_level = .8)
 
 # ecosystem_type Estimate 10.0 % 90.0 %
 #    Terrestrial   -1.466  -2.18 -0.775
 #    Freshwater    -0.303  -1.20  0.543


avg_comparisons(best_model_species,
                newdata = datagrid(ecosystem_type=unique),
                variables="ecosystem_type",
                re_formula=NA,,
                conf_level = .8)

# Estimate 10.0 % 90.0 %
#   1.15  0.432   1.82
# 
# Term: ecosystem_type
# Type:  response 
# Comparison: Freshwater - Terrestrial

avg_predictions(best_model_traits, 
                newdata = datagrid(ecosystem_type=unique),
                variables="ecosystem_type", 
                re_formula=NA,
                conf_level = .8)

 # ecosystem_type Estimate 10.0 % 90.0 %
 #    Terrestrial    0.615 -0.228   1.34
 #    Freshwater     0.227 -0.662   1.02
# 
# Type: response

avg_comparisons(best_model_traits, 
                newdata = datagrid(ecosystem_type=unique),
                variables="ecosystem_type",
                re_formula=NA,
                conf_level = .8)

 # Estimate 10.0 % 90.0 %
 #   -0.385 -0.805 0.0563

# 
# Term: ecosystem_type
# Type:  response 
# Comparison: Freshwater - Terrestrial

avg_predictions(best_model_species, 
                newdata = datagrid(main_land_use_type=unique),
                re_formula=NA,
                variables = "main_land_use_type", conf_level = .8)

 # main_land_use_type Estimate 10.0 % 90.0 %
 #        Multiple      -1.082 -1.802 -0.264
 #        Agriculture   -0.859 -1.629 -0.075
 #        Forest        -1.466 -2.184 -0.775
 #        Urban          2.793  0.998  4.528

avg_comparisons(
  best_model_species,
  newdata = datagrid(main_land_use_type=unique),
  variables = list(main_land_use_type = "pairwise"),
  re_formula=NA,
  conf_level = .8
) 

 #               Contrast Estimate 10.0 % 90.0 %
 # Agriculture - Multiple    0.221 -0.584 0.9625
 # Forest - Agriculture     -0.600 -1.281 0.0724
 # Forest - Multiple        -0.378 -1.059 0.3249
 # Urban - Agriculture       3.641  1.926 5.5071
 # Urban - Forest            4.216  2.375 6.1659
 # Urban - Multiple          3.834  2.013 5.7462


avg_predictions(best_model_traits, 
                newdata = datagrid(main_land_use_type=unique),
                variables = "main_land_use_type", 
                re_formula=NA,
                conf_level = .8)

 # main_land_use_type Estimate 10.0 % 90.0 %
 #        Multiple       0.294 -0.545  1.062
 #        Agriculture    0.318 -0.432  1.163
 #        Forest         0.615 -0.228  1.343
 #        Urban         -0.575 -1.684  0.576


avg_comparisons(
  best_model_traits,
  newdata = datagrid(main_land_use_type=unique),
  re_formula=NA,
  variables = list(main_land_use_type = "pairwise"),
  conf_level = .8
) 

 #               Contrast Estimate  10.0 %  90.0 %
 # Agriculture - Multiple   0.0352 -0.4305  0.4524
 # Forest - Agriculture     0.2821 -0.1382  0.6808
 # Forest - Multiple        0.3184 -0.0485  0.6951
 # Urban - Agriculture     -0.8993 -1.8385  0.0330
 # Urban - Forest          -1.1819 -2.1698 -0.1815
 # Urban - Multiple        -0.8683 -1.8537  0.0728
# 
# Term: main_land_use_type
# Type:  response 



#Exploratory figures -------------------------------------------------------------
#Taxonomic turnover

# Predictions ecosystem_type ----


species_draws_eco <- predictions(
  best_model_species,
  newdata = datagrid(ecosystem_type=unique),
  re_formula = NA,
  type = "response"
) |>
  posterior_draws() |>
  mutate(facet = "Taxonomic turnover")

trait_draws_eco <- predictions(
  best_model_traits,
  newdata = datagrid(ecosystem_type=unique),
  re_formula = NA,
  type = "response"
) |>
  posterior_draws() |>
  mutate(facet = "Functional turnover")

pred_eco_type <- bind_rows(species_draws_eco, trait_draws_eco) |>
  select(draw, ecosystem_type, facet) |>
  mutate(
    facet = factor(facet,
                   levels = c("Taxonomic turnover",
                              "Functional turnover")),
    ecosystem_type = factor(ecosystem_type, levels = c("Freshwater","Terrestrial"))
  )|> 
  mutate(level="Ecosystem type")


#Land use values 


species_draws_lu <- predictions(
  best_model_species,
  newdata = datagrid(main_land_use_type=unique),
  re_formula = NA,
  type = "response"
) |>
  posterior_draws() |>
  mutate(facet = "Taxonomic turnover")

trait_draws_lu <- predictions(
  best_model_traits,
  newdata = datagrid(main_land_use_type=unique),
  re_formula = NA,
  type = "response"
) |>
  posterior_draws() |>
  mutate(facet = "Functional turnover")

pred_main_lu <- bind_rows(species_draws_lu, trait_draws_lu) |>
  select(draw, main_land_use_type, facet) |>
  mutate(
    facet = factor(facet,
                   levels = c("Taxonomic turnover",
                              "Functional turnover")),
    main_land_use_type = factor(main_land_use_type, levels = rev(c("Agriculture", "Forest", "Urban", "Multiple")))
  ) |> 
  mutate(level="Main land use")


# Define colors outside the individual figures
blue_green <- c(Freshwater = "#118ab2", Terrestrial = "#06d6a0")
land_use_colors <- c(Agriculture = "#a36627", Forest = "#848c04", Urban = "#1c1c0c", Multiple = "#dc7c5c")


annotations_df <- data.frame(
  label = c(
    "<span>&larr; Dissimilarity</span>",
    "<span>Similarity &rarr;</span>",
    "<span>&larr; Dissimilarity</span>",
    "<span>Similarity &rarr;</span>"
  ),
  x = c(-0.03, 0.03, -7, -7),
  y = c(6.5, 6.5, -0.2, 0.2),
  angle = c(0, 0, 90, 90),
  hjust = c(1, 0, 1, 0)
)

p1<-model_data |>
  select(dataset, facet, direction_harrel_davis) |>
  pivot_wider(names_from = facet, values_from = direction_harrel_davis) |>
  ggplot(aes(y = Species, x = Traits)) +
  geom_richtext(
    data = annotations_df,
    aes(x = x, y = y, label = label, angle = angle, hjust = hjust), 
    size=convert_size(12),
    label.color = NA,
    inherit.aes=FALSE
  ) +
  geom_point(size = 2, alpha = .7) +
  coord_cartesian(xlim = c(-7, NA), ylim = c(-10, 7),clip="on") +
  geom_hline(yintercept = 0,
             linetype = 3,
             linewidth = .1) +
  geom_vline(xintercept = 0,
             linetype = 3,
             linewidth = .1) +
  ggplot2::theme_void(base_family = "sans", base_size = 12) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    plot.background = ggplot2::element_rect(fill = "white", color = NA),
    strip.background = ggplot2::element_rect(fill = NA, color = NA),
    strip.text = ggplot2::element_text(
      face = "bold",
      margin = margin(b = 5)
    ),
    axis.title.x = ggtext::element_markdown(
      family = "sans",
      hjust = 0.5,
      margin = margin(t = 5)
    ),
    axis.title.y = ggtext::element_markdown(
      family = "sans",
      angle = 90,
      margin = margin(r = 4)
    ),
    axis.text.x = element_markdown(
      family = "sans",
      color = "grey30",
      margin = margin(t = 5)
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
    plot.tag = element_text(face="bold",size=14),
    legend.position = "none",
  ) +
  labs(
    x   = "Standardised difference in model fit<br><br><b>Functional turnover</b>",
    y   = "<b>Taxonomic turnover</b><br><br>Standardised difference in model fit",
    tag = "a"
  ) +
  scale_x_continuous(breaks=c(-6.5, -4.5, -2.5, 0, 2.5))+
  guides(
    x = guide_axis(cap = TRUE),
    y = guide_axis(cap = TRUE)
  )

annotations_p2 <- data.frame(
  label = c(
    "<span>&larr; Dissimilarity</span>",
    "<span>Similarity &rarr;</span>"
  ),
  x      = rep(c(-0.03,  0.03),2),   # tweak these to your actual x‐limits
  y      = c( 2.5,     2.5,2.5,2.5),   # y position above your points
  angle  = c( 0,     0,0,0   ),
  hjust  = c( 1,     0,1,0   ),
  facet  = c("Taxonomic turnover", "Taxonomic turnover", "Functional turnover","Functional turnover")  # or repeat for "Trait replacement"
) |>  mutate(facet = factor(facet,levels=c("Taxonomic turnover","Functional turnover")))


p2<- ggplot(pred_eco_type, aes(y = ecosystem_type, x = draw, colour = ecosystem_type)) +
  stat_pointinterval(.width = c(.8, .9, .95)) +
  geom_richtext(
    data = annotations_p2,
    aes(x = x, y = y, label = label, angle = angle, hjust = hjust), 
    size=convert_size(10),
    label.color = NA,
    inherit.aes=FALSE
  ) +
  geom_vline(xintercept = 0, linetype = 3) +
  facet_grid(~facet)+
  scale_colour_manual(values = c(blue_green, land_use_colors)) +
  ggplot2::theme_void(base_family = "sans", base_size = 14) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    plot.background = ggplot2::element_rect(fill = "white", color = NA),
    strip.background = ggplot2::element_rect(fill = NA, color = NA),
    strip.text = ggtext::element_markdown(family = "sans", face = "bold", hjust = 0.5, size = 12, margin=margin(t=2)),
    axis.title.x = ggtext::element_markdown(family = "sans", hjust = 0.5, size = 12, margin=margin(t=2)),
    axis.text.y = ggtext::element_markdown(family = "sans", size = 12, margin=margin(t=5),hjust=1),
    axis.title.y = element_blank(),
    axis.text.x = ggtext::element_markdown(family = "sans", size = 12, margin=margin(t=5)),
    panel.spacing.x = unit(0.2, "lines"),
    panel.spacing.y = unit(0.2, "lines"),
    axis.line.x = element_line(color = "grey20", linewidth=0.1),
    axis.ticks.x = element_line(color = "grey20", linewidth=0.1),
    axis.ticks.length.x = unit(0.3, "lines"),
    plot.margin = margin(3, 0, 3, 0, unit="mm"),
    legend.position="none",
    plot.tag = element_text(face="bold",size=14)
  )+
  scale_y_discrete(expand=c(.1,1), labels=c("Freshwater (n=30)","Terrestrial (n=131)"))+
  labs(x="",tag="b")+
  guides(
    x = guide_axis(cap = TRUE)
  )+
  coord_cartesian(xlim=c(-3,3),clip="off")

p3<- ggplot(pred_main_lu, aes(y = main_land_use_type, x = draw, colour = main_land_use_type)) +
  stat_pointinterval(.width = c(.8, .9, .95)) +
  geom_vline(xintercept = 0, linetype = 3) +
  facet_grid(~facet)+
  scale_colour_manual(values = c(blue_green, land_use_colors)) +
  ggplot2::theme_void(base_family = "sans", base_size = 14) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    plot.background = ggplot2::element_rect(fill = "white", color = NA),
    strip.background = ggplot2::element_rect(fill = NA, color = NA),
    strip.text = element_blank(),
    axis.title.x = ggtext::element_markdown(family = "sans", hjust = 0.5, size = 12, margin=margin(t=5)),
    axis.text.y = ggtext::element_markdown(family = "sans", size = 12, margin=margin(t=2),hjust=1),
    axis.title.y = element_blank(),
    axis.text.x = ggtext::element_markdown(family = "sans", size = 12, margin=margin(t=5)),
    panel.spacing.x = unit(0.2, "lines"),
    panel.spacing.y = unit(0.2, "lines"),
    axis.line.x = element_line(color = "grey20", linewidth=0.1),
    axis.ticks.x = element_line(color = "grey20", linewidth=0.1),
    axis.ticks.length.x = unit(0.3, "lines"),
    plot.margin = margin(3, 0, 3, 0, unit="mm"),
    legend.position="none",
    plot.tag = element_text(face="bold",size=14)
  )+
  scale_y_discrete(labels=rev(c("Agriculture (n=38)","Forest (n=56)", "Urban (n=15)","Multiple (n=52)")))+
  labs(x="Standardised difference in model fit<br><br><b>Posterior predictions</b>",tag="c")+
  guides(
    x = guide_axis(cap = TRUE)
  )+
  coord_cartesian(xlim=c(-5,8),clip="off")

Figure_2 <-  p1|(p2/p3)

ggsave(
  plot = Figure_2,
  filename = here::here("S7_Model_outputs_figures_and_tables", "main_figures", "Figure_2.pdf"),
  device = cairo_pdf(),
  width=12,
  height=5,
  units="in"
)

#Extended Data

direction_figure_data <- generate_figure_data(best_model_species, best_model_traits, "direction", draw_min = -3, draw_max =8 )

library(dplyr)
library(ggplot2)
library(ggh4x)

# 1. Tag each row as either a draw or a text annotation
draws_df <- direction_figure_data %>%
  mutate(layer = " ")

text_df <- direction_figure_data %>%
  distinct(parameter, facet, text, draw_min, draw_max) %>%
  # position text slightly beyond the max draw
  mutate(posterior = draw_max * 1.05,
         layer     = "  ")

plot_df <- bind_rows(draws_df, text_df) %>%
  mutate(layer = factor(layer, levels = c(" ", "  ")))

# 2. Compute x‐limits
draw_limits <- range(plot_df$draw_min, plot_df$draw_max, na.rm = TRUE)
# give the text panel a bit of extra room to the right
text_limits <- c(draw_limits[1], draw_limits[2] * 1.15)
# 3. Build the plot
Extended_Figure_3 <- ggplot(plot_df) +
  facet_nested(
    ~facet+layer,
    scales      = "free_x",
    space       = "free_x",
    strip       = strip_nested(by_layer_y = TRUE)    # ← drop the Posterior/Text strip
  ) +
  # posterior draws + intervals
  stat_pointinterval(
    data    = ~ .x |>  filter(layer == " "),
    aes(x = posterior, y = parameter, colour = significant),
    .width = c(.8, .89, .95)
  ) +
  # text labels
  geom_text(
    data    =~ .x |> filter(layer == "  "),
    aes(x = posterior, y = parameter, label = text),
    hjust   = 0,
    size    = 3
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  # nested facets: rows = Posterior/Text, cols = your facet (Species vs Trait)
  # per‐strip x scales
  facetted_pos_scales(
    x = list(
      # for the Posterior strip
      layer == " " ~ scale_x_continuous(
        name   = "Posterior draws",
        expand = expansion(mult = 0.02)
      ),
      # for the Text strip
      layer == "  " ~ scale_x_continuous(
        limits = c(8,13),
        breaks=NULL,
        name   = NULL,
        expand = expansion(mult = c(0, 0.02))
      )
    )
  ) +
  scale_colour_manual(values = c("gray60", "black")) +
  theme_void(base_family = "sans", base_size = 12) +
  theme(
    panel.grid.minor   = element_blank(),
    plot.background    = element_blank(),
    panel.background =  element_blank(),
    strip.text.y = element_text(face="bold",angle=90,size=12),
    strip.text.x = element_text(face="bold",size=12),,
    axis.title.x       = ggtext::element_markdown(face = "bold", size = 12),
    axis.text.y        = ggtext::element_markdown(size = 10, hjust = 1),
    axis.text.x        = ggtext::element_markdown(size = 12),
    axis.title.y       = element_blank(),
    axis.line.x        = element_line(color = "grey20", linewidth = 0.1),
    axis.ticks.x       = element_line(color = "grey20", linewidth = 0.1),
    axis.ticks.length  = unit(0.1, "lines"),
    legend.position    = "none",
    plot.tag           = element_text(face = "bold", size = 14),
    panel.spacing.x = unit(0,"lines"),
    plot.margin = margin(2,2,2,2, unit="mm"),
    plot.title= element_text(face="bold",size=12)
  ) +
  guides(
    x = guide_axis(cap = TRUE)
  ) +
  scale_x_continuous(breaks=seq(-2.5,7.5,by=2.5))

ggsave(filename = here::here("S7_Model_outputs_figures_and_tables", "extended_data", "Extended_Figure_3.pdf"),
       plot = Extended_Figure_3,
       device = cairo_pdf,
       width=10.5,height=5,units="in")






