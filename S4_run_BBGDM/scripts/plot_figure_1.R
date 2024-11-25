# --------------------------------------------------------------------------------------------------------------
# Script Name: plot_figure_1.R
# Purpose: Generate Figure 1 for manuscript submission, including two subplots:
#          - Figure 1A: Effect size and predictor relationships across species/trait replacement.
#          - Figure 1B: Parlament plot summarizing predictor effects on biotic groups.
# Dependencies:
#          - Data: Processed outputs from BBGDM analyses and dataset metadata.
#          - Libraries: tidyverse, ggbrace, ggpattern, ggnewscale, cowplot, etc.
# Outputs:
#          - PDF file: S8_Model_outputs_figures_and_tables/main_figures/Figure_1.pdf
# Author: Caio Graco-Roza
# Date: 2024-11-24
# --------------------------------------------------------------------------------------------------------------

library(tidyverse)
library(ggbrace) #remotes::install_github("nicolash2/ggbrace")
library(ggpattern) #remotes::install_github("coolbutuseless/ggpattern", force=FALSE)
library(ggnewscale)
library(colorspace) #for manipulating colours in figures
library(ggtext) # for formatted titles and subtitles
library(showtext) #for using alternative fonts in ggplot2
library(ggimage) # add sillhouette
library(ggh4x) #Extra features for plots
library(cowplot)#make plot panels

# extrafont::font_import(prompt = FALSE)  # This might take a while
#   extrafont::loadfonts() 
# 
# showtextdb::font_install(showtextdb::google_fonts("sans"))
# showtextdb::font_install(showtextdb::google_fonts("sans"))
showtext::showtext_auto(TRUE)

source("S4_run_BBGDM/functions/helper_functions_S4.R") #get silly functions to help
# Combine data for parlament plot  --------------------------------------------------------------------------------

files <- tools::file_path_sans_ext(list.files("S4_run_BBGDM/bbgdm_output"))
data_summary <- readxl::read_xlsx("S1_Preprocessing/Miscellaneous/dataset_info_all.xlsx") 
bbgdms_output <- purrr::map(files, ~ read_rds(glue::glue("S4_run_BBGDM/bbgdm_output/{.x}.rds"))) |> 
  set_names(gsub("_bbgdm", "", files)) |> #set names to bbgdm results
  discard_at(c("N67TTP", "N78TTP")) #Remove datasets that do not cross land use gradient
# SUMMARY ##### ------------------------------------------------------------------------------------------

direction_colors <-  c("#ff7f32", "#9f5cc0")
ecosystem_type_colors <- c("#118ab2", "#06d6a0")

#source("S7_Synthesis_model/functions/helper_functions_S7.R", echo=FALSE)

my_theme <- function() {
  ggplot2::theme_void(base_family = "sans", base_size=10) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   plot.background = ggplot2::element_rect(fill = "white", color = NA),
                   strip.background = ggplot2::element_rect(fill = NA, color = NA),
                   axis.title.x = ggtext::element_markdown(family= "sans", face="bold", hjust = 0, margin = margin(t=10)),
                   axis.title.y = ggtext::element_markdown(family= "sans", face= "bold", hjust = 0),
                   legend.title = ggplot2::element_text(family="sans", face = "bold"),
                   axis.text.x = element_markdown(family = "sans", color = "grey30", margin = margin(t = 5), hjust=.9),
                   axis.text.y = element_markdown(family = "sans", color = "grey30", margin = margin(r = 10), hjust=1),  # Added right margin
                   panel.spacing.x = unit(0.5, "lines"),  # Reduced facet spacing (horizontal)
                   panel.spacing.y = unit(0.5, "lines"),  # Reduced facet spacing (vertical)
                   axis.line.x = element_line(color = "grey70"),
                   axis.ticks.x = element_line(color = "grey70"),
                   axis.ticks.length.x = unit(.2, "lines"),
                   axis.line.y = element_line(color = "grey70"),
                   axis.ticks.y = element_line(color = "grey70"),
                   axis.ticks.length.y = unit(.2, "lines"),
                   plot.title = element_text(family = "sans", size = 84, 
                                             color = "grey10", hjust = 0, 
                                             margin = margin(15, 0, 30, 0)),
                   plot.subtitle = element_markdown(family = "sans", size = 25,
                                                    color = "grey30", lineheight = 1.2,
                                                    hjust = 0, margin = margin(0, 0, 30, 0)),
                   plot.caption = element_markdown(family = "sans", color = "grey30", 
                                                   size = 19, face = "italic", 
                                                   hjust = .5, margin = margin(60, 0, 0, 0)),
                   plot.title.position = "plot",
                   plot.caption.position = "plot",
                   plot.margin = margin(60, 120, 45, 90)
    )
}

# Taxonomic ----------------------------------------------------------------------------------------------

# The `bbgdms_output` object is a nested list with the following structure: 
#lobstr::tree(data, index_unnamed = TRUE, max_length=22, max_depth = 4, show_attributes =FALSE)
# <list>
#   ├─barbaro_birds_1: <list> #dataset modelled
#   │ ├─Taxonomic: <list> #facet modelled
#   │ │ ├─differentiation: <list> #direction modelled
#   │ │ │ ├─intercept<dbl [1,000]>: 0.397258592889227, 0.33688530664247, 0.365121075080379, 0.371143120760553, 0.367574324886847, 0.371114461789452, 0.390882060046496, 0.397754441978083, 0.377890645208343, 0.339788719233044, ......
#   │ │ ├─coefs: <list>...
#   │ │ └─r2<dbl [1,000]>: 1.36060207276438, 5.18530504900768, 2.26887496189099, 4.96445929073578, 4.29128938977275, 2.73097996924366, 2.2565718987008, 0.719167076708205, 3.00024053196958, 1.46016377001666, ......
#   │ └─homogenization: <list>
#   │ │   ├─intercept<dbl [1,000]>: 1.06242410234085, 1.04944756338073, 1.06011393379697, 1.09410130085229, 1.06238094378768, 1.04738111754974, 1.05079011515717, 1.06302036510654, 1.06510158769613, 1.05888047919068, ......
#   │   ├─coefs: <list>...
#   │   └─r2<dbl [1,000]>: 2.37587178081908, 0.427849662757995, 0.656967913409756, 1.05227157755122, 0.0598293761075008, 0.178305884868513, 0.904675500354601, 1.57421716917264, 0.356187424210441, 0.0954000741912653, ......
#   └─Functional: <list>
#   │   ├─differentiation: <list>
#   │   │ ├─intercept<dbl [1,000]>: 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ......
#   │   │ ├─coefs: <list>...
#   │   │ └─r2<dbl [1,000]>: 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ......
#   │   └─homogenization: <list>
#   │     ├─intercept<dbl [1,000]>: 2.14199263003067, 2.17879728521932, 2.16997847243116, 2.20789267743642, 2.14874187230854, 2.21100489569246, 2.23766446704733, 2.17344263394147, 2.16000198007066, 2.21029066761949, ......
#   │     ├─coefs: <list>...
#   │     └─r2<dbl [1,000]>: 2.98890245823764, 2.87378131197802, 1.29120540775047, 1.80288251717702, 1.34186027579942, 0.786796283140234, 0.741647228164843, 1.3396641254916, 1.73101022381718, 1.47797480421219, ......
#   ├─barbaro_birds_2: <list>
#   │ ├─Taxonomic: <list>
#   ... 
# Here we use the `map()` function from the `purrr` package to apply a series
# of data manipulations to each nested list element

## Begin example explanation ----
## The figure below shows the relationship between
## the mean effect size (magnitude) and the mean effect size without weighting (magnitude_2)
## for each predictor, grouped by facet and direction. 
## The magnitude represents the average effect size calculated by 
## summing the coefficients (coefficient.1, coefficient.2, and coefficient.3) 
## and multiplying them by the number of observations (n) within each group. 
## On the other hand, the magnitude_2 is the average effect size without considering the number of valid bootstraps.
# 
# bbgdms_output |> 
#   enframe(name = "dataset", value = "dataset_value") |> 
#   mutate(dataset_value = map(dataset_value, ~ enframe(.x, name = "facet", value = "facet_value"))) |> 
#   unnest(dataset_value, keep_empty = TRUE) |> 
#   mutate(facet_value = map(facet_value, ~ enframe(.x, name = "direction", value = "direction_value"))) |> 
#   unnest(facet_value, keep_empty = TRUE) |> 
#   mutate(coef = map(direction_value, ~ .x |>  
#                       pluck("coefs") |>  
#                       bind_rows() |>  
#                       rownames_to_column("predictor")),.keep="unused") |>
#   unnest(coef) |>
#   mutate(predictor = str_remove(predictor,"\\.{3}.*")) |> 
#   add_count(dataset,facet,direction,predictor) |> 
#   mutate(effect_size = (coefficient.1+coefficient.2+coefficient.3)*n,
#          effect_size_2 = (coefficient.1+coefficient.2+coefficient.3), .keep="unused") |>
#   group_by(dataset,facet,direction,predictor) |> 
#   summarise(magnitude = mean(effect_size), magnitude_2=mean(effect_size_2)) |> 
#   ggplot(aes(x=magnitude,y=magnitude_2)) + 
#   geom_point() + 
#   facet_wrap(facet~direction)
#
# The points in the plot represent the relationship between the two measures. 
# As you can see, there are cases where the mean effect size (magnitude) is lower than the mean effect size without weighting (magnitude_2).
# This difference arises due to the presence of datasets with a lower number of valid bootstraps.
# 
# When functional differentiation is analyzed without considering the number of valid bootstraps, 
# datasets with a few valid bootstraps may have coefficients that appear higher simply due to chance. 
# On the other hand, datasets with a larger number of valid bootstraps provide more reliable coefficient estimates. 
# By weighting the coefficients based on the number of valid bootstraps, 
# we can account for this variation in reliability and obtain a more accurate measure of the mean effect size.
# 
# Therefore, the weighting by the number of valid bootstraps is essential to ensure that datasets
# with a greater number of valid bootstraps have a larger influence on the mean effect size estimation.
# It helps mitigate the impact of datasets with a small number of valid bootstraps, 
# ensuring a more robust and reliable estimate of the mean effect size for each predictor, particularly
# in the context of functional differentiation.

mean_effect_size <- bbgdms_output |> 
  enframe(name = "dataset", value = "dataset_value") |> # Convert list to tibble with dataset names as a column
  mutate(dataset_value = purrr::map(dataset_value, ~ enframe(.x, name = "facet", value = "facet_value"))) |> # Map over each dataset, creating nested tibbles with facets
  unnest(dataset_value, keep_empty = TRUE) |> # Expand the nested tibble for each dataset
  mutate(facet_value = purrr::map(facet_value, ~ enframe(.x, name = "direction", value = "direction_value"))) |> # Map over each facet, creating nested tibbles with directions
  unnest(facet_value, keep_empty = TRUE) |> # Expand the nested tibble for each facet
  mutate(coef = purrr::map(direction_value, ~ .x |> pluck("coefs") |> bind_rows() |> rownames_to_column("predictor")),.keep="unused") |> # Extract coefficients and bind them into a single dataframe
  unnest(coef) |> # Expand the coefficients into individual rows
  mutate(predictor = str_remove(predictor,"\\.{3}.*")) |> # Clean predictor names by removing suffixes
  mutate(effect_size = (coefficient.1 + coefficient.2 + coefficient.3), .keep="unused") |> # Calculate total effect size from coefficients
  group_by(dataset, facet, direction, predictor) |> # Group by dataset, facet, direction, and predictor
  dplyr::summarise(magnitude = median(effect_size) * (n() / 1000)) |> # Summarize magnitude as scaled median of effect size
  separate(predictor, into = c("predictor", "buffer")) # Separate predictor names into predictor and buffer columns


# iterate through the nested lists and extract the mean r2 values
# Calculate the mean r2 value for each combination of dataset, facet, and direction

mean_r2 <- bbgdms_output |> 
  enframe(name = "dataset", value = "dataset_value") |> # Convert list to tibble with dataset names as a column
  mutate(dataset_value = purrr::map(dataset_value, ~ enframe(.x, name = "facet", value = "facet_value"))) |> # Map over each dataset, creating nested tibbles with facets
  unnest(dataset_value, keep_empty = TRUE) |> # Expand the nested tibble for each dataset
  mutate(facet_value = purrr::map(facet_value, ~ enframe(.x, name = "direction", value = "direction_value"))) |> # Map over each facet, creating nested tibbles with directions
  unnest(facet_value, keep_empty = TRUE) |> # Expand the nested tibble for each facet
  dplyr::mutate(r2 = purrr::map(direction_value, ~ .x |> pluck("r2")), .keep="unused") |> # Extract R² values from each direction
  unnest(r2) |> # Expand the R² values into individual rows
  mutate(across(c(r2), ~ ifelse(.x <= 0, NA, .x / 100))) |> # Normalize R² values and replace non-positive values with NA
  group_by(dataset, facet) |> # Group by dataset and facet
  dplyr::summarise( # Summarize the R² difference between directions
    direction_r2 = gammaEffectSize(
      y = r2[direction == "homogenization"], # R² values for homogenization direction
      x = r2[direction == "differentiation"], # R² values for differentiation direction
      prob = .5 # Probability level for the effect size calculation
    )
  )


# mean_r2 |>
#   mutate(dir = ifelse(sign(direction_r2) > 0 , "Homogenization","differentiation")) |>
#   janitor::tabyl(facet,dir) |>
#   knitr::kable()

# |facet      | differentiation| Homogenization|
# |:----------|---------------:|--------------:|
# |Functional |              89|             71|
# |Taxonomic  |             118|             42|

# Summarizing R² and effect size results for each dataset and direction, aligning with predictors for analysis or visualization.
bbgdm_summary <- mean_r2 |> 
  mutate(direction = ifelse(sign(direction_r2) > 0 , "homogenization","differentiation")) |> # Determine the direction based on the sign of the R² effect size
  left_join(mean_effect_size, by = c("dataset", "facet", "direction")) |> # Combine R² and effect size data for each dataset, facet, and direction
  separate(predictor, into = c("predictor"), extra = "drop") |> # Clean up predictor column by dropping additional information
  mutate(predictor = factor(predictor, levels = rev(c("hfp", "modis", "het", "Geographic", "Temp", "Prec")))) |>  # Order predictors for plotting or further analysis
  mutate(facet = str_replace_all(facet, c("Taxonomic" = "Species replacement", "Functional" = "Trait replacement"))) |> # Rename facet levels for readability
  group_by(dataset, direction) %>% # Group by dataset and direction for calculations
  mutate(effect_size = magnitude / sum(magnitude)) %>% # Calculate effect size as a proportion of total magnitude
  ungroup() %>% # Ungroup for operations on the full dataset
  arrange(dataset, predictor) 

# Creates a named list of formatted HTML strings for facet labels, each including:
# - A bold facet title ("Species replacement" or "Trait replacement").
# - Counts of homogenization (in purple) and differentiation (in orange) directions dynamically computed using `count_directions`.
# - Styled for use in plots or reports to enhance interpretability and aesthetics.
facet_labels <- c(
  "Species replacement" = glue::glue("<span style='font-size:10pt;'><b>Species replacement</b></span><br>
                                     <span style = 'color:#9f5cc0; font-size:8pt;'>Homogenisation = {count_directions(bbgdm_summary,'Species replacement','homogenization')}</span><br>
                                     <span style='color:#ff7f32; font-size:8pt;'>Differentiation = {count_directions(bbgdm_summary,'Species replacement','differentiation')}</span>"), 
  "Trait replacement" = glue::glue("<span style='font-size:10pt;'><b>Trait replacement</b></span><br>
                                   <span style = 'color:#9f5cc0; font-size:8pt;'>Homogenisation = {count_directions(bbgdm_summary,'Trait replacement','homogenization')}</span><br>
                                   <span style='color:#ff7f32; font-size:8pt;'>Differentiation = {count_directions(bbgdm_summary,'Trait replacement','differentiation')}</span>")
)

Figure_1A <- bbgdm_summary %>%
  ggplot(aes(x = effect_size, y = predictor, fill = direction)) + # Set up ggplot with effect size and predictor, colored by direction
  facet_wrap(~facet, labeller = as_labeller(facet_labels), ncol = 2) + # Add facets for species and trait replacement
  geom_jitter(aes(colour = direction, group = direction), alpha = 0.05, shape = 19, size = 1, position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0)) + # Add jittered points for distribution
  stat_summary(aes(group = direction), fun.data = "mean_cl_boot", geom = "linerange", linewidth = .7, position = position_dodge(0.9)) + # Add confidence intervals for group means
  stat_summary(aes(group = direction, fill = direction), fun = "mean", geom = "point", size = 1.5, shape = 21, position = position_dodge(0.9), stroke = .7) + # Add group mean points
  scale_y_discrete(labels = rev(c("<b>Human Footprint (HFP)</b>", "<b>Human Land Cover (HLC)</b>", "<b>Habitat Heterogeneity (HH)</b>", "Spatial Distance (SD)", "Temperature (T)", "Precipitation (PPT)"))) + # Set custom labels for predictors
  scale_fill_manual(values = direction_colors) + # Define manual fill colors for directions
  scale_colour_manual(values = direction_colors) + # Define manual line colors for directions
  labs(
    x = "Magnitude of effect (standardised)", # Label for x-axis
    y = "" # Remove label for y-axis
  ) +
  my_theme() + # Apply custom theme
  theme(
    axis.text = element_text(family = "sans", size = 8), # Customize axis text appearance
    axis.title.x = element_markdown(size = 10, hjust = 0.5, margin = margin(t = 10)), # Customize x-axis title
    axis.text.y = element_markdown(), # Use markdown for y-axis text
    strip.text = ggtext::element_markdown(), # Use markdown for facet strip text
    strip.text.x = element_markdown(size = 10, hjust = 0), # Customize facet strip title appearance
    plot.title.position = "plot", # Align plot title with the plot
    legend.position = "none", # Remove legend
    plot.tag = element_text(face = "bold"), # Customize plot tag appearance,
    plot.margin = unit(c(t = 0, l = 2, b = 1, r = 2), "lines"), panel.spacing = unit(0.5, "cm")
  ) +
  guides(x = guide_axis_truncated(trunc_lower = 0, trunc_upper = 1), # Add truncated axis guide for x-axis
         y = guide_axis_truncated(trunc_lower = 1, trunc_upper = 6)) + # Add truncated axis guide for y-axis
  coord_cartesian(xlim = c(0, 1), clip = "off") + # Set Cartesian limits for x-axis and disable clipping
  scale_x_continuous(breaks = seq(0, 1, by = .25), expand = c(.1, 0), limits = c(0, 1)) # Customize x-axis breaks and limits

#Summary statistics of HFP using bootstrap for the text
# bbgdm_summary %>%
#   filter(predictor == "hfp") %>% # Filter for the predictor "hfp"
#   select(facet, direction, effect_size) %>% # Select relevant columns
#   group_by(facet) %>% # Group by facet and direction
#   reframe({
#     stats <- smean.cl.boot(effect_size, conf.int = .95, B = 1000, na.rm = TRUE) # Calculate stats once
#     tibble::tibble(
#       Mean = stats[1],
#       Lower = stats[2],
#       Upper = stats[3]
#     )
#   }) |>
#   knitr::kable()
# |facet               |      Mean|     Lower|    Upper|
# |:-------------------|---------:|---------:|--------:|
# |Species replacement | 0.2145276| 0.1810296| 0.248589|
# |Trait replacement   | 0.2056555| 0.1755219| 0.240798|


# PARLAMENT ##### ----------------------------------------------------------------------------------------

# Taxonomic ----------------------------------------------------------------------------------------------

#clean_dataset for plotting 
bbgdm_data_parsed <- bbgdm_summary %>% 
  left_join(data_summary %>% select(dataset_name, system, taxa, group, disturbance), by = c("dataset" = "dataset_name")) %>% # Add metadata columns from `data_summary` based on dataset name
  mutate(direction = factor(direction, levels = c("homogenization", "differentiation"))) %>% # Set the order of levels for the direction factor
  mutate(system = factor(system, levels = c("aquatic", "terrestrial"))) %>% # Set the order of levels for the system factor
  mutate(predictor = factor(predictor, levels = rev(c("hfp", "modis", "het", "Geographic", "Temp", "Prec")))) %>% # Set the order of levels for the predictor factor
  mutate(disturbance = factor(disturbance, levels = c("agriculture", "forest", "urban", "multiple"))) %>% # Set the order of levels for the disturbance factor
  mutate(taxa = str_to_title(taxa)) |> # Convert taxa names to title case
  mutate(taxa = case_when( # Adjust taxa names 
    taxa %in% c("Fish", "Fungi") ~ taxa,
    TRUE ~ paste0(taxa, "s")
  )) |> 
  mutate(taxa = gsub("Waterinsects", "Aquatic\ninvertebrates", taxa)) |> # Replace "Waterinsects" with "Aquatic invertebrates"
  mutate(taxa = gsub("Insects", "Arthropods", taxa)) |> # Replace "Insects" with "Arthropods"
  mutate(taxa = factor(taxa, levels = c("Microorganisms", "Fish", "Amphibians", "Aquatic\ninvertebrates", "Fungi", "Plants", "Gastropods", "Arthropods", "Birds", "Bats"))) # Set the order of levels for the taxa factor


bbgdm_data_ordered <- bbgdm_data_parsed %>% 
  filter(predictor %in% "hfp") %>% 
  mutate(ordering = case_when(direction %in% "differentiation" ~ magnitude,
                              direction %in% "homogenization" ~ -magnitude)) %>% 
  arrange(facet, system, taxa, disturbance, ordering) %>% distinct(dataset) %>%  pull(dataset)

# bbgdm_data_parla <- bbgdm_data_parsed %>% 
#     mutate(dataset = factor(dataset, levels= bbgdm_data_ordered)) %>% 
#     arrange(dataset) %>% 
#     mutate(data_num = as.numeric(dataset),
#            var_num = as.numeric(predictor),
#            dist_num = as.numeric(disturbance)) %>%  
#     group_by(system, taxa) %>% 
#     mutate(brace_middle = ceiling(mean(c(min(data_num),max(data_num)))),
#            fig_center = 10) %>% 
#     mutate(image = case_when(
#       #these are codes for the phylopic sillhouettes to each taxonomic group included
#       taxa == "Aquatic\ninvertebrates" ~ "e6ba2bd6-b459-4a20-bd7a-4c4bfc6f9eda" 
#       , taxa == "Fish" ~ "84c7e672-2593-44a6-a807-cffbd3156cc5"
#       , taxa == "Microorganisms" ~ "ef65efc6-f6a6-4a9d-adf8-c6d9b56ac119"
#       , taxa == "Amphibians" ~ "7e9286b9-0e49-46c2-88f6-218291fb06a4"
#       , taxa == "Arthropods" ~ "724b30b5-0155-40ef-9d11-5ea201323bfa"
#       , taxa == "Birds" ~ "dfdfb59e-8126-44e1-a7a9-1bf698113e1c"
#       , taxa == "Gastropods" ~ "75d1f2b1-53e2-4896-958d-d5204d81ced9"
#       , taxa == "Plants" ~ "13a5a7e9-be33-4d40-8c90-9ffcf367d7cb"
#       , taxa == "Fungi" ~ "8cff2d66-6549-44d2-8304-d2dfecf53d78"
#       , taxa == "Bats" ~ "a2ef021b-2000-49f0-a36e-916840755480"
#         ))

bbgdm_data_parla <- bbgdm_data_parsed %>% 
  mutate(dataset = factor(dataset, levels = bbgdm_data_ordered)) %>% # Reorder the dataset factor based on a pre-defined order
  arrange(dataset) %>% # Sort the data by dataset
  mutate(
    data_num = as.numeric(dataset), # Convert the dataset factor to numeric for plotting
    var_num = as.numeric(predictor), # Convert the predictor factor to numeric for plotting
    dist_num = as.numeric(disturbance) # Convert the disturbance factor to numeric for plotting
  ) %>%  
  group_by(system, taxa) %>% # Group data by system and taxa for further calculations
  mutate(
    brace_middle = ceiling(mean(c(min(data_num), max(data_num)))), # Calculate the middle position for braces in the plot
    fig_center = 10 # Set a fixed center for figures (used for consistent alignment)
  ) %>% 
  mutate(image = case_when( # Assign unique PhyloPic silhouette codes to each taxonomic group
    taxa == "Aquatic\ninvertebrates" ~ "e6ba2bd6-b459-4a20-bd7a-4c4bfc6f9eda", # Code for Aquatic invertebrates
    taxa == "Fish" ~ "84c7e672-2593-44a6-a807-cffbd3156cc5", # Code for Fish
    taxa == "Microorganisms" ~ "ef65efc6-f6a6-4a9d-adf8-c6d9b56ac119", # Code for Microorganisms
    taxa == "Amphibians" ~ "7e9286b9-0e49-46c2-88f6-218291fb06a4", # Code for Amphibians
    taxa == "Arthropods" ~ "724b30b5-0155-40ef-9d11-5ea201323bfa", # Code for Arthropods
    taxa == "Birds" ~ "dfdfb59e-8126-44e1-a7a9-1bf698113e1c", # Code for Birds
    taxa == "Gastropods" ~ "75d1f2b1-53e2-4896-958d-d5204d81ced9", # Code for Gastropods
    taxa == "Plants" ~ "13a5a7e9-be33-4d40-8c90-9ffcf367d7cb", # Code for Plants
    taxa == "Fungi" ~ "8cff2d66-6549-44d2-8304-d2dfecf53d78", # Code for Fungi
    taxa == "Bats" ~ "a2ef021b-2000-49f0-a36e-916840755480" # Code for Bats
  )) # Assign the silhouette code based on the taxa group

Figure_1B_draft <- 
  bbgdm_data_parla %>% 
  mutate(facet = str_replace_all(facet, c("Taxonomic" = "Species replacement", "Functional" = "Trait replacement"))) %>% # Rename facets for better interpretability
  mutate(direction = gsub("z", "s", str_to_title(direction))) %>% # Format direction values for consistency
  group_by(dataset, facet) %>% # Group by dataset and facet for calculations
  mutate(es_rank = magnitude / sum(magnitude)) %>% # Calculate rank-normalized effect size
  mutate(direction = factor(direction, levels = c("Homogenisation", "Differentiation"))) %>% # Convert direction to factor with specific order
  select(dataset, facet, predictor, direction, magnitude, es_rank, data_num, var_num, taxa, system, disturbance, brace_middle, fig_center, image) %>% # Select relevant columns
  ggplot(aes(x = data_num, y = var_num)) + # Define ggplot base aesthetics
  facet_wrap(~facet) + # Create facets for species and trait replacement
  coord_polar(clip = "off") + # Use polar coordinates for circular visualization
  scale_x_continuous(limits = c(-5, 169), breaks = seq(1, 163, 1)) + # Set x-axis scale
  scale_y_continuous(limits = c(-3, 12), expand = expansion(mult = c(0, -0.3))) + # Set y-axis scale with expansion
  
  # Add braces for each biotic group
  stat_brace(aes(group = taxa), outerstart = 8, bending = 0.4, width = 0.5, linewidth = 0.3) +
  
  # Add blue and green tiles for ecosystem type
  geom_tile(aes(x = data_num, y = var_num, fill = system), alpha = 0.3, colour = "gray70") +
  scale_fill_manual(values = ecosystem_type_colors, labels = c("Freshwater", "Terrestrial"), name = "Ecosystem type") +
  
  # Add labels for predictors
  annotate("text", x = 167, y = 1, label = "PPT", size = convert_size(8), family = "sans", vjust = 0.5) + # Precipitation
  annotate("text", x = 167, y = 2, label = "T", size = convert_size(8), family = "sans") + # Temperature
  annotate("text", x = 167, y = 3, label = "SD", size = convert_size(8), family = "sans") + # Spatial Distance
  annotate("text", x = 167, y = 4, label = "HH", size = convert_size(8), family = "sans", fontface = "bold") + # Habitat Heterogeneity
  annotate("text", x = 167, y = 5, label = "HLC", size = convert_size(8), family = "sans", fontface = "bold") + # Human Land Cover
  annotate("text", x = 167, y = 6, label = "HFP", size = convert_size(8), family = "sans", fontface = "bold") + # Human Footprint
  
  # Add tiles for disturbance types
  ggnewscale::new_scale_fill() + 
  geom_tile(aes(y = 7, x = data_num, fill = disturbance), colour = NA) +
  scale_fill_manual(name = "Land use type", 
                    labels = c("Agriculture", "Forest", "Urban", "Multiple"),
                    values = c("#a36627", "#848c04", "#1c1c0c", "#dc7c5c"),
                    guide = guide_legend(direction = "horizontal", title.position = "top", nrow = 2, byrow = TRUE)) +
  
  # Add dots for effect sizes of predictors on species/trait replacement
  ggnewscale::new_scale_fill() +
  geom_hline(yintercept = seq(1, 6, by = 1), colour = "gray80", linewidth = 0.1, linetype = "dotted") + # Horizontal lines for grid
  geom_point(aes(x = data_num, y = var_num, size = es_rank, fill = direction), stroke = 0.1, shape = 21) + # Points for effect sizes
  scale_fill_manual(values = rev(direction_colors), labels = c("Homogenisation", "Differentiation"), 
                    na.translate = FALSE, name = "Direction", guide = guide_legend(override.aes = list(size = 5))) +
  scale_size_continuous("Magnitude", range = c(0, 2), breaks = c(0, 0.25, 0.5, 1), 
                        guide = guide_legend(override.aes = list(fill = "black"), direction = "horizontal", title.position = "top", nrow = 2, byrow = TRUE)) +
  
  # Add silhouettes for taxa
  geom_phylopic(data = . %>% distinct(brace_middle, .keep_all = TRUE), aes(x = brace_middle, y = fig_center, image = image), 
                size = 0.06, inherit.aes = FALSE) +
  
  # Customize plot theme
  labs(size = "Magnitude") + 
  my_theme() +
  theme(
    panel.grid = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    text = element_text(size = 8, colour = "black"),
    axis.text = element_text(size = 8, colour = "black"),
    axis.text.y =  element_blank(),
    axis.title.y=  element_blank(),
    axis.line.x =  element_blank(),
    axis.line.y=   element_blank(),
    axis.title.x=  element_blank(),
    axis.text.x =  element_blank(),
    axis.ticks.y = element_blank(),
    line = element_line(colour = "black", linewidth = 1),
    strip.text = element_text(colour = "black",size=10, hjust=.5,face="bold",family="sans"),
    ,plot.title.position = "plot"
    ,plot.background = element_blank()
    ,panel.spacing.x = unit(2, "lines")  # Reduced facet spacing (horizontal)
    ,legend.position = 'bottom'
    ,legend.box = 'horizontal'
    ,legend.box.just = 'top'
    ,legend.direction = "vertical"
    ,legend.justification = "center"
    ,legend.key.size = unit(.5, 'cm')
    ,legend.text = element_text(size=8)
    ,legend.title = element_text(size=8, face="bold")
    ,plot.margin = unit(c(0,2,0,2),"lines")
    ,panel.spacing = unit(1.3, "cm")
    ,plot.tag=element_text(face="bold", size=14)
  )   

  # Define nudging adjustments for text labels
  nudge_x <- c(2, 1, 1.5, 2, 2.3, 4, 3, -1, -3, -2) # Horizontal adjustments for text positions
  nudge_y <- c(1, 1, 0.5, 0, 0, 0, 0, 1, 1, -0.3)   # Vertical adjustments for text positions
  
  # Extract and adjust data for text annotations
  Figure_1B_text_data <- ggplot2::layer_data(Figure_1B_draft, 12) |> # Extract data for layer 12 of the draft figure
    left_join(bbgdm_data_parla |> ungroup() |> distinct(taxa, image)) |> # Join with taxa and image information
    mutate(
      label = ifelse(PANEL == 1, as.character(taxa), ""), # Add labels for species replacement facet only
      facet = ifelse(PANEL == 1, "Species replacement", "Trait replacement") # Assign facets based on PANEL
    ) |> 
    distinct(PANEL, image, .keep_all = TRUE) |> # Retain unique PANEL and image combinations
    group_by(PANEL) |> # Group by PANEL to apply nudging adjustments
    mutate(
      x = x + nudge_x, # Apply horizontal nudging
      y = y + nudge_y  # Apply vertical nudging
    )
  
  # Add text annotations to Figure_1B
  Figure_1B <- Figure_1B_draft + 
    geom_text(
      data = Figure_1B_text_data, # Use prepared text data
      aes(y = y, x = x, label = label), # Map adjusted positions and labels
      hjust = 0, # Align text to the left
      size = convert_size(8), # Adjust text size
      family = "sans" # Use sans-serif font
    )
  
  # Combine Figure_1A and Figure_1B into a single figure
  Figure_1 <- plot_grid(
    Figure_1A + theme(plot.margin = margin(t = 10, b = 10)), # Add margin adjustments to Figure_1A
    Figure_1B + theme(plot.margin = margin(t = 20, b = 10)), # Add margin adjustments to Figure_1B
    ncol = 1, # Arrange plots vertically
    rel_heights = c(0.3, 0.7), # Set relative heights for the two plots
    labels = c("a", "b"), # Add subplot labels
    label_size = 14, # Set size for subplot labels
    label_x = c(0, 0), # Align both labels to the far left
    label_y = c(1, 1) # Align both labels vertically
  )
  
  # Save the final figure as a PDF
  ggsave(
    here::here("S8_Model_outputs_figures_and_tables", "main_figures", "Figure_1.pdf"), # Define save path
    device = cairo_pdf, # Use Cairo PDF device for high-quality output
    Figure_1, # Specify the figure to save
    height = 190, width = 247, units = "mm" # Define dimensions in millimeters
  )
