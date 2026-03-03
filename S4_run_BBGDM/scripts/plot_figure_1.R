# --------------------------------------------------------------------------------------------------------------
# Script Name: plot_figure_1.R
# Purpose: Generate Figure 1 for manuscript submission, including two subplots:
#          - Figure 1A: Effect size and predictor relationships across species/Functional turnover.
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


showtext::showtext_auto(TRUE)

source("S4_run_BBGDM/functions/helper_functions_S4.R") #get silly functions to help

beta_type = "Podani_abun"
# Combine data for parlament plot  --------------------------------------------------------------------------------

files <- tools::file_path_sans_ext(list.files(glue::glue("S4_run_BBGDM/bbgdm_output/{beta_type}")))
data_summary <- readxl::read_xlsx("S1_Preprocessing/Miscellaneous/dataset_info_all.xlsx") 
bbgdms_output <- purrr::map(files, ~ read_rds(glue::glue("S4_run_BBGDM/bbgdm_output/{beta_type}/{.x}.rds"))) |> 
  set_names(gsub("_bbgdm", "", files)) |> 
  discard_at(c("N67TTP", "N78TTP")) #Remove datasets that do not cross land use gradient
# SUMMARY ##### ------------------------------------------------------------------------------------------

direction_colors <-  c("#ff7f32", "#9f5cc0")
ecosystem_type_colors <- c("#118ab2", "#06d6a0")

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

# Define a function to map beta_type to the desired plot title
getPlotTitle <- function(beta_type) {
  # Ensure beta_type is treated as a character string
  beta_type <- as.character(beta_type)
  
  # Use switch to assign the title based on beta_type
  title <- switch(beta_type,
                  "Podani_pa"   = "Presence-absence based  replacement (Podani partition)",
                  "Podani_abun" = "Abundance based  replacement (Podani partition)",
                  "Baselga_pa"  = "Presence-absence based turnover (Baselga partition)",
                  "Baselga_abun"= "Abundance based turnover (Baselga partition)",
                  "Unknown beta_type"  # default if no match is found
  )
  return(title)
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

mean_r2 <- bbgdms_output |> # Read BBGDM results for each file
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


mean_r2 |>
  mutate(dir = ifelse(sign(harrel_davis) > 0 , "Homogenisation","differentiation")) |>
  janitor::tabyl(facet,dir) |>
  knitr::kable()

# |facet      | differentiation| Homogenisation|
# |:----------|---------------:|--------------:|
# |Functional |              47|            115|
# |Taxonomic  |             110|             52|

# Summarizing R² and effect size results for each dataset and direction, aligning with predictors for analysis or visualization.
bbgdm_summary <- mean_r2 |> 
  mutate(direction = ifelse(sign(harrel_davis) > 0 , "homogenization","differentiation")) |> # Determine the direction based on the sign of the R² effect size
  left_join(mean_effect_size, by = c("dataset", "facet", "direction")) |> # Combine R² and effect size data for each dataset, facet, and direction
  separate(predictor, into = c("predictor"), extra = "drop") |> # Clean up predictor column by dropping additional information
  mutate(predictor = factor(predictor, levels = rev(c("hfp", "modis", "het", "Geographic", "Temp", "Prec")))) |>  # Order predictors for plotting or further analysis
  mutate(facet = str_replace_all(facet, c("Taxonomic" = "Taxonomic turnover", "Functional" = "Functional turnover"))) |> # Rename facet levels for readability
  group_by(dataset, direction) %>% # Group by dataset and direction for calculations
  mutate(effect_size = magnitude / sum(magnitude)) %>% # Calculate effect size as a proportion of total magnitude
  ungroup() %>% # Ungroup for operations on the full dataset
  arrange(dataset, predictor) |> 
  mutate(facet = factor(facet,levels=c("Taxonomic turnover","Functional turnover")))



# Creates a named list of formatted HTML strings for facet labels, each including:
# - A bold facet title ("Taxonomic turnover" or "Functional turnover").
# - Counts of homogenization (in purple) and differentiation (in orange) directions dynamically computed using `count_directions`.
# - Styled for use in plots or reports to enhance interpretability and aesthetics.

facet_labels <- c(
  "Taxonomic turnover" = glue::glue("<span style='font-size:10pt;'><b>Taxonomic turnover</b></span><br>
                                     <span style = 'color:#9f5cc0; font-size:10pt;'>Homogenisation = {count_directions(bbgdm_summary,'Taxonomic turnover','homogenization')}</span><br>
                                     <span style='color:#ff7f32; font-size:10pt;'>Differentiation = {count_directions(bbgdm_summary,'Taxonomic turnover','differentiation')}</span>"), 
  "Functional turnover" = glue::glue("<span style='font-size:10pt;'><b>Functional turnover</b></span><br>
                                   <span style = 'color:#9f5cc0; font-size:10pt;'>Homogenisation = {count_directions(bbgdm_summary,'Functional turnover','homogenization')}</span><br>
                                   <span style='color:#ff7f32; font-size:10pt;'>Differentiation = {count_directions(bbgdm_summary,'Functional turnover','differentiation')}</span>")
)

bbgdm_summary |>
  group_by(dataset,facet) |>
  slice_max(magnitude) |>
  janitor::tabyl(predictor,facet)

 #  predictor Taxonomic turnover Functional turnover
 #       Prec                  20                33
 #       Temp                  35                22
 # Geographic                  18                12
 #        het                  30                45
 #      modis                  12                14
 #        hfp                  47                36

  bbgdm_summary |>
  group_by(facet, predictor) |>
  summarise(mean = Hmisc::smean.cl.boot(magnitude)[1],
            low = Hmisc::smean.cl.boot(magnitude)[2],
            high = Hmisc::smean.cl.boot(magnitude)[3])
  
#   facet               predictor    mean    low   high
#    <chr>               <fct>       <dbl>  <dbl>  <dbl>
#  1 Taxonomic turnover Prec       0.135  0.108  0.167 
#  2 Taxonomic turnover Temp       0.160  0.134  0.190 
#  3 Taxonomic turnover Geographic 0.0898 0.0717 0.111 
#  4 Taxonomic turnover het        0.168  0.141  0.202 
#  5 Taxonomic turnover modis      0.0862 0.0642 0.110 
#  6 Taxonomic turnover hfp        0.199  0.158  0.250 
#  7 Functional turnover   Prec       0.142  0.110  0.180 
#  8 Functional turnover   Temp       0.0878 0.0680 0.108 
#  9 Functional turnover   Geographic 0.0636 0.0403 0.0886
# 10 Functional turnover   het        0.219  0.162  0.277 
# 11 Functional turnover   modis      0.0710 0.0493 0.0959
# 12 Functional turnover   hfp        0.172  0.135  0.215 


#   bbgdm_summary |>
#     group_by(dataset,facet) |>
#     mutate(total_effect=sum(magnitude)) |> 
#     mutate(rel_effect=magnitude/total_effect) |> 
#     group_by(facet, predictor) |>
#     summarise(mean = Hmisc::smean.cl.boot(rel_effect)[1],
#               low = Hmisc::smean.cl.boot(rel_effect)[2],
#               high = Hmisc::smean.cl.boot(rel_effect)[3])
#   
#      facet               predictor    mean    low  high
#    <chr>               <fct>       <dbl>  <dbl> <dbl>
#  1 Taxonomic turnover Prec       0.158  0.134  0.185
#  2 Taxonomic turnover Temp       0.200  0.170  0.231
#  3 Taxonomic turnover Geographic 0.113  0.0918 0.138
#  4 Taxonomic turnover het        0.204  0.181  0.231
#  5 Taxonomic turnover modis      0.0969 0.0789 0.116
#  6 Taxonomic turnover hfp        0.228  0.198  0.258
#  7 Functional turnover   Prec       0.182  0.149  0.218
#  8 Functional turnover   Temp       0.142  0.116  0.170
#  9 Functional turnover   Geographic 0.0949 0.0724 0.120
# 10 Functional turnover   het        0.271  0.236  0.308
# 11 Functional turnover   modis      0.0974 0.0756 0.121
# 12 Functional turnover   hfp        0.213  0.181  0.249


 summary_stats <- bbgdm_summary %>%
    group_by(dataset,facet) %>%
    mutate(rel_ef = magnitude / sum(magnitude)) %>%
   group_by(facet, predictor) %>%
    summarise(
      mean_ef = smean.cl.boot(rel_ef, B = 5000)[1],
      lower_ci = smean.cl.boot(rel_ef, B = 5000)[2],
      upper_ci = smean.cl.boot(rel_ef, B = 5000)[3],
      .groups = "drop"
    )
  
  
Figure_1A <- 
  bbgdm_summary %>%
  group_by(dataset,facet) |> 
  mutate(rel_ef = magnitude/sum(magnitude)) |> 
  ggplot(aes(x = rel_ef, y = predictor)) + # Set up ggplot with effect size and predictor, colored by direction
  facet_wrap(~facet, ncol = 2) + # Add facets for species and Functional turnover
  geom_jitter(alpha = 0.05, shape = 19, size = 1) + # Add jittered points for distribution
  
  stat_summary(fun.data = "mean_cl_boot", geom = "linerange", fun.args = list(B = 5000)) + # Add confidence intervals for group means
  stat_summary(fun.data = "mean_cl_boot", geom = "point", size = 2, shape = 19, stroke = .7) + # Add group mean points
  
    # Add labels for mean and CI
    geom_text(
      data = summary_stats,
      aes(
        x = 0.5,  # Fixed x-position for labels
        y = predictor,
        label = sprintf(
          "Mean: %.2f, 95%% CI: [%.2f, %.2f]",  # \n for line break,  # Forces 2 decimal places
          mean_ef, 
          lower_ci, 
          upper_ci
        )
      ),
      hjust = 0,
      size = 3.5,
      color = "black",
      family = "sans"
    )+
    
  scale_y_discrete(labels = rev(c("<b>Human Footprint</b>", "<b>Human Land Cover</b>", "<b>Habitat Heterogeneity</b>", "Spatial Distance", "Temperature", "Precipitation"))) + # Set custom labels for predictors
  labs(
    x = "Cumulative fitted turnover (standardised)", # Label for x-axis
    y = "" # Remove label for y-axis
  ) +
    theme_void(base_family = "sans", base_size = 14) +
    theme(
      panel.grid.minor   = element_blank(),
      plot.background    = element_blank(),
      panel.background =  element_blank(),
      strip.text = element_text(face="bold", size=14),
      axis.title.x       = ggtext::element_markdown(face = "bold", size = 12),
      axis.text.y        = element_markdown(size = 12, hjust = 1),
      axis.text.x        = ggtext::element_markdown(size = 12),
      axis.title.y       = element_blank(),
      axis.line.x        = element_line(color = "grey20", linewidth = 0.1),
      axis.ticks.x       = element_line(color = "grey20", linewidth = 0.1),
      axis.ticks.length  = unit(0.2, "lines"),
      legend.position    = "none",
      plot.tag           = element_text(face = "bold", size = 14),
      panel.spacing.x = unit(2,"lines"),
      plot.margin = margin(0,2,1,2)
    ) +
  guides(
    x = guide_axis(cap = TRUE)
  ) 
  #scale_x_continuous(breaks = seq(0, 1, by = .25), expand = c(.1, .1) , limits = c(0, 1.1)) # Customize x-axis breaks and limits

# PARLAMENT ##### ----------------------------------------------------------------------------------------

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
    fig_center = 9.5 # Set a fixed center for figures (used for consistent alignment)
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
  )) |> # Assign the silhouette code based on the taxa group
  mutate(facet = factor(facet,levels=c("Taxonomic turnover","Functional turnover"))) 
  
Figure_1B_draft <- 
  bbgdm_data_parla %>% 
  mutate(direction = gsub("z", "s", str_to_title(direction))) %>% # Format direction values for consistency
  group_by(dataset, facet) %>% # Group by dataset and facet for calculations
  mutate(es_rank = magnitude / sum(magnitude)) %>% # Calculate rank-normalized effect size
  select(dataset, facet, predictor, direction, magnitude, es_rank, data_num, var_num, taxa, system, disturbance, brace_middle, fig_center, image) %>% # Select relevant columns
  mutate(facet = factor(facet,levels=c("Taxonomic turnover","Functional turnover"))) |> 
  ggplot(aes(x = data_num, y = var_num)) + # Define ggplot base aesthetics
  facet_wrap(~facet)+  # Create facets for species and Functional turnover
  coord_radial(start=0.2 * pi, end=-0.2 * pi,expand=TRUE) + # Use polar coordinates for circular visualization
  scale_y_continuous(limits = c(-4, 12), expand = expansion(mult = c(0, -.3)), labels= c('Precipitation','Temperature','Spatial distance','<b>Habitat heterogeneity</b>','<b>Human Land Cover</b>','<b>Human Footprint</b>'), breaks=1:6)+ # Set y-axis scale with expansion  # Add braces for each biotic group
  scale_x_continuous(expand=c(0,.1))+
    theme_void(base_size=14)+
    theme(
      axis.text.y =  element_markdown(size=12,hjust=.99,angle=35),
      axis.text.x =  element_blank(),
      axis.title.x = element_blank(),
      axis.title.y= element_blank(),
      strip.background = element_blank(),
      panel.grid.major.y =   element_line(colour = "gray90", linewidth = .3,linetype="11"),
      strip.text = element_text(colour = "black",size=14, hjust=.5,face="bold",family="sans"),
      ,plot.title.position = "plot"
      ,plot.background = element_blank()
      ,panel.spacing.x = unit(2, "lines")  # Reduced facet spacing (horizontal)
      ,legend.position = 'bottom'
      ,legend.box = 'horizontal'
      ,legend.box.just = 'top'
      ,legend.direction = "vertical"
      ,legend.justification = "center"
      ,legend.key.size = unit(.5, 'cm')
      ,legend.text = element_text(size=14)
      ,legend.title = element_text(size=14, face="bold")
      ,plot.margin = unit(c(0,2,0,2),"lines")
      ,panel.spacing = unit(5, "cm")
      ,plot.tag=element_text(face="bold", size=14)
    )   +
    
  labs(size = "Relative contribution") +
  stat_brace(aes(group = taxa), outerstart = 8, bending = 0.9, width = 0.5, linewidth = 0.3) +
  
  # Add blue and green tiles for ecosystem type
  geom_tile(aes(x = data_num, y = var_num, fill = system), alpha = 0.3, colour = "gray70") +
  scale_fill_manual(values = ecosystem_type_colors, labels = c("Freshwater", "Terrestrial"), name = "Ecosystem type") +
  
  # Add tiles for disturbance types
  ggnewscale::new_scale_fill() +
  geom_tile(aes(y = 7, x = data_num, fill = disturbance), colour = NA) +
  scale_fill_manual(name = "Land use type", 
                    labels = c("Agriculture", "Forest", "Urban", "Multiple"),
                    values = c("#a36627", "#848c04", "#1c1c0c", "#dc7c5c"),
                    guide = guide_legend(direction = "horizontal", title.position = "top", nrow = 2, byrow = TRUE)) +
  
  # Add dots for effect sizes of predictors on species/Functional turnover
  ggnewscale::new_scale_fill() +
  geom_point(aes(x = data_num, y = var_num, size = es_rank), stroke = 0.1, shape = 21, fill="black") + # Points for effect sizes
  scale_size_continuous("Contribution to turnover", range = c(0, 3), breaks = c(0, 0.25, 0.5, 1), 
                        guide = guide_legend(override.aes = list(fill = "black"), direction = "horizontal", title.position = "top", nrow = 2, byrow = TRUE))+
  
  # Add silhouettes for taxa
   geom_phylopic(data = . %>% distinct(brace_middle, .keep_all = TRUE), aes(x = brace_middle, y = fig_center, image = image),
                 size = 0.06, inherit.aes = FALSE)
  #
  # # Define nudging adjustments for text labels
  nudge_x <- c(2, 1, 1.5, 2, 2.3, 4, 3, -1, -3, -3) # Horizontal adjustments for text positions
  nudge_y <- c(1, 1, 0.5, 0, 0, 0, 0, 1, 1, -0.2)   # Vertical adjustments for text positions
  #
  # # Extract and adjust data for text annotations
  Figure_1B_text_data <- ggplot2::layer_data(Figure_1B_draft, 5) |> # Extract data for layer 12 of the draft figure
    left_join(bbgdm_data_parla |> ungroup() |> distinct(taxa, image)) |> # Join with taxa and image information
    mutate(
      label = ifelse(PANEL == 1, as.character(taxa), ""), # Add labels for Taxonomic turnover facet only
      facet = ifelse(PANEL == 1, "Taxonomic turnover", "Functional turnover") # Assign facets based on PANEL
    ) |>
    distinct(PANEL, image, .keep_all = TRUE) |> # Retain unique PANEL and image combinations
    group_by(PANEL) |> # Group by PANEL to apply nudging adjustments
    mutate(
      x = x + nudge_x, # Apply horizontal nudging
      y = y + nudge_y  # Apply vertical nudging
    ) |>
    mutate(facet = factor(facet,levels=c("Taxonomic turnover","Functional turnover")))

  # Add text annotations to Figure_1B
  Figure_1B <- Figure_1B_draft +
    geom_text(
      data = Figure_1B_text_data  , # Use prepared text data
      aes(y = y, x = x, label = label), # Map adjusted positions and labels
      hjust = 0, # Align text to the left
      size = convert_size(10), # Adjust text size
      family = "sans" # Use sans-serif font
    )

  # Create your subplots as before
  Figure_1 <- plot_grid(
    Figure_1A + theme(plot.margin = margin(t = 10, b = 10)),
    Figure_1B_draft + theme(plot.margin = margin(t = 20, b = 10)),
    ncol = 1,
    rel_heights = c(0.3, 0.7),
    labels = c("a", "b"),
    label_size = 14,
    label_x = c(0, 0),
    label_y = c(1, 1)
  )

  # Define the save folder based on beta_type within the S7 folder
  if(beta_type == "Podani_abun"){
    save_folder <- here::here("S7_Model_outputs_figures_and_tables", "main_figures")
  } else {
    save_folder <- here::here("S7_Model_outputs_figures_and_tables", "supplementary_material", beta_type)
    if(!dir.exists(save_folder)){
      dir.create(save_folder, recursive = TRUE)
    }
  }
  
  # Save the final figure as a PDF
  ggsave(
    filename = file.path(save_folder, "Figure_1.pdf"),  # Save path with file name
    plot = Figure_1,         # Specify the figure to save
    device = cairo_pdf,      # Use Cairo PDF device for high-quality output
    height = 190*1.3, width = 247*1.3, units = "mm" # Define dimensions in millimeters
  )

    
  
bbgdm_data_parsed |> 
   group_by(facet,predictor,taxa) |> 
   summarise(mean = Hmisc::smean.cl.boot(effect_size)[1],
             lower = Hmisc::smean.cl.boot(effect_size)[2],
             upper = Hmisc::smean.cl.boot(effect_size)[3]) |>  print(n=120)
# A tibble: 120 × 6
# Groups:   facet, predictor [12]
#     facet               predictor  taxa                        mean    lower   upper
#     <chr>               <fct>      <fct>                      <dbl>    <dbl>   <dbl>
#   1 Taxonomic turnover Prec       "Microorganisms"         0.160    0.0777   0.252 
#   2 Taxonomic turnover Prec       "Fish"                   0.0923   0.0210   0.176 
#   3 Taxonomic turnover Prec       "Amphibians"             0       NA       NA     
#   4 Taxonomic turnover Prec       "Aquatic\ninvertebrates" 0.104    0.0674   0.140 
#   5 Taxonomic turnover Prec       "Fungi"                  0.256   NA       NA     
#   6 Taxonomic turnover Prec       "Plants"                 0.144    0.0738   0.241 
#   7 Taxonomic turnover Prec       "Gastropods"             0.00959 NA       NA     
#   8 Taxonomic turnover Prec       "Arthropods"             0.142    0.109    0.183 
#   9 Taxonomic turnover Prec       "Birds"                  0.120    0.0912   0.153 
#  10 Taxonomic turnover Prec       "Bats"                   0.0665  NA       NA     
#  11 Taxonomic turnover Temp       "Microorganisms"         0.159    0.0737   0.241 
#  12 Taxonomic turnover Temp       "Fish"                   0.244    0.0183   0.470 
#  13 Taxonomic turnover Temp       "Amphibians"             0.0396  NA       NA     
#  14 Taxonomic turnover Temp       "Aquatic\ninvertebrates" 0.147    0.0922   0.208 
#  15 Taxonomic turnover Temp       "Fungi"                  0.264   NA       NA     
#  16 Taxonomic turnover Temp       "Plants"                 0.143    0.0769   0.213 
#  17 Taxonomic turnover Temp       "Gastropods"             0.0347  NA       NA     
#  18 Taxonomic turnover Temp       "Arthropods"             0.128    0.0935   0.163 
#  19 Taxonomic turnover Temp       "Birds"                  0.198    0.160    0.245 
#  20 Taxonomic turnover Temp       "Bats"                   0.275   NA       NA     
#  21 Taxonomic turnover Geographic "Microorganisms"         0.102    0.0286   0.174 
#  22 Taxonomic turnover Geographic "Fish"                   0.114    0.0370   0.172 
#  23 Taxonomic turnover Geographic "Amphibians"             0.514   NA       NA     
#  24 Taxonomic turnover Geographic "Aquatic\ninvertebrates" 0.0634   0.0320   0.0997
#  25 Taxonomic turnover Geographic "Fungi"                  0.330   NA       NA     
#  26 Taxonomic turnover Geographic "Plants"                 0.0971   0.0549   0.141 
#  27 Taxonomic turnover Geographic "Gastropods"             0.109   NA       NA     
#  28 Taxonomic turnover Geographic "Arthropods"             0.101    0.0731   0.130 
#  29 Taxonomic turnover Geographic "Birds"                  0.104    0.0694   0.141 
#  30 Taxonomic turnover Geographic "Bats"                   0.493   NA       NA     
#  31 Taxonomic turnover het        "Microorganisms"         0.165    0.0988   0.246 
#  32 Taxonomic turnover het        "Fish"                   0.145    0.0600   0.202 
#  33 Taxonomic turnover het        "Amphibians"             0.0159  NA       NA     
#  34 Taxonomic turnover het        "Aquatic\ninvertebrates" 0.182    0.113    0.257 
#  35 Taxonomic turnover het        "Fungi"                  0.0944  NA       NA     
#  36 Taxonomic turnover het        "Plants"                 0.217    0.130    0.326 
#  37 Taxonomic turnover het        "Gastropods"             0.0600  NA       NA     
#  38 Taxonomic turnover het        "Arthropods"             0.157    0.122    0.198 
#  39 Taxonomic turnover het        "Birds"                  0.177    0.141    0.215 
#  40 Taxonomic turnover het        "Bats"                   0.0860  NA       NA     
#  41 Taxonomic turnover modis      "Microorganisms"         0.130    0.0540   0.216 
#  42 Taxonomic turnover modis      "Fish"                   0.0640   0.0202   0.106 
#  43 Taxonomic turnover modis      "Amphibians"             0       NA       NA     
#  44 Taxonomic turnover modis      "Aquatic\ninvertebrates" 0.0934   0.0448   0.151 
#  45 Taxonomic turnover modis      "Fungi"                  0.0555  NA       NA     
#  46 Taxonomic turnover modis      "Plants"                 0.0742   0.0407   0.109 
#  47 Taxonomic turnover modis      "Gastropods"             0.00940 NA       NA     
#  48 Taxonomic turnover modis      "Arthropods"             0.117    0.0827   0.155 
#  49 Taxonomic turnover modis      "Birds"                  0.0640   0.0400   0.0899
#  50 Taxonomic turnover modis      "Bats"                   0       NA       NA     
#  51 Taxonomic turnover hfp        "Microorganisms"         0.139    0.0722   0.212 
#  52 Taxonomic turnover hfp        "Fish"                   0.242    0.0740   0.410 
#  53 Taxonomic turnover hfp        "Amphibians"             0.431   NA       NA     
#  54 Taxonomic turnover hfp        "Aquatic\ninvertebrates" 0.236    0.149    0.330 
#  55 Taxonomic turnover hfp        "Fungi"                  0       NA       NA     
#  56 Taxonomic turnover hfp        "Plants"                 0.133    0.0811   0.182 
#  57 Taxonomic turnover hfp        "Gastropods"             0.216   NA       NA     
#  58 Taxonomic turnover hfp        "Arthropods"             0.172    0.143    0.204 
#  59 Taxonomic turnover hfp        "Birds"                  0.197    0.155    0.241 
#  60 Taxonomic turnover hfp        "Bats"                   0.0801  NA       NA     
#  61 Functional turnover   Prec       "Microorganisms"         0.178    0.0413   0.345 
#  62 Functional turnover   Prec       "Fish"                   0.264    0.0790   0.400 
#  63 Functional turnover   Prec       "Amphibians"             0       NA       NA     
#  64 Functional turnover   Prec       "Aquatic\ninvertebrates" 0.0900   0.0292   0.185 
#  65 Functional turnover   Prec       "Fungi"                  0.283   NA       NA     
#  66 Functional turnover   Prec       "Plants"                 0.163    0.0599   0.289 
#  67 Functional turnover   Prec       "Gastropods"             0.133   NA       NA     
#  68 Functional turnover   Prec       "Arthropods"             0.132    0.0766   0.189 
#  69 Functional turnover   Prec       "Birds"                  0.145    0.0970   0.193 
#  70 Functional turnover   Prec       "Bats"                   0.523   NA       NA     
#  71 Functional turnover   Temp       "Microorganisms"         0.0681   0.0396   0.0987
#  72 Functional turnover   Temp       "Fish"                   0.0911   0.0365   0.146 
#  73 Functional turnover   Temp       "Amphibians"             0.276   NA       NA     
#  74 Functional turnover   Temp       "Aquatic\ninvertebrates" 0.0530   0.0145   0.101 
#  75 Functional turnover   Temp       "Fungi"                  0.104   NA       NA     
#  76 Functional turnover   Temp       "Plants"                 0.114    0.0338   0.216 
#  77 Functional turnover   Temp       "Gastropods"             0.0184  NA       NA     
#  78 Functional turnover   Temp       "Arthropods"             0.0736   0.0406   0.113 
#  79 Functional turnover   Temp       "Birds"                  0.127    0.0871   0.170 
#  80 Functional turnover   Temp       "Bats"                   0       NA       NA     
#  81 Functional turnover   Geographic "Microorganisms"         0.0616   0.0133   0.124 
#  82 Functional turnover   Geographic "Fish"                   0.00590  0        0.0177
#  83 Functional turnover   Geographic "Amphibians"             0       NA       NA     
#  84 Functional turnover   Geographic "Aquatic\ninvertebrates" 0.0626   0.00942  0.135 
#  85 Functional turnover   Geographic "Fungi"                  0       NA       NA     
#  86 Functional turnover   Geographic "Plants"                 0.0199   0.00633  0.0398
#  87 Functional turnover   Geographic "Gastropods"             0.135   NA       NA     
#  88 Functional turnover   Geographic "Arthropods"             0.0693   0.0294   0.122 
#  89 Functional turnover   Geographic "Birds"                  0.0860   0.0526   0.125 
#  90 Functional turnover   Geographic "Bats"                   0       NA       NA     
#  91 Functional turnover   het        "Microorganisms"         0.298    0.123    0.522 
#  92 Functional turnover   het        "Fish"                   0.315    0.233    0.470 
#  93 Functional turnover   het        "Amphibians"             0       NA       NA     
#  94 Functional turnover   het        "Aquatic\ninvertebrates" 0.227    0.150    0.313 
#  95 Functional turnover   het        "Fungi"                  0.293   NA       NA     
#  96 Functional turnover   het        "Plants"                 0.136    0.0591   0.233 
#  97 Functional turnover   het        "Gastropods"             0.0182  NA       NA     
#  98 Functional turnover   het        "Arthropods"             0.124    0.0778   0.185 
#  99 Functional turnover   het        "Birds"                  0.210    0.158    0.264 
# 100 Functional turnover   het        "Bats"                   0       NA       NA     
# 101 Functional turnover   modis      "Microorganisms"         0.0442   0.0134   0.0794
# 102 Functional turnover   modis      "Fish"                   0.0238   0        0.0648
# 103 Functional turnover   modis      "Amphibians"             0.182   NA       NA     
# 104 Functional turnover   modis      "Aquatic\ninvertebrates" 0.0820   0.0404   0.125 
# 105 Functional turnover   modis      "Fungi"                  0.0234  NA       NA     
# 106 Functional turnover   modis      "Plants"                 0.0675   0.0338   0.114 
# 107 Functional turnover   modis      "Gastropods"             0.0873  NA       NA     
# 108 Functional turnover   modis      "Arthropods"             0.0630   0.0314   0.104 
# 109 Functional turnover   modis      "Birds"                  0.0549   0.0291   0.0850
# 110 Functional turnover   modis      "Bats"                   0.477   NA       NA     
# 111 Functional turnover   hfp        "Microorganisms"         0.120    0.0542   0.188 
# 112 Functional turnover   hfp        "Fish"                   0.148    0        0.301 
# 113 Functional turnover   hfp        "Amphibians"             0.543   NA       NA     
# 114 Functional turnover   hfp        "Aquatic\ninvertebrates" 0.190    0.105    0.283 
# 115 Functional turnover   hfp        "Fungi"                  0.296   NA       NA     
# 116 Functional turnover   hfp        "Plants"                 0.0807   0.0416   0.126 
# 117 Functional turnover   hfp        "Gastropods"             0.169   NA       NA     
# 118 Functional turnover   hfp        "Arthropods"             0.221    0.161    0.300 
# 119 Functional turnover   hfp        "Birds"                  0.164    0.119    0.216 
# 120 Functional turnover   hfp        "Bats"                   0       NA       NA     







#Effect of spatial extent on magnitude of effect 


metadata <- readxl::read_excel("S6_Synthesis_model/data/synthesis_data_complete.xlsx") |> 
 select(dataset,spatial.extent) |> 
 distinct(dataset, .keep_all=TRUE)

eco_breaks <- c(0, 1, 100, 1e4, 1e6, Inf)  # Example: patch (<1 km²), landscape (1-100 km²), regional (100-10k km²), continental (>10k km²)
eco_labels <- c("Micro (<1 km²)", "Local (1-100 km²)", "Landscape (100-10k km²)", "Regional (10k-1M km²)", "Continental (>1M km²)")
color_pred <- c(
  "hfp" = "#C2185B",  # Bold magenta-red
  "modis" = "#E53935", 
  "het" = "#C20291",  # Darker green
  "Geographic" = "#8BC34A",
  "Prec" = "#1976D2", 
  "Temp" = "#90CAF9"  # Very light blue
)

predictor_rank_plot <- bbgdm_summary |> 
  select(dataset, facet, predictor, effect_size) |> 
  left_join(metadata, by = "dataset") |>  # Added explicit join column
  mutate(
    spatial_bin = cut(
      spatial.extent, 
      breaks = eco_breaks, 
      labels = eco_labels, 
      include.lowest = TRUE
    )
  ) |> 
  group_by(dataset, facet) |> 
  mutate(pred_rank = rank(-effect_size, ties.method = "random")) |> 
  ungroup() |>  # Important: ungroup before plotting
  # Plot
  ggplot(aes(x = spatial_bin, y = pred_rank, colour = predictor, group = predictor)) +
  geom_jitter(alpha = 0.1, width = 0.1, height = 0) +  # Added height=0 to prevent vertical jitter
  stat_summary(
    fun = mean,  # Changed from fun.data for line geometry
    geom = "line",
    linewidth = 1
  ) + 
  stat_summary(
    fun.data = "mean_cl_boot",
    geom = "pointrange", 
    size = 0.5,
    fun.args = list(B = 5000)  # Added bootstrap replicates
  ) +
  facet_wrap(~facet) +
  theme_bw() +
  scale_y_reverse(breaks = 1:6) +  # Explicit breaks for ranks
  labs(
    x = "Spatial Extent", 
    y = "Predictor Rank (1 = Strongest Effect)", 
    colour = "Predictor"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  ) +
  scale_colour_manual(
    values=color_pred,
    labels = c("Precipitation", "Temperature", "Spatial distance",
               "Habitat heterogeneity", "Human land cover", "Human footprint")
  ) +
  guides(color = guide_legend(nrow = 2))  # Wrap legend items




# Save the final figure as a PDF
ggsave(
  filename = file.path("S7_Model_outputs_figures_and_tables/extended_data", "Extended_Figure_2.pdf"),  # Save path with file name
  plot = predictor_rank_plot,         # Specify the figure to save
  device = cairo_pdf,      # Use Cairo PDF device for high-quality output
  height = 5, width = 6, units = "in" # Define dimensions in millimeters
)



















rank_summary <- bbgdm_summary |> 
  select(dataset, facet, predictor, effect_size) |> 
  left_join(metadata, by = "dataset") |> 
  mutate(
    spatial_bin = cut(
      spatial.extent, 
      breaks = eco_breaks, 
      labels = eco_labels, 
      include.lowest = TRUE
    )
  ) |> 
  group_by(dataset, facet) |> 
  mutate(pred_rank = rank(-effect_size, ties.method = "random")) |> 
  ungroup() |> 
  # Calculate mean rank and bootstrap CIs by predictor, facet, and spatial bin
  group_by(facet, predictor, spatial_bin) |> 
  summarise(
    mean_rank = mean(pred_rank),
    lower_ci = smean.cl.boot(pred_rank, B = 5000)[["Lower"]],
    upper_ci = smean.cl.boot(pred_rank, B = 5000)[["Upper"]],
    n_datasets = n(),
    .groups = "drop"
  ) |> 
  arrange(facet, spatial_bin, mean_rank)  # Sort for readability

rank_summary |>  print(n=60)

# A tibble: 60 × 7
#    facet               predictor  spatial_bin             mean_rank lower_ci upper_ci n_datasets
#    <chr>               <fct>      <fct>                       <dbl>    <dbl>    <dbl>      <int>
#  1 Taxonomic turnover het        Micro (<1 km²)               2.5      1.67     3.67          6
#  2 Taxonomic turnover Geographic Micro (<1 km²)               3.17     1.83     4.5           6
#  3 Taxonomic turnover Prec       Micro (<1 km²)               3.33     2.5      4.17          6
#  4 Taxonomic turnover hfp        Micro (<1 km²)               3.5      1.83     5.17          6
#  5 Taxonomic turnover Temp       Micro (<1 km²)               3.83     2.33     5.33          6
#  6 Taxonomic turnover modis      Micro (<1 km²)               4.67     3.83     5.33          6
#  7 Taxonomic turnover hfp        Local (1-100 km²)            2.63     2.09     3.23         35
#  8 Taxonomic turnover Temp       Local (1-100 km²)            3.06     2.51     3.63         35
#  9 Taxonomic turnover het        Local (1-100 km²)            3.17     2.71     3.63         35
# 10 Taxonomic turnover Prec       Local (1-100 km²)            3.63     3.09     4.14         35
# 11 Taxonomic turnover Geographic Local (1-100 km²)            4.2      3.66     4.71         35
# 12 Taxonomic turnover modis      Local (1-100 km²)            4.31     3.77     4.83         35
# 13 Taxonomic turnover het        Landscape (100-10k km²)      2.85     2.47     3.24         66
# 14 Taxonomic turnover hfp        Landscape (100-10k km²)      3.08     2.65     3.48         66
# 15 Taxonomic turnover Temp       Landscape (100-10k km²)      3.26     2.86     3.65         66
# 16 Taxonomic turnover Prec       Landscape (100-10k km²)      3.70     3.29     4.11         66
# 17 Taxonomic turnover Geographic Landscape (100-10k km²)      3.91     3.52     4.29         66
# 18 Taxonomic turnover modis      Landscape (100-10k km²)      4.21     3.82     4.62         66
# 19 Taxonomic turnover hfp        Regional (10k-1M km²)        3.04     2.6      3.5          50
# 20 Taxonomic turnover Temp       Regional (10k-1M km²)        3.18     2.7      3.68         50
# 21 Taxonomic turnover het        Regional (10k-1M km²)        3.38     2.98     3.78         50
# 22 Taxonomic turnover Prec       Regional (10k-1M km²)        3.46     2.96     3.96         50
# 23 Taxonomic turnover Geographic Regional (10k-1M km²)        3.82     3.36     4.3          50
# 24 Taxonomic turnover modis      Regional (10k-1M km²)        4.12     3.66     4.54         50
# 25 Taxonomic turnover Temp       Continental (>1M km²)        2.4      1.6      3.2           5
# 26 Taxonomic turnover Prec       Continental (>1M km²)        3        1.8      4.2           5
# 27 Taxonomic turnover Geographic Continental (>1M km²)        3        1.2      4.8           5
# 28 Taxonomic turnover hfp        Continental (>1M km²)        4        2.4      5.6           5
# 29 Taxonomic turnover het        Continental (>1M km²)        4.2      3.6      4.8           5
# 30 Taxonomic turnover modis      Continental (>1M km²)        4.4      3        5.8           5
# 31 Functional turnover   het        Micro (<1 km²)               1.83     1        3.17          6
# 32 Functional turnover   Temp       Micro (<1 km²)               2.83     1.83     4             6
# 33 Functional turnover   Prec       Micro (<1 km²)               3.33     2.5      4.17          6
# 34 Functional turnover   Geographic Micro (<1 km²)               3.67     2.83     4.67          6
# 35 Functional turnover   hfp        Micro (<1 km²)               4.5      2.83     5.83          6
# 36 Functional turnover   modis      Micro (<1 km²)               4.83     4        5.67          6
# 37 Functional turnover   hfp        Local (1-100 km²)            2.91     2.43     3.4          35
# 38 Functional turnover   het        Local (1-100 km²)            3.2      2.63     3.8          35
# 39 Functional turnover   Prec       Local (1-100 km²)            3.43     2.86     4.03         35
# 40 Functional turnover   Geographic Local (1-100 km²)            3.54     3        4.11         35
# 41 Functional turnover   Temp       Local (1-100 km²)            3.71     3.11     4.31         35
# 42 Functional turnover   modis      Local (1-100 km²)            4.2      3.74     4.63         35
# 43 Functional turnover   hfp        Landscape (100-10k km²)      2.79     2.45     3.14         66
# 44 Functional turnover   het        Landscape (100-10k km²)      3.06     2.62     3.5          66
# 45 Functional turnover   Prec       Landscape (100-10k km²)      3.32     2.92     3.71         66
# 46 Functional turnover   Temp       Landscape (100-10k km²)      3.36     2.98     3.76         66
# 47 Functional turnover   Geographic Landscape (100-10k km²)      4.11     3.71     4.5          66
# 48 Functional turnover   modis      Landscape (100-10k km²)      4.36     3.97     4.74         66
# 49 Functional turnover   het        Regional (10k-1M km²)        2.52     2.14     2.92         50
# 50 Functional turnover   hfp        Regional (10k-1M km²)        2.76     2.32     3.2          50
# 51 Functional turnover   modis      Regional (10k-1M km²)        3.48     3.06     3.9          50
# 52 Functional turnover   Prec       Regional (10k-1M km²)        3.66     3.2      4.1          50
# 53 Functional turnover   Temp       Regional (10k-1M km²)        3.88     3.42     4.36         50
# 54 Functional turnover   Geographic Regional (10k-1M km²)        4.7      4.28     5.08         50
# 55 Functional turnover   Prec       Continental (>1M km²)        1.6      1        2.4           5
# 56 Functional turnover   hfp        Continental (>1M km²)        2.8      1.8      4             5
# 57 Functional turnover   het        Continental (>1M km²)        3        2.2      3.8           5
# 58 Functional turnover   modis      Continental (>1M km²)        3.4      2        4.8           5
# 59 Functional turnover   Temp       Continental (>1M km²)        5        4.2      5.8           5
# 60 Functional turnover   Geographic Continental (>1M km²)        5.2      4        6             5