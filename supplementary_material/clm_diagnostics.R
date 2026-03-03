#'#######################################
#'Climate, Spatial distance and Spatial extent 
#'#######################################
#'
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

metadata <- readxl::read_excel("S6_Synthesis_model/data/synthesis_data_revision.xlsx") |> 
  mutate(metric_type = ifelse(metric_type == "pa","Presence/Absence","Abundance")) |> 
  mutate(facet = ifelse(facet == "Taxonomic","Species replacement","Trait replacement")) |> 
  rename(framework = "beta_type")

source("S4_run_BBGDM/functions/helper_functions_S4.R") #get silly functions to help

beta_type = "Podani_abun"
# Combine data for parlament plot  --------------------------------------------------------------------------------

files <- tools::file_path_sans_ext(list.files(glue::glue("S4_run_BBGDM/bbgdm_output/{beta_type}")))
data_summary <- readxl::read_xlsx("S1_Preprocessing/Miscellaneous/dataset_info_all.xlsx") 
bbgdms_output <- purrr::map(files, ~ read_rds(glue::glue("S4_run_BBGDM/bbgdm_output/{beta_type}/{.x}.rds"))) |> 
  set_names(gsub("_bbgdm", "", files)) |> 
  discard_at(c("N67TTP", "N78TTP"))

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


bbgdm_data_parsed <- bbgdm_summary |> 
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


ranked_predictors <- bbgdm_data_parsed |> 
  left_join(metadata |> distinct(dataset,spatial.extent)) |> 
  select(dataset,facet,predictor,magnitude,spatial.extent) |> 
  group_by(dataset,facet) |> 
  arrange(dataset,facet,predictor) |> 
  mutate(ranked_effect = rank(magnitude, ties.method="random")) 

selected_preds <- ranked_predictors %>%
  mutate(extent_bin = cut(spatial.extent,
                          breaks = c(0, 100, 1000, 10000, 100000, Inf),
                          labels = c("0–100 km²", 
                                     "100–1,000 km²", 
                                     "1,000–10,000 km²", 
                                     "10,000–100,000 km²",
                                     "&lt;100,000 km²"),
                          include.lowest = TRUE))


# Replace predictor codes with full names
selected_preds <- selected_preds %>%
  mutate(predictor = recode(predictor,
                            "Prec" = "Precipitation",
                            "Temp" = "Temperature",
                            "Geographic" = "Geographic distance",
                            "modis" = "Land-cover (MODIS)",
                            "het" = "Habitat heterogeneity",
                            "hfp" = "Human Footprint")) |> 
  mutate(predictor = factor(predictor, levels = c(
    "Human Footprint",
    "Habitat heterogeneity",
    "Land-cover (MODIS)",
    "Geographic distance",
    "Temperature",
    "Precipitation"
  )))


# Make sure ranked_effect is an ordered factor
selected_preds <- selected_preds |> 
  mutate(ranked_label = factor(ranked_effect,
                               levels = 1:6,
                               labels = c(
                                 "Least important", 
                                 "5th most important", 
                                 "4th most important", 
                                 "3rd most important", 
                                 "2nd most important", 
                                 "Most important"
                               )
  )) |>  drop_na()


#Models ----

# Load necessary libraries
library(ordinal)
library(performance)
library(marginaleffects)
library(dplyr)
library(ggplot2)

# Filter data for the species replacement facet
df_species <- filter(selected_preds, facet == "Species replacement") |> 
  mutate(spatial.extent = log10(spatial.extent)) |> 
  drop_na()

df_traits <- filter(selected_preds, facet == "Trait replacement") |> 
  mutate(spatial.extent = log10(spatial.extent)) |> 
  drop_na()

# Fit the original model
mod <- clm(ranked_label ~ spatial.extent * predictor, data = df)

# 1. Test the proportional odds assumption
mod_nominal <- clm(ranked_label ~ spatial.extent * predictor,
                   nominal = ~ predictor,
                   data = df)
cat("Proportional odds assumption test:\n")
print(anova(mod, mod_nominal))

# 2. Assess model fit
cat("\nPseudo R²:\n")
print(r2_nagelkerke(mod))

cat("\nModel Summary:\n")
print(summary(mod))

# 3. Check predicted vs. observed values
pred <- predict(mod, type = "class")
tab <- table(Predicted = pred$fit, Observed = df$ranked_label)
cat("\nPrediction vs. Observed:\n")
print(tab)

# 5. Check collinearity with VIF (using lm as proxy)
vif_model <- lm(as.numeric(ranked_label) ~ spatial.extent * predictor, data = df)
cat("\nVIF for predictors:\n")
print(car::vif(vif_model,type = 'predictor'))


library(gtsummary)

gtsummary::tbl_merge(
list(tbl_regression(
  clm(ranked_label ~ spatial.extent * predictor, data = df_species),
  intercept = FALSE,
  estimate_fun = purrr::partial(style_ratio, digits = 3),
  pvalue_fun = purrr::partial(style_sigfig, digits = 3)),
tbl_regression(clm(ranked_label ~ spatial.extent * predictor, data = df_traits),
intercept = FALSE,
estimate_fun = purrr::partial(style_ratio, digits = 3),
pvalue_fun = purrr::partial(style_sigfig, digits = 3))),
tab_spanner=c("**Species replacement**","**Trait replacement**"))


selected_preds |> 
  group_by(predictor) |> 
  slice_max(spatial.extent) |> 
  arrange(facet)
# Filter data for the species replacement facet
df <- filter(selected_preds, facet == "Trait replacement") |> 
  mutate(spatial.extent = log10(spatial.extent)) |> 
  drop_na()

# Fit the original model
mod <- clm(ranked_label ~ spatial.extent * predictor, data = df)

# 1. Test the proportional odds assumption
mod_nominal <- clm(ranked_label ~ spatial.extent * predictor,
                   nominal = ~ predictor,
                   data = df)
cat("Proportional odds assumption test:\n")
print(anova(mod, mod_nominal))

# 2. Assess model fit
cat("\nPseudo R²:\n")
print(r2_nagelkerke(mod))

cat("\nModel Summary:\n")
print(summary(mod))

# 3. Check predicted vs. observed values
pred <- predict(mod, type = "class")
tab <- table(Predicted = pred$fit, Observed = df$ranked_label)
cat("\nPrediction vs. Observed:\n")
print(tab)

# 5. Check collinearity with VIF (using lm as proxy)
vif_model <- lm(as.numeric(ranked_label) ~ spatial.extent * predictor, data = df)
cat("\nVIF for predictors:\n")
print(car::vif(vif_model,type = 'predictor'))

