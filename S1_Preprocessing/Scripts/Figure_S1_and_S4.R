###############################################################################
# SCRIPT NAME: Figure_S1_and_S4.R
#
# DESCRIPTION:
#   Generates map and dataset-feature figures used in supplementary outputs
#   (Figure S1 and Figure S4).
#
# USAGE:
#   Run after Step 1 preprocessing and synthesis data compilation are available.
#
# INPUTS:
#   - `S1_Preprocessing/Processed/*.rds`
#   - `S1_Preprocessing/Miscellaneous/dataset_info_all.xlsx`
#   - `S6_Synthesis_model/data/synthesis_data_new.xlsx`
#
# OUTPUTS:
#   - Publication-ready figure objects saved/exported within this script.
#
# AUTHOR: Caio Graco-Roza
# LAST UPDATED: 2026-03-03
###############################################################################
pacman::p_load(ggmap,
               tidyverse,
               ggthemes,
               here,
               showtext,
               readxl,
               ggplot2,
               ggtext,
               patchwork,
               colorspace,
               showtext,
               ggh4x)
#registering my google token for using ggmap ----

google_key <- Sys.getenv("GGMAP_GOOGLE_API_KEY", unset = "")
if (nzchar(google_key)) {
  ggmap::register_google(key = google_key)
}

# install fonts  ----

showtext::showtext_auto(TRUE)

#Load dataset_features
dataset_features <- "S1_Preprocessing/Processed" |>
  list.files() |>
  purrr::map(~ here("S1_Preprocessing/Processed", .x) |>
               readr::read_rds() |>
               pluck("predictors")) |>
  set_names(tools::file_path_sans_ext(list.files("S1_Preprocessing/Processed"))) |>
  bind_rows(.id = "dataset") |>
  mutate(dataset = str_remove(dataset, "_clean")) |>
  left_join(
    "S1_Preprocessing/Miscellaneous/dataset_info_all.xlsx" |>
      readxl::read_xlsx(),
    by = c("dataset" = "dataset_name")
  ) |> 
  drop_na(disturbance)

# Set theme parameters
my_base_theme <- function() {
  ggplot2::theme_void(base_family = "sans", base_size = 5) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      strip.background = ggplot2::element_rect(fill = NA, color = NA),
      strip.text = ggplot2::element_text(face = "bold", margin=margin(b=2)),
      axis.title.x = ggtext::element_markdown(family = "sans", face = "bold", hjust = 0.5, margin=margin(t=5)),
      plot.title = ggtext::element_markdown(family = "sans", size=5, face = "bold", hjust = 0, margin=margin(b=2)),
      axis.text.x = element_markdown(family = "sans", color = "grey30", margin = margin(t = 2)),
      axis.text.y = element_blank(),
      panel.spacing.x = unit(0.2, "lines"),
      panel.spacing.y = unit(0.2, "lines"),
      axis.line.x = element_line(color = "grey70", linewidth=0.1),
      axis.ticks.x = element_line(color = "grey70", linewidth=0.1),
      axis.ticks.length.x = unit(0.1, "lines"),
      plot.margin = margin(5, 5, 5, 5),
      legend.position="none",
      plot.tag = element_text(face="bold")
    )
}

additional_features <- 
  "S6_Synthesis_model/data/synthesis_data_new.xlsx" |>
  readxl::read_xlsx() |>
  filter(!dataset %in% c("N67TTP","N78TTP","S47TTP")) |> 
  dplyr::mutate(
    ,facet = base::factor(facet, levels=c("Taxonomic","Functional"))
    ,realm = gsub("aquatic","freshwater",realm),
    ,realm = base::factor(realm, levels = c("terrestrial","freshwater"))
    ,biotic.group = base::factor(biotic.group, levels = c("invertebrate", "vertebrate", "plant", "microorganism"))
    ,disturbance = base::factor(disturbance,levels = c("multiple","agriculture", "forest", "urban"))
    ,species.number =base::log10(species.number+1)
    ,spatial.min =base::log10(spatial.min+1)
    ,hfp.range=hfp.max-hfp.min
    ,latitude.mean=abs(latitude.mean)
    ,direction = gsub("z","s",direction),
    ,direction= base::factor(direction, levels=c("differentiation", "homogenisation"))
  )

# Figure 1 (Map)
blue_green <- c("#118ab2", "#06d6a0")

world_map <- map_data("world") |>
  filter(!long > 180)

map <- world_map |>
  distinct(region) |>
  ggplot(aes(map_id = region)) +
  geom_map(map = world_map, fill = "gray40") +
  expand_limits(x = world_map$long, y = world_map$lat) +
  coord_map("moll") +
  theme_void()+
  theme(axis.text.x  = element_blank(),
        axis.text.y  = element_blank(),
        axis.title.x  = element_blank(),
        axis.title.y  = element_blank(),
        plot.title = element_text(size=5, face="bold"),
        plot.margin = margin(5, 5, 5, 5))

map_aquatic <- map + 
  geom_point(
    data = dataset_features %>% filter(system %in% "aquatic"),
    aes(x = x, y = y),
    fill = blue_green[1],
    inherit.aes = FALSE,
    shape = 21,
    size = 1,
    alpha = 1
  ) + 
  labs(title = "Freshwater sites")

map_terrestrial <- map + 
  geom_point(
    data = dataset_features %>% filter(system %in% "terrestrial"),
    aes(x = x, y = y),
    fill = blue_green[2],
    inherit.aes = FALSE,
    shape = 21,
    size = 1,
    alpha = 1
  ) + 
  labs(title = "Terrestrial sites") 

spatial_extent <- additional_features |> 
  ggplot(aes(x = log10(spatial.extent), fill = realm)) +
  geom_density(position = "identity", alpha = .5, colour ="white", linewidth=0.2) +
  scale_fill_manual(values = blue_green) +
  labs(x = " Spatial extent (Log<sub>10</sub> km)", y = "") +
  my_base_theme() +
 # guide_axis(cap = TRUE)+
  scale_x_continuous(breaks=seq(-1,7,1)) 

species_richness <- additional_features |> 
  ggplot(aes(x = 10^species.number, fill = realm)) +
  geom_density(position = "identity", alpha = .5, colour ="white", linewidth=0.2) +
  scale_fill_manual(values = blue_green) +
  labs(x = "Species number", y = "") +
  my_base_theme() +
  scale_x_continuous(breaks=seq(0,1150,230)) +
  guides(x = guide_axis_truncated(trunc_lower = 0, trunc_upper = 1150)) 

# Latitude

latitude <- dataset_features |>
  mutate(system = factor(system, levels=c("aquatic", "terrestrial"))) |>
  group_by(dataset, system) |>
  ggplot(aes(x = abs(y), fill = system)) +
  geom_density(position = "identity", alpha = .5, colour ="white", linewidth=0.2) +
  scale_fill_manual(values = blue_green) +
  labs(x = "Absolute latitude (degrees)", y = "") +
  my_base_theme() +
  scale_x_continuous(breaks=seq(0,70,10)) +
  guides(x = guide_axis_truncated(trunc_lower = 0, trunc_upper = 70)) 

# Human footprint
HFP <- dataset_features |>
  ggplot(aes(x = hfp_2000, fill = system)) +
  geom_density(position = "identity", alpha = .5, colour ="white", linewidth=0.2) +
  scale_fill_manual(values = blue_green) +
  labs(x = "Human footprint index", y = "") +
  my_base_theme() +
  scale_x_continuous(breaks=seq(0,50,10), limits=c(0,50)) +
  guides(x = guide_axis_truncated(trunc_lower = 0, trunc_upper = 50)) 

# Habitat Heterogeneity
HH <- dataset_features |>
  group_by(dataset, system) |>
  ggplot(aes(x = het_2000, fill = system)) +
  geom_density(position = "identity", alpha = .5, colour ="white", linewidth=0.2) +
  scale_fill_manual(values = blue_green) +
  labs(x = "Habitat heterogeneity", y = "") +
  my_base_theme() +
  scale_x_continuous(breaks=seq(0,2,.5)) +
  guides(x = guide_axis_truncated(trunc_lower = 0, trunc_upper = 2)) 

# Human Land Cover
HLC <- dataset_features |>
  group_by(dataset, system) |>
  ggplot(aes(x = modis_1500, fill = system)) +
  geom_density(position = "identity", alpha = .5, colour ="white", linewidth=0.2) +
  scale_fill_manual(values = blue_green) +
  labs(x = "Human land cover", y = "") +
  my_base_theme() +
  scale_x_continuous(breaks=seq(0,1,.2)) +
  guides(x = guide_axis_truncated(trunc_lower = 0, trunc_upper = 1)) 


# Main taxa
biotic_groups <- dataset_features |>
  arrange(group) |>
  group_by(system) |>
  count(group) |>
  mutate(group = stringr::str_to_title(group)) |>
  mutate(group = factor(group, levels = rev(c("Vertebrate", "Plant", "Invertebrate", "Microorganism")))) |>
  ggplot(aes(y = group, x = n, fill = system)) +
  geom_bar(stat = "identity", position = position_dodge2(preserve = "single"), colour="white",linewidth=.1) +
  scale_fill_manual(values = blue_green) +
  labs(title = "Biotic groups", x = "Number of sites", y="") +
  my_base_theme() +
  theme(axis.title.y = ggtext::element_markdown(family = "sans", face = "bold", hjust = 0.5,  margin=margin(t=2), angle=90),
        axis.text.y = element_markdown(family = "sans", color = "grey30", margin = margin(r = 2), hjust=1),
        axis.line.y = element_line(color = "grey70", linewidth=0.1),
        axis.ticks.y = element_line(color = "grey70", linewidth=0.1),
        axis.ticks.length.y = unit(0.1, "lines"))+
  scale_x_continuous(breaks=seq(0,6000,1000), limits=c(0,6000)) +
  guides(x = guide_axis_truncated(trunc_lower = 0, trunc_upper = 6000),
         y = guide_axis_truncated(trunc_lower = function(x) x-0.2, trunc_upper = function(x) x+0.2)) 

# Realm taxa
realm <- dataset_features |>
  mutate(system = gsub("aquatic", "freshwater", system)) |>
  count(system) |>
  mutate(system = str_to_title(system)) |>
  ggplot(aes(y = system, x = n, fill = system, colour = after_scale(darken(fill, .1)))) +
  geom_bar(stat = "identity", position = position_dodge2(), colour="white",linewidth=.1) +
  scale_fill_manual(values = blue_green) +
  labs(title = "Ecosystem type", x = "Number of sites", y="") +
my_base_theme() +
  theme(axis.title.y = ggtext::element_markdown(family = "sans", face = "bold", hjust = 0.5,  margin=margin(t=2), angle=90),
        axis.text.y = element_markdown(family = "sans", color = "grey30", margin = margin(r = 2), hjust=1),
        axis.line.y = element_line(color = "grey70", linewidth=0.1),
        axis.ticks.y = element_line(color = "grey70", linewidth=0.1),
        axis.ticks.length.y = unit(0.1, "lines"))+
  scale_x_continuous(breaks=seq(0,14000,2000), limits=c(0,14000)) +
  guides(x = guide_axis_truncated(trunc_lower = 0, trunc_upper = 12000),
         y = guide_axis_truncated(trunc_lower = function(x) x-0.2, trunc_upper = function(x) x+0.2))

# Disturbance
disturbance <- dataset_features |>
  group_by(system) |>
  count(disturbance) |>
  mutate(disturbance = stringr::str_to_title(disturbance)) |>
  mutate(disturbance = factor(disturbance, levels = rev(c("Agriculture","Forest","Urban","Multiple")))) |>
  ggplot(aes(y = disturbance, x = n, fill = system)) +
  geom_bar(stat = "identity", position = position_dodge2(preserve = "single"), colour="white",linewidth=.1) +
  scale_fill_manual(values = blue_green) +
  labs(title = "Main land use type", x = "Number of sites", y="") +
  my_base_theme() +
  theme(axis.title.y = ggtext::element_markdown(family = "sans", face = "bold", hjust = 0.5,  margin=margin(t=2), angle=90),
        axis.text.y = element_markdown(family = "sans", color = "grey30", margin = margin(r = 2), hjust=1),
        axis.line.y = element_line(color = "grey70", linewidth=0.1),
        axis.ticks.y = element_line(color = "grey70", linewidth=0.1),
        axis.ticks.length.y = unit(0.1, "lines"))+
  guides(x = guide_axis_truncated(trunc_lower = 0, trunc_upper = 4000),
         y = guide_axis_truncated(trunc_lower = function(x) x-0.2, trunc_upper = function(x) x+0.2))

# Combined figure
  Extended_figure_1 <-
    ({map_terrestrial | map_aquatic} /
    {realm | biotic_groups | disturbance} /
    {(latitude | HFP)/
    (HH | HLC)}/
    {spatial_extent| species_richness})+ 
    plot_layout(heights=c(.3,.2,.3,.1)) +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(face = "bold", size=8, hjust=1))
  
 
library(dplyr)
library(tidyr)
library(stringr)
  
land_use_colors <- c(Agriculture = "#a36627", Forest = "#848c04", Urban = "#1c1c0c", Multiple = "#dc7c5c")
  
  # Step 1: Get unique buffer used per dataset and predictor
predictor_buffers<-read.csv("S1_Preprocessing/Miscellaneous/dataset_selected_buffer.csv") |> 
  select(-var_name)

  # Step 2: Pivot dataset_features to long format for matching
long_features <- dataset_features %>%
    select(dataset,site,disturbance, hfp_1000:het_2000) |> 
    pivot_longer(cols = matches("^(hfp|modis|het)_\\d+$"),
                 names_to = "var_name", values_to = "value") |> 
  separate(var_name, into=c("predictor","buffer")) |> 
  mutate(buffer = as.numeric(buffer))
  
Extended_figure_6<-long_features |> 
  filter(predictor != "modis") |> 
  mutate(predictor = stringr::str_replace_all(predictor, c("hfp" = "Human Footprint", "het" = "Habitat Heterogeneity"))) |> 
  mutate(predictor = factor(predictor, levels=c("Human Footprint","Habitat Heterogeneity"))) |> 
  mutate(disturbance = factor(stringr::str_to_title(disturbance), levels = rev(c("Agriculture","Forest","Urban","Multiple")))) |> 
  ggplot(aes(y=disturbance,x=value, fill = disturbance)) +
  ggridges::geom_density_ridges(colour="white", scale=1.2)+
  facet_grid(~predictor, scales="free") + 
  theme_bw()+
  scale_fill_manual(values=land_use_colors)+
  labs(x="Human pressure value", y = "Main land use type")+theme(legend.position="none")
  
ggsave(here("S7_Model_outputs_figures_and_tables", "extended_data", "Extended_Figure_6.pdf"), Extended_figure_6, device = cairo_pdf, height = 150, width = 183, units="mm")

