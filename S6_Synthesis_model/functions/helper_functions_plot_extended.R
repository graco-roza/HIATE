require(janitor)

my_theme_extended_models <- function() {
  ggplot2::theme_void(base_family = "sans", base_size = 3) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      strip.background = ggplot2::element_rect(fill = NA, color = NA),
      strip.text = ggplot2::element_text(size = 3, face = "bold", margin = margin(b = 2)),
      strip.text.y = ggplot2::element_text(size = 3, face = "bold", margin = margin(r = 2), angle=90),
      axis.title.x = ggtext::element_markdown(family = "sans", face = "bold", hjust = 0.5, size = 3, margin = margin(t = 2)),
      axis.text.x = element_markdown(family = "sans", color = "grey30", size = 2.5, margin = margin(t = 2)),
      axis.text.y = element_markdown(family = "sans", color = "grey30", size = 2.5, hjust=1, margin = margin(t = 2)),
      panel.spacing = unit(0.2, "lines"),
      axis.line.x = element_line(color = "grey70", linewidth = 0.1),
      axis.ticks.x = element_line(color = "grey70", linewidth = 0.1),
      axis.ticks.length.x = unit(0.1, "lines"),
      plot.margin = margin(3, 3, 3, 3),
      legend.position = "bottom",
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      ,legend.title = element_text(size = 2.5, face="bold")
      ,legend.text = element_text(size = 2.5)
      ,legend.key.size = unit(1, 'mm')
      ,legend.box="horizontal"
      ,legend.margin = margin(0,0,0,0)
      ,legend.box.margin = margin(0,0,0,0)
    )
}

# Helper function for annotations, applies only for direction figures
add_annotations <- function(add_annotations) {
  if (isTRUE(add_annotations)) {
    list(
      ggplot2::annotate("text", family = "sans", label = "Differentiation", 
                        size = convert_size(2.5), x = -.25, y = 13.5, hjust = 1, col = "gray30"),
      ggplot2::annotate("text", family = "sans", label = "Homogenisation", 
                        size = convert_size(2.5), x = 0.25, y = 13.5, hjust = 0, col = "gray30"),
      ggplot2::annotate("segment", linewidth = 0.1, x = -2.4, xend = -2.8, y = 13.5, yend = 13.5, 
                        arrow = arrow(type = "closed", length = unit(0.02, "inches")), col = "gray30"),
      ggplot2::annotate("segment", linewidth = 0.1, x = 2.8, xend = 3.2, y = 13.5, yend = 13.5, 
                        arrow = arrow(type = "closed", length = unit(0.02, "inches")), col = "gray30")
    )
  } else {
    return(NULL)  # No specific annotations for magnitude or shape
  }
}

# Helper function for parameter cleaning and renaming
clean_parameters <- function(parameter) {
  parameter <- tolower(parameter)
  parameter <- gsub("\\_", " ", parameter)
  parameter <- gsub("main land use type|buffer", "", parameter)
  parameter <- gsub("human pressure baseline", "Human Pressure [Baseline]", parameter)
  parameter <- gsub("human pressure range", "Human Pressure [Range]", parameter)
  parameter <- gsub("spatial distance min", "Spatial Distance [Minimum]", parameter)
  parameter <- gsub("absolute latitude mean", "Absolute Latitude [Mean]", parameter)
  parameter <- stringr::str_replace_all(parameter, c(
    "agriculture" = "Agriculture – Multiple",
    "forest" = "Forest – Multiple",
    "urban" = "Urban – Multiple",
    "1500" = "1.5km – 1km",
    "2000" = "2km – 1km"
  ))
  parameter <- case_when(
    grepl("^ecosystem typefreshwater", parameter) ~ "Freshwater – Terrestrial",
    grepl("^direction", parameter) ~ "Homogenisation – Differentiation",
    TRUE ~ parameter
  )
  stringr::str_to_title(parameter)
}


# Function to get support and generate text annotations for figures
get_support <- function(data) {
  cis <- c(.95, .89, .80)  # Define confidence intervals
  foo <- map_dfr(cis, ~ data.frame(bayestestR::sexit(data, ci = .x))) |>
    rowwise() |>
    filter(!between(0, round_half_up(CI_low, 2), round_half_up(CI_high, 2)) ||
             between(0, round_half_up(CI_low, 2), round_half_up(CI_high, 2)) & CI == 0.80) |>
    group_by(Parameter) |>
    slice_max(CI) |>
    mutate(across(c(Median, CI_low, CI_high), ~ round_half_up(.x, 2)),
           Direction = round_half_up(Direction, 5),
           CI = ifelse(between(0, CI_low, CI_high), "<80", as.character(CI * 100))) |>
    ungroup() |> 
    select(Parameter, Median, CI, CI_low, CI_high, Direction) |>
    mutate(text = glue::glue("{format(Median, digits = 2)}, {CI}% [{format(CI_low, digits = 2)}, {format(CI_high, digits = 2)}]"))
  
  return(foo)
}

# Generalized function to extract posterior draws and generate figure data
generate_figure_data <- function(species_model, trait_model, type, draw_min, draw_max) {

  
  parameter_order <- rev(c("Species Number",
                           "Spatial Distance [Minimum]",
                           "Spatial Extent",
                           "Absolute Latitude [Mean]",
                           "<i>Homogenisation – Differentiation</i>",
                           "<i>Freshwater – Terrestrial</i>",
                           "<i>Agriculture – Multiple</i>",
                           "<i>Forest – Multiple</i>",
                           "<i>Urban – Multiple</i>",
                           "Human Pressure [Baseline]",
                           "Human Pressure [Range]",
                           "<i>1.5km – 1km</i>",
                           "<i>2km – 1km</i>"))
  
  # Remove specific parameters for direction models
  if (type == "direction") {
    custom_parameter_order <- parameter_order[!parameter_order %in% c("<i>Homogenisation – Differentiation</i>")]
    var_prefix <- "b_"
    names <- c("coef","parameter")
  } else if (type == "magnitude") {
    custom_parameter_order <- parameter_order
    names <- c("coef","parameter")
    var_prefix <- "b_"
  } else if (type == "shape") {
    custom_parameter_order <- parameter_order
    var_prefix <- "b_mu"  # For shape, use `b_mu` prefix
    names <- c("coef","shape","parameter")
  }
  
  # Get support for taxonomic and functional datasets
  text_species<- get_support(as_draws_df(species_model) |> select(contains(var_prefix))) |>
    separate_wider_delim(Parameter, delim = "_", names = names, too_many="merge") |>
    mutate(text = glue::glue("{format(Median, digits = 2)}, {CI}% [{format(CI_low, digits = 2)}, {format(CI_high, digits = 2)}]")) |>
    select(any_of(c("coef","shape", "parameter", "text"))) |> 
    mutate(text = as.character(text))
  
  text_trait <- get_support(as_draws_df(trait_model) |> select(contains(var_prefix))) |>
    separate_wider_delim(Parameter, delim = "_",  names = names , too_many="merge") |>
    mutate(text = glue::glue("{format(Median, digits = 2)}, {CI}% [{format(CI_low, digits = 2)}, {format(CI_high, digits = 2)}]")) |>
    select(any_of(c("coef","shape", "parameter", "text"))) |> 
     mutate(text = as.character(text))
  
  # # Extract posterior draws and process them for taxonomic and functional datasets
  draws_species <- species_model |>
    as_draws_df(reserved = FALSE, variable = var_prefix, regex = TRUE) |> # Extract taxonomic draws
    as_tibble() |>
    pivot_longer(cols = !c(.chain, .iteration, .draw), names_to = "variable", values_to = "posterior") |>
    separate(variable, into = names, extra = "merge") |>
    mutate(facet = "Species replacement") |>
    left_join(text_species)

  draws_trait <- trait_model |>
    as_draws_df(reserved = FALSE, variable = var_prefix, regex = TRUE) |> # Extract functional draws
    as_tibble() |>
    pivot_longer(cols = !c(.chain, .iteration, .draw), names_to = "variable", values_to = "posterior") |>
    separate(variable, into = names, extra = "merge") |>
    mutate(facet = "Trait replacement") |>
    left_join(text_trait)
  
  # # Bind taxonomic and functional datasets
  figure_data <- bind_rows(draws_species, draws_trait) |>
    mutate(facet = factor(facet, levels = c("Species replacement", "Trait replacement"))) |> # Create a factor for applying facet
    select(any_of(c("shape","parameter","posterior","facet","text"))) |>
    group_by(across(where(is.character))) |>
    mutate(mean_hdi(posterior, .width = .8)) |>
    mutate(parameter = clean_parameters(parameter)) |>  # Clean and rename parameters
    filter(!parameter %in% "Intercept") |>
    mutate(parameter = case_when(grepl("–",parameter) ~ glue::glue("<i>{parameter}</i>"),
                                 TRUE ~ parameter)) |> 
    ungroup() |>
    mutate(parameter = factor(parameter, levels = custom_parameter_order)) |>
    #arrange(facet, parameter) |>
    mutate(grid_fill = as.numeric(parameter) %% 2 == 0,  # Create gray/white tire pattern
           draw_min = draw_min,  # Adjust based on figure type
           draw_max = draw_max) |> # Adjust based on figure type
    rowwise() |>
    mutate(significant = !between(0, ymin, ymax)) |>
    ungroup() |>
    mutate(parameter = factor(parameter, levels = custom_parameter_order)) 
  
  if (type == "shape") {
    figure_data<- figure_data |>  
      mutate(shape = gsub("mu","",shape))  |> 
      mutate(shape = gsub("Revlog","Reverse logistic",shape)) 
  }
   
   return(figure_data)
}


# Example usage:
# direction_figure_data <- generate_figure_data(fit_direction_tax, fit_direction_fun, "direction")
# magnitude_figure_data <- generate_figure_data(fit_magnitude_tax, fit_magnitude_fun, "magnitude")
# shape_figure_data <- generate_figure_data(fit_shape_tax, fit_shape_fun, "shape")

# Generalized plotting function with text
plot_posterior_draws <- function(data,  add_annotations = FALSE,
                                   facet_var = "facet", 
                                   trunc_upper = NULL, trunc_lower = NULL, 
                                   tick_factor = NULL, xlims = NULL) {
   # data=magnitude_figure_data
   # trunc_upper = 6
   # trunc_lower = -2
   # tick_factor = 2
   # xlims = c(-3,7)
   # add_annotations = FALSE
   # facet_var = "facet"
  # 
    # Define limits and truncation based on data if not provided
    posterior_range <- round_half_up(range(data$posterior, na.rm = TRUE),0)
    xlims <- xlims %||% c(round_half_up(posterior_range[1],0) - 0.1, round_half_up(posterior_range[2],0) + 0.1)
    trunc_lower <- trunc_lower %||% round_half_up(posterior_range[1],0)
    trunc_upper <- trunc_upper %||% round_half_up(posterior_range[2],0)
    tick_factor <- tick_factor %||% round_half_up(abs(diff(posterior_range)) / 5,0)
  
  plot_data <- data |>
    ggplot() +
    aes(x = posterior, y = as.numeric(parameter), group = parameter) +
    my_theme_extended_models()  # Use the new base theme
  
  if (length(facet_var) == 2) {
    plot_data <- plot_data + 
      ggh4x::facet_nested(cols = vars(!!sym(facet_var[1])), rows = vars(!!sym(facet_var[2])))
  } else { 
    plot_data <- plot_data + 
      facet_wrap(facet_var)
  }
  
   plot_data <- 
      plot_data + scale_y_continuous(breaks = seq_along(levels(data$parameter)), labels = levels(data$parameter)) +
    coord_cartesian(xlim = xlims) +
    geom_rect(
      data = ~ .x |> filter(grid_fill == TRUE),
      aes(
        xmin = draw_min,
        xmax = draw_max,
        ymin = as.numeric(parameter) - .3,
        ymax = as.numeric(parameter) + .7
      ),
      fill = "gray95"
    ) +
    stat_slab(
      aes(
        fill = significant,
        fill_ramp = after_stat(level)
      ),
      .width = c(.8, 0.89, .95),
      point_interval = median_hdci,
      interval_size_range = c(.5, 1),
      fatten_point = 1.2,
      height = 1,
      slab_size = .1,
      slab_colour = "gray70",
      normalize = "groups")+
    geom_text(
      data = data |>  distinct(across(any_of(c("response","shape", "parameter", "facet"))), .keep_all = TRUE), # Add text annotations from get_support
      aes(y = as.numeric(parameter) + 0.1, x = draw_max, label = text),
      hjust = 1, size = 1
    ) +
    scale_fill_manual(values = rev(c("#bc5090", "gray50"))) +
    geom_vline(aes(xintercept = 0), lty = "11", linewidth = 0.1) +
    scale_fill_ramp_discrete(na.translate = FALSE) +
    labs(
      fill_ramp = "Credible intervals",
      x = "Posterior draws",
      y = "Predictors",
      fill = "Supported at 80% HDI"
    ) +
  #  add_annotations(add_annotations) +  # Add annotations only for direction
      guides(x = guide_axis_truncated(trunc_lower = trunc_lower,
                                      trunc_upper = trunc_upper),
             fill_ramp = guide_legend(legend.box = "horizontal"),
             fill = guide_legend(legend.box = "horizontal")) +
      scale_x_continuous(breaks = seq(trunc_lower, trunc_upper, by = tick_factor))
    
    return(plot_data)
}

# Example usage:
# plot_figure(direction_figure_data, "direction", "Extended_Figure_2.pdf")
# plot_figure(magnitude_figure_data, "magnitude", "Extended_Figure_3.pdf")
# plot_figure(shape_figure_data, "shape", "Extended_Figure_4.pdf")


generate_ses_figure_data <- function(fit_model, model_type, draw_min = -7, draw_max = 12) {

  # Set up custom parameter order excluding specific terms for each model
  parameter_order <- rev(c("Species Number",
                           "Spatial Distance [Minimum]",
                           "Spatial Extent",
                           "Absolute Latitude [Mean]",
                           "<i>Homogenisation – Differentiation</i>",
                           "<i>Freshwater – Terrestrial</i>",
                           "<i>Agriculture – Multiple</i>",
                           "<i>Forest – Multiple</i>",
                           "<i>Urban – Multiple</i>",
                           "Human Pressure [Baseline]",
                           "Human Pressure [Range]",
                           "<i>1.5km – 1km</i>",
                           "<i>2km – 1km</i>"))
  
  # Remove model-specific terms from parameter order
  custom_parameter_order <- if (model_type == "direction") {
    parameter_order[!parameter_order %in% c("<i>Homogenisation – Differentiation</i>")]
  } else {
    parameter_order
  }
  
  # Extract model draws and compute posterior statistics
  text_fun <- get_support(as_draws_df(fit_model) |> 
                            select(contains("b_mu"))) 
  if(model_type == "shape"){
    text_fun<-text_fun |> separate_wider_delim(Parameter, delim = "_", names = c("coef", "response","shape", "parameter"), too_many = "merge")
  } else {
    text_fun<-text_fun |>separate_wider_delim(Parameter, delim = "_", names = c("coef", "response", "parameter"), too_many = "merge")
  }
  text_fun <- text_fun |> mutate(response = gsub("mu", "", response)) |>
    mutate(text = glue::glue("{format(Median, digits = 2)}, {CI}% [{format(CI_low, digits = 2)}, {format(CI_high, digits = 2)}]")) |> 
    select(any_of(c("coef", "response", "parameter", "text", "shape"))) |>
    mutate(text = as.character(text))
  
  # Process model draws
  model_draws <- fit_model |> 
    as_draws_df(reserved = FALSE, variable = "^b_mu", regex = TRUE) |> 
    as_tibble() |>
    pivot_longer(cols = !c(.chain, .iteration, .draw), names_to = "variable", values_to = "posterior") 
    
    if(model_type == "shape"){
      model_draws<-model_draws |>  separate(variable, into = c("coef", "response", "shape","parameter"), extra = "merge") 
    } else {
      model_draws<-model_draws |> separate(variable, into = c("coef", "response","parameter"), extra = "merge") 
    }
  
  model_draws<-model_draws |> 
    mutate(response = gsub("mu", "", response)) |> 
    mutate(facet = "Traits") |> 
    left_join(text_fun)
  
  # Final data wrangling for plotting
  processed_data <- model_draws |>
    select(any_of(c("response","shape", "parameter", "posterior", "facet", "text"))) |>
    mutate(parameter = clean_parameters(parameter)) |> 
    filter(parameter != "Intercept") |>
    mutate(parameter = case_when(
      grepl("–", parameter) ~ glue::glue("<i>{parameter}</i>"),
      TRUE ~ parameter
    )) |> 
    mutate(parameter = factor(parameter, levels = custom_parameter_order)) |> 
    mutate(grid_fill = as.numeric(parameter) %% 2 == 0,  # Create alternating fill pattern
           draw_min = draw_min, 
           draw_max = draw_max,
           significant = !grepl("<", text),
           response = str_to_title(response))
  
  return(processed_data)
}

# For direction model
# direction_data <- generate_ses_figure_data(fit_model = fit_divergence_direction, model_type = "direction")
# 
# # For magnitude model
# magnitude_data <- generate_ses_figure_data(fit_model = fit_magnitude_model, model_type = "magnitude")
# 
# # For shape model
# shape_data <- generate_ses_figure_data(fit_model = fit_shape_model, model_type = "shape")


extract_model_ses_data <- function(ses_data, model_data, model_type) {
  
  # Conditionally process based on model type
  if (model_type == "direction") {
    processed_data <-  
      model_data |> 
      filter(facet == "Traits") |> 
      left_join(ses_data |>  
                  select(dataset,direction, direction_SES,direction_pvalue) |>  
                  mutate(direction=str_to_title(gsub("z","s",direction)))) |> 
      drop_na() |> 
      mutate(direction_ses = factor(case_when(direction_pvalue <= 0.05 & direction_SES < 0 ~ "Smaller",
                                              direction_pvalue <= 0.05 & direction_SES > 0 ~ "Higher",
                                              TRUE ~ "Random"), levels=c("Random","Smaller","Higher")), .keep = "unused") |> 
      mutate(response = case_when(direction_harrel_davis < 0 & direction_ses == "Smaller"~"homogenisation",
                                  direction_harrel_davis < 0 & direction_ses == "Higher"~"differentiation",
                                  direction_harrel_davis > 0 & direction_ses == "Smaller"~"homogenisation",
                                  direction_harrel_davis > 0 & direction_ses == "Higher"~"differentiation", TRUE ~"Random"))  |> 
      mutate(response=factor(response))
    
  } else if (model_type == "magnitude") {
    processed_data <- 
      model_data |> 
      filter(facet == "Traits") |> 
      left_join(ses_data |>  
                  select(dataset,direction, magnitude_SES,magnitude_pvalue) |>
                  mutate(direction=gsub("z","s",direction))) |> 
      drop_na() |> 
      mutate(response = factor(case_when(magnitude_pvalue <= 0.05 & magnitude_SES < 0 ~ "Weaker",
                                         magnitude_pvalue <= 0.05 & magnitude_SES > 0 ~ "Stronger",
                                         TRUE ~ "Random"), levels=c("Random","Weaker","Stronger")))
    
    
  } else if (model_type == "shape") {
    processed_data <- 
      ses_data |> 
      select(dataset, direction, Absent_ses,Absent_pvalue,Revlog_ses,Revlog_pvalue,Saturating_ses,Saturating_pvalue, Exponential_ses, Exponential_pvalue) |> 
      mutate(direction=gsub("z","s",direction)) |> 
      ungroup() |> 
      pivot_longer(cols = contains("_ses") | contains("_pvalue"),
                   names_to = c(".value", "measure"),
                   names_pattern = "(.*)_([^_]+)") |>
      pivot_longer(Absent:Exponential, values_to = "value") |>
      pivot_wider(names_from = measure, values_from = value) |> 
      mutate(response = factor(case_when(pvalue <= 0.05 & sign(ses) == -1 ~ "Smaller",
                                         pvalue <= 0.05 & sign(ses) ==  1 ~ "Higher",
                                         TRUE ~ "Random"), levels=c("Random","Smaller","Higher"))) |> 
      select(dataset,direction,name,response) |> 
      pivot_wider(names_from="name", values_from="response") |> 
      right_join(model_data |> 
                   filter(facet == "Traits") |> 
                   select(-c(absent:revlog))) |> 
      drop_na()
  } else {
    stop("Invalid model_type. Choose from 'direction', 'magnitude', or 'shape'.")
  }
  
  return(processed_data)
}

run_and_save_model <- function(data, formula) {
  # Get formula name as a string to use in the file path
  formula_name <- deparse(substitute(formula))
  brms::brm(
    formula = formula,
    data = data,
    chains = CHAINS,
    family = brms::categorical(refcat = "Random"),
    #control = list(max_treedepth = 12, adapt_delta = .99),
    save_pars = save_pars(all = TRUE),
    seed = BAYES_SEED,
    file = paste0("S7_Model_outputs_figures_and_tables/model/convergence/", gsub("formula_", "", formula_name)),
    file_refit = "on_change"
  )
}
