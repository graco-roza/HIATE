library(tidybayes)
library(ggplot2)

# Extract random intercept estimates for taxa_coarse.
# Adjust the parameter name if your model uses a different naming convention.
ranef_taxa <- md_species_base %>%
  tidybayes::spread_draws(r_taxa_coarse[taxa, `(Intercept)`]) |> 
  mutate(facet = "Species replacement") |> 
  bind_rows(
    md_trait_base %>%
      tidybayes::spread_draws(r_taxa_coarse[taxa, `(Intercept)`]) |> 
      mutate(facet = "Traits replacement")
  )

# Plot the random intercept estimates for each taxa (ordered by their median value)
ggplot(ranef_taxa, aes(x = r_taxa_coarse, y = taxa)) +
  stat_halfeye(
    aes(fill_ramp = after_stat(x)),
    point_interval = median_qi,
    .width = c(0.8, 0.9, 0.95),
    alpha = 1
  ) +
    labs(x = "Random Intercept Estimate",
       y = "Coarse biotic group",
       title = "Random Effects for relationship direction") +
  theme_minimal(base_size = 14)+facet_wrap(~facet)+
  coord_cartesian(xlim=c(-.5,.5),clip="off")


library(ggplot2)
library(tidybayes)

# Example: random intercept draws for a grouping factor
# Replace with your own data + parameter name
ranef_taxa <- best_model_traits %>%
  spread_draws(r_taxa_coarse[taxa, "(Intercept)"])

ggplot(ranef_taxa, aes(x = `(Intercept)`, y = reorder(taxa, `(Intercept)`))) +
  # Use fill = after_stat(x) to color by the x-value (i.e. the intercept)
  stat_halfeye(
    aes(fill = after_stat(x)),
    point_interval = median_qi,
    .width = c(0.8, 0.9, 0.95),
    alpha = 0.9
  ) +
  scale_fill_gradient2(
    midpoint = 0,  # or another midpoint if needed
    low = "blue",
    mid = "white",
    high = "red"
  ) +
  labs(
    x = "Random Intercept Estimate",
    y = "Taxa coarse",
    fill = "Intercept Value"
  ) +
  theme_minimal(base_size = 14)



