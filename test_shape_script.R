
  
df <- het_species_land |> 
  filter(Predictor == "Multiple") |> 
  select(.draw,.category,direction,.epred) |> 
  rename(draw=.draw, shape = .category, probability = .epred)

regime_test <- function(df, width = 0.80){
  
  library(dplyr)
  library(tidyr)
  library(tidybayes)
  
  # ─────────────────────────────────────────────────────────────
  # helper: HDI separation test
  # ─────────────────────────────────────────────────────────────
  separated_hdi <- function(x, y, width){
    delta <- x - y
    h <- tidybayes::hdi(delta, .width = width)
    (h[,1] > 0) | (h[,2] < 0)
  }
  
  # ─────────────────────────────────────────────────────────────
  # build per-draw regime probabilities
  # ─────────────────────────────────────────────────────────────
  regime_draws <- df |>
    group_by(draw, shape, direction) |>
    summarise(p = mean(probability), .groups = "drop") |>
    mutate(
      shape = as.character(shape),
      direction = as.character(direction)
    )
  
  regime_wide <- regime_draws |>
    unite(regime, shape, direction, remove = FALSE) |>
    select(draw, regime, p) |>
    pivot_wider(names_from = regime, values_from = p)
  # ─────────────────────────────────────────────────────────────
  # STAGE 1 — resolve direction WITHIN each shape
  # one representative regime per shape
  # ─────────────────────────────────────────────────────────────
  shapes <- unique(regime_draws$shape)
  
  shape_candidates <- lapply(shapes, function(s){
    
    tmp <- regime_draws |> filter(shape == s)
    
    dirs <- unique(tmp$direction)
    
    # if only one direction exists
    if(length(dirs) == 1){
      draws <- tmp$p
      return(tibble(
        shape = s,
        direction = dirs,
        shape_draws = list(draws),
        regime_probability = median(draws)
      ))
    }
    
    # get per-draw vectors
    d1 <- tmp |> filter(direction == dirs[1]) |> pull(p)
    d2 <- tmp |> filter(direction == dirs[2]) |> pull(p)
    
    # direction contrast
    dir_sep <- separated_hdi(d1, d2, width)
    
    if(dir_sep){
      # dominant direction
      if(median(d1) >= median(d2)){
        dom_draws <- d1
        dom_dir   <- dirs[1]
      } else {
        dom_draws <- d2
        dom_dir   <- dirs[2]
      }
    } else {
      # direction inconclusive → keep best direction draws
      if(median(d1) >= median(d2)){
        dom_draws <- d1
      } else {
        dom_draws <- d2
      }
      dom_dir <- "inconclusive"
    }
    
    tibble(
      shape = s,
      direction = dom_dir,
      shape_draws = list(dom_draws),
      regime_probability = median(dom_draws)
    )
    
  }) |> bind_rows()
  
  # rank shapes by probability
  shape_candidates <- shape_candidates |>
    arrange(desc(regime_probability))
  
  
  # ─────────────────────────────────────────────────────────────
  # STAGE 2 — hierarchical shape dominance
  # stop at first separation
  # ─────────────────────────────────────────────────────────────
  kept <- list()
  kept[[1]] <- shape_candidates[1,]
  
  for(i in 1:(nrow(shape_candidates)-1)){
    i=1
    s1 <- shape_candidates$shape[i]
    s2 <- shape_candidates$shape[i+1]
    
    d1 <- shape_candidates$shape_draws[[i]]
    d2 <- shape_candidates$shape_draws[[i+1]]

    sep <- separated_hdi(d1, d2, width)
    
    if(sep){
      # dominance boundary found → stop
      break
    } else {
      kept[[i+1]] <- shape_candidates[i+1,]
    }
  }
  
  bind_rows(kept) |>
    select(shape, direction, regime_probability) |>
    arrange(desc(regime_probability))
}

