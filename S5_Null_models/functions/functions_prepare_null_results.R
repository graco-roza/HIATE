require(Hmisc)

pooled <- function(x, y, FUN) {
  nx <- length(x)
  ny <- length(y)
  sqrt(((nx - 1) * FUN(x) ^ 2 + (ny - 1) * FUN(y) ^ 2) / (nx + ny - 2))
}
hdmedian <- function(x) as.numeric(hdquantile(x, 0.5))
hdmad <- function(x) 1.4826 * hdmedian(abs(x - hdmedian(x)))
phdmad <- function(x, y) pooled(x, y, hdmad)
gammaEffectSize <- function(x, y, prob){
  #if(length(na.exclude(y)) < length(na.exclude(x)))  y<-c(na.exclude(y), rep(0, length(na.exclude(x))  - length(na.exclude(y))))
  if(length(na.exclude(y)) < 300)  y<-c(na.exclude(y), rep(0, length(na.exclude(x))  - length(na.exclude(y))))
  #if(length(na.exclude(x)) < length(na.exclude(y))) x<-c(na.exclude(x), rep(0, length(na.exclude(y))  - length(na.exclude(x))))
  if(length(na.exclude(x)) < 300)  x<-c(na.exclude(x), rep(0, length(na.exclude(y))  - length(na.exclude(x))))
  res<-as.numeric((hdquantile(na.exclude(y), prob) - hdquantile(na.exclude(x), prob)) / phdmad(na.exclude(x), na.exclude(y)))
  return(res)
}


estimate_ses <- function (obs, est, param = TRUE, p = TRUE) 
{
  if (param) {
    est_sd <- sd(est, na.rm = TRUE)
    if (is.na(est_sd)) est_sd <- 0
    res = (obs - mean(est, na.rm=TRUE))/(est_sd)
    res = ifelse(is.nan(res),0,res)
    pval = pnorm(res)
  }
  else {
    est = c(obs, est)
    pval = (sum(est < obs, na.rm=TRUE) + sum(est == obs, na.rm=TRUE)/2)/length(na.exclude(est))
    res = qnorm(pval)
  }
  if (p) {
    pval = pval * 2
    if (pval > 1) 
      pval = 1 - (pval - 1)
    res = c(res, pval)
    names(res) = c("ses","pvalue")
  }
  return(res)
}

estimate_ses_weighted <- function(obs, est, weights, param = TRUE, p = TRUE) {
  if (param) {
    weighted_mean = sum(est * weights) / sum(weights)
    weighted_sd = sqrt(sum((est - weighted_mean)^2 * weights) / sum(weights))
    if (is.na(weighted_sd)) weighted_sd <- 0
    res = (obs - weighted_mean) / weighted_sd
    res = ifelse(is.nan(res),0,res)
    pval = pnorm(res)
  } else {
    est = c(obs, est)
    pval = (sum(est < obs) + sum(est == obs)/2) / length(est)
    res = qnorm(pval)
  }
  if (p) {
    pval = pval * 2
    if (pval > 1) 
      pval = 1 - (pval - 1)
    res = c(res, pval)
    names(res) = c("SES", "pvalue")
  }
  return(res)
}


# Function to calculate SES for observed and null models
# Arguments:
#   null_model: Null model data
#   observed_model: Observed model data
#   param: Boolean indicating whether to use parametric estimation (default: TRUE)
#   p: Boolean indicating whether to calculate p-values (default: TRUE)
# Returns:
#   SES results including SES, p-values, and other relevant information
estimate_SES_direction <-  function(observed_model=NULL, null_model=NULL,param=TRUE) {
  
  # Calculate null model output
  output_null <- null_model %>%
map_depth(2, ~ .x |>  pluck("r2")) |>
map_df(~ tibble(
ID = seq_along(.x$differentiation),
differentiation = .x$differentiation,
homogenization =  .x$homogenization
), .id = "source") |>
    pivot_longer(cols=!c(source,ID), names_to = "direction", values_to = "r2") |> 
    dplyr::mutate(across(c(r2), ~ ifelse(.x <= 0, NA, .x/100))) |>
    group_by(source) |> 
    summarise(direction_r2 = gammaEffectSize(y=r2[direction != "differentiation"],x=r2[direction == "differentiation"],prob=.5)
    )

# Calculate observed model output
  output_observed <- observed_model |> 
    pluck("Functional") |> 
    map_depth(1, ~ .x |>  pluck("r2"))  |> 
    bind_rows(.id = "direction") |> 
    mutate(differentiation = differentiation,
           homogenization = homogenization,
           ID = seq_along(differentiation),
           source = "observed", .keep = "unused") |> 
    pivot_longer(cols=!c(source,ID), names_to = "direction", values_to = "r2") |> 
    dplyr::mutate(across(c(r2), ~ ifelse(.x <= 0, NA, .x/100))) |>
    group_by(source) |> 
    summarise(direction_r2 = gammaEffectSize(y=r2[direction != "differentiation"],x=r2[direction == "differentiation"],prob=.5)
    )
# Combine observed and null model output
  combined_output <- bind_rows(output_observed, output_null)

  # Calculate SES results

    SES_r2 <- combined_output %>%
    summarise(ses_results_diff = list(estimate_ses(direction_r2[source == "observed"], direction_r2[source != "observed"], param = param, p = TRUE))) |> 
    unnest_wider(ses_results_diff) |> 
    set_names(c("direction_SES","direction_pvalue"))
  
  return(SES_r2)
}

# Function to estimate the Standardized Effect Size (SES) of the strength of human footprint effect on traits
estimate_SES_magnitude <- function(observed_model=NULL, null_model=NULL,param=TRUE) {
  require(tidyverse)
    
  # Process null_model data
   output_null<-  null_model |> 
     enframe(name = "source", value = "iteration_value") |> 
     mutate(iteration_value = map(iteration_value, ~ enframe(.x, name = "direction", value = "direction_value"))) |> 
     unnest(iteration_value, keep_empty = TRUE) |>
     mutate(coef = map(direction_value, ~ .x |>  
                         pluck("coefs") |>  
                         bind_rows() |>  
                         rownames_to_column("predictor") |> 
                         filter(grepl("^hfp_", predictor))),.keep="unused") |> unnest(coef) |> 
     group_by(source, direction) |> 
     summarise(magnitude = mean(coefficient.1+coefficient.2+coefficient.3)*(n()/1000)) |> 
     mutate(direction = factor(direction, levels=c("differentiation","homogenization"))) |> 
     complete(direction) |> 
     mutate(across(c(magnitude),~ ifelse(is.na(.x),0,.x)))
  
  

  output_observed <- observed_model |> 
    pluck("Functional") |> 
    enframe(name = "direction", value = "direction_value") |> 
    mutate(coef = map(direction_value, ~ .x |>  
                        pluck("coefs") |>  
                        bind_rows() |>  
                        rownames_to_column("predictor") |> 
                        filter(grepl("^hfp_", predictor))),.keep="unused") |> unnest(coef) |> 
    select(-predictor) |> 
    group_by(direction) |> 
    summarise(magnitude = mean(coefficient.1+coefficient.2+coefficient.3)*(n()/1000)) |> 
    mutate(direction = factor(direction, levels=c("differentiation","homogenization"))) |> 
    complete(direction) |> 
    mutate(across(c(magnitude),~ ifelse(is.na(.x),0,.x))) |> 
    mutate(source = "observed") |> 
    ungroup()
  
  
  combined_output <- bind_rows(output_observed, output_null)
  
  # Calculate SES results
  SES_coef <- combined_output %>%
    group_by(direction) |> 
    summarise(magnitude_ses = list(estimate_ses(magnitude[source == "observed"], magnitude[source != "observed"], param = param, p = TRUE))) |> 
    unnest_wider(magnitude_ses) |> 
    set_names(c("direction","magnitude_SES","magnitude_pvalue"))
  
  return(SES_coef)
  
}



# Function to estimate the SES of the probability of a given shape between dissimilarities and human footprint
estimate_SES_shape <-  function(observed_model=NULL, null_model=NULL,param=TRUE) {
  require(tidyverse)

  output_null <- null_model |> 
    enframe(name = "source", value = "iteration_value") |> 
    mutate(iteration_value = map(iteration_value, ~ enframe(.x, name = "direction", value = "direction_value")),.keep="unused") |> 
    unnest(iteration_value, keep_empty = TRUE) |> 
    mutate(coef = map(direction_value, ~ .x |>  
                        pluck("coefs") |>  
                        bind_rows() |>  
                        rownames_to_column("predictor") |> 
                        filter(grepl("^hfp_", predictor)))) |> unnest(coef) |> 
    dplyr::mutate(shape = dplyr::case_when(coefficient.1 + coefficient.2 + coefficient.3 == 0 ~ "Absent",
                                           coefficient.1 <= coefficient.2 & coefficient.2 < coefficient.3 ~ "Exponential",
                                           coefficient.1 > coefficient.2 & coefficient.2 > coefficient.3 ~ "Saturating",
                                           coefficient.1 < coefficient.2 & coefficient.2 > coefficient.3 ~ "Saturating",
                                           coefficient.1 < coefficient.2 & coefficient.2 < coefficient.3 ~ "Saturating",
                                           coefficient.1 > coefficient.2 & coefficient.2 == coefficient.3 ~ "Saturating",
                                           coefficient.1 > coefficient.2 & coefficient.2 < coefficient.3 & coefficient.3 > 0 ~ "Revlog"))  |> 
    group_by(source, direction) |> 
    mutate(shape = factor(shape, levels =c("Absent", "Saturating","Exponential","Revlog"))) |> 
    mutate(nrow = row_number()) |> 
    count(shape) |>  
    complete(shape) |> 
    mutate(n = ifelse(is.na(n),0,n)) |> 
    pivot_wider(names_from = shape, values_from = n, values_fill = 0) |> 
    mutate(across(Absent:Revlog, ~ .x / (Absent + Saturating + Exponential + Revlog))) 
  # Process observed data
  output_observed <- observed_model |> 
    pluck("Functional") |> 
    enframe(name = "direction", value = "direction_value") |> 
    mutate(coef = map(direction_value, ~ .x |>  
                        pluck("coefs") |>  
                        bind_rows() |>  
                        rownames_to_column("predictor") |> 
                        filter(grepl("^hfp_", predictor))),.keep="unused") |> unnest(coef) |> 
    mutate(direction = factor(direction, levels=c("differentiation","homogenization"))) |> 
    complete(direction) |> 
    mutate(across(c(coefficient.1,coefficient.2,coefficient.3),~ ifelse(is.na(.x),0,.x))) |> 
    dplyr::mutate(shape = dplyr::case_when(coefficient.1 + coefficient.2 + coefficient.3 == 0 ~ "Absent",
                                           coefficient.1 <= coefficient.2 & coefficient.2 < coefficient.3 ~ "Exponential",
                                           coefficient.1 > coefficient.2 & coefficient.2 > coefficient.3 ~ "Saturating",
                                           coefficient.1 < coefficient.2 & coefficient.2 > coefficient.3 ~ "Saturating",
                                           coefficient.1 < coefficient.2 & coefficient.2 < coefficient.3 ~ "Saturating",
                                           coefficient.1 > coefficient.2 & coefficient.2 == coefficient.3 ~ "Saturating",
                                           coefficient.1 > coefficient.2 & coefficient.2 < coefficient.3 & coefficient.3 > 0 ~ "Revlog"))  |> 


    group_by(direction) |> 
    mutate(shape = factor(shape, levels =c("Absent", "Saturating","Exponential","Revlog"))) |> 
    count(shape) |>  
    complete(shape) |> 
    mutate(n = ifelse(is.na(n),0,n)) |> 
    pivot_wider(names_from = shape, values_from = n, values_fill = 0) |> 
    mutate(across(Absent:Revlog, ~ .x / (Absent + Saturating + Exponential + Revlog))) |> 
    mutate(source = "observed") 
  
  # Combine null and observed outputs
  combined_output <- bind_rows(output_observed, output_null)
  
  # Calculate SES results
  SES_shape <- combined_output %>%
    group_by(direction) |> 
    summarise(across(c(Absent, Exponential, Saturating, Revlog), ~ list(estimate_ses(.x[source == "observed"], .x[source != "observed"], param = param, p = TRUE)))) |> 
    unnest_wider(c(Absent, Exponential, Saturating, Revlog), names_sep = "_")
  
  return(SES_shape)
}


# output_null <- null_model |> 
#   enframe(name = "source", value = "iteration_value") |> 
#   mutate(iteration_value = map(iteration_value, ~ enframe(.x, name = "direction", value = "direction_value")),.keep="unused") |> 
#   unnest(iteration_value, keep_empty = TRUE) |> 
#   mutate(coef = map(direction_value, ~ .x |>  
#                       pluck("coefs") |>  
#                       bind_rows() |>  
#                       rownames_to_column("predictor") |> 
#                       filter(grepl("^hfp_", predictor)))) |> unnest(coef) |> 
#   dplyr::mutate(shape = dplyr::case_when(coefficient.1 + coefficient.2 + coefficient.3 == 0 ~ "Absent",
#                                          coefficient.1 <= coefficient.2 & coefficient.2 < coefficient.3 ~ "Exponential",
#                                          coefficient.1 > coefficient.2 & coefficient.2 > coefficient.3 ~ "Saturating",
#                                          coefficient.1 < coefficient.2 & coefficient.2 > coefficient.3 ~ "Saturating",
#                                          coefficient.1 < coefficient.2 & coefficient.2 < coefficient.3 ~ "Saturating",
#                                          coefficient.1 > coefficient.2 & coefficient.2 == coefficient.3 ~ "Saturating",
#                                          coefficient.1 > coefficient.2 & coefficient.2 < coefficient.3 & coefficient.3 > 0 ~ "Revlog"))  |> 
#   group_by(source, direction) |> 
#   mutate(shape = factor(shape, levels =c("Absent", "Saturating","Exponential","Revlog"))) |> 
#   mutate(nrow = row_number()) |> 
#   count(shape) |>  
#   complete(shape) |> 
#   mutate(n = ifelse(is.na(n),0,n)) |> 
#   pivot_wider(names_from = shape, values_from = n, values_fill = 0) |> 
#   mutate(across(Absent:Revlog, ~ .x / (Absent + Saturating + Exponential + Revlog))) 
# 
# # Process observed data
# output_observed <- observed_model |> 
#   pluck("Functional") |> 
#   enframe(name = "direction", value = "direction_value") |> 
#   mutate(coef = map(direction_value, ~ .x |>  
#                       pluck("coefs") |>  
#                       bind_rows() |>  
#                       rownames_to_column("predictor") |> 
#                       filter(grepl("^hfp_", predictor))),.keep="unused") |> unnest(coef) |> 
#   dplyr::mutate(shape = dplyr::case_when(coefficient.1 + coefficient.2 + coefficient.3 == 0 ~ "Absent",
#                                          coefficient.1 <= coefficient.2 & coefficient.2 < coefficient.3 ~ "Exponential",
#                                          coefficient.1 > coefficient.2 & coefficient.2 > coefficient.3 ~ "Saturating",
#                                          coefficient.1 < coefficient.2 & coefficient.2 > coefficient.3 ~ "Saturating",
#                                          coefficient.1 < coefficient.2 & coefficient.2 < coefficient.3 ~ "Saturating",
#                                          coefficient.1 > coefficient.2 & coefficient.2 == coefficient.3 ~ "Saturating",
#                                          coefficient.1 > coefficient.2 & coefficient.2 < coefficient.3 & coefficient.3 > 0 ~ "Revlog"))  |> 
#   group_by(direction) |> 
#   mutate(shape = factor(shape, levels =c("Absent", "Saturating","Exponential","Revlog"))) |> 
#   mutate(nrow = row_number()) |> 
#   count(shape) |>  
#   complete(shape) |> 
#   mutate(n = ifelse(is.na(n),0,n)) |> 
#   pivot_wider(names_from = shape, values_from = n, values_fill = 0) |> 
#   mutate(across(Absent:Revlog, ~ .x / (Absent + Saturating + Exponential + Revlog))) |> 
#   mutate(source = "observed") 
