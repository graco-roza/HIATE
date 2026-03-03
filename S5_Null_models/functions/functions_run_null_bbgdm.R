# ==============================================================================
# Script: functions_run_null_bbgdm.R
# Author: Caio Graco-Roza
# Last update: 2024-11-24
# Description: 
#   This script contains functions and workflows for analyzing beta diversity 
#   using hypervolume models, calculating standardized effect sizes (SES),
#   and estimating the relationship between observed and null models.
# 
# Workflow:
# 1. Append nested lists recursively for combining data.
# 2. Compute alpha diversity using hypervolume Gaussian models.
# 3. Calculate pairwise beta diversity components (total, replacement, richness).
# 4. Run Bayesian Bootstrap Generalized Dissimilarity Modeling (BBGDM).
# 5. Estimate SES for direction, magnitude, and shape of predictor effects.
# 
# Key Functions:
# - appendList: Recursively appends elements from one list to another.
# - run_null: Executes null model calculations for beta diversity components.
# - pooled: Calculates pooled variance or standard deviation.
# - hdmedian, hdmad: High-density median and MAD calculations.
# - gammaEffectSize: Computes gamma effect size for predictors.
# - estimate_ses, estimate_ses_weighted: Standardized Effect Size calculations.
# - estimate_SES_direction, estimate_SES_magnitude, estimate_SES_shape:
#   Functions to calculate SES for direction, magnitude, and shape of predictors.
# 
# Dependencies:
# - Libraries: tidyverse, future, snow, doSNOW, hypervolume, and custom functions.
# - External Functions: S4_run_BBGDM/functions/helper_functions_S4.R.
#
# Notes:
# - Parallel computation is used for efficiency in high-performance environments.
# - Ensure adequate memory and cores for large datasets.
# ==============================================================================

#combine the results from beta diversite into organized lists 
source("S4_run_BBGDM/functions/helper_functions_S4.R") #borrow some functions from the previous step

#' Append Nested Lists
#' 
#' Recursively appends elements of one list to another while handling nested structures.
#' 
#' @param x A list to which elements will be appended.
#' @param val A list of values to append to `x`.
#' @return A combined list with appended elements.
#' @examples
#' list1 <- list(a = 1, b = list(c = 2))
#' list2 <- list(b = list(d = 3), e = 4)
#' appendList(list1, list2)
appendList <- function (x, val) {
  stopifnot(is.list(x), is.list(val))
  xnames <- names(x)
  for (v in names(val)) {
    x[[v]] <- if (v %in% xnames && is.list(x[[v]]) && is.list(val[[v]])) 
      appendList(x[[v]], val[[v]])
    else c(x[[v]], val[[v]])
  }
  x
}
# function to run null beta models --------------------------------------------------------------------------------

#' Run Null Model for Beta Diversity
#' 
#' Computes alpha diversity using hypervolume models and calculates beta diversity components (total, replacement, and richness).
#' 
#' @param comm A community matrix where rows represent sites and columns represent species.
#' @param traits A matrix of traits associated with species in `comm`.
#' @param pred Predictor variables for BBGDM modeling.
#' @param component A character string specifying the beta diversity component to use ("Btotal", "Brepl", or "Brich").
#' @return A list containing the results of the BBGDM model for the specified beta diversity component.
#' @examples
#' comm <- matrix(sample(0:1, 20, replace = TRUE), nrow = 5)
#' traits <- matrix(rnorm(15), nrow = 5)
#' pred <- list(Functional = matrix(rnorm(10), nrow = 5))
#' run_null(comm, traits, pred, "Btotal")
run_null <- function(comm,traits,pred,component){
  
  cores <- future::availableCores()-1
  
  # 3a. Alpha hypervolumes – Abundance-based
  
  cl <- makeSOCKcluster(cores)
  registerDoSNOW(cl)
  withDoRNG(
    alpha.FD_abun <- foreach(
      i = 1:nrow(comm),
      .combine = hypervolume::hypervolume_join,
      .multicombine = TRUE,
      .errorhandling = 'remove',
      .options.future = list(seed = TRUE)
    ) %dopar%  {
      .libPaths(c("/projappl/project_2004932/project_rpackages_4.4.0", .libPaths()))
      set.seed(123)
      present <- which(comm[i, ] > 0)
      subset <- traits[names(comm[i, present]), 1:3]
      abun <- comm[i, present]
      hypervolume::hypervolume_gaussian(
        data = subset,
        verbose = FALSE,
        kde.bandwidth = hypervolume::estimate_bandwidth(subset, method="silverman-1d"),
        name = rownames(comm)[i],
        weight = as.numeric(abun / sum(abun))
      )
    })
  stopCluster(cl)
  
  # 3b. Alpha hypervolumes – Presence-absence (equal weights)
  cl <- makeSOCKcluster(cores)
  registerDoSNOW(cl)
  withDoRNG(
    alpha.FD_pa <- foreach(
      i = 1:nrow(comm),
      .combine = hypervolume::hypervolume_join,
      .multicombine = TRUE,
      .errorhandling = 'remove',
      .options.future = list(seed = TRUE)
    ) %dopar%  {
      set.seed(123)
      .libPaths(c("/projappl/project_2004932/project_rpackages_4.4.0", .libPaths()))
      present <- which(comm[i, ] > 0)
      subset <- traits[names(comm[i, present]), 1:3]
      abun <- comm[i, present]
      hypervolume::hypervolume_gaussian(
        data = subset,
        verbose = FALSE,
        kde.bandwidth = hypervolume::estimate_bandwidth(subset, method="silverman-1d"),
        name = rownames(comm)[i]
      )
    })
  stopCluster(cl)
  
  # Helper function to compute both Podani and Baselga partitions from two hypervolumes.
  # It returns a list with elements "Podani" and "Baselga".
  compute_partitions <- function(hv1, hv2) {
    hvset <- hypervolume::hypervolume_set(hv1, hv2,
                                          check.memory = FALSE,
                                          verbose = FALSE)
    shared   <- hvset[[3]]@Volume
    union   <- hvset[[4]]@Volume
    unique1 <- hvset[[5]]@Volume
    unique2 <- hvset[[6]]@Volume
    
    #sorensen total
    union_adj <- 2 * shared + unique1 + unique2
    
    
    # Podani partitions:
    podani_Btotal <- (unique1 + unique2) / union_adj
    podani_Brepl  <- 2 * min(unique1, unique2) / union_adj
    podani_Brich  <- abs(unique1 - unique2) / union_adj
    
    # Baselga partitions:
    baselga_Btotal <- (unique1 + unique2) / union_adj
    baselga_Bturnover <- min(unique1, unique2) / (shared + min(unique1, unique2))
    baselga_Bnestedness <- baselga_Btotal - baselga_Bturnover
    
    list(Podani = list(Btotal = round(podani_Btotal,3),
                       Brepl  = round(podani_Brepl,3),
                       Brich  = round(podani_Brich,3)),
         Baselga = list(Btotal = round(baselga_Btotal,3),
                        Bturn = round(baselga_Bturnover,3),
                        Bnes = round(baselga_Bnestedness,3)))
  }
  
  # 3c. Pairwise comparisons for Functional Beta Diversity in one loop.
  # We assume both hypervolume lists have the same length.
  alpha_n <- length(alpha.FD_abun@HVList)
  name_sites <- sapply(seq_along(alpha.FD_abun@HVList), function(i) alpha.FD_abun@HVList[[i]]@Name)
  
  cl <- makeSOCKcluster(cores)
  registerDoSNOW(cl)
    pairwise_all <- foreach(i = 1:alpha_n, .combine = appendList) %:%
      foreach(j = i:alpha_n, .combine = appendList) %dopar% {
        .libPaths(c("/projappl/project_2004932/project_rpackages_4.4.0", .libPaths()))
        res_abun <- compute_partitions(alpha.FD_abun@HVList[[i]], alpha.FD_abun@HVList[[j]])
        res_pa   <- compute_partitions(alpha.FD_pa@HVList[[i]], alpha.FD_pa@HVList[[j]])
        list(abundance = res_abun, pa = res_pa)
      }
  stopCluster(cl)
  
  # Now, pairwise_all contains two main components:
  #   - pairwise_all$abundance: results for abundance-based hypervolumes (with Podani and Baselga partitions)
  #   - pairwise_all$pa: results for presence-absence-based hypervolumes.
  
  # Assemble distance matrices for Functional Beta Diversity:
  
  # For abundance-based hypervolumes:
  # Podani partitions
  output_beta_Podani_abun <- list(
    Btotal = matrix(NA, alpha_n, alpha_n),
    Brepl  = matrix(NA, alpha_n, alpha_n),
    Brich  = matrix(NA, alpha_n, alpha_n)
  )
  output_beta_Podani_abun$Btotal[lower.tri(output_beta_Podani_abun$Btotal, diag = TRUE)] <- 
    round(pairwise_all$abundance$Podani$Btotal, 3)
  output_beta_Podani_abun$Brepl[lower.tri(output_beta_Podani_abun$Brepl, diag = TRUE)] <- 
    round(pairwise_all$abundance$Podani$Brepl, 3)
  output_beta_Podani_abun$Brich[lower.tri(output_beta_Podani_abun$Brich, diag = TRUE)] <- 
    round(pairwise_all$abundance$Podani$Brich, 3)
  Functional_Podani_abun <- lapply(output_beta_Podani_abun, as.dist)
  Functional_Podani_abun <- lapply(Functional_Podani_abun, function(x) set_names(x, name_sites[1:alpha_n]))
  
  # Baselga partitions for abundance-based hypervolumes:
  output_beta_Baselga_abun <- list(
    Btotal = matrix(NA, alpha_n, alpha_n),
    Bturn = matrix(NA, alpha_n, alpha_n),
    Bnes = matrix(NA, alpha_n, alpha_n)
  )
  output_beta_Baselga_abun$Btotal[lower.tri(output_beta_Baselga_abun$Btotal, diag = TRUE)] <- 
    round(pairwise_all$abundance$Baselga$Btotal, 3)
  output_beta_Baselga_abun$Bturn[lower.tri(output_beta_Baselga_abun$Bturn, diag = TRUE)] <- 
    round(pairwise_all$abundance$Baselga$Bturn, 3)
  output_beta_Baselga_abun$Bnes[lower.tri(output_beta_Baselga_abun$Bnes, diag = TRUE)] <- 
    round(pairwise_all$abundance$Baselga$Bnes, 3)
  Functional_Baselga_abun <- lapply(output_beta_Baselga_abun, as.dist)
  Functional_Baselga_abun <- lapply(Functional_Baselga_abun, function(x) set_names(x, name_sites[1:alpha_n]))
  
  # For presence-absence-based hypervolumes:
  # Podani partitions
  output_beta_Podani_pa <- list(
    Btotal = matrix(NA, alpha_n, alpha_n),
    Brepl  = matrix(NA, alpha_n, alpha_n),
    Brich  = matrix(NA, alpha_n, alpha_n)
  )
  output_beta_Podani_pa$Btotal[lower.tri(output_beta_Podani_pa$Btotal, diag = TRUE)] <- 
    round(pairwise_all$pa$Podani$Btotal, 3)
  output_beta_Podani_pa$Brepl[lower.tri(output_beta_Podani_pa$Brepl, diag = TRUE)] <- 
    round(pairwise_all$pa$Podani$Brepl, 3)
  output_beta_Podani_pa$Brich[lower.tri(output_beta_Podani_pa$Brich, diag = TRUE)] <- 
    round(pairwise_all$pa$Podani$Brich, 3)
  Functional_Podani_pa <- lapply(output_beta_Podani_pa, as.dist)
  Functional_Podani_pa <- lapply(Functional_Podani_pa, function(x) set_names(x, name_sites[1:alpha_n]))
  
  # Baselga partitions for presence-absence-based hypervolumes:
  output_beta_Baselga_pa <- list(
    Btotal = matrix(NA, alpha_n, alpha_n),
    Bturn = matrix(NA, alpha_n, alpha_n),
    Bnes = matrix(NA, alpha_n, alpha_n)
  )
  output_beta_Baselga_pa$Btotal[lower.tri(output_beta_Baselga_pa$Btotal, diag = TRUE)] <- 
    round(pairwise_all$pa$Baselga$Btotal, 3)
  output_beta_Baselga_pa$Bturn[lower.tri(output_beta_Baselga_pa$Bturn, diag = TRUE)] <- 
    round(pairwise_all$pa$Baselga$Bturn, 3)
  output_beta_Baselga_pa$Bnes[lower.tri(output_beta_Baselga_pa$Bnes, diag = TRUE)] <- 
    round(pairwise_all$pa$Baselga$Bnes, 3)
  Functional_Baselga_pa <- lapply(output_beta_Baselga_pa, as.dist)
  Functional_Baselga_pa <- lapply(Functional_Baselga_pa, function(x) set_names(x, name_sites[1:alpha_n]))
  
  # -------------------------------------------------------------------------
  # 4. SAVE THE RESULTS
  # -------------------------------------------------------------------------
  # Each output is a list with two elements: Taxonomic and Functional.
  # Save each to its corresponding folder.


# Run BBGDM for the specified component.
bbgdm_pd_abun <- run_bbgdm(Functional_Podani_abun %>% pluck(2), pred$Podani_abun)
bbgdm_pd_pa <- run_bbgdm(Functional_Podani_pa %>% pluck(2), pred$Podani_pa)
bbgdm_bs_abun <- run_bbgdm(Functional_Baselga_abun %>% pluck(2), pred$Baselga_abun)
bbgdm_bs_pa <- run_bbgdm(Functional_Baselga_pa %>% pluck(2), pred$Baselga_pa)

fun_bbgdm<-list(Podani_abun = bbgdm_pd_abun,
                Podani_pa = bbgdm_pd_pa,
                Baselga_abun = bbgdm_bs_abun,
                Baselga_pa = bbgdm_bs_pa)
return(fun_bbgdm)
}

#' #' Pooled Variance Estimator
#' #'
#' #' Calculates a pooled variance or standard deviation from two vectors using a specified function (e.g., `sd` or `var`).
#' #'
#' #' @param x Numeric vector.
#' #' @param y Numeric vector.
#' #' @param FUN Function to compute variance or standard deviation (e.g., `sd` or `var`).
#' #' @return A single numeric value representing the pooled result.
#' #' @examples
#' #' x <- rnorm(10)
#' #' y <- rnorm(15)
#' #' pooled(x, y, sd)
#' pooled <- function(x, y, FUN) {
#'   nx <- length(x)
#'   ny <- length(y)
#'   sqrt(((nx - 1) * FUN(x) ^ 2 + (ny - 1) * FUN(y) ^ 2) / (nx + ny - 2))
#' }
#' 
#' #' High-Density Median
#' #'
#' #' Computes the high-density median using high-density quantile estimation.
#' #'
#' #' @param x Numeric vector.
#' #' @return Numeric value representing the high-density median.
#' #' @examples
#' #' x <- rnorm(10)
#' #' hdmedian(x)
#' hdmedian <- function(x) as.numeric(hdquantile(x, 0.5))
#' 
#' #' High-Density Median Absolute Deviation
#' #'
#' #' Computes the high-density median absolute deviation (MAD).
#' #'
#' #' @param x Numeric vector.
#' #' @return Numeric value representing the high-density MAD.
#' #' @examples
#' #' x <- rnorm(10)
#' #' hdmad(x)
#' hdmad <- function(x) 1.4826 * hdmedian(abs(x - hdmedian(x)))
#' 
#' #' Pooled High-Density MAD
#' #'
#' #' Calculates the pooled high-density MAD for two numeric vectors.
#' #'
#' #' @param x Numeric vector.
#' #' @param y Numeric vector.
#' #' @return Numeric value representing the pooled MAD.
#' #' @examples
#' #' x <- rnorm(10)
#' #' y <- rnorm(15)
#' #' phdmad(x, y)
#' phdmad <- function(x, y) pooled(x, y, hdmad)
#' 
#' 
#' #' Gamma Effect Size
#' #'
#' #' Computes a gamma effect size between two numeric vectors at a specified quantile.
#' #'
#' #' @param x Numeric vector.
#' #' @param y Numeric vector.
#' #' @param prob Quantile probability (e.g., 0.5 for the median).
#' #' @return Numeric value representing the gamma effect size.
#' #' @examples
#' #' x <- rnorm(10)
#' #' y <- rnorm(15)
#' #' gammaEffectSize(x, y, 0.5)
#' gammaEffectSize <- function(x, y, prob){
#'   #if(length(na.exclude(y)) < length(na.exclude(x)))  y<-c(na.exclude(y), rep(0, length(na.exclude(x))  - length(na.exclude(y))))
#'   if(length(na.exclude(y)) < 300)  y<-c(na.exclude(y), rep(0, length(na.exclude(x))  - length(na.exclude(y))))
#'   #if(length(na.exclude(x)) < length(na.exclude(y))) x<-c(na.exclude(x), rep(0, length(na.exclude(y))  - length(na.exclude(x))))
#'   if(length(na.exclude(x)) < 300)  x<-c(na.exclude(x), rep(0, length(na.exclude(y))  - length(na.exclude(x))))
#'   res<-as.numeric((hdquantile(na.exclude(y), prob) - hdquantile(na.exclude(x), prob)) / phdmad(na.exclude(x), na.exclude(y)))
#'   return(res)
#' }
#' 
#' 
#' #' Estimate Standardized Effect Size (SES)
#' #'
#' #' Computes SES using observed and estimated null model values.
#' #'
#' #' @param obs Observed value.
#' #' @param est Numeric vector of null model values.
#' #' @param param Logical, if `TRUE`, parametric estimation is used. Defaults to `TRUE`.
#' #' @param p Logical, if `TRUE`, p-values are computed. Defaults to `TRUE`.
#' #' @return Named numeric vector containing SES and p-value.
#' #' @examples
#' #' obs <- 0.5
#' #' est <- rnorm(100)
#' #' estimate_ses(obs, est)
#' estimate_ses <- function (obs, est, param = TRUE, p = TRUE) {
#'   if (param) {
#'     est_sd <- sd(est, na.rm = TRUE)
#'     if (is.na(est_sd)) est_sd <- 0
#'     res = (obs - mean(est, na.rm=TRUE))/(est_sd)
#'     res = ifelse(is.nan(res),0,res)
#'     pval = pnorm(res)
#'   }
#'   else {
#'     est = c(obs, est)
#'     pval = (sum(est < obs, na.rm=TRUE) + sum(est == obs, na.rm=TRUE)/2)/length(na.exclude(est))
#'     res = qnorm(pval)
#'   }
#'   if (p) {
#'     pval = pval * 2
#'     if (pval > 1) 
#'       pval = 1 - (pval - 1)
#'     res = c(res, pval)
#'     names(res) = c("ses","pvalue")
#'   }
#'   return(res)
#' }
#' 
#' #' Estimate Weighted Standardized Effect Size (SES)
#' #'
#' #' Computes SES with weighted null model values.
#' #'
#' #' @param obs Observed value.
#' #' @param est Numeric vector of null model values.
#' #' @param weights Numeric vector of weights.
#' #' @param param Logical, if `TRUE`, parametric estimation is used. Defaults to `TRUE`.
#' #' @param p Logical, if `TRUE`, p-values are computed. Defaults to `TRUE`.
#' #' @return Named numeric vector containing SES and p-value.
#' #' @examples
#' #' obs <- 0.5
#' #' est <- rnorm(100)
#' #' weights <- runif(100)
#' #' estimate_ses_weighted(obs, est, weights)
#' estimate_ses_weighted <- function(obs, est, weights, param = TRUE, p = TRUE) {
#'   if (param) {
#'     weighted_mean = sum(est * weights) / sum(weights)
#'     weighted_sd = sqrt(sum((est - weighted_mean)^2 * weights) / sum(weights))
#'     if (is.na(weighted_sd)) weighted_sd <- 0
#'     res = (obs - weighted_mean) / weighted_sd
#'     res = ifelse(is.nan(res),0,res)
#'     pval = pnorm(res)
#'   } else {
#'     est = c(obs, est)
#'     pval = (sum(est < obs) + sum(est == obs)/2) / length(est)
#'     res = qnorm(pval)
#'   }
#'   if (p) {
#'     pval = pval * 2
#'     if (pval > 1) 
#'       pval = 1 - (pval - 1)
#'     res = c(res, pval)
#'     names(res) = c("SES", "pvalue")
#'   }
#'   return(res)
#' }
#' 
#' 
#' #' Estimate SES for Direction
#' #'
#' #' Computes SES for observed and null models in terms of directional R² values.
#' #'
#' #' @param observed_model Observed model data.
#' #' @param null_model Null model data.
#' #' @param param Logical, if `TRUE`, parametric estimation is used. Defaults to `TRUE`.
#' #' @return Data frame containing SES results and p-values.
#' #' @examples
#' #' observed_model <- list(...)
#' #' null_model <- list(...)
#' #' estimate_SES_direction(observed_model, null_model)
#' estimate_SES_direction <-  function(observed_model=NULL, null_model=NULL,param=TRUE) {
#'   
#'   # Calculate null model output
#'   output_null <- null_model %>%
#'     map_depth(2, ~ .x |>  pluck("r2")) |>
#'     map_df(~ tibble(
#'       ID = seq_along(.x$differentiation),
#'       differentiation = .x$differentiation,
#'       homogenization =  .x$homogenization
#'     ), .id = "source") |>
#'     pivot_longer(cols=!c(source,ID), names_to = "direction", values_to = "r2") |> 
#'     dplyr::mutate(across(c(r2), ~ ifelse(.x <= 0, NA, .x/100))) |>
#'     group_by(source) |> 
#'     summarise(direction_r2 = gammaEffectSize(y=r2[direction != "differentiation"],x=r2[direction == "differentiation"],prob=.5)
#'     )
#'   
#'   # Calculate observed model output
#'   output_observed <- observed_model |> 
#'     pluck("Functional") |> 
#'     map_depth(1, ~ .x |>  pluck("r2"))  |> 
#'     bind_rows(.id = "direction") |> 
#'     mutate(differentiation = differentiation,
#'            homogenization = homogenization,
#'            ID = seq_along(differentiation),
#'            source = "observed", .keep = "unused") |> 
#'     pivot_longer(cols=!c(source,ID), names_to = "direction", values_to = "r2") |> 
#'     dplyr::mutate(across(c(r2), ~ ifelse(.x <= 0, NA, .x/100))) |>
#'     group_by(source) |> 
#'     summarise(direction_r2 = gammaEffectSize(y=r2[direction != "differentiation"],x=r2[direction == "differentiation"],prob=.5)
#'     )
#'   # Combine observed and null model output
#'   combined_output <- bind_rows(output_observed, output_null)
#'   
#'   # Calculate SES results
#'   
#'   SES_r2 <- combined_output %>%
#'     summarise(ses_results_diff = list(estimate_ses(direction_r2[source == "observed"], direction_r2[source != "observed"], param = param, p = TRUE))) |> 
#'     unnest_wider(ses_results_diff) |> 
#'     set_names(c("direction_SES","direction_pvalue"))
#'   
#'   return(SES_r2)
#' }
#' 
#' #' Estimate SES for Magnitude
#' #'
#' #' Computes SES for the strength of human footprint effects on traits.
#' #'
#' #' @param observed_model Observed model data.
#' #' @param null_model Null model data.
#' #' @param param Logical, if `TRUE`, parametric estimation is used. Defaults to `TRUE`.
#' #' @return Data frame containing SES results for magnitude and p-values.
#' #' @examples
#' #' observed_model <- list(...)
#' #' null_model <- list(...)
#' #' estimate_SES_magnitude(observed_model, null_model)
#' estimate_SES_magnitude <- function(observed_model=NULL, null_model=NULL,param=TRUE) {
#'   require(tidyverse)
#'   
#'   # Process null_model data
#'   output_null<-  null_model |> 
#'     enframe(name = "source", value = "iteration_value") |> 
#'     mutate(iteration_value = map(iteration_value, ~ enframe(.x, name = "direction", value = "direction_value"))) |> 
#'     unnest(iteration_value, keep_empty = TRUE) |>
#'     mutate(coef = map(direction_value, ~ .x |>  
#'                         pluck("coefs") |>  
#'                         bind_rows() |>  
#'                         rownames_to_column("predictor") |> 
#'                         filter(grepl("^hfp_", predictor))),.keep="unused") |> unnest(coef) |> 
#'     group_by(source, direction) |> 
#'     summarise(magnitude = mean(coefficient.1+coefficient.2+coefficient.3)*(n()/1000)) |> 
#'     mutate(direction = factor(direction, levels=c("differentiation","homogenization"))) |> 
#'     complete(direction) |> 
#'     mutate(across(c(magnitude),~ ifelse(is.na(.x),0,.x)))
#'   
#'   
#'   
#'   output_observed <- observed_model |> 
#'     pluck("Functional") |> 
#'     enframe(name = "direction", value = "direction_value") |> 
#'     mutate(coef = map(direction_value, ~ .x |>  
#'                         pluck("coefs") |>  
#'                         bind_rows() |>  
#'                         rownames_to_column("predictor") |> 
#'                         filter(grepl("^hfp_", predictor))),.keep="unused") |> unnest(coef) |> 
#'     select(-predictor) |> 
#'     group_by(direction) |> 
#'     summarise(magnitude = mean(coefficient.1+coefficient.2+coefficient.3)*(n()/1000)) |> 
#'     mutate(direction = factor(direction, levels=c("differentiation","homogenization"))) |> 
#'     complete(direction) |> 
#'     mutate(across(c(magnitude),~ ifelse(is.na(.x),0,.x))) |> 
#'     mutate(source = "observed") |> 
#'     ungroup()
#'   
#'   
#'   combined_output <- bind_rows(output_observed, output_null)
#'   
#'   # Calculate SES results
#'   SES_coef <- combined_output %>%
#'     group_by(direction) |> 
#'     summarise(magnitude_ses = list(estimate_ses(magnitude[source == "observed"], magnitude[source != "observed"], param = param, p = TRUE))) |> 
#'     unnest_wider(magnitude_ses) |> 
#'     set_names(c("direction","magnitude_SES","magnitude_pvalue"))
#'   
#'   return(SES_coef)
#'   
#' }
#' 
#' 
#' 
#' #' Estimate SES for Shape
#' #'
#' #' Computes SES for the probability of different shapes between dissimilarities and human footprint effects.
#' #'
#' #' @param observed_model Observed model data.
#' #' @param null_model Null model data.
#' #' @param param Logical, if `TRUE`, parametric estimation is used. Defaults to `TRUE`.
#' #' @return Data frame containing SES results for shape and p-values.
#' #' @examples
#' #' observed_model <- list(...)
#' #' null_model <- list(...)
#' #' estimate_SES_shape(observed_model, null_model)
#' estimate_SES_shape <-  function(observed_model=NULL, null_model=NULL,param=TRUE) {
#'   require(tidyverse)
#'   
#'   output_null <- null_model |> 
#'     enframe(name = "source", value = "iteration_value") |> 
#'     mutate(iteration_value = map(iteration_value, ~ enframe(.x, name = "direction", value = "direction_value")),.keep="unused") |> 
#'     unnest(iteration_value, keep_empty = TRUE) |> 
#'     mutate(coef = map(direction_value, ~ .x |>  
#'                         pluck("coefs") |>  
#'                         bind_rows() |>  
#'                         rownames_to_column("predictor") |> 
#'                         filter(grepl("^hfp_", predictor)))) |> unnest(coef) |> 
#'     dplyr::mutate(shape = dplyr::case_when(coefficient.1 + coefficient.2 + coefficient.3 == 0 ~ "Absent",
#'                                            coefficient.1 <= coefficient.2 & coefficient.2 < coefficient.3 ~ "Exponential",
#'                                            coefficient.1 > coefficient.2 & coefficient.2 > coefficient.3 ~ "Saturating",
#'                                            coefficient.1 < coefficient.2 & coefficient.2 > coefficient.3 ~ "Saturating",
#'                                            coefficient.1 < coefficient.2 & coefficient.2 < coefficient.3 ~ "Saturating",
#'                                            coefficient.1 > coefficient.2 & coefficient.2 == coefficient.3 ~ "Saturating",
#'                                            coefficient.1 > coefficient.2 & coefficient.2 < coefficient.3 & coefficient.3 > 0 ~ "Revlog"))  |> 
#'     group_by(source, direction) |> 
#'     mutate(shape = factor(shape, levels =c("Absent", "Saturating","Exponential","Revlog"))) |> 
#'     mutate(nrow = row_number()) |> 
#'     count(shape) |>  
#'     complete(shape) |> 
#'     mutate(n = ifelse(is.na(n),0,n)) |> 
#'     pivot_wider(names_from = shape, values_from = n, values_fill = 0) |> 
#'     mutate(across(Absent:Revlog, ~ .x / (Absent + Saturating + Exponential + Revlog))) 
#'   # Process observed data
#'   output_observed <- observed_model |> 
#'     pluck("Functional") |> 
#'     enframe(name = "direction", value = "direction_value") |> 
#'     mutate(coef = map(direction_value, ~ .x |>  
#'                         pluck("coefs") |>  
#'                         bind_rows() |>  
#'                         rownames_to_column("predictor") |> 
#'                         filter(grepl("^hfp_", predictor))),.keep="unused") |> unnest(coef) |> 
#'     mutate(direction = factor(direction, levels=c("differentiation","homogenization"))) |> 
#'     complete(direction) |> 
#'     mutate(across(c(coefficient.1,coefficient.2,coefficient.3),~ ifelse(is.na(.x),0,.x))) |> 
#'     dplyr::mutate(shape = dplyr::case_when(coefficient.1 + coefficient.2 + coefficient.3 == 0 ~ "Absent",
#'                                            coefficient.1 <= coefficient.2 & coefficient.2 < coefficient.3 ~ "Exponential",
#'                                            coefficient.1 > coefficient.2 & coefficient.2 > coefficient.3 ~ "Saturating",
#'                                            coefficient.1 < coefficient.2 & coefficient.2 > coefficient.3 ~ "Saturating",
#'                                            coefficient.1 < coefficient.2 & coefficient.2 < coefficient.3 ~ "Saturating",
#'                                            coefficient.1 > coefficient.2 & coefficient.2 == coefficient.3 ~ "Saturating",
#'                                            coefficient.1 > coefficient.2 & coefficient.2 < coefficient.3 & coefficient.3 > 0 ~ "Revlog"))  |> 
#'     
#'     
#'     group_by(direction) |> 
#'     mutate(shape = factor(shape, levels =c("Absent", "Saturating","Exponential","Revlog"))) |> 
#'     count(shape) |>  
#'     complete(shape) |> 
#'     mutate(n = ifelse(is.na(n),0,n)) |> 
#'     pivot_wider(names_from = shape, values_from = n, values_fill = 0) |> 
#'     mutate(across(Absent:Revlog, ~ .x / (Absent + Saturating + Exponential + Revlog))) |> 
#'     mutate(source = "observed") 
#'   
#'   # Combine null and observed outputs
#'   combined_output <- bind_rows(output_observed, output_null)
#'   
#'   # Calculate SES results
#'   SES_shape <- combined_output %>%
#'     group_by(direction) |> 
#'     summarise(across(c(Absent, Exponential, Saturating, Revlog), ~ list(estimate_ses(.x[source == "observed"], .x[source != "observed"], param = param, p = TRUE)))) |> 
#'     unnest_wider(c(Absent, Exponential, Saturating, Revlog), names_sep = "_")
#'   
#'   return(SES_shape)
#' }
#' 
#' 
#' 
