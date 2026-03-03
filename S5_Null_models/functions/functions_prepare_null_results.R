require(Hmisc)
require(gdm)
pooled <- function(x, y, FUN) {
  nx <- length(x)
  ny <- length(y)
  sqrt(((nx - 1) * FUN(x) ^ 2 + (ny - 1) * FUN(y) ^ 2) / (nx + ny - 2))
}
hdmedian <- function(x) as.numeric(hdquantile(x, 0.5))
hdmad <- function(x) 1.4826 * hdmedian(abs(x - hdmedian(x)))
phdmad <- function(x, y) pooled(x, y, hdmad)
median_overlap <- function(homog, diffe, dens_n = 512) {
  n <- 1000
  if (length(homog) != n || length(diffe) != n) 
    stop("homog and diffe must both have length 1000")
  
  ok_h <- !is.na(homog)
  ok_d <- !is.na(diffe)
  n_h <- sum(ok_h)
  n_d <- sum(ok_d)
  prop_h <- n_h / n
  prop_d <- n_d / n
  
  med_h <- if (n_h > 0) median(homog[ok_h]) else 0
  med_d <- if (n_d > 0) median(diffe[ok_d]) else 0
  
  if (n_h == 0 && n_d == 0) {
    return(list(
      med_h = 0, med_d = 0,
      prop_h = 0, prop_d = 0,
      overlap = 1,
      hybrid = 0
    ))
  }
  
  # Continuous overlap only if both have >= 2 successful draws
  if (n_h >= 2 && n_d >= 2) {
    rng <- range(c(homog[ok_h], diffe[ok_d]))
    dx <- diff(rng) / (dens_n - 1)
    fx <- stats::density(homog[ok_h], from = rng[1], to = rng[2], n = dens_n)$y
    fy <- stats::density(diffe[ok_d], from = rng[1], to = rng[2], n = dens_n)$y
    cont_ovl <- sum(pmin(fx, fy)) * dx
  } else {
    cont_ovl <- 0
  }
  
  # NA-overlap = shared failure rate
  fail_overlap <- min(1 - prop_h, 1 - prop_d)
  overlap <- cont_ovl + fail_overlap
  overlap <- min(max(overlap, 0), 1)
  
  corr_h <- med_h * prop_h
  corr_d <- med_d * prop_d
  
  hybrid <- (corr_h - corr_d) * (1 - overlap)
  
  list(
    r2_diff = (corr_h - corr_d),
    hybrid  = hybrid
  )
}

harrel_davis <- function(homog, diffe, prob = 0.5, n = 1000) {
  ok_h <- !is.na(homog)
  ok_d <- !is.na(diffe)
  n_h  <- sum(ok_h)
  n_d  <- sum(ok_d)
  
  prop_h <- n_h / n
  prop_d <- n_d / n
  
  if (n_h == 0 && n_d == 0) return(0)
  
  # Use fallback: if only one value, take it; else try hdquantile
  med_h <- if (n_h == 1) homog[ok_h] else if (n_h > 1) hdquantile(homog[ok_h], prob) else 0
  med_d <- if (n_d == 1) diffe[ok_d] else if (n_d > 1) hdquantile(diffe[ok_d], prob) else 0
  
  adj_h <- med_h * prop_h
  adj_d <- med_d * prop_d
  
  # Safe pooled MAD
  safe_hd_mad <- function(x, y) {
    x0 <- na.omit(x)
    y0 <- na.omit(y)
    
    if (length(x0) < 2 && length(y0) < 2) return(NA_real_)
    if (length(x0) < 2) return(mad(y0))
    if (length(y0) < 2) return(mad(x0))
    
    phdmad(x0, y0)
  }
  
  pooled_mad <- safe_hd_mad(homog[ok_h], diffe[ok_d])
  
  # Handle edge cases
  if (is.na(pooled_mad) || pooled_mad == 0) {
    raw_diff <- adj_h - adj_d
    if (raw_diff > 0) return(1)
    if (raw_diff < 0) return(-1)
    return(0)
  }
  
  return(as.numeric((adj_h - adj_d) / pooled_mad))
}
#’ Build a single I‐spline curve for one predictor (BBGDM)
#’
#’ Calls the internal GDM C routine to get the fitted spline values for
#’ a three-knot I-spline, then returns the (x,y) curve.
#’
#’ @param k1,k2,k3 Numeric. The three knot positions (min, median, max).
#’ @param c1,c2,c3 Numeric. The corresponding I-spline coefficients.
#’ @param PSAMPLE Integer. Number of equally-spaced x points (default 200).
#’
#’ @return A tibble with columns:
#’   \item{x}{Predictor values, from k1→k3.}
#’   \item{y}{Cumulative turnover (I-spline sum at each x).}
#’
#’ @examples
#’ \dontrun{
#’   curve <- ispline_curve_row(0, 10, 20, 0.3, 0.5, 0.2)
#’   plot(curve$x, curve$y, type="l")
#’ }
#’ @export
ispline_curve_row <- function(k1, k2, k3,
                              c1, c2, c3,
                              PSAMPLE = 200) {
  preddata <- numeric(PSAMPLE)
  out <- .C("GetPredictorPlotData",
            pdata      = as.double(preddata),
            as.integer(PSAMPLE),
            as.double(c(c1, c2, c3)),
            as.double(c(k1, k2, k3)),
            as.integer(3),
            PACKAGE    = "gdm")
  x <- seq(k1, k3, length.out = PSAMPLE)
  y <- out$pdata
  tibble::tibble(x = x, y = y)
}


#’ Classify an I-spline curve’s shape via two-segment slope logic
#’
#’ Splits the first derivative into three equal sections, computes their
#’ means, then exhaustively maps the pair of adjacent‐leg comparisons
#’ (down, flat, up) to one of five shapes.
#’
#’ @param c1,c2,c3 Numeric coefficients for the I-spline basis functions.
#’ @param k1,k2,k3 Numeric knot positions (min, median, max) for the predictor.
#’ @param PSAMPLE Integer grid resolution (default 200).
#’ @param rel_tol Numeric ∈ [0,1], relative tolerance on slope differences.
#’ @param tail_cut Numeric ∈ (0,1], fraction of derivative to keep (1 = no trim).
#’
#’ @return A character scalar:  
#’   “Absent”, “Linear”, “Exponential”, “Saturating”, “Revlog”, or “Uncertain”.
#’
#’ @examples
#’ \dontrun{
#’   classify_shapes(0.3, 0.5, 0.2, 0, 10, 20)
#’ }
#’ @export
classify_shapes <- function(c1, c2, c3,
                            k1, k2, k3,
                            PSAMPLE  = 200,
                            rel_tol  = 0.05,
                            tail_cut = 1) {
  if (c1 == 0 && c2 == 0 && c3 == 0) return("Absent")
  pd <- ispline_curve_row(k1, k2, k3, c1, c2, c3, PSAMPLE)
  d  <- diff(pd$y) / diff(pd$x)
  L  <- length(d)
  if (tail_cut < 1) {
    drop <- floor(((1 - tail_cut) / 2) * L)
    d    <- d[(drop + 1):(L - drop)]
    L    <- length(d)
  }
  t3  <- floor(L / 3)
  s1  <- mean(d[       1:t3       ])
  s2  <- mean(d[(t3+1):(2*t3)     ])
  s3  <- mean(d[(2*t3+1):     L   ])
  tol <- rel_tol * max(abs(s1), abs(s2), abs(s3))
  
  eq12 <- abs(s1 - s2) <= tol
  eq23 <- abs(s2 - s3) <= tol
  lt12 <- (s2 - s1)  > tol
  gt12 <- (s1 - s2)  > tol
  lt23 <- (s3 - s2)  > tol
  gt23 <- (s2 - s3)  > tol
  
  shape <- dplyr::case_when(
    eq12 & eq23             ~ "Linear",
    lt12 & lt23             ~ "Exponential",
    gt12 & gt23             ~ "Saturating",
    gt12 & lt23             ~ "Revlog",
    eq12 & lt23             ~ "Exponential",
    lt12 & eq23             ~ "Exponential",
    eq12 & gt23             ~ "Saturating",
    gt12 & eq23             ~ "Saturating",
    lt12 & gt23             ~ "Saturating",
    TRUE                     ~ "Uncertain"
  )
  shape
}


#’ Helper: choose knots based on predictor name
#’
#’ @param predictor  Character. “hfp…” or “het…” prefix.
#’ @param hfp_min    Numeric. Minimum of HFP range.
#’ @param hfp_med    Numeric. Median of HFP range.
#’ @param hfp_max    Numeric. Maximum of HFP range.
#’ @param het_min    Numeric. Minimum of HET range.
#’ @param het_med    Numeric. Median of HET range.
#’ @param het_max    Numeric. Maximum of HET range.
#’ @return Numeric vector c(k1,k2,k3).
#’ @noRd
get_knots <- function(predictor,
                      hfp_min, hfp_med, hfp_max,
                      het_min, het_med, het_max) {
  if (stringr::str_detect(unique(predictor), "^het")) {
    c(het_min, het_med, het_max)
  } else {
    c(hfp_min, hfp_med, hfp_max)
  }
}

#’ Classify one bootstrap by predictor (hfp or het)
#’
#’ @inheritParams classify_shapes
#’ @param predictor   Character. “hfp_…” or “het_…”.
#’ @param hfp_min, hfp_median, hfp_max  Numeric HFP knots.
#’ @param het_min, het_median, het_max  Numeric HET knots.
#’ @return Character. One of “Absent”, “Linear”, “Exponential”, “Saturating”, “Revlog”, “Uncertain”.
#’ @export
classify_predictor_row <- function(predictor,
                                   c1, c2, c3,
                                   hfp_min, hfp_median, hfp_max,
                                   het_min, het_median, het_max,
                                   PSAMPLE  = 200,
                                   rel_tol  = 0.05,
                                   tail_cut = 1) {
  knots <- get_knots(predictor,
                     hfp_min, hfp_median, hfp_max,
                     het_min, het_median, het_max)
  classify_shapes(c1, c2, c3,
                  knots[1], knots[2], knots[3],
                  PSAMPLE  = PSAMPLE,
                  rel_tol  = rel_tol,
                  tail_cut = tail_cut)
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

# estimate_ses_weighted <- function(obs, est, weights, param = TRUE, p = TRUE) {
#   if (param) {
#     weighted_mean = sum(est * weights) / sum(weights)
#     weighted_sd = sqrt(sum((est - weighted_mean)^2 * weights) / sum(weights))
#     if (is.na(weighted_sd)) weighted_sd <- 0
#     res = (obs - weighted_mean) / weighted_sd
#     res = ifelse(is.nan(res),0,res)
#     pval = pnorm(res)
#   } else {
#     est = c(obs, est)
#     pval = (sum(est < obs) + sum(est == obs)/2) / length(est)
#     res = qnorm(pval)
#   }
#   if (p) {
#     pval = pval * 2
#     if (pval > 1) 
#       pval = 1 - (pval - 1)
#     res = c(res, pval)
#     names(res) = c("SES", "pvalue")
#   }
#   return(res)
# }


# Function to calculate SES for observed and null models
# Arguments:
#   null_model: Null model data
#   observed_model: Observed model data
#   param: Boolean indicating whether to use parametric estimation (default: TRUE)
#   p: Boolean indicating whether to calculate p-values (default: TRUE)
# Returns:
#   SES results including SES, p-values, and other relevant information
estimate_SES_direction <-  function(observed_model=NULL, null_model=NULL,param=FALSE) {
  
  # Calculate null model output
output_null <- null_model %>%
map_depth(2, ~ .x |>  pluck("r2")) |>
map_df(~ tibble(
ID = seq_along(.x$differentiation),
differentiation = .x$differentiation,
homogenisation =  .x$homogenization
), .id = "source") |> 
  mutate(across(differentiation:homogenisation, ~ifelse(.x<=0,NA,.x/100))) |> 
  #  pivot_longer(cols=!c(source,ID), names_to = "direction", values_to = "r2") |> 
  #  dplyr::mutate(across(c(r2), ~ ifelse(.x <= 0, NA, .x/100))) |>
    group_by(source) |> 
    summarise(direction_r2 = harrel_davis(homogenisation,differentiation))


# Calculate observed model output
  output_observed <- observed_model |> 
    pluck("Functional") |> 
    map_depth(1, ~ .x |>  pluck("r2"))  |> 
    bind_rows(.id = "direction") |> 
    mutate(differentiation = differentiation,
           homogenisation = homogenization,
           ID = seq_along(differentiation),
           source = "observed", .keep = "unused") |> 
    mutate(across(differentiation:homogenisation, ~ifelse(.x<=0,NA,.x/100))) |> 
    group_by(source) |> 
    summarise(direction_r2 = harrel_davis(homogenisation,differentiation))
    
# Combine observed and null model output
  combined_output <- bind_rows(output_observed, output_null)

  # Calculate SES results
  
  SES_r2 <- combined_output %>%
    summarise(
      res = list(
        tryCatch(
          estimate_ses(
            direction_r2[source == "observed"],
            direction_r2[source != "observed"],
            param = param,
            p = TRUE
          ),
          error = function(e) list(ses = NA, pval = NA)
        )
      )
    ) %>%
    unnest_wider(res) %>%
    rename(direction_SES = ses, direction_pvalue = pvalue)

    # SES_r2 <- combined_output %>%
    # summarise(ses_results_diff = list(estimate_ses(direction_r2[source == "observed"], direction_r2[source != "observed"], param = param, p = TRUE))) |> 
    # unnest_wider(ses_results_diff) |> 
    # set_names(c("direction_SES","direction_pvalue"))

  return(SES_r2)
}

# Function to estimate the Standardized Effect Size (SES) of the strength of human footprint effect on traits
estimate_SES_magnitude <- function(observed_model = NULL, null_model = NULL, param = FALSE) {
  require(tidyverse)
  
  # --- helpers ---------------------------------------------------------------
  normalize_direction <- function(x) {
    x <- tolower(as.character(x))
    x <- gsub("z", "s", x) # homogenization -> homogenisation
    case_when(
      grepl("homogenis", x) ~ "homogenisation",
      grepl("differen",  x) ~ "differentiation",
      TRUE ~ x
    )
  }
  
  sum_coeffs <- function(df) {
    # robustly sum all columns named coefficient.*
    coeff_cols <- grep("^coefficient\\.", names(df), value = TRUE)
    if (length(coeff_cols) == 0) return(NA_real_)
    rowSums(df[, coeff_cols, drop = FALSE], na.rm = TRUE)
  }
  
  safe_estimate_ses <- function(obs, nulls, param = TRUE) {
    obs  <- obs[!is.na(obs)]
    nulls <- nulls[!is.na(nulls)]
    if (length(obs) == 0 || length(nulls) == 0) {
      return(tibble(ses = NA_real_, pvalue = NA_real_))
    }
    res <- tryCatch(estimate_ses(obs, nulls, param = param, p = TRUE),
                    error = function(e) NULL)
    if (is.null(res)) return(tibble(ses = NA_real_, pvalue = NA_real_))
    # standardise possible name variants from estimate_ses()
    get_first <- function(lst, keys) {
      for (k in keys) if (!is.null(lst[[k]])) return(as.numeric(lst[[k]]))
      if (is.numeric(lst) && length(lst) >= 1) return(as.numeric(lst[1]))
      return(NA_real_)
    }
    tibble(
      ses    = if (is.list(res)) get_first(res, c("SES","ses","z","stat")) else as.numeric(res[1]),
      pvalue = if (is.list(res)) get_first(res, c("pval","pvalue","p"))  else if (length(res) >= 2) as.numeric(res[2]) else NA_real_
    )
  }
  
  # --- extract observed magnitudes ------------------------------------------
  obs_df <- observed_model %>%
    pluck("Functional") %>%
    enframe(name = "direction", value = "obj") %>%
    mutate(direction = normalize_direction(direction)) %>%
    mutate(coefs = map(obj, ~ .x %>% pluck("coefs") %>% bind_rows() %>%
                         rownames_to_column("predictor"))) %>%
    select(direction, coefs) %>% unnest(coefs) %>%
    mutate(gradient = case_when(
      str_starts(predictor, "het_") ~ "het",
      str_starts(predictor, "hfp_") ~ "hfp",
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(gradient)) %>%
    mutate(coef_sum = sum_coeffs(cur_data_all())) %>%
    group_by(direction, gradient) %>%
    summarise(magnitude = mean(coef_sum), .groups = "drop") %>%
    complete(direction = c("differentiation","homogenisation"),
             gradient  = c("hfp","het"),
             fill = list(magnitude = NA_real_)) %>%
    mutate(source = "observed")
  
  # --- extract null magnitudes ----------------------------------------------
  null_df <- enframe(null_model, name = "iter", value = "val") %>%
    mutate(val = map(val, ~ enframe(.x, name = "direction", value = "obj"))) %>%
    unnest(val, keep_empty = TRUE) %>%
    mutate(direction = normalize_direction(direction)) %>%
    mutate(coefs = map(obj, ~ .x %>% pluck("coefs") %>% bind_rows() %>%
                         rownames_to_column("predictor"))) %>%
    select(iter, direction, coefs) %>% unnest(coefs) %>%
    mutate(gradient = case_when(
      str_starts(predictor, "het_") ~ "het",
      str_starts(predictor, "hfp_") ~ "hfp",
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(gradient)) %>%
    mutate(coef_sum = sum_coeffs(cur_data_all())) %>%
    group_by(iter, direction, gradient) %>%
    summarise(magnitude = mean(coef_sum), .groups = "drop") %>%
    complete(iter, direction = c("differentiation","homogenisation"),
             gradient = c("hfp","het"),
             fill = list(magnitude = NA_real_)) %>%
    mutate(source = "null")
  
  combined <- bind_rows(obs_df, null_df)

  # --- SES per direction × gradient -----------------------------------------
  SES_tbl <- combined %>%
    group_by(direction, gradient) %>%
    summarise(
      res = list(
        safe_estimate_ses(
          obs   = magnitude[source == "observed"],
          nulls = magnitude[source == "null"],
          param = param
        )
      ),
      .groups = "drop"
    ) %>%
    unnest(res)
  
  # --- pivot wider to requested columns --------------------------------------
  out <- SES_tbl %>%
    mutate(direction = as.character(direction)) %>%
    pivot_wider(
      id_cols = direction,
      names_from = gradient,
      values_from = c(ses, pvalue),
      names_glue = "magnitude_{gradient}_{.value}"
    ) %>%
    # ensure all expected columns exist
    complete(direction = c("differentiation","homogenisation")) %>%
    mutate(
      magnitude_hfp_ses     = as.numeric(magnitude_hfp_ses),
      magnitude_hfp_pvalue  = as.numeric(magnitude_hfp_pvalue),
      magnitude_het_ses     = as.numeric(magnitude_het_ses),
      magnitude_het_pvalue  = as.numeric(magnitude_het_pvalue)
    )
  
  return(out)
}


# Function to estimate the SES of the probability of a given shape between dissimilarities and human footprint
estimate_SES_shape <-  function(observed_model=NULL, null_model=NULL,param=FALSE) {
  require(tidyverse)

 output_null <-
    null_model |> 
    enframe(name = "source", value = "iteration_value") |> 
    mutate(iteration_value = map(iteration_value, ~ enframe(.x, name = "direction", value = "direction_value")),.keep="unused") |> 
    unnest(iteration_value, keep_empty = TRUE) |> 
    mutate(coef = map(direction_value, ~ .x |>  
                        pluck("coefs") |>  
                        bind_rows() |>  
                        rownames_to_column("predictor") |> 
                        filter(grepl("^hfp_|^het_", predictor)))) |> 
    select(-direction_value) |> 
    unnest(coef) |> 
    mutate(dataset = filename) |> 
    mutate(direction=gsub("z","s",direction)) |> 
    mutate(predictor = substr(predictor, 1, 3)) |> 
    left_join(metadata |>  
                select(dataset,direction,
                       hfp_min,hfp_median,hfp_max,
                       het_min,het_median,het_max)) |>
    drop_na() |> 
    group_by(source) |> 
    mutate(iter = row_number()) |>
    group_by(source,iter) |> 
    dplyr::mutate(shape=classify_predictor_row(predictor,c1=coefficient.1,c2=coefficient.2,c3=coefficient.3,
                                               hfp_min,hfp_median,hfp_max,
                                               het_min,het_median,het_max)) |> 
    group_by(source, predictor,direction) |> 
    count(shape) |> 
    complete(shape) |>      
    mutate(n = ifelse(is.na(n),0,n)) |> 
    pivot_wider(names_from = c(predictor, shape), values_from = n, values_fill = 0, names_sep = "_") |> 
      mutate(
        hfp_trials = rowSums(across(starts_with("hfp_"))),
        het_trials = rowSums(across(starts_with("het_"))))
  # Process observed data
  output_observed <- observed_model |> 
    pluck("Functional") |> 
    enframe(name = "direction", value = "direction_value") |> 
    mutate(coef = map(direction_value, ~ .x |>  
                        pluck("coefs") |>  
                        bind_rows() |>  
                        rownames_to_column("predictor") |> 
                        filter(grepl("^hfp_|^het_", predictor)))) |> 
    select(-direction_value) |> 
    unnest(coef) |> 
    mutate(dataset = filename) |> 
    mutate(direction=gsub("z","s",direction)) |> 
    mutate(predictor = substr(predictor, 1, 3)) |> 
    left_join(metadata |>  
                select(dataset,direction,
                       hfp_min,hfp_median,hfp_max,
                       het_min,het_median,het_max)) |>
    mutate(iter = row_number()) |> 
    drop_na() |> 
    group_by(iter) |> 
    dplyr::mutate(shape=classify_predictor_row(predictor,c1=coefficient.1,c2=coefficient.2,c3=coefficient.3,
                                               hfp_min,hfp_median,hfp_max,
                                               het_min,het_median,het_max)) |> 
    group_by(predictor,direction) |> 
    count(shape) |> 
    complete(shape) |> 
    mutate(n = ifelse(is.na(n),0,n)) |> 
    pivot_wider(names_from = c(predictor, shape),, values_from = n, values_fill = 0, names_sep = "_") |> 
    mutate(
      direction = first(direction),
      hfp_trials = rowSums(across(starts_with("hfp_"))),
      het_trials = rowSums(across(starts_with("het_")))) |> 
    mutate(source = "observed") 
  # Combine null and observed outputs

  het_cols <- paste0("het_", c("Absent", "Exponential", "Saturating", "Revlog"))
  hfp_cols <- paste0("hfp_", c("Absent", "Exponential", "Saturating", "Revlog"))
  
  combined_output <- bind_rows(output_observed, output_null) %>%
    # Add missing columns with 0
    add_column(!!!setNames(
      rep(list(0), sum(!(het_cols %in% names(.)))),
      het_cols[!(het_cols %in% names(.))]
    )) %>%
    add_column(!!!setNames(
      rep(list(0), sum(!(hfp_cols %in% names(.)))),
      hfp_cols[!(hfp_cols %in% names(.))]
    )) %>%
    # Replace any NA with 0
    mutate(across(everything(), ~replace_na(., 0))) %>%
    # Do division safely now that all cols exist
    mutate(across(any_of(het_cols), ~ . / het_trials),
           across(any_of(hfp_cols), ~ . / hfp_trials))
  
  # Calculate SES results
  SES_shape <- combined_output %>%
    group_by(direction) |> 
    summarise(across(any_of(het_cols), ~ list(estimate_ses(.x[source == "observed"], .x[source != "observed"], param = param, p = TRUE))),
              across(any_of(hfp_cols), ~ list(estimate_ses(.x[source == "observed"], .x[source != "observed"], param = param, p = TRUE)))) |> 
    unnest_wider(any_of(c(het_cols,hfp_cols)), names_sep = "_")
  
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
