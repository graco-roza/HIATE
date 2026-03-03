# Install if needed
if (!require(effsize)) install.packages("effsize")
if (!require(DescTools)) install.packages("DescTools")
if (!require(WRS2)) install.packages("WRS2")
if (!require(Hmisc)) install.packages("Hmisc")

library(effsize)
library(DescTools) # For HodgesLehmann estimator
library(WRS2)      # For robust D / Yuen
library(Hmisc)


#' Robust Harrell–Davis Standardized Effect Size (Self-contained)
#'
#' Computes a robust standardized effect size between two numeric vectors,
#' using Harrell–Davis quantile estimators for the median and a pooled
#' robust MAD. Suitable for highly skewed, non-normal, or heteroscedastic data.
#'
#' @param x Numeric vector for homogenisation.
#' @param y Numeric vector for differentiation.
#' @param n Optional. Total number of samples expected (default is 1000).
#' @param prob Quantile probability to use (default is 0.5, i.e., median).
#' @return Numeric: the robust standardized effect size.
#' @export
harrel_davis <- function(x, y, n = 1000, prob = 0.5) {
  if (!requireNamespace("Hmisc", quietly = TRUE)) {
    stop("Please install the Hmisc package for hdquantile.")
  }
  hdmedian <- function(z) as.numeric(Hmisc::hdquantile(z, 0.5))
  hdmad <- function(z) {
    med <- hdmedian(z)
    1.4826 * hdmedian(abs(z - med))
  }
  pooled <- function(a, b, FUN) {
    na <- length(a)
    nb <- length(b)
    sqrt(((na - 1) * FUN(a)^2 + (nb - 1) * FUN(b)^2) / (na + nb - 2))
  }
  
  ok_x <- !is.na(x)
  ok_y <- !is.na(y)
  n_x  <- sum(ok_x)
  n_y  <- sum(ok_y)
  prop_x <- n_x / n
  prop_y <- n_y / n
  if (n_x == 0 && n_y == 0) return(0)
  med_x <- if (n_x == 1) x[ok_x] else if (n_x > 1) Hmisc::hdquantile(x[ok_x], prob) else 0
  med_y <- if (n_y == 1) y[ok_y] else if (n_y > 1) Hmisc::hdquantile(y[ok_y], prob) else 0
  adj_x <- med_x * prop_x
  adj_y <- med_y * prop_y
  pooled_mad <- {
    x0 <- x[ok_x]
    y0 <- y[ok_y]
    if (length(x0) < 2 && length(y0) < 2) NA_real_
    else if (length(x0) < 2) mad(y0)
    else if (length(y0) < 2) mad(x0)
    else pooled(x0, y0, hdmad)
  }
  if (is.na(pooled_mad) || pooled_mad == 0) {
    raw_diff <- adj_x - adj_y
    if (raw_diff > 0) return(1)
    if (raw_diff < 0) return(-1)
    return(0)
  }
  as.numeric((adj_x - adj_y) / pooled_mad)
}

median_overlap <- function(homog, diffe, dens_n = 512, n = 1000) {

  
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

hl <- function(x, y = NULL) {
  if (is.null(y)) {
    walsh <- outer(x, x, "+") / 2
    median(walsh[lower.tri(walsh, diag = TRUE)])
  } else {
    median(outer(x, y, "-"),na.rm=TRUE)
  }
}

#' Robust Cohen's d (median/MAD)
#'
#' Computes a robust version of Cohen's d using the difference of medians
#' divided by the mean of the MADs for two independent samples.
#'
#' @param x Numeric vector (group 1)
#' @param y Numeric vector (group 2)
#' @return Numeric value: robust Cohen's d
#' @export
cohen_d_robust <- function(x, y) {
  mdiff <- (median(y, na.rm = TRUE)*(length(na.exclude(y))/1000)) - (median(x, na.rm = TRUE)*(length(na.exclude(y))/1000)) 
  mad_x <- mad(x,  na.rm = TRUE)  # classic MAD (median absolute deviation from median)
  mad_y <- mad(y,  na.rm = TRUE)
  denom <- mean(c(mad_x, mad_y),na.rm=TRUE)
  # Avoid division by zero
  if (denom == 0) {
    if (mdiff > 0) return(Inf)
    if (mdiff < 0) return(-Inf)
    return(0)
  }
  mdiff / denom
}



# set.seed(123)
# 
# # Scenario 1: Minimal separation, unequal sizes
# setA1 <- runif(999, 0.1, 0.2)    # 10 values
# setB1 <- runif(999, 0.3, 0.4)    # 30 values
# 
# # Scenario 2: Maximal separation, unequal sizes
# setA2 <- runif(999, 0.1, 0.2)    # 10 values
# setB2 <- runif(999, 0.8, 0.9)    # 30 values
# 
# hl(y=setA1,x=setB1) #0.1995391
# hl(y=setA2,x=setB2) #0.6989917
# 
# cliff.delta(setB1,setA1)$estimate # 1
# cliff.delta(setB2,setA2)$estimate # 1
# 
# cohen_d_robust(setA1,setB1) #7.965724
# cohen_d_robust(setA2,setB2) #27.29989
# 
# harrel_davis(y=setA1,x=setB1) #2.687952
# harrel_davis(y=setA2,x=setB2) #9.218537
# 
# 
# # Scenario 1: Minimal separation, unequal sizes
# set2A1 <- runif(999, 0.1, 0.2)    # 10 values
# set2B1 <- runif(200, 0.1, 0.4)    # 30 values
# 
# # Scenario 2: Maximal separation, unequal sizes
# set2A2 <- runif(999, 0.1, 0.2)    # 10 values
# set2B2 <- runif(200, 0.8, 0.9)    # 30 values
# 
# hl(y=set2A1,x=set2B1) #0.1987752
# hl(y=set2A2,x=set2B2) #0.6987098
# 
# cliff.delta(set2B1,set2A1) # 1
# cliff.delta(set2B2,set2A2) # 1
# 
# cohen_d_robust(set2A1,set2B1) #8.305365
# cohen_d_robust(set2A2,set2B2) #27.60784
# 
# harrel_davis(y=set2A1,x=set2B1) #-2.337019
# harrel_davis(y=set2A2,x=set2B2) #0.5182828

#1. Discover all BBGDM output files and strip extensions for dataset names
files_all <- list.files(
  "S4_run_BBGDM/bbgdm_output", 
  recursive = TRUE,
  full.names = TRUE
)
files <- tools::file_path_sans_ext(files_all)

# 2. Read predictor metadata for each processed dataset
dataset_features <- list.files("S1_Preprocessing/Processed") %>%
  purrr::map(~ here("S1_Preprocessing/Processed", .x) %>%
               read_rds() %>%
               pluck("predictors")) %>%
  set_names(tools::file_path_sans_ext(list.files("S1_Preprocessing/Processed"))) %>%
  bind_rows(.id = "dataset") %>%
  mutate(dataset = str_remove(dataset, "_clean")) %>%
  left_join(
    read_xlsx("S1_Preprocessing/Miscellaneous/dataset_info_all.xlsx"),
    by = c("dataset" = "dataset_name")
  )

# 3. Extract dataset identifiers from file paths
basic_structure <- tibble(fullpath = files) %>%
  mutate(
    dataset     = basename(fullpath),                # raw dataset name
    beta_folder = basename(dirname(fullpath)),       # e.g. "Baselga_abun"
    beta_type   = str_extract(beta_folder, "^(Baselga|Podani)"),
    metric_type = str_extract(beta_folder, "(abun|pa)$")
  ) %>%
  select(dataset, beta_type, metric_type)

# --- Helper function: read BBGDM outputs and unwrap list structure ------------
read_bbgdm <- function(files, .pattern) {
  # .pattern = "r2" or "coefs" or etc.
  purrr::map(files, ~ read_rds(glue::glue("{.x}.rds"))) %>%
    set_names(gsub("_bbgdm$", "", files)) %>%
    enframe(name = "dataset", value = "data_list") %>%
    mutate(
      data_list = purrr::map(data_list, ~ enframe(.x, name = "facet", value = "facet_list")),
      facet_list = purrr::map(data_list, ~ unnest(.x, cols = "facet_list")),
      dir_list   = purrr::map(facet_list, ~ unnest(enframe(.x$facet_list, name="direction", value="direction_list")))
    ) %>%
    unnest(dir_list)
}

# 4A. Compute R² metrics per dataset/facet/direction ---------------------------
mean_r2 <- purrr::map(files, ~ readr::read_rds(glue::glue("{.x}.rds"))) |> # Read BBGDM results for each file
  set_names(gsub("_bbgdm", "", files)) |> 
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
  filter(grepl("Podani_abun", dataset)) |> 
  group_by(dataset,facet) |> 
  summarise(
    cliff_delta = tryCatch(cliff.delta(homogenization,differentiation)$estimate, error = function(e) NA),
    cohen_d = tryCatch(cohen_d_robust(differentiation,homogenization), error = function(e) NA),
    hodges_lehmann = hl(homogenization,differentiation),
    harrel_davis =  harrel_davis(homogenization, differentiation),
    median_overlap     = median_overlap(homogenization, differentiation)$hybrid,
    median_diff =  median_overlap(homogenization, differentiation)$r2_diff
  )


cliff_d<- mean_r2 |> 
 ggplot(aes(x=harrel_davis, cliff_delta)) + 
  geom_point()+
  geom_smooth(method="loess")+
  theme_bw()

c_d<- mean_r2 |> 
  ggplot(aes(x=harrel_davis, cohen_d)) + 
  geom_point()+
  geom_smooth(method="loess")+
  theme_bw()

hodges_l<- mean_r2 |> 
  ggplot(aes(x=harrel_davis, hodges_lehmann)) + 
  geom_point()+
  geom_smooth(method="loess")+
  theme_bw()

med_diff<-mean_r2 |> 
  ggplot(aes(x=harrel_davis, median_diff)) + 
  geom_point()+
  geom_smooth(method="loess")+
  theme_bw()

med_ovl<-mean_r2 |> 
  ggplot(aes(x=harrel_davis, median_overlap)) + 
  geom_point()+
  geom_smooth(method="loess")+
  theme_bw()

library(patchwork)

cliff_d+c_d+hodges_l+med_diff+med_ovl

direction_sensitivity <- mean_r2 |> 
 mutate(across(cliff_delta:median_diff, ~ ifelse(.x >0,"homogenisation","differentiation")))


 direction_sensitivity %>%
    filter(facet == "Functional") |> 
  janitor::tabyl(cohen_d)
  rename(Harrell_Davis = harrel_davis, Cohen_D = cohen_d, Freq = n)
which(direction_sensitivity$cohen_d != direction_sensitivity$harrel_davis)
mean_r2$dataset[12]


mean_r2
# Print table (optional)
print(conf_mat)

# Plot as a heatmap
ggplot(conf_mat, aes(x = Harrell_Davis, y = Cohen_D, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Freq), color = "black", size = 5) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_minimal() +
  labs(title = "Confusion Matrix: Harrell-Davis vs Cohen's d",
       x = "Harrell-Davis Classification",
       y = "Cohen's d Classification",
       fill = "Count")


test<-purrr::map(files, ~ readr::read_rds(glue::glue("{.x}.rds"))) |> # Read BBGDM results for each file
  set_names(gsub("_bbgdm", "", files)) |> 
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
  filter(grepl("Podani_abun", dataset)) |> 
  group_by(dataset,facet) 


foo<-test |> 
  filter(dataset == "S4_run_BBGDM/bbgdm_output/Podani_abun/N29FMI", facet== "Taxonomic") |> 
  pull(homogenization) |> na.exclude()


foo2<-test |> 
  filter(dataset == "S4_run_BBGDM/bbgdm_output/Podani_abun/N29FMI", facet== "Taxonomic") |> 
  pull(differentiation) |> na.exclude() 



par(mfrow=c(2,1))
plot(density(foo), col="red", lwd=2, xlim=c(0,1),, main="homogenization [red] / differentiation [blue]") 
lines(density(foo2), col="blue",lwd=2,lty=2)

plot(density(foo*(length(foo)/1000)), col="red", lwd=2, xlim=c(0,1), ylim=c(0,60),main="Weighted by bootstrap success") 
lines(density(foo2*(length(foo2)/1000)), col="blue",lwd=2,lty=2)


foo<-test |> 
  filter(dataset == "S4_run_BBGDM/bbgdm_output/Podani_abun/N29FMI", facet== "Functional") |> 
  pull(homogenization) |> na.exclude()


foo2<-test |> 
  filter(dataset == "S4_run_BBGDM/bbgdm_output/Podani_abun/N29FMI", facet== "Functional") |> 
  pull(differentiation) |> na.exclude() 


par(mfrow=c(2,1))
plot(density(foo), col="red", lwd=2, xlim=c(0,1),, main="homogenization [red] / differentiation [blue]") 
lines(density(foo2), col="blue",lwd=2,lty=2)

plot(density(foo*(length(foo)/1000)), col="red", lwd=2, xlim=c(0,1), main="Weighted by bootstrap success") 
lines(density(c(foo2)*(length(foo2)/1000)), col="blue",lwd=2,lty=2)
