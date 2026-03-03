###############################################################################
# 1) Join model outputs (models_collected) to your metadata (combined_data)
###############################################################################
joined_df <- right_join(
  models_collected,
  combined_data[ , c("dataset","facet","beta_type","metric_type","direction",
                     "harrel_davis",
                     "hfp_min","hfp_median","hfp_max",
                     "het_min","het_median","het_max")],
  by = c("dataset","facet","beta_type","metric_type","direction")
)

###############################################################################
# 2) Clean up the predictor column: split off any “_buffer” suffix
###############################################################################
# Extract only the part before the underscore
joined_df$predictor <- sub("^(.*)_.*$", "\\1", joined_df$predictor)

###############################################################################
# 3) Drop any rows with missing critical info
###############################################################################
# Here we drop rows where predictor or coefficients are NA
keep_rows <- complete.cases(
  joined_df[ , c("predictor","coefficient.1","coefficient.2","coefficient.3")]
)
joined_df <- joined_df[keep_rows, ]

###############################################################################
# 4) Prepare an index of each unique group we’ll loop over
###############################################################################
group_cols <- c("dataset","facet","beta_type","metric_type","direction")
group_keys <- unique(joined_df[ , group_cols])

###############################################################################
# 5) Initialize an empty list to hold each group’s results
###############################################################################
result_list <- list()
counter <- 1L

###############################################################################
# 6) Outer loop: for each unique (dataset,facet,beta_type,metric_type,direction)
###############################################################################
for (i in seq_len(nrow(group_keys))) {

  # Pull out the i-th group key
  this_key <- group_keys[i, ]
  ds       <- this_key$dataset
  fc       <- this_key$facet
  bt       <- this_key$beta_type
  mt       <- this_key$metric_type
  dir      <- this_key$direction
  
  # Subset joined_df to just this group
  sel <- with(joined_df,
              dataset == ds &
                facet == fc &
                beta_type == bt &
                metric_type == mt &
                direction == dir)
  sub_g <- joined_df[sel, ]
  
  # If no rows, skip
  if (nrow(sub_g) == 0) next
  
  #############################################################################
  # 6a) Row‐wise loop: compute shape, q_25/50/75, b1–b4 for each bootstrap row
  #############################################################################
  # Pre‐allocate new columns
  sub_g$shape <- NA_character_
  sub_g$q_25  <- NA_real_; sub_g$q_50  <- NA_real_; sub_g$q_75  <- NA_real_
  sub_g$b1    <- NA_real_; sub_g$b2    <- NA_real_; sub_g$b3    <- NA_real_; sub_g$b4    <- NA_real_
  
  for (j in seq_len(nrow(sub_g))) {
    # pull the j-th row
    rowj <- sub_g[j, ]
    
    # classify shape
    shp <- classify_predictor_row(
      rowj$predictor,
      rowj$coefficient.1, rowj$coefficient.2, rowj$coefficient.3,
      rowj$hfp_min, rowj$hfp_median, rowj$hfp_max,
      rowj$het_min, rowj$het_median, rowj$het_max
    )
    
    # compute cumulative % at quartiles
    props_vec <- pct_by_props_predictor(
      rowj$predictor,
      rowj$coefficient.1, rowj$coefficient.2, rowj$coefficient.3,
      rowj$hfp_min, rowj$hfp_median, rowj$hfp_max,
      rowj$het_min, rowj$het_median, rowj$het_max
    )
    
    # compute incremental % gains by bands
    bands_vec <- pct_by_bands_predictor(
      rowj$predictor,
      rowj$coefficient.1, rowj$coefficient.2, rowj$coefficient.3,
      rowj$hfp_min, rowj$hfp_median, rowj$hfp_max,
      rowj$het_min, rowj$het_median, rowj$het_max
    )
    
    # write back into sub_g
    sub_g$shape[j] <- shp
    sub_g$q_25[j]  <- props_vec["q_25"]
    sub_g$q_50[j]  <- props_vec["q_50"]
    sub_g$q_75[j]  <- props_vec["q_75"]
    sub_g$b1[j]    <- bands_vec["b1"]
    sub_g$b2[j]    <- bands_vec["b2"]
    sub_g$b3[j]    <- bands_vec["b3"]
    sub_g$b4[j]    <- bands_vec["b4"]
  }
  
  #############################################################################
  # 6b) Summarise over all bootstraps *by predictor* within this group
  #############################################################################
  preds <- unique(sub_g$predictor)
  
  for (pred in preds) {
    # subset to this predictor
    sel2 <- sub_g$predictor == pred
    sub2 <- sub_g[sel2, ]
    
    #  median quartiles
    q25 <- median(sub2$q_25, na.rm=TRUE)
    q50 <- median(sub2$q_50, na.rm=TRUE)
    q75 <- median(sub2$q_75, na.rm=TRUE)
    #  median band gains
    b1m <- median(sub2$b1, na.rm=TRUE)
    b2m <- median(sub2$b2, na.rm=TRUE)
    b3m <- median(sub2$b3, na.rm=TRUE)
    b4m <- median(sub2$b4, na.rm=TRUE)
    #  counts for multinomial
    nA <- sum(sub2$shape=="Absent")
    nL <- sum(sub2$shape=="Linear")
    nE <- sum(sub2$shape=="Exponential")
    nS <- sum(sub2$shape=="Saturating")
    nR <- sum(sub2$shape=="Revlog")
    nU <- sum(sub2$shape=="Uncertain")
    #  modal shape via integer index
    counts <- c(nA,nL,nE,nS,nR,nU)
    idx    <- which.max(counts)
    modal  <- switch(idx,
                     "Absent","Linear","Exponential","Saturating","Revlog","Uncertain"
    )
    
    # Build a named list for this predictor×group
    out <- list(
      dataset         = ds,
      facet           = fc,
      beta_type       = bt,
      metric_type     = mt,
      direction       = dir,
      predictor       = pred,
      q_25            = q25,
      q_50            = q50,
      q_75            = q75,
      b1              = b1m,
      b2              = b2m,
      b3              = b3m,
      b4              = b4m,
      n_Absent        = nA,
      n_Linear        = nL,
      n_Exponential   = nE,
      n_Saturating    = nS,
      n_Revlog        = nR,
      n_Uncertain     = nU,
      modal_shape     = modal
    )
    
    # Store in result_list
    result_list[[counter]] <- out
    counter <- counter + 1L
  }
}

###############################################################################
# 7) Bind all predictor×group summaries into one data.frame
###############################################################################
result_df <- bind_rows(result_list)

###############################################################################
# 8) Pivot so that each predictor gets its own suffix on all summary columns
###############################################################################
library(tidyr)
summary_df <- pivot_wider(
  result_df,
  names_from   = predictor,
  values_from  = c(q_25, q_50, q_75,
                   b1, b2, b3, b4,
                   n_Absent, n_Linear, n_Exponential,
                   n_Saturating, n_Revlog, n_Uncertain,
                   modal_shape),
  names_glue   = "{predictor}_{.value}"
)

###############################################################################
# 9) Re‐order: IDs first, then hfp_*, then het_*
###############################################################################
summary_df <- relocate(
  summary_df,
  dataset, facet, direction, beta_type, metric_type
)
summary_df <- relocate(
  summary_df,
  matches("^hfp_"),
  .after = metric_type
)
summary_df <- relocate(
  summary_df,
  matches("^het_"),
  .after = last_col()
)

# `summary_df` is now your final wide table, easy to debug at each step
