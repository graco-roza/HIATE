###############################################################################
# SCRIPT NAME: Beta Diversity Estimation (All Methods Combined in One Loop)
#
# DESCRIPTION:
#   This script calculates taxonomic and functional beta diversity for the dataset,
#   simultaneously producing:
#     - Podani (abundance and presence–absence) beta diversity,
#     - Baselga (abundance and presence–absence) beta diversity.
#
#   For taxonomic beta, it uses:
#     - BAT::beta (for Podani partitions)
#     - betapart::beta.pair.abund / betapart::beta.pair (for Baselga partitions)
#
#   For functional beta, it calculates two sets of alpha hypervolumes (one
#   abundance-based and one presence–absence-based) and then computes pairwise
#   comparisons with both Podani and Baselga equations in a single parallel loop.
#
# OUTPUTS:
#   Results are saved in the following folders:
#     - Podani_abun, Podani_pa, Baselga_abun, Baselga_pa
#
###############################################################################

.libPaths(c("/projappl/project_2004932/project_rpackages_4.4.0", .libPaths()))

# Load required packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  BAT,          # for taxonomic beta (Podani)
  betapart,     # for taxonomic beta (Baselga)
  hypervolume,  # for functional beta (hypervolumes)
  tidyverse,    # for data manipulation
  magrittr,     # for piping and set_names
  glue,         # for string interpolation
  readxl,       # reading files if needed
  future,       # for parallel processing
  doSNOW,       # for parallel processing
  doFuture
)

# Utility function to combine nested list results from foreach loops
appendList <- function(x, val) {
  stopifnot(is.list(x), is.list(val))
  xnames <- names(x)
  for (v in names(val)) {
    x[[v]] <- if (v %in% xnames && is.list(x[[v]]) && is.list(val[[v]]))
      appendList(x[[v]], val[[v]])
    else c(x[[v]], val[[v]])
  }
  x
}

# -------------------------------------------------------------------------
# 1. Data Input & Preparation
# -------------------------------------------------------------------------

# HPC cluster: get dataset index from environment variable (if applicable)
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
ii <- as.numeric(slurm_arrayid)
# Get the focal dataset name and read the pre-processed file
files <- tools::file_path_sans_ext(list.files("S2_get_beta_diversity/betadiv_input"))
focal_dataset <- gsub("_beta_Input", "", files[ii])
write(focal_dataset, stderr())
species_diff <- glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds") %>% read_rds()

comm <- species_diff %>% pluck("comm")         # community data
traits <- species_diff %>% pluck("trait_syndrome")  # trait (PCoA axes) data

# Restrict community to species with trait data
comm <- comm %>% select(any_of(rownames(traits)))

# Create two versions of the community matrix:
comm_abun <- comm                              # abundance‐based
comm_pa   <- ifelse(comm > 0, 1, 0)              # presence–absence

# -------------------------------------------------------------------------
# 2. TAXONOMIC BETA DIVERSITY
# -------------------------------------------------------------------------

# 2a. Podani partitions using BAT::beta
Podani_tax_abun <- BAT::beta(comm_abun, abund = TRUE, func = "sorensen")
Podani_tax_pa   <- BAT::beta(comm_pa, abund = FALSE, func = "sorensen")

# 2b. Baselga partitions using betapart package
Baselga_tax_abun <- betapart::beta.pair.abund(comm_abun, index.family = "bray")
Baselga_tax_pa   <- betapart::beta.pair(comm_pa, index.family = "sorensen")
# (Adjust 'index.family' if needed)

Baselga_tax_abun <- Baselga_tax_abun[c("beta.bray", "beta.bray.bal", "beta.bray.gra")]
names(Baselga_tax_abun) <- c("Btotal", "Bturn", "Bnes")


Baselga_tax_pa <- Baselga_tax_pa[c("beta.sor", "beta.sim", "beta.sne")]
names(Baselga_tax_pa) <- c("Btotal", "Bturn", "Bnes")
betapart::functional.beta.pair(comm,)
# -------------------------------------------------------------------------
# 3. FUNCTIONAL BETA DIVERSITY (Using Hypervolumes)
# -------------------------------------------------------------------------

# Compute two sets of alpha hypervolumes: abundance-based and presence–absence-based.
cores <- parallelly::availableCores()-1

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
pairwise_all <- foreach(i = 1:alpha_n, .combine = appendList,   .options.future = list(seed = TRUE)) %:%
  foreach(j = i:alpha_n, .combine = appendList) %dopar% {
    set.seed(123)
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

# Podani - abundance based
betadiv_Podani_abun <- list(
  Taxonomic = Podani_tax_abun,
  Functional = Functional_Podani_abun
)
saveRDS(object = betadiv_Podani_abun,
        file = glue::glue("/scratch/project_2004932/HIATE/S2_get_beta_diversity/betadiv_output/Podani_abun/{focal_dataset}_beta_Output.rds"))

# Podani - presence-absence
betadiv_Podani_pa <- list(
  Taxonomic = Podani_tax_pa,
  Functional = Functional_Podani_pa
)
saveRDS(object = betadiv_Podani_pa,
        file = glue::glue("/scratch/project_2004932/HIATE/S2_get_beta_diversity/betadiv_output/Podani_pa/{focal_dataset}_beta_Output.rds"))

# Baselga - abundance based
betadiv_Baselga_abun <- list(
  Taxonomic = Baselga_tax_abun,
  Functional = Functional_Baselga_abun
)
saveRDS(object = betadiv_Baselga_abun,
        file = glue::glue("/scratch/project_2004932/HIATE/S2_get_beta_diversity/betadiv_output/Baselga_abun/{focal_dataset}_beta_Output.rds"))

# Baselga - presence-absence
betadiv_Baselga_pa <- list(
  Taxonomic = Baselga_tax_pa,
  Functional = Functional_Baselga_pa
)
saveRDS(object = betadiv_Baselga_pa,
        file = glue::glue("/scratch/project_2004932/HIATE/S2_get_beta_diversity/betadiv_output/Baselga_pa/{focal_dataset}_beta_Output.rds"))
