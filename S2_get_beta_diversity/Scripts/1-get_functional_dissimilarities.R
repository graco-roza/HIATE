###############################################################################
# SCRIPT NAME: Functional Dissimilarities Estimation
#
# DESCRIPTION:
#   This script calculates functional dissimilarities between species using
#   the `gawdis` function, based on selected functional traits. 
#
# USAGE:
#   - This script is designed to process individual datasets included in the
#     analysis. Each dataset corresponds to a specific set of functional traits
#     and species.
#   - To run the analysis for a specific dataset, uncomment the relevant lines
#     and execute the script.
#
# INPUTS:
#   - Functional trait data for species
#   - Dataset files corresponding to the selected traits
#
# OUTPUTS:
#   - Functional dissimilarity matrices for each processed dataset
#
# AUTHOR: Caio Graco-Roza
# LAST UPDATED: 2024-11-20

# NOTES:
#   - This script is modular, allowing users to process datasets one at a time
#     by uncommenting specific sections.
#   - Ensure all required libraries are installed, and the input files are
#     correctly formatted.
###############################################################################

# Packages .............................................................................................................
if(!require("pacman")) {install.packages("pacman")}
pacman::p_load(
  tidyverse,
  magrittr,
  glue,
  gawdis,
  readxl,
  vegan,
  FD)

#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ barbaro_birds_1 ####
#' #'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "barbaro_birds_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%
  dplyr::select(mass,diet,forag,range) %>%  #mass - diet - forag - range
  mutate(across(.fns=ordered)) %>%  #traits are all ordinal
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

trait_dissimilarity <-  gawdis::gawdis(trait_subset)
trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ barbaro_birds_2 ####
#' #'-----------------------------------------------------------------------------------------------------------------
#'
focal_dataset <- "barbaro_birds_2"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%
  dplyr::select(mass,diet,forag,move) %>%  #mass - diet - forag - move
  mutate(across(.fns=ordered)) %>%  #traits are all ordinal
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty
trait_dissimilarity <-  gawdis::gawdis(trait_subset)
trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ barber_beetles ####
#' #'-----------------------------------------------------------------------------------------------------------------
#'

focal_dataset <- "barber_beetles"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%
  dplyr::select(SIZE_mean.body.length
               ,DIET_diet
               ,HABITAT_grasslands, HABITAT_cultivated.fields, HABITAT_open.forests.and.other.open.habitats
               ,MOBILITY_wing.morphology) %>%
  mutate(across(c(DIET_diet,MOBILITY_wing.morphology), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty


trait_dissimilarity <-  gawdis::gawdis(trait_subset,
                                       groups = c(1,
                                                  2,
                                                  3,3,3,
                                                  4), w.type = "optimized")
trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ bdm_birds ####
#' #'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "bdm_birds"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

trait_dissimilarity <-  gawdis(trait_subset,groups=c(1
                                                             ,rep(2,10)
                                                             ,rep(3,5)
                                                             ,4)
                               , fuzzy=c(2,3)
                               ,w.type = "optimized")
trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm %>% dplyr::select(any_of(rownames(trait_syndrome))),#remove species with no traits
               trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'
#'
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ bdm_butterflies ####
#' #'-----------------------------------------------------------------------------------------------------------------
#'
focal_dataset <- "bdm_butterflies"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%
mutate(across(where(is.character),factor)) %>%  #convert character to factor
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

trait_dissimilarity <- gawdis(trait_subset, groups = c(1,
                                                 2,
                                                 rep(3, 24),
                                                 rep(4, 4)),
 w.type = "optimized")
names(trait_subset)
trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm %>% dplyr::select(any_of(rownames(trait_syndrome))),#remove species with no traits
               trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())


#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ bettergardens_bowltrap_apidae ####
#' #'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "bettergardens_bowltrap_apidae"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%
  dplyr::select(ITD_mean_ffww,feeding_spec,habitat_generalism,feeddist_km_ffww) %>%
  mutate(across(where(is.character),factor)) %>% #convert character to factor
  mutate(across(c(feeding_spec,habitat_generalism), ordered)) %>%
  slice(which(apply(., 1, function(x) all(!is.na(x)))))  #remove species where all traits are empty

trait_dissimilarity <- gawdis(trait_subset, w.type="analytic")
trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm %>% dplyr::select(any_of(rownames(trait_syndrome))),#remove species with no traits
               trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'
#'
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ bettergardens_bowltrap_carabidae ####   Only 2 sites
#' #'-----------------------------------------------------------------------------------------------------------------
# focal_dataset <- "bettergardens_bowltrap_carabidae"
# data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
# data_trait <- data_clean %>%  pluck("trait")
# trait_subset<- data_trait %>%
#   dplyr::select(bodylength,feeding_guild,habitatgeneralism,wing_form) %>%
#   mutate(across(where(is.character),factor)) %>% #convert character to factor
# slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

# trait_dissimilarity <- gawdis(trait_subset, w.type="optimized")
#
# trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#
# output <- list(comm = data_clean$comm %>% dplyr::select(any_of(rownames(trait_syndrome))),#remove species with no traits
#                trait_syndrome = trait_syndrome)
#
# saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
# rm(list=ls())
#'
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ bettergardens_bowltrap_heteroptera ####
#' #'-----------------------------------------------------------------------------------------------------------------
#'
# focal_dataset <- "bettergardens_bowltrap_heteroptera"
# data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
# data_trait <- data_clean %>%  pluck("trait")
# trait_subset<- data_trait %>%
#   dplyr::select(bodylength,Feeding_guild,Stratum,Wing_form) %>%
#   mutate(across(where(is.character),factor)) %>%
#   slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

#' trait_dissimilarity <- gawdis(trait_subset, w.type="optimized")
#'
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#' plot(trait_syndrome)
#' output <- list(comm = data_clean$comm %>% dplyr::select(any_of(rownames(trait_syndrome))),#remove species with no traits
#'                trait_syndrome = trait_syndrome)
#'
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#'
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ bettergardens_bowltrap_homoptera ####
#' #'-----------------------------------------------------------------------------------------------------------------
#'
# focal_dataset <- "bettergardens_bowltrap_homoptera"
# data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
# data_trait <- data_clean %>%  pluck("trait")
# trait_subset<- data_trait %>%
#   dplyr::select(Size_mean,Diet_width,Stratum,Wing_form) %>%
#   mutate(across(where(is.character),factor)) %>%
#   slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty   #convert to factor


#' trait_dissimilarity <- gawdis(trait_subset)
#'
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#'
#' output <- list(comm = data_clean$comm %>% dplyr::select(any_of(rownames(trait_syndrome))),#remove species with no traits
#'                trait_syndrome = trait_syndrome)
#'
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#'
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ bettergardens_bowltrap_syrphidae ####
#' #'-----------------------------------------------------------------------------------------------------------------
#'
# focal_dataset <- "bettergardens_bowltrap_syrphidae"
# data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
# data_trait <- data_clean %>%  pluck("trait")
# trait_subset<- data_trait %>%
#   dplyr::select(ITD,food_specialization_adult,saproxylic_larvae,WingLoading) %>%
#   mutate(across(food_specialization_adult,ordered)) %>%
#   slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

#'
#' trait_dissimilarity <- gawdis(trait_subset)
#'
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#'
#' output <- list(comm = data_clean$comm %>% dplyr::select(any_of(rownames(trait_syndrome))),#remove species with no traits
#'                trait_syndrome = trait_syndrome)
#'
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ bettergardens_gastropoda ####
#' #'-----------------------------------------------------------------------------------------------------------------
#'
# focal_dataset <- "bettergardens_gastropoda"
# data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
# data_trait <- data_clean %>%  pluck("trait")
# trait_subset<- data_trait %>%
#   dplyr::select(shell_size,feeding_guild,humidity,snail_mobility_index) %>%
#   #convert all to numeric
#   mutate(across(.fns= as.numeric)) %>%
#   slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

#' trait_dissimilarity <- gawdis(trait_subset, w.type="optimized")
#'
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#'
#' output <- list(comm = data_clean$comm %>% dplyr::select(any_of(rownames(trait_syndrome))),#remove species with no traits
#'                trait_syndrome = trait_syndrome)
#'
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#'
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ bettergardens_pitfall_araneae ####
#' #'-----------------------------------------------------------------------------------------------------------------
# focal_dataset <- "bettergardens_pitfall_araneae"
# data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
# data_trait <- data_clean %>%  pluck("trait")
# trait_subset<- data_trait %>%
#    dplyr::select(bodylength,food_preference,habitat_width,mobility_mode) %>%
#   #convert to factor
#   mutate(across(where(is.character),factor)) %>%
#   slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

#'
#' trait_dissimilarity <- gawdis(trait_subset)
#'
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#'
#' output <- list(comm = data_clean$comm %>% dplyr::select(any_of(rownames(trait_syndrome))),#remove species with no traits
#'                trait_syndrome = trait_syndrome)
#'
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ bettergardens_pitfall_carabidae ####
#' #'-----------------------------------------------------------------------------------------------------------------
# focal_dataset <- "bettergardens_pitfall_carabidae"
# data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
# data_trait <- data_clean %>%  pluck("trait")
# trait_subset<- data_trait %>%
#   dplyr::select(bodylength,feeding_guild,habitatgeneralism,wing_form) %>%
#   #convert to factor
#   mutate(across(where(is.character),factor)) %>%
#   slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

#'
#' trait_dissimilarity <- gawdis(trait_subset)
#'
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#'
#' output <- list(comm = data_clean$comm %>% dplyr::select(any_of(rownames(trait_syndrome))),#remove species with no traits
#'                trait_syndrome = trait_syndrome)
#'
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ bettergardens_pitfall_heteroptera ####
#' #'-----------------------------------------------------------------------------------------------------------------
# focal_dataset <- "bettergardens_pitfall_heteroptera"
# data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
# data_trait <- data_clean %>%  pluck("trait")
# trait_subset<- data_trait %>%
#   dplyr::select(bodylength,Feeding_guild,Stratum,Wing_form) %>%
#   #convert to factor
#   mutate(across(where(is.character),factor)) %>%
#   slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

#'
#' trait_dissimilarity <- gawdis(trait_subset)
#'
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#'
#' output <- list(comm = data_clean$comm %>% dplyr::select(any_of(rownames(trait_syndrome))),#remove species with no traits
#'                trait_syndrome = trait_syndrome)
#'
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ bettergardens_pitfall_homoptera ####
#' #'-----------------------------------------------------------------------------------------------------------------
# focal_dataset <- "bettergardens_pitfall_homoptera"
# data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
# data_trait <- data_clean %>%  pluck("trait")
# trait_subset<- data_trait %>%
#   dplyr::select(Size_mean,Diet_width,Stratum,Wing_form) %>%
#   #convert to factor
#   mutate(across(where(is.character),factor)) %>%
#   mutate(across(Diet_width,ordered)) %>%
#   slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

#'
#' trait_dissimilarity <- gawdis(trait_subset)
#'
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#'
#' output <- list(comm = data_clean$comm %>% dplyr::select(any_of(rownames(trait_syndrome))),#remove species with no traits
#'                trait_syndrome = trait_syndrome)
#'
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ bexis_forest_birds ####
#' #'-----------------------------------------------------------------------------------------------------------------
# focal_dataset <- "bexis_forest_birds"
# data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
# data_trait <- data_clean %>%  pluck("trait")
# trait_subset<- data_trait %>%
#   dplyr::select(size_bodymass,diet,contains("habitat_"),contains("mobility")) %>%
#   #convert to factor
#   mutate(across(diet,factor)) %>%
#   #convert to numeric
#   mutate(across(!diet,as.numeric)) %>%
#   slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

#' trait_dissimilarity <-  gawdis(trait_subset, groups = c(1, # size
#'                                                                 2, # diet
#'                                                                 rep(3,7), # habitat
#'                                                                 4), fuzzy=c(3), w.type="optimized")  # mobility
#'
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#'
#' output <- list(comm = data_clean$comm %>% dplyr::select(any_of(rownames(trait_syndrome))),#remove species with no traits
#'                trait_syndrome = trait_syndrome)
#'
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ bexis_forest_flightinterceptioncanopy_beetles ####
#' #'-----------------------------------------------------------------------------------------------------------------
# focal_dataset <- "bexis_forest_flightinterceptioncanopy_beetles"
# data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
# data_trait <- data_clean %>%  pluck("trait")
# trait_subset<- data_trait %>%
#   dplyr::select(Mean_BodySize,Feeding_guild_short,Stratum_use_short,Dispersal_ability) %>%
#   #convert to factor
#   mutate(across(Dispersal_ability,as.numeric)) %>%
#   mutate(across(where(is.character),factor)) %>%
#   slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

#' trait_dissimilarity <- gawdis(trait_subset)
#'
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#'
#' output <- list(comm = data_clean$comm %>% dplyr::select(any_of(rownames(trait_syndrome))),#remove species with no traits
#'                trait_syndrome = trait_syndrome)
#'
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ bexis_forest_flightinterceptionunderstory_beetles ####
#' #'-----------------------------------------------------------------------------------------------------------------
# focal_dataset <- "bexis_forest_flightinterceptionunderstory_beetles"
# data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
# data_trait <- data_clean %>%  pluck("trait")
# trait_subset<- data_trait %>%
#   dplyr::select(Mean_BodySize,Feeding_guild_short,Stratum_use_short,Dispersal_ability) %>%
#   #convert to factor
#   mutate(across(Dispersal_ability,as.numeric)) %>%
#   mutate(across(where(is.character),factor)) %>%
#   slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

#' trait_dissimilarity <- gawdis(trait_subset)
#'
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#'
#' output <- list(comm = data_clean$comm %>% dplyr::select(any_of(rownames(trait_syndrome))),#remove species with no traits
#'                trait_syndrome = trait_syndrome)
#'
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ bexis_forest_hemiptera ####
#' #'-----------------------------------------------------------------------------------------------------------------
# focal_dataset <- "bexis_forest_hemiptera"
# data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
# data_trait <- data_clean %>%  pluck("trait")
# trait_subset<- data_trait %>%
#   dplyr::select(Mean_BodySize,Feeding_guild_short,Stratum_use_short,Dispersal_ability) %>%
#   #convert to factor
#   mutate(across(Dispersal_ability,ordered)) %>%
#   mutate(across(where(is.character),factor)) %>%
#   slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

#' trait_dissimilarity <- gawdis(trait_subset)
#'
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#'
#' output <- list(comm = data_clean$comm %>% dplyr::select(any_of(rownames(trait_syndrome))),#remove species with no traits
#'                trait_syndrome = trait_syndrome)
#'
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ bexis_forest_pitfall_beetles ####
#' #'-----------------------------------------------------------------------------------------------------------------
# focal_dataset <- "bexis_forest_pitfall_beetles"
# data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
# data_trait <- data_clean %>%  pluck("trait")
# trait_subset<- data_trait %>%
#   dplyr::select(Mean_BodySize,Feeding_guild_short,Stratum_use_short,Dispersal_ability) %>%
#   #convert to factor
#   mutate(across(Dispersal_ability,ordered)) %>%
#   mutate(across(where(is.character),factor)) %>%
#   slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

#' trait_dissimilarity <- gawdis(trait_subset)
#'
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#'
#' output <- list(comm = data_clean$comm %>% dplyr::select(any_of(rownames(trait_syndrome))),#remove species with no traits
#'                trait_syndrome = trait_syndrome)
#'
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ bexis_forest_plants ####
#' #'-----------------------------------------------------------------------------------------------------------------
# focal_dataset <- "bexis_forest_plants"
# data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
# data_trait <- data_clean %>%  pluck("trait")
# trait_subset<- data_trait %>%
#   dplyr::select(Height,SLA_all,T:N,Seed_mass) %>%
#   slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

#' names_remove<- FD::gowdis(trait_subset2) %>%
#'   as.matrix %>%
#'   as.data.frame.table %>%
#'   slice(which(is.na(.$Freq))) %>%
#'   select(Var1,Var2) %>%
#'   t() %>%
#'   c() %>%
#'   unique()
#'
#' trait_subset2<-trait_subset %>%
#'   rownames_to_column("species") %>%
#'   filter(!species %in% names_remove) %>%
#'   column_to_rownames("species")
#'
#' trait_dissimilarity <-  gawdis(trait_subset2, groups = c(1
#'                                                   ,2
#'                                                   ,rep(3,5)
#'                                                   ,4), w.type="optimized")
#'
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#'
#' output <- list(comm = data_clean$comm %>% dplyr::select(any_of(rownames(trait_syndrome))),#remove species with no traits
#'                trait_syndrome = trait_syndrome)
#'
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ bexis_forest_spiders ####
#' #'-----------------------------------------------------------------------------------------------------------------
# focal_dataset <- "bexis_forest_spiders"
# data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
# data_trait <- data_clean %>%  pluck("trait")
# trait_subset<- data_trait %>%
#   dplyr::select(Mean_BodySize,contains("diet"),Stratum_use_short,Dispersal_ability) %>%
#   mutate(across(where(is.character),as.factor)) %>%
#   slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

#' trait_dissimilarity <-  gawdis(trait_subset, w.type="optimized")
#'
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#'
#' output <- list(comm = data_clean$comm %>% dplyr::select(any_of(rownames(trait_syndrome))),#remove species with no traits
#'                trait_syndrome = trait_syndrome)
#'
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ bexis_forest_waterinsects ####
#' #'-----------------------------------------------------------------------------------------------------------------
# focal_dataset <- "bexis_forest_waterinsects"
# data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
# data_trait <- data_clean %>%  pluck("trait")
# trait_subset<- data_trait %>%
#   mutate(across(where(is.character),as.factor)) %>%
#   slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

#' trait_dissimilarity <-  gawdis(trait_subset, w.type="optimized")
#'
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#'
#' output <- list(comm = data_clean$comm %>% dplyr::select(any_of(rownames(trait_syndrome))),#remove species with no traits
#'                trait_syndrome = trait_syndrome)
#'
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ bexis_grassland_birds ####
#' #'-----------------------------------------------------------------------------------------------------------------
# focal_dataset <- "bexis_grassland_birds"
# data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
# data_trait <- data_clean %>%  pluck("trait")
# trait_subset<- data_trait %>%
#   dplyr::select(size_bodymass,diet,contains("habitat"),contains("mobility")) %>%
#   mutate(across(where(is.character),as.factor)) %>%
#   slice(which(apply(., 1, function(x) !all(is.na(x))))) #remove species where all traits are empty

#' trait_dissimilarity <-  gawdis(trait_subset, groups = c(1, # size
#'                                                   2, # diet
#'                                                   3,3,3,3,3,3,3, # habitat
#'                                                   4),# mobility
#'                                fuzzy = c(3), w.type="optimized")
#'
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#'
#' output <- list(comm = data_clean$comm %>% dplyr::select(any_of(rownames(trait_syndrome))),#remove species with no traits
#'                trait_syndrome = trait_syndrome)
#'
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ bexis_grassland_pitfall_beetles ####
#' #'-----------------------------------------------------------------------------------------------------------------
# focal_dataset <- "bexis_grassland_pitfall_beetles"
# data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
# data_trait <- data_clean %>%  pluck("trait")
# trait_subset<- data_trait %>%
#   mutate(across(where(is.character),as.factor)) %>%
#   mutate(Wingclass = ordered(Wingclass)) %>%
#   slice(which(apply(., 1, function(x) !all(is.na(x))))) #remove species where all traits are empty

#' trait_dissimilarity <-  gawdis(trait_subset)
#'
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#'
#' output <- list(comm = data_clean$comm %>% dplyr::select(any_of(rownames(trait_syndrome))),#remove species with no traits
#'                trait_syndrome = trait_syndrome)
#'
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#'
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ bexis_grassland_pitfall_spiders ####
#' #'-----------------------------------------------------------------------------------------------------------------
# focal_dataset <- "bexis_grassland_pitfall_spiders"
# data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
# data_trait <- data_clean %>%  pluck("trait")
# trait_subset<- data_trait %>%
#   select(Size_mean,contains("diet"),Stratum,Mobility) %>%
#   mutate(across(where(is.character),as.factor)) %>%
#   mutate(Mobility = ordered(Mobility)) %>%
#   slice(which(apply(., 1, function(x) !all(is.na(x))))) #remove species where all traits are empty

#' trait_dissimilarity <-  gawdis(trait_subset)
#'
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#'
#' output <- list(comm = data_clean$comm %>% dplyr::select(any_of(rownames(trait_syndrome))),#remove species with no traits
#'                trait_syndrome = trait_syndrome)
#'
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#'
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ bexis_grassland_plants ####
#' #'-----------------------------------------------------------------------------------------------------------------
#' focal_dataset <- "bexis_grassland_plants"
#' data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
#' data_trait <- data_clean %>%  pluck("trait")
#'
#' trait_subset<- data_trait %>%
#'   mutate(across(where(is.character),as.factor)) %>%
#'   slice(which(apply(., 1, function(x) !all(is.na(x))))) #remove species where all traits are empty
#'
#' #there was some issue with a couple of species that had NA values in trait dissimilarity, so we remove those
#' names_remove<- FD::gowdis(trait_subset) %>%
#'   as.matrix %>%
#'   as.data.frame.table %>%
#'     filter(is.na(value)) %>% distinct(Var2) %>% pull(Var2) %>%  as.character
#'
#' trait_subset2<-trait_subset %>%
#'   rownames_to_column("species") %>%
#'   filter(!species %in% names_remove) %>%
#'   column_to_rownames("species")
#'
#' trait_dissimilarity <-  gawdis(trait_subset2, groups=c(1,2,rep(3,5),4,3), groups.weight = TRUE, w.type="optimized")
#'
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#'
#'
#' output <- list(comm = data_clean$comm %>% dplyr::select(any_of(rownames(trait_syndrome))),#remove species with no traits
#'                trait_syndrome = trait_syndrome)
#'
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#'
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ bexis_grassland_sweepnet_beetles ####
#' #'-----------------------------------------------------------------------------------------------------------------
#' focal_dataset <- "bexis_grassland_sweepnet_beetles"
#' data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
#' data_trait <- data_clean %>%  pluck("trait")
#' trait_subset<- data_trait %>%
#'   mutate(across(where(is.character),as.factor)) %>%
#'   slice(which(apply(., 1, function(x) !all(is.na(x))))) #remove species where all traits are empty
#'
#' trait_dissimilarity <-  gawdis(trait_subset, w.type="optimized")
#'
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#'
#' output <- list(comm = data_clean$comm %>% dplyr::select(any_of(rownames(trait_syndrome))),#remove species with no traits
#'                trait_syndrome = trait_syndrome)
#'
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#'
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ bexis_grassland_sweepnet_heteroptera ####
#' #'-----------------------------------------------------------------------------------------------------------------
#' focal_dataset <- "bexis_grassland_sweepnet_heteroptera"
#' data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
#' data_trait <- data_clean %>%  pluck("trait")
#' trait_subset<- data_trait %>%
#'   dplyr::select(Mean_BodySize,Feeding_guild_short, Stratum_use_short, Dispersal_ability) %>%
#'   mutate(across(where(is.character),as.factor)) %>%
#'   mutate(Dispersal_ability = ordered(Dispersal_ability)) %>%
#'   slice(which(apply(., 1, function(x) !all(is.na(x))))) #remove species where all traits are empty
#'
#' trait_dissimilarity <-  gawdis(trait_subset, w.type="optimized")
#'
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#'
#' output <- list(comm = data_clean$comm %>% dplyr::select(any_of(rownames(trait_syndrome))),#remove species with no traits
#'                trait_syndrome = trait_syndrome)
#'
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#'
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ bexis_grassland_sweepnet_homoptera ####
#' #'-----------------------------------------------------------------------------------------------------------------
#' focal_dataset <- "bexis_grassland_sweepnet_homoptera"
#' data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
#' data_trait <- data_clean %>%  pluck("trait")
#' trait_subset<- data_trait %>%
#'   dplyr::select(Mean_BodySize,diet_feedingtyp, Stratum_use_short, Dispersal_ability) %>%
#'   mutate(across(where(is.character),as.factor)) %>%
#'   mutate(Dispersal_ability = ordered(Dispersal_ability)) %>%
#'   slice(which(apply(., 1, function(x) !all(is.na(x))))) #remove species where all traits are empty
#'
#' trait_dissimilarity <-  gawdis(trait_subset, w.type="optimized")
#'
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#'
#' output <- list(comm = data_clean$comm %>% dplyr::select(any_of(rownames(trait_syndrome))),#remove species with no traits
#'                trait_syndrome = trait_syndrome)
#'
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#'
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ bexis_grassland_sweepnet_orthoptera ####
#' # we had to remove Stratum_use_short because all species had same values
#' #'-----------------------------------------------------------------------------------------------------------------
#' focal_dataset <- "bexis_grassland_sweepnet_orthoptera"
#' data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
#' data_trait <- data_clean %>%  pluck("trait")
#' trait_subset<- data_trait %>%
#'   dplyr::select(Mean_BodySize,Feeding_guild_short, Dispersal_ability) %>%
#'   mutate(across(where(is.character),as.factor)) %>%
#'   mutate(Dispersal_ability = ordered(Dispersal_ability)) %>%
#'   slice(which(apply(., 1, function(x) !all(is.na(x))))) #remove species where all traits are empty
#'
#' trait_dissimilarity <-  gawdis(trait_subset, w.type="optimized")
#'
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#'
#' output <- list(comm = data_clean$comm %>% dplyr::select(any_of(rownames(trait_syndrome))),#remove species with no traits
#'                trait_syndrome = trait_syndrome)
#'
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ bexis_grassland_sweepnet_spiders ####
#' #'-----------------------------------------------------------------------------------------------------------------
#' focal_dataset <- "bexis_grassland_sweepnet_spiders"
#' data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
#' data_trait <- data_clean %>%  pluck("trait")
#' trait_subset<- data_trait %>%
#'   dplyr::select(Mean_BodySize,contains("diet"), Stratum_use_short, Dispersal_ability) %>%
#'   mutate(across(where(is.character),as.factor)) %>%
#'   mutate(Dispersal_ability = ordered(Dispersal_ability)) %>%
#'   slice(which(apply(., 1, function(x) !all(is.na(x))))) #remove species where all traits are empty
#'
#' trait_dissimilarity <-  gawdis(trait_subset, w.type="optimized")
#'
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#'
#' output <- list(comm = data_clean$comm %>% dplyr::select(any_of(rownames(trait_syndrome))),#remove species with no traits
#'                trait_syndrome = trait_syndrome)
#'
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ bonada_waterinsects ####
#' #'-----------------------------------------------------------------------------------------------------------------
#' focal_dataset <- "bonada_waterinsects"
#' data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
#' data_trait <- data_clean %>%  pluck("trait")
#' trait_subset<- data_trait %>%
#'   dplyr::select(a1:a7 #size
#'                ,h1:h9 #diet
#'                ,u1:u8 #habitat
#'                ,f1:f4 #mobility
#'                ) %>%
#'   slice(which(apply(., 1, function(x) !all(is.na(x))))) #remove species where all traits are empty
#'
#' trait_dissimilarity <-  gawdis(trait_subset, groups =c(rep(1,7)
#'                                                        ,rep(2,9)
#'                                                        ,rep(3,8)
#'                                                        ,rep(4,4))
#'                                , fuzzy=1:4
#'                                ,w.type="optimized")
#'
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#'
#' output <- list(comm = data_clean$comm %>% dplyr::select(any_of(rownames(trait_syndrome))),#remove species with no traits
#'                trait_syndrome = trait_syndrome)
#'
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ brindamour_fish ####
#' #'-----------------------------------------------------------------------------------------------------------------
#' focal_dataset <- "brindamour_fish"
#' data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
#' data_trait <- data_clean %>%  pluck("trait")
#' trait_subset<- data_trait %>%
#'   dplyr::select(max_length_fishbase, contains("Diet"),contains("Feed"),contains("migr")) %>%
#'   mutate(across(where(is.character),as.factor)) %>%
#'   slice(which(apply(., 1, function(x) !all(is.na(x))))) #remove species where all traits are empty
#'
#' trait_dissimilarity <-  gawdis(trait_subset,groups=c(1
#'                                                      ,rep(2,4)
#'                                                      ,rep(3,4)
#'                                                      ,rep(4,2))
#'                                ,fuzzy=c(2:4)
#'                                ,w.type='optimized')
#'
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#'
#' output <- list(comm = data_clean$comm %>% dplyr::select(any_of(rownames(trait_syndrome))),#remove species with no traits
#'                trait_syndrome = trait_syndrome)
#'
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ brodersen_fish ####
#' #'-----------------------------------------------------------------------------------------------------------------
# focal_dataset <- "brodersen_fish"
# data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
# data_trait <- data_clean %>%  pluck("trait")
# trait_subset<- data_trait %>%
#   dplyr::select(Weight_g, # size
#                 contains("diet"), # diet
#                 freshwater,brackish,marine,demersal,benthopelagic,pelagic,anadromous,potamodromous,catadromous, # habitat
#                 aspect.ratio.of.caudal.fin) %>% #mobility
#   mutate(across(.fns=as.numeric))
#   slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty
#   names(trait_subset)
#

#' trait_dissimilarity <-  gawdis::gawdis(trait_subset
#'                                        , groups = c(1 # size
#'                                                     , rep(2, 7) # diet
#'                                                     , rep(3, 9) # habitat
#'                                                     , 4) #mobility
#'                                        , fuzzy = 2:3
#'                                        , w.type="optimized")
#'
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#'
#' output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)
#'
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ carvalho_fish ####
#' #'-----------------------------------------------------------------------------------------------------------------
# focal_dataset <- "carvalho_fish"
# data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
# data_trait <- data_clean %>%  pluck("trait")
# trait_subset<- data_trait %>%
#   dplyr::select(BM # size
#                 ,paste0("troph",c(1:3,5:6)) # diet (troph 4 and 7 only zero)
#                 ,contains("watpos") # habitat
#                 ,contains("migra")) %>% #mobility
#   mutate(across(.fns=as.numeric)) %>%
# slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty
#

#' trait_dissimilarity <-  gawdis::gawdis(trait_subset
#'                                        , groups = c(1 # size
#'                                                     , rep(2, 5) # diet
#'                                                     , rep(3, 3) # habitat
#'                                                     , rep(4,2)) #mobility
#'                                        , fuzzy = 2:4
#'                                        , w.type="optimized")
#'
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#'
#' output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)
#'
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ chapman_birds ####
#' #'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "chapman_birds"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%
  dplyr::select(pca.bodysize # size
                ,pca.trophic.shape # diet (troph 4 and 7 only zero)
                ,habitat # habitat
                ,pca.dispersal.shape) %>% #mobility
  mutate(across(habitat, .fns=ordered)) %>%
slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty


#' 
#' trait_dissimilarity <-  gawdis::gawdis(trait_subset
#'                                        , w.type="optimized")
#' 
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#' 
#' output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)
#' 
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ decastro_waterinsects ####
#' #'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "decastro_waterinsects"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%
  dplyr::select(contains("size") # size
                ,contains("diet") # diet (troph 4 and 7 only zero)
                ,contains("habitat") # habitat
                ,contains("mobility")) %>% #mobility
  mutate(across(.fns=as.numeric)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

#' trait_dissimilarity <-  gawdis::gawdis(trait_subset
#'                                        , groups = c(rep(1,6) # size
#'                                                     , rep(2, 5) # diet
#'                                                     , rep(3, 3) # habitat
#'                                                     , rep(4, 5)) #mobility
#'                                        , fuzzy = 2:4
#'                                        , w.type="optimized")
#' 
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#' 
#' output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)
#' 
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ desouzaqueiroz_amphibians ####
#' #'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "desouzaqueiroz_amphibians"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%
  dplyr::select(size_max # size
                ,Feeding_behaviour # diet (troph 4 and 7 only zero)
                ,Neustonic:Nektonic # habitat
                ,WTM.TL) %>% #mobility
  mutate(across(Feeding_behaviour, .fns=factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty
# Example output you want to copy
output<-summarytools::dfSummary(trait_subset, varnumbers = FALSE, labels.col=TRUE,na.col = FALSE, style="grid", graph.col=FALSE)
write.table(output, file = pipe("pbcopy"), sep = "\t", row.names = FALSE, col.names = TRUE)
#' trait_dissimilarity <-  gawdis::gawdis(trait_subset
#'                                        , groups = c(rep(1,1) # size
#'                                                     , rep(2, 1) # diet
#'                                                     , rep(3, 3) # habitat
#'                                                     , rep(4, 1)) #mobility
#'                                        , fuzzy = c(3)
#'                                        , w.type="equal")
#' 
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#' 
#' output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)
#' 
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ drose_ants ####
#' #'-----------------------------------------------------------------------------------------------------------------
# focal_dataset <- "drose_ants"
# data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
# data_trait <- data_clean %>%  pluck("trait")
# trait_subset<- data_trait %>%
#   dplyr::select(contains("SIZE") # size
#                 ,contains("DIET") # diet (troph 4 and 7 only zero)
#                 ,contains("HABITAT") # habitat
#                 ,contains("MOBILITY")) %>% #mobility
#   slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty


#' trait_dissimilarity <-  gawdis::gawdis(trait_subset
#'                                        , w.type="optimized")
#' 
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#' 
#' output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)
#' 
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ gallardo_waterinsects ####
#' #'-----------------------------------------------------------------------------------------------------------------
# focal_dataset <- "gallardo_waterinsects"
# data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
# data_trait <- data_clean %>%  pluck("trait")
# trait_subset<- data_trait %>%
#   dplyr::select(a1:a7 # size
#                 ,k1:k9 # diet (troph 4 and 7 only zero)
#                 ,j1:j9 # habitat
#                 ,h1:h8) %>% #mobility
#   mutate(across(.fns=as.numeric)) %>%
#   slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

#' 
#' trait_dissimilarity <-  gawdis::gawdis(trait_subset
#'                                        , groups = c(rep(1,7) # size
#'                                                     , rep(2, 9) # diet
#'                                                     , rep(3, 9) # habitat
#'                                                     , rep(4, 8)) #mobility
#'                                        , fuzzy = 1:4
#'                                        , w.type="analytic")
#' 
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#' 
#' output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)
#' 
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ gracoroza_phytoplankton_1 ####
#' #'-----------------------------------------------------------------------------------------------------------------
# focal_dataset <- "gracoroza_phytoplankton_1"
# data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
# data_trait <- data_clean %>%  pluck("trait")
# trait_subset<- data_trait %>%
#   dplyr::select(MLD # size
#                 ,Mucilage:Aerotopes # diet (troph 4 and 7 only zero)
#                 ,Heterocyte:Silica # habitat
#                 ,Flagella) %>% #mobility
#   mutate(across(.fns=as.numeric)) %>%
#   slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

#' trait_dissimilarity <-  gawdis::gawdis(trait_subset
#'                                        , groups = c(rep(1,1) # size
#'                                                     , rep(2, 2) # diet
#'                                                     , rep(3, 2) # habitat
#'                                                     , rep(4, 1)) #mobility
#'                                                                       , w.type="optimized")
#' 
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#' 
#' output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)
#' 
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ gracoroza_phytoplankton_2 ####
#' #'-----------------------------------------------------------------------------------------------------------------
# focal_dataset <- "gracoroza_phytoplankton_2"
# data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
# data_trait <- data_clean %>%  pluck("trait")
# trait_subset<- data_trait %>%
#   dplyr::select(MLD # size
#                 ,Mucilage:Aerotopes # diet (troph 4 and 7 only zero)
#                 ,Heterocyte:Silica # habitat
#                 ,Flagella) %>% #mobility
#   mutate(across(.fns=as.numeric)) %>%
#   slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

#' trait_dissimilarity <-  gawdis::gawdis(trait_subset
#'                                        , groups = c(rep(1,1) # size
#'                                                     , rep(2, 2) # diet
#'                                                     , rep(3, 2) # habitat
#'                                                     , rep(4, 1)) #mobility
#'                                        , w.type="optimized")
#' 
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#' 
#' output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)
#' 
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ gracoroza_phytoplankton_3 ####
#' #'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "gracoroza_phytoplankton_3"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%
  dplyr::select(MLD,sv # size
                ,Mucilage:Aerotopes # diet (troph 4 and 7 only zero)
                ,Silica # habitat
                ,Flagella) %>% #mobility
  mutate(across(.fns=as.numeric)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

#' trait_dissimilarity <-  gawdis::gawdis(trait_subset
#'                                        , groups = c(rep(1,2) # size
#'                                                     , rep(2, 2) # diet
#'                                                     , rep(3, 1) # habitat
#'                                                     , rep(4, 1)) #mobility
#'                                        , w.type="optimized")
#' 
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#' 
#' output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)
#' 
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ jeliazkov_waterinsects ####
#' #'-----------------------------------------------------------------------------------------------------------------
#' 
# focal_dataset <- "jeliazkov_waterinsects"
# data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
# data_trait <- data_clean %>%  pluck("trait")
# trait_subset<- data_trait %>%
#   dplyr::select(contains("Size") # size
#                 ,contains("Food") # diet (troph 4 and 7 only zero)
#                 ,contains("Sub") # habitat
#                 ,contains("Loco")) %>% #mobility
#   mutate(across(.fns=as.numeric)) %>%
#   slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

#' trait_dissimilarity <-  gawdis::gawdis(trait_subset
#'                                        , groups = c(rep(1,7) # size
#'                                                     , rep(2, 9) # diet
#'                                                     , rep(3, 9) # habitat
#'                                                     , rep(4, 8)) #mobility
#'                                        , w.type="analytic")
#' 
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#' 
#' output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)
#' 
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#' 
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ lee_ants ####
#' #'-----------------------------------------------------------------------------------------------------------------
#' 
# focal_dataset <- "lee_ants"
# data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
# data_trait <- data_clean %>%  pluck("trait")
# trait_subset<- data_trait %>%
#   dplyr::select(WL # size
#                 ,ML # diet (troph 4 and 7 only zero)
#                 ,HW # habitat
#                 ,LL) %>% #mobility
#   mutate(across(.fns=as.numeric)) %>%
#   slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

#' trait_dissimilarity <-  gawdis::gawdis(trait_subset
#'                                           , w.type="analytic")
#' 
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#' 
#' output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)
#' 
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#' 
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ lehikoinen_birds ####
#' #'-----------------------------------------------------------------------------------------------------------------
# 
# focal_dataset <- "lehikoinen_birds"
# data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
# data_trait <- data_clean %>%  pluck("trait")
# trait_subset<- data_trait %>%
#   dplyr::select(Mass # size
#                 ,Diet.Inv:Diet.PlantO # diet (troph 4 and 7 only zero)
#                 ,contains("ForStrat") # habitat
#                 ,contains("Mig")) %>% #mobility
#   mutate(across(Mig, .fns=factor)) %>%
#   slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty
# 

#' trait_dissimilarity <-  gawdis::gawdis(trait_subset
#'                                        , groups = c(rep(1,1) # size
#'                                                     , rep(2, 10) # diet
#'                                                     , rep(3, 7) # habitat
#'                                                     , rep(4, 1)) #mobility
#'                                        ,fuzzy=c(2,3)
#'                                        , w.type="optimized")
#' 
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#' 
#' output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)
#' 
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#' 
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ lowe_spiders ####
#' #'-----------------------------------------------------------------------------------------------------------------
# focal_dataset <- "lowe_spiders"
# data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
# data_trait <- data_clean %>%  pluck("trait")
# trait_subset<- data_trait %>%
#   dplyr::select(contains("SIZE") # size
#                 ,contains("DIET") # diet (troph 4 and 7 only zero)
#                 ,contains("HABITAT") # habitat
#                 ,contains("MOBILITY")) %>% #mobility
#   mutate(across(c(contains("DIET")), .fns=factor)) %>%
#   mutate(SIZE_Size = ordered(SIZE_Size, levels=c("small","medium","large"))) %>%
#   slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty


#' trait_dissimilarity <-  gawdis::gawdis(trait_subset
#'                                        , w.type="optimized")
#' 
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#' 
#' output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)
#' 
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ luoto_butterflies ####
#' #'-----------------------------------------------------------------------------------------------------------------
# focal_dataset <- "luoto_butterflies"
# data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
# data_trait <- data_clean %>%  pluck("trait")
# trait_subset<- data_trait %>%  
#   dplyr::select(Mass # size
#                 ,Diet.Inv:Diet.PlantO # diet (troph 4 and 7 only zero)
#                 ,contains("ForStrat") # habitat
#                 ,contains("Mig")) %>% #mobility
#   mutate(across(Mig, .fns=factor)) %>% 
#   slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 
# 
# trait_dissimilarity <-  gawdis::gawdis(trait_subset
#                                        , groups = c(rep(1,1) # size
#                                                     , rep(2, 10) # diet
#                                                     , rep(3, 7) # habitat
#                                                     , rep(4, 1)) #mobility
#                                        ,fuzzy=c(2,3)
#                                        , w.type="optimized")
# 
# trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
# 
# output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)
# 
# saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
# rm(list=ls())
#' 
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ meffert_birds ####
#' #'-----------------------------------------------------------------------------------------------------------------
#' 
# focal_dataset <- "meffert_birds"
# data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
# data_trait <- data_clean %>%  pluck("trait")
# trait_subset<- data_trait %>%
#   dplyr::select(body.mass # size
#                 ,food # diet (troph 4 and 7 only zero)
#                 ,foraging.technique # habitat
#                 ,mig1.strategy) %>% #mobility
#   mutate(across(where(is.character), .fns=factor)) %>%
#   mutate(mig1.strategy = ordered(mig1.strategy, c("sedentary","short","long"))) %>%
#   slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty
#' trait_dissimilarity <-  gawdis::gawdis(trait_subset
#'                                        , w.type="optimized")
#' 
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#' 
#' output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)
#' 
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#' 
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ N21FBD_2 ####
#' #'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "N21FBD_2"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%
  dplyr::select(contains("SIZE") # size
                ,contains("DIET") # diet (troph 4 and 7 only zero)
                ,contains("HABITAT") # habitat
                ,contains("MOBILITY")) %>% #mobility
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

#' trait_dissimilarity <-  gawdis::gawdis(trait_subset
#'                                        , groups = c(rep(1,1) # size
#'                                                     , rep(2, 1) # diet
#'                                                     , rep(3, 4) # habitat
#'                                                     , rep(4, 1)) #mobility
#'                                                                     , w.type="optimized")
#' 
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#' 
#' output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)
#' 
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#' 
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ N21FMI ####
#' #'-----------------------------------------------------------------------------------------------------------------
# focal_dataset <- "N21FMI"
# data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
# data_trait <- data_clean %>%  pluck("trait")
# trait_subset<- data_trait %>%
#   dplyr::select(contains("SIZE") # size
#                 ,contains("DIET") # diet (troph 4 and 7 only zero)
#                 ,contains("HABITAT") # habitat
#                 ,contains("MOBILITY")) %>% #mobility
#   mutate(across(.fns=ordered)) %>%
#   slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty
# 

 
#' trait_dissimilarity <-  gawdis::gawdis(trait_subset
#'                                        ,w.type="optimized")
#' 
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#' 
#' output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)
#' 
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ N25FBD #### DIET were all the same across species 
#' #'-----------------------------------------------------------------------------------------------------------------
#' focal_dataset <- "N25FBD"
#' data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
#' data_trait <- data_clean %>%  pluck("trait")
#' trait_subset<- data_trait %>%  
#'   dplyr::select(contains("SIZE") # size
#'                 ,High.profile.guild:Planktonic # diet
#'                 ,contains("HABITAT") # habitat
#'                 ,contains("MOBILITY")) %>% #mobility
#'   mutate(SIZE_size.class = ordered(SIZE_size.class)) %>% 
#'   slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 
#' 
#' trait_dissimilarity <-  gawdis::gawdis(trait_subset
#'                                        , groups = c(rep(1,1) # size
#'                                                    , rep(2, 4) # diet
#'                                                     , rep(3, 4) # habitat
#'                                                     , rep(4, 1)) #mobility
#'                                        , w.type="optimized")
#' 
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#' 
#' output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)
#' 
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ N25FMI ####
#' #'-----------------------------------------------------------------------------------------------------------------
# focal_dataset <- "N25FMI"
# data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
# data_trait <- data_clean %>%  pluck("trait")
# trait_subset<- data_trait %>%
#   dplyr::select(contains("SIZE") # size
#                 ,contains("DIET") # die
#                 ,contains("HABITAT") # habitat
#                 ,contains("MOBILITY")) %>% #mobility
#   mutate(across(.fns=ordered)) %>%
#   slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

#' trait_dissimilarity <-  gawdis::gawdis(trait_subset
#'                                        , w.type="optimized")
#' 
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#' 
#' output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)
#' 
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ N29FBD #### check if needs the guild as diet
#' #'-----------------------------------------------------------------------------------------------------------------
#' focal_dataset <- "N29FBD"
#' data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
#' data_trait <- data_clean %>%  pluck("trait")
#' trait_subset<- data_trait %>%  
#'   dplyr::select(contains("SIZE") # size
#'                 ,contains("DIET") # diet (troph 4 and 7 only zero)
#'                 ,contains("HABITAT") # habitat
#'                 ,contains("MOBILITY")) %>% #mobility
#'   slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 
#' 
#' trait_dissimilarity <-  gawdis::gawdis(trait_subset
#'                                        , groups = c(rep(1,1) # size
#'                                                     , rep(2, 1) # diet
#'                                                     , rep(3, 4) # habitat
#'                                                     , rep(4, 1)) #mobility
#'                                        , w.type="optimized")
#' 
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#' 
#' output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)
#' 
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ N29FMI ####
#' #'-----------------------------------------------------------------------------------------------------------------
#' focal_dataset <- "N29FMI"
#' data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
#' data_trait <- data_clean %>%  pluck("trait")
#' trait_subset<- data_trait %>%  
#'   dplyr::select(contains("SIZE") # size
#'                 ,contains("DIET") # diet (troph 4 and 7 only zero)
#'                 ,contains("HABITAT") # habitat
#'                 ,contains("MOBILITY")) %>% #mobility
#'   slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 
#' 
#' trait_dissimilarity <-  gawdis::gawdis(trait_subset
#'                                        , w.type="optimized")
#' 
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#' 
#' output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)
#' 
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ N34FBD ####
#' #'-----------------------------------------------------------------------------------------------------------------
#' focal_dataset <- "N34FBD"
#' data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
#' data_trait <- data_clean %>%  pluck("trait")
#' trait_subset<- data_trait %>%  
#'   dplyr::select(contains("SIZE") # size
#'                 ,contains("DIET") # diet (troph 4 and 7 only zero)
#'                 ,contains("HABITAT") # habitat
#'                 ,contains("MOBILITY")) %>% #mobility
#'   slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 
#' 
#' trait_dissimilarity <-  gawdis::gawdis(trait_subset
#'                                        , groups = c(rep(1,1) # size
#'                                                     , rep(2, 1) # diet
#'                                                     , rep(3, 4) # habitat
#'                                                     , rep(4, 1)) #mobility
#'                                        , w.type="optimized")
#' 
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#' 
#' output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)
#' 
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ N34FMI ####
#' #'-----------------------------------------------------------------------------------------------------------------
#' #'
#' focal_dataset <- "N34FMI"
#' data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
#' data_trait <- data_clean %>%  pluck("trait")
#' trait_subset<- data_trait %>%  
#'   dplyr::select(contains("SIZE") # size
#'                 ,contains("DIET") # diet (troph 4 and 7 only zero)
#'                 ,contains("HABITAT") # habitat
#'                 ,contains("MOBILITY")) %>% #mobility
#'   slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 
#' 
#' trait_dissimilarity <-  gawdis::gawdis(trait_subset
#'                                        , w.type="optimized")
#' 
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#' 
#' output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)
#' 
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ N38TTP ####
#' #'-----------------------------------------------------------------------------------------------------------------
# focal_dataset <- "N38TTP"
# data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
# data_trait <- data_clean %>%  pluck("trait")
# trait_subset<- data_trait %>%
#   dplyr::select(contains("SIZE") # size
#                 ,contains("DIET") # diet (troph 4 and 7 only zero)
#                 ,T:N # habitat
#                 ,contains("MOBILITY")) %>% #mobility
#   slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty
# 
# names_remove<- FD::gowdis(trait_subset) %>%
#   as.matrix %>%
#   as.data.frame.table %>%
#   slice(which(is.na(.$Freq))) %>%
#   select(Var2) %>%
#   t() %>%
#   c() %>%
#   unique()
# 
# trait_subset2<-trait_subset %>%
#   rownames_to_column("species") %>%
#   filter(!species %in% names_remove) %>%
#   column_to_rownames("species")
# 
# 
# 
# 

#' trait_dissimilarity <-  gawdis::gawdis(trait_subset2
#'                                        , groups = c(rep(1,1) # size
#'                                                     , rep(2, 1) # diet
#'                                                     , rep(3, 5) # habitat
#'                                                     , rep(4, 1)) #mobility
#'                                        , w.type="optimized")
#' 
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#' 
#' output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)
#' 
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ N39FMI_2 ####
#' #'-----------------------------------------------------------------------------------------------------------------
# focal_dataset <- "N39FMI_2"
# data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
# data_trait <- data_clean %>%  pluck("trait")
# trait_subset<- data_trait %>%
#   dplyr::select(contains("SIZE") # size
#                 ,contains("food") # diet (troph 4 and 7 only zero)
#                 ,contains("feedhabit") # habitat
#                 ,contains("loc")) %>% #mobility
#   slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty
# 
# 

#' trait_dissimilarity <-  gawdis::gawdis(trait_subset
#'                                        , groups = c(rep(1,2) # size
#'                                                     , rep(2, 2) # diet
#'                                                     , rep(3, 2) # habitat
#'                                                     , rep(4, 2)) #mobility
#'                                        , w.type="optimized")
#' 
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#' 
#' output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)
#' 
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ N39FMI ####
#' #'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "N39FMI"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%
  dplyr::select(contains("SIZE") # size
                ,contains("food") # diet (troph 4 and 7 only zero)
                ,contains("feedhabit") # habitat
                ,contains("loc")) %>% #mobility
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

output<-summarytools::dfSummary(trait_subset, varnumbers = FALSE, labels.col=TRUE,na.col = FALSE, style="grid", graph.col=FALSE)
write.table(output, file = pipe("pbcopy"), sep = "\t", row.names = FALSE, col.names = TRUE) 

#' trait_dissimilarity <-  gawdis::gawdis(trait_subset
#'                                        , groups = c(rep(1,2) # size
#'                                                     , rep(2, 2) # diet
#'                                                     , rep(3, 2) # habitat
#'                                                     , rep(4, 2)) #mobility
#'                                        , w.type="optimized")
#' 
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#' 
#' output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)
#' 
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ N39TTP ####
#' #'-----------------------------------------------------------------------------------------------------------------
#' #'
focal_dataset <- "N39TTP"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%
  dplyr::select(contains("SIZE") # size
                ,contains("DIET") # diet (troph 4 and 7 only zero)
                ,T:N # habitat
                ,contains("MOBILITY")) %>% #mobility
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

names_remove<- FD::gowdis(trait_subset) %>%
  as.matrix %>%
  as.data.frame.table %>%
  slice(which(is.na(.$Freq))) %>%
  select(Var2) %>%
  t() %>%
  c() %>%
  unique()

trait_subset2<-trait_subset %>%
  rownames_to_column("species") %>%
  filter(!species %in% names_remove) %>%
  column_to_rownames("species")
#' 

output<-summarytools::dfSummary(trait_subset, varnumbers = FALSE, labels.col=TRUE,na.col = FALSE, style="grid", graph.col=FALSE)
write.table(output, file = pipe("pbcopy"), sep = "\t", row.names = FALSE, col.names = TRUE)  
#' trait_dissimilarity <-  gawdis::gawdis(trait_subset2
#'                                        , groups = c(rep(1,1) # size
#'                                                     , rep(2, 1) # diet
#'                                                     , rep(3, 5) # habitat
#'                                                     , rep(4, 1)) #mobility
#'                                        , w.type="optimized")
#' 
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#' 
#' output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)
#' 
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ N40FMI ####
#' #'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "N40FMI"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%
  rename(Loc_sec = Log_sec) %>%
  dplyr::select(contains("SIZE") # size
                ,contains("Food") # diet (troph 4 and 7 only zero)
                ,contains("Feed") # habitat
                ,contains("Loc")) %>% #mobility
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty


output<-summarytools::dfSummary(trait_subset, varnumbers = FALSE, labels.col=TRUE,na.col = FALSE, style="grid", graph.col=FALSE)
write.table(output, file = pipe("pbcopy"), sep = "\t", row.names = FALSE, col.names = TRUE)  
#' trait_dissimilarity <-  gawdis::gawdis(trait_subset
#'                                        , groups = c(rep(1,2) # size
#'                                                     , rep(2, 2) # diet
#'                                                     , rep(3, 2) # habitat
#'                                                     , rep(4, 2)) #mobility
#'                                        , w.type="optimized")
#' 
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#' 
#' output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)
#' 
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ N41TTP ####
#' #'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "N41TTP"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%
  dplyr::select(contains("SIZE") # size
                ,contains("DIET") # diet (troph 4 and 7 only zero)
                ,T:N # habitat
                ,contains("MOBILITY")) %>% #mobility
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty


names_remove<- FD::gowdis(trait_subset) %>%
  as.matrix %>%
  as.data.frame.table %>%
  slice(which(is.na(.$Freq))) %>%
  select(Var2) %>%
  t() %>%
  c() %>%
  unique()

trait_subset2 <- trait_subset %>%
  rownames_to_column("species") %>%
  filter(!species %in% names_remove) %>%
  column_to_rownames("species")


output<-summarytools::dfSummary(trait_subset, varnumbers = FALSE, labels.col=TRUE,na.col = FALSE, style="grid", graph.col=FALSE)
write.table(output, file = pipe("pbcopy"), sep = "\t", row.names = FALSE, col.names = TRUE)  
#' trait_dissimilarity <-  gawdis::gawdis(trait_subset2
#'                                        , groups = c(rep(1,1) # size
#'                                                     , rep(2, 1) # diet
#'                                                     , rep(3, 5) # habitat
#'                                                     , rep(4, 1)) #mobility
#'                                        , w.type="optimized")
#' 
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#' 
#' output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)
#' 
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ N42TTP_2 ####
#' #'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "N42TTP_2"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%
  dplyr::select(contains("SIZE") # size
                ,contains("DIET") # diet (troph 4 and 7 only zero)
                ,T:N # habitat
                ,contains("MOBILITY")) %>% #mobility
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

names_remove<- FD::gowdis(trait_subset) %>%
  as.matrix %>%
  as.data.frame.table %>%
  slice(which(is.na(.$Freq))) %>%
  select(Var2) %>%
  t() %>%
  c() %>%
  unique()

trait_subset2<-trait_subset %>%
  rownames_to_column("species") %>%
  filter(!species %in% names_remove) %>%
  column_to_rownames("species")


output<-summarytools::dfSummary(trait_subset, varnumbers = FALSE, labels.col=TRUE,na.col = FALSE, style="grid", graph.col=FALSE)
write.table(output, file = pipe("pbcopy"), sep = "\t", row.names = FALSE, col.names = TRUE)  
#' trait_dissimilarity <-  gawdis::gawdis(trait_subset2
#'                                        , groups = c(rep(1,1) # size
#'                                                     , rep(2, 1) # diet
#'                                                     , rep(3, 5) # habitat
#'                                                     , rep(4, 1)) #mobility
#'                                        , w.type="optimized")
#' 
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#' 
#' output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)
#' 
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ N42TTP_3 ####
#' #'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "N42TTP_3"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%
  dplyr::select(contains("SIZE") # size
                ,contains("DIET") # diet (troph 4 and 7 only zero)
                ,T:N # habitat
                ,contains("MOBILITY")) %>% #mobility
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

names_remove<- FD::gowdis(trait_subset) %>%
  as.matrix %>%
  as.data.frame.table %>%
  slice(which(is.na(.$Freq))) %>%
  select(Var2) %>%
  t() %>%
  c() %>%
  unique()

trait_subset2<-trait_subset %>%
  rownames_to_column("species") %>%
  filter(!species %in% names_remove) %>%
  column_to_rownames("species")

output<-summarytools::dfSummary(trait_subset, varnumbers = FALSE, labels.col=TRUE,na.col = FALSE, style="grid", graph.col=FALSE)
write.table(output, file = pipe("pbcopy"), sep = "\t", row.names = FALSE, col.names = TRUE)  
#' trait_dissimilarity <-  gawdis::gawdis(trait_subset2
#'                                        , groups = c(rep(1,1) # size
#'                                                     , rep(2, 1) # diet
#'                                                     , rep(3, 5) # habitat
#'                                                     , rep(4, 1)) #mobility
#'                                        , w.type="optimized")
#' 
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  magrittr::set_colnames(paste0("trait_",1:4))
#' 
#' output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)
#' 
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ N42TTP ####
#' #'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "N42TTP"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%
  dplyr::select(contains("SIZE") # size
                ,contains("DIET") # diet (troph 4 and 7 only zero)
                ,T:N # habitat
                ,contains("MOBILITY")) %>% #mobility
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

names_remove<- FD::gowdis(trait_subset) %>%
  as.matrix %>%
  as.data.frame.table %>%
  slice(which(is.na(.$Freq))) %>%
  select(Var2) %>%
  t() %>%
  c() %>%
  unique()

trait_subset2<-trait_subset %>%
  rownames_to_column("species") %>%
  filter(!species %in% names_remove) %>%
  column_to_rownames("species")


output<-summarytools::dfSummary(trait_subset, varnumbers = FALSE, labels.col=TRUE,na.col = FALSE, style="grid", graph.col=FALSE)
write.table(output, file = pipe("pbcopy"), sep = "\t", row.names = FALSE, col.names = TRUE)  
#' trait_dissimilarity <-  gawdis::gawdis(trait_subset2
#'                                        , groups = c(rep(1,1) # size
#'                                                     , rep(2, 1) # diet
#'                                                     , rep(3, 5) # habitat
#'                                                     , rep(4, 1)) #mobility
#'                                        , w.type="optimized")
#' 
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#' 
#' output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)
#' 
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ N43TTP_1 ####
#' #'-----------------------------------------------------------------------------------------------------------------
#' focal_dataset <- "N43TTP_1"
#' data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
#' data_trait <- data_clean %>%  pluck("trait")
#' trait_subset<- data_trait %>%  
#'   dplyr::select(contains("SIZE") # size
#'                 ,contains("DIET") # diet (troph 4 and 7 only zero)
#'                 ,T:N # habitat
#'                 ,contains("MOBILITY")) %>% #mobility
#'   slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 
#' 
#' names_remove<- FD::gowdis(trait_subset) %>%
#'   as.matrix %>%
#'   as.data.frame.table %>%
#'   slice(which(is.na(.$Freq))) %>%
#'   select(Var2) %>%
#'   t() %>%
#'   c() %>%
#'   unique()
#' 
#' trait_subset2<-trait_subset %>%
#'   rownames_to_column("species") %>%
#'   filter(!species %in% names_remove) %>%
#'   column_to_rownames("species")
#' 
#' trait_dissimilarity <-  gawdis::gawdis(trait_subset2
#'                                        , groups = c(rep(1,1) # size
#'                                                     , rep(2, 1) # diet
#'                                                     , rep(3, 5) # habitat
#'                                                     , rep(4, 1)) #mobility
#'                                        , w.type="optimized")
#' 
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#' 
#' output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)
#' 
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ N43TTP_2 ####
#' #'-----------------------------------------------------------------------------------------------------------------
#' focal_dataset <- "N43TTP_2"
#' data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
#' data_trait <- data_clean %>%  pluck("trait")
#' trait_subset<- data_trait %>%  
#'   dplyr::select(contains("SIZE") # size
#'                 ,contains("DIET") # diet (troph 4 and 7 only zero)
#'                 ,T:N # habitat
#'                 ,contains("MOBILITY")) %>% #mobility
#'   slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 
#' 
#' names_remove<- FD::gowdis(trait_subset) %>%
#'   as.matrix %>%
#'   as.data.frame.table %>%
#'   slice(which(is.na(.$Freq))) %>%
#'   select(Var2) %>%
#'   t() %>%
#'   c() %>%
#'   unique()
#' 
#' trait_subset2<-trait_subset %>%
#'   rownames_to_column("species") %>%
#'   filter(!species %in% names_remove) %>%
#'   column_to_rownames("species")
#' 
#' trait_dissimilarity <-  gawdis::gawdis(trait_subset2
#'                                        , groups = c(rep(1,1) # size
#'                                                     , rep(2, 1) # diet
#'                                                     , rep(3, 5) # habitat
#'                                                     , rep(4, 1)) #mobility
#'                                        , w.type="optimized")
#' 
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#' 
#' output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)
#' 
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ N43TTP_3 ####
#' #'-----------------------------------------------------------------------------------------------------------------
#' focal_dataset <- "N43TTP_3"
#' data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
#' data_trait <- data_clean %>%  pluck("trait")
#' trait_subset<- data_trait %>%  
#'   dplyr::select(contains("SIZE") # size
#'                 ,contains("DIET") # diet (troph 4 and 7 only zero)
#'                 ,T:N # habitat
#'                 ,contains("MOBILITY")) %>% #mobility
#'   slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 
#' 
#' names_remove<- FD::gowdis(trait_subset) %>%
#'   as.matrix %>%
#'   as.data.frame.table %>%
#'   slice(which(is.na(.$Freq))) %>%
#'   select(Var2) %>%
#'   t() %>%
#'   c() %>%
#'   unique()
#' 
#' trait_subset2<-trait_subset %>%
#'   rownames_to_column("species") %>%
#'   filter(!species %in% names_remove) %>%
#'   column_to_rownames("species")
#' 
#' trait_dissimilarity <-  gawdis::gawdis(trait_subset2
#'                                        , groups = c(rep(1,1) # size
#'                                                     , rep(2, 1) # diet
#'                                                     , rep(3, 5) # habitat
#'                                                     , rep(4, 1)) #mobility
#'                                        , w.type="optimized")
#' 
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#' 
#' output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)
#' 
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ N44TBT ####
#' #'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "N44TBT"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%
  dplyr::select(contains("SIZE") # size
                ,contains("DIET") # diet (troph 4 and 7 only zero)
                ,contains("HABITAT") # habitat
                ,contains("MOBILITY")) %>% #mobility
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty


output<-summarytools::dfSummary(trait_subset, varnumbers = FALSE, labels.col=TRUE,na.col = FALSE, style="grid", graph.col=FALSE)
write.table(output, file = pipe("pbcopy"), sep = "\t", row.names = FALSE, col.names = TRUE) 
#' trait_dissimilarity <-  gawdis::gawdis(trait_subset
#'                                        , w.type="optimized")
#' 
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#' 
#' output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)
#' 
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ N44TTP ####
#' #'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "N44TTP"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%
  dplyr::select(PlantHeight_mean # size
                ,SLA_mean # diet (troph 4 and 7 only zero)
                ,e_light:e_temp # habitat
                ,SeedMass_mean) %>% #mobility
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty


names_remove<- FD::gowdis(trait_subset) %>%
  as.matrix %>%
  as.data.frame.table %>%
  slice(which(is.na(.$Freq))) %>%
  select(Var2) %>%
  t() %>%
  c() %>%
  unique()

trait_subset2<-trait_subset %>%
  rownames_to_column("species") %>%
  filter(!species %in% names_remove) %>%
  column_to_rownames("species")


output<-summarytools::dfSummary(trait_subset, varnumbers = FALSE, labels.col=TRUE,na.col = FALSE, style="grid", graph.col=FALSE)
write.table(output, file = pipe("pbcopy"), sep = "\t", row.names = FALSE, col.names = TRUE)
#' trait_dissimilarity <-  gawdis::gawdis(trait_subset2
#'                                        , groups = c(rep(1,1) # size
#'                                                     , rep(2, 1) # diet
#'                                                     , rep(3, 6) # habitat
#'                                                     , rep(4, 1)) #mobility
#'                                        , w.type="optimized")
#' 
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#' 
#' output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)
#' 
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ N47FFI ####
#' #'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "N47FFI"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%
  dplyr::select(contains("SIZE_") # size
                ,contains("DIET") # diet (troph 4 and 7 only zero)
                ,contains("HABITAT") # habitat
                ,contains("MOBILITY")) %>% #mobility
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

output<-summarytools::dfSummary(trait_subset, varnumbers = FALSE, labels.col=TRUE,na.col = FALSE, style="grid", graph.col=FALSE)
write.table(output, file = pipe("pbcopy"), sep = "\t", row.names = FALSE, col.names = TRUE)
#' 
#' trait_dissimilarity <-  gawdis::gawdis(trait_subset
#'                                        , groups = c(rep(1,1) # size
#'                                                     , rep(2, 6) # diet
#'                                                     , rep(3, 3) # habitat
#'                                                     , rep(4, 7)) #mobility
#'                                        , fuzzy=c(2:4)
#'                                        , w.type="optimized")
#' 
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#' 
#' output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)
#' 
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#' #'-----------------------------------------------------------------------------------------------------------------
#' # @ N47TTP_2 ####
#' #'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "N47TTP_2"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%
  dplyr::select(contains("SIZE") # size
                ,contains("DIET") # diet (troph 4 and 7 only zero)
                ,T:N # habitat
                ,contains("MOBILITY")) %>% #mobility
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

names_remove<- FD::gowdis(trait_subset) %>%
  as.matrix %>%
  as.data.frame.table %>%
  slice(which(is.na(.$Freq))) %>%
  select(Var2) %>%
  t() %>%
  c() %>%
  unique()

trait_subset2<-trait_subset %>%
  rownames_to_column("species") %>%
  filter(!species %in% names_remove) %>%
  column_to_rownames("species")


output<-summarytools::dfSummary(trait_subset, varnumbers = FALSE, labels.col=TRUE,na.col = FALSE, style="grid", graph.col=FALSE)
write.table(output, file = pipe("pbcopy"), sep = "\t", row.names = FALSE, col.names = TRUE)
#' trait_dissimilarity <-  gawdis::gawdis(trait_subset2
#'                                        , groups = c(rep(1,1) # size
#'                                                     , rep(2, 1) # diet
#'                                                     , rep(3, 5) # habitat
#'                                                     , rep(4, 1)) #mobility
#'                                        , w.type="optimized")
#' 
#' trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#' 
#' output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)
#' 
#' saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
#' rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ N47TTP ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "N47TTP"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%
  dplyr::select(contains("SIZE") # size
                ,contains("DIET") # diet (troph 4 and 7 only zero)
                ,T:N # habitat
                ,contains("MOBILITY")) %>% #mobility
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

names_remove<- FD::gowdis(trait_subset) %>%
  as.matrix %>%
  as.data.frame.table %>%
  slice(which(is.na(.$Freq))) %>%
  select(Var2) %>%
  t() %>%
  c() %>%
  unique()

trait_subset2<-trait_subset %>%
  rownames_to_column("species") %>%
  filter(!species %in% names_remove) %>%
  column_to_rownames("species")



trait_dissimilarity <-  gawdis::gawdis(trait_subset2
                                       , groups = c(rep(1,1) # size
                                                    , rep(2, 1) # diet
                                                    , rep(3, 5) # habitat
                                                    , rep(4, 1)) #mobility
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ N48FBD ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "N48FBD"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%
  dplyr::select(contains("SIZE") # size
                ,contains("DIET") # diet (troph 4 and 7 only zero)
                ,contains("HABITAT") # habitat
                ,contains("MOBILITY")) %>% #mobility
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , groups = c(rep(1,1) # size
                                                    , rep(2, 1) # diet
                                                    , rep(3, 4) # habitat
                                                    , rep(4, 1)) #mobility
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ N48FMI ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "N48FMI"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%
  dplyr::select(contains("SIZE") # size
                ,contains("DIET") # diet (troph 4 and 7 only zero)
                ,contains("HABITAT") # habitat
                ,contains("MOBILITY")) %>% #mobility
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ N49TBI ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "N49TBI"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%
  dplyr::select(contains("SIZE") # size
                ,contains("DIET") # diet (troph 4 and 7 only zero)
                ,contains("HABITAT") # habitat
                ,contains("MOBILITY")) %>% #mobility
  mutate(across(where(is.character),factor)) %>%
  mutate(MOBILITY_home = ordered(MOBILITY_home, levels=c("small","mid","large"))) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

foo()

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type=ifelse(all(!is.na(trait_subset)),"analytic","optimized"))

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ N49TBT ####
#'-----------------------------------------------------------------------------------------------------------------


focal_dataset <- "N49TBT"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%
  dplyr::select(contains("SIZE") # size
                ,contains("DIET") # diet (troph 4 and 7 only zero)
                ,contains("HABITAT") # habitat
                ,contains("MOBILITY")) %>% #mobility
  mutate(across(where(is.character),factor)) %>%
  mutate(MOBILITY_home = ordered(MOBILITY_home, levels=c("small","mid","large"))) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

foo()

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ N50TFG ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "N50TFG"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%
  dplyr::select(contains("SIZE") # size
                ,contains("DIET") # diet (troph 4 and 7 only zero)
                ,contains("HABITAT") # habitat
                ,contains("MOBILITY")) %>% #mobility
  mutate(across(where(is.character),factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty


foo()


trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ N51TTP ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "N51TTP"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%
  dplyr::select(contains("SIZE") # size
                ,contains("DIET") # diet (troph 4 and 7 only zero)
                ,T:N # habitat
                ,contains("MOBILITY")) %>% #mobility
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty


foo() 

names_remove<- FD::gowdis(trait_subset) %>%
  as.matrix %>%
  as.data.frame.table %>%
  slice(which(is.na(.$Freq))) %>%
  select(Var2) %>%
  t() %>%
  c() %>%
  unique()

trait_subset2<-trait_subset %>%
  rownames_to_column("species") %>%
  filter(!species %in% names_remove) %>%
  column_to_rownames("species")

trait_dissimilarity <-  gawdis::gawdis(trait_subset2
                                       , groups = c(rep(1,1) # size
                                                    , rep(2, 1) # diet
                                                    , rep(3, 5) # habitat
                                                    , rep(4, 1)) #mobility
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ N55TAT ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "N55TAT"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%
  dplyr::select(contains("SIZE") # size
                ,contains("DIET") # diet (troph 4 and 7 only zero)
                ,contains("HABITAT") # habitat
                ,contains("MOBILITY")) %>% #mobility
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

foo() 

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ N55TTP ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "N55TTP"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%
  dplyr::select(contains("SIZE") # size
                ,contains("DIET") # diet (troph 4 and 7 only zero)
                ,T:N # habitat
                ,contains("MOBILITY")) %>% #mobility
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty



names_remove<- FD::gowdis(trait_subset) %>%
  as.matrix %>%
  as.data.frame.table %>%
  slice(which(is.na(.$Freq))) %>%
  select(Var2) %>%
  t() %>%
  c() %>%
  unique()

trait_subset2<-trait_subset %>%
  rownames_to_column("species") %>%
  filter(!species %in% names_remove) %>%
  column_to_rownames("species")

trait_dissimilarity <-  gawdis::gawdis(trait_subset2
                                       , groups = c(rep(1,1) # size
                                                    , rep(2, 1) # diet
                                                    , rep(3, 5) # habitat
                                                    , rep(4, 1)) #mobility
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ N62FMI ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "N62FMI"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%
  dplyr::select(contains("SIZE") # size
                ,contains("DIET") # diet (troph 4 and 7 only zero)
                ,contains("HABITAT") # habitat
                ,contains("MOBILITY")) %>% #mobility
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

foo()

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , groups = c(rep(1,1) # size
                                                    , rep(2, 1) # diet
                                                    , rep(3, 1) # habitat
                                                    , rep(4, 8)) #mobility
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ N67TTP ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "N67TTP"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%
  dplyr::select(Height # size
                ,SLA # diet (troph 4 and 7 only zero)
                ,T:N # habitat
                ,Seed_mass) %>% #mobility
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

foo()

names_remove<- FD::gowdis(trait_subset) %>%
  as.matrix %>%
  as.data.frame.table %>%
  slice(which(is.na(.$Freq))) %>%
  select(Var2) %>%
  t() %>%
  c() %>%
  unique()

trait_subset2<-trait_subset %>%
  rownames_to_column("species") %>%
  filter(!species %in% names_remove) %>%
  column_to_rownames("species")

trait_dissimilarity <-  gawdis::gawdis(trait_subset2
                                       , groups = c(rep(1,1) # size
                                                    , rep(2, 1) # diet
                                                    , rep(3, 5) # habitat
                                                    , rep(4, 1)) #mobility
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ N69FMI ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "N69FMI"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%
  dplyr::select(contains("SIZE") # size
                ,contains("DIET") # diet (troph 4 and 7 only zero)
                ,contains("HABITAT") # habitat
                ,contains("MOBILITY")) %>% #mobility
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

foo()
trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , groups = c(rep(1,1) # size
                                                    , rep(2, 1) # diet
                                                    , rep(3, 1) # habitat
                                                    , rep(4, 4)) #mobility
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ N78TTP ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "N78TTP"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%
  dplyr::select(Height # size
                ,SLA # diet (troph 4 and 7 only zero)
                ,T:N # habitat
                ,Seed_mass) %>% #mobility
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

names_remove<- FD::gowdis(trait_subset) %>%
  as.matrix %>%
  as.data.frame.table %>%
  slice(which(is.na(.$Freq))) %>%
  select(Var2) %>%
  t() %>%
  c() %>%
  unique()

trait_subset2<-trait_subset %>%
  rownames_to_column("species") %>%
  filter(!species %in% names_remove) %>%
  column_to_rownames("species")

trait_dissimilarity <-  gawdis::gawdis(trait_subset2
                                       , groups = c(rep(1,1) # size
                                                    , rep(2, 1) # diet
                                                    , rep(3, 5) # habitat
                                                    , rep(4, 1)) #mobility
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ ossola_ants ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "ossola_ants"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%
  dplyr::select(contains("SIZE"),  # size
                ,contains("DIET_") # diet (troph 4 and 7 only zero)
                ,contains("HABITAT") # habitat
                ,contains("MOBILITY")) %>% #mobility
  dplyr::select(-DIET_Diet) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty


foo()

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , groups = c(rep(1,1) # size
                                                    , rep(2, 6) # diet
                                                    , rep(3, 1) # habitat
                                                    , rep(4, 1)) #mobility
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ penone_orthoptera ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "penone_orthoptera"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%
  dplyr::select(contains("SIZE") # size
                ,contains("DIET") # diet (troph 4 and 7 only zero)
                ,contains("HABITAT") # habitat
                ,contains("dispersal")) %>% #mobility
  select(-diet.from.pyrgus.de, - habitat.from.pyrgus.de) %>%  #remove some traits that are not relevant
  mutate(across(where(is.character),factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty


foo()
trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , groups = c(rep(1,3) # size
                                                    , rep(2, 5) # diet
                                                    , rep(3, 1) # habitat
                                                    , rep(4, 1)) #mobility
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ piano_beetles ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "piano_beetles"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%
  dplyr::select(contains("SIZE") # size
                ,contains("DIET") # diet (troph 4 and 7 only zero)
                ,contains("HABITAT") # habitat
                ,contains("MOBILITY")) %>% #mobility
  mutate(across(where(is.character),factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty


foo()


trait_dissimilarity <-  gawdis::gawdis(trait_subset
                               , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ piano_spiders ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "piano_spiders"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%
  dplyr::select(contains("SIZE") # size
                ,contains("DIET") # diet (troph 4 and 7 only zero)
                ,contains("HABITAT") # habitat
                ,contains("MOBILITY")) %>% #mobility
  mutate(across(where(is.character),factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

foo()

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_cc1_2013_waite_1 ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "Predicts_cc1_2013_waite_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%
  dplyr::select(Body_mass_g # size
                ,Trophic_level # diet (troph 4 and 7 only zero)
                ,Habitat_breadth_IUCN # habitat
                ,Hand.Wing.Index) %>% #mobility
  mutate(across(where(is.character),factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

foo()
trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_db1_2010_dures_1 ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "Predicts_db1_2010_dures_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%
  dplyr::select(Body_mass_g # size
                ,Trophic_level # diet (troph 4 and 7 only zero)
                ,Habitat_breadth_IUCN # habitat
                ,Hand.Wing.Index) %>% #mobility
  mutate(across(where(is.character),factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_di1_2004_naidoo_1 ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "Predicts_di1_2004_naidoo_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%
  dplyr::select(Body_mass_g # size
                ,Trophic_level # diet (troph 4 and 7 only zero)
                ,Habitat_breadth_IUCN # habitat
                ,Hand.Wing.Index) %>% #mobility
  mutate(across(where(is.character),factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_di1_2010_milder_2 ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "Predicts_di1_2010_milder_2"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%
  dplyr::select(Body_mass_g # size
                ,Trophic_level # diet (troph 4 and 7 only zero)
                ,Habitat_breadth_IUCN # habitat
                ,Hand.Wing.Index) %>% #mobility
  mutate(across(where(is.character),factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_di1_2011_dawson_1 ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "Predicts_di1_2011_dawson_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%
  dplyr::select(Body_mass_g # size
                ,Trophic_level # diet (troph 4 and 7 only zero)
                ,Habitat_breadth_IUCN # habitat
                ,Hand.Wing.Index) %>% #mobility
  mutate(across(where(is.character),factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_di1_2011_neuschulz_1 ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "Predicts_di1_2011_neuschulz_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%
  dplyr::select(Body_mass_g # size
                ,Trophic_level # diet (troph 4 and 7 only zero)
                ,Habitat_breadth_IUCN # habitat
                ,Hand.Wing.Index) %>% #mobility
  mutate(across(where(is.character),factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_di1_2012_reid_1 ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "Predicts_di1_2012_reid_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%
  dplyr::select(Body_mass_g # size
                ,Trophic_level # diet (troph 4 and 7 only zero)
                ,Habitat_breadth_IUCN # habitat
                ,Hand.Wing.Index) %>% #mobility
  mutate(across(where(is.character),factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_di1_2013_azhar_1 ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "Predicts_di1_2013_azhar_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%
  dplyr::select(Body_mass_g # size
                ,Trophic_level # diet (troph 4 and 7 only zero)
                ,Habitat_breadth_IUCN # habitat
                ,Hand.Wing.Index) %>% #mobility
  mutate(across(where(is.character),factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_di1_2013_de_lima_1 ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "Predicts_di1_2013_de_lima_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%
  dplyr::select(Body_mass_g # size
                ,Trophic_level # diet (troph 4 and 7 only zero)
                ,Habitat_breadth_IUCN # habitat
                ,Hand.Wing.Index) %>% #mobility
  mutate(across(where(is.character),factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_dl1_2009_woinarski_1 ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "Predicts_dl1_2009_woinarski_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%
  dplyr::select(Body_mass_g # size
                ,Trophic_level # diet (troph 4 and 7 only zero)
                ,Habitat_breadth_IUCN # habitat
                ,Hand.Wing.Index) %>% #mobility
  mutate(across(where(is.character),factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_dl1_2010_proenca_2 ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "Predicts_dl1_2010_proenca_2"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%
  dplyr::select(Body_mass_g # size
                ,Trophic_level # diet (troph 4 and 7 only zero)
                ,Habitat_breadth_IUCN # habitat
                ,Hand.Wing.Index) %>% #mobility
  mutate(across(where(is.character),factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_dl1_2011_mallari_1 ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "Predicts_dl1_2011_mallari_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%
  dplyr::select(Body_mass_g # size
                ,Trophic_level # diet (troph 4 and 7 only zero)
                ,Habitat_breadth_IUCN # habitat
                ,Hand.Wing.Index) %>% #mobility
  mutate(across(where(is.character),factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_dl1_2011_moreno_mateos_1 ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "Predicts_dl1_2011_moreno_mateos_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%
  dplyr::select(Body_mass_g # size
                ,Trophic_level # diet (troph 4 and 7 only zero)
                ,Habitat_breadth_IUCN # habitat
                ,Hand.Wing.Index) %>% #mobility
  mutate(across(where(is.character),factor)) %>%
  select(-Trophic_level)  %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_dl1_2012_dallimer_1 ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "Predicts_dl1_2012_dallimer_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%
  dplyr::select(Body_mass_g # size
                ,Trophic_level # diet (troph 4 and 7 only zero)
                ,Habitat_breadth_IUCN # habitat
                ,Hand.Wing.Index) %>% #mobility
  mutate(across(where(is.character),factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_dl1_2013_bartolommei_1 ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "Predicts_dl1_2013_bartolommei_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%
  dplyr::select(Body_mass_g # size
                ,Trophic_level # diet (troph 4 and 7 only zero)
                ,Habitat_breadth_IUCN # habitat
                ,Hand.Wing.Index) %>% #mobility
  mutate(across(where(is.character),factor)) %>%
  dplyr::select(-Trophic_level) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_dl1_2013_de_thoisy_1 ####
#'-----------------------------------------------------------------------------------------------------------------
# focal_dataset <- "Predicts_dl1_2013_de_thoisy_1"
# data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
# data_trait <- data_clean %>%  pluck("trait")
# trait_subset<- data_trait %>%
#   dplyr::select(Body_mass_g # size
#                 ,Trophic_level # diet (troph 4 and 7 only zero)
#                 ,Habitat_breadth_IUCN # habitat
#                 ,Hand.Wing.Index) %>% #mobility
#   mutate(across(where(is.character),factor)) %>%
#   slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty
#
# trait_dissimilarity <-  gawdis::gawdis(trait_subset
#                                        , w.type="optimized")
#
# trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
#
# output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)
#
# saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
# rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_gp1_2007_kutt_2 ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "Predicts_gp1_2007_kutt_2"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%
  dplyr::select(Body_mass_g # size
                ,Trophic_level # diet (troph 4 and 7 only zero)
                ,Habitat_breadth_IUCN # habitat
                ,Hand.Wing.Index) %>% #mobility
  mutate(across(where(is.character),factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_hb1_2001_aumann_1 ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "Predicts_hb1_2001_aumann_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%
  dplyr::select(Body_mass_g # size
                ,Trophic_level # diet (troph 4 and 7 only zero)
                ,Habitat_breadth_IUCN # habitat
                ,Hand.Wing.Index) %>% #mobility
  mutate(across(where(is.character),factor)) %>%
  dplyr::select(-Trophic_level) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_hp1_2006_wunderle_1 ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "Predicts_hp1_2006_wunderle_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%
  dplyr::select(Body_mass_g # size
                ,Trophic_level # diet (troph 4 and 7 only zero)
                ,Habitat_breadth_IUCN # habitat
                ,Hand.Wing.Index) %>% #mobility
  mutate(across(where(is.character),factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_hp1_2007_borges_1 ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "Predicts_hp1_2007_borges_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%
  dplyr::select(Body_mass_g # size
                ,Trophic_level # diet (troph 4 and 7 only zero)
                ,Habitat_breadth_IUCN # habitat
                ,Hand.Wing.Index) %>% #mobility
  mutate(across(where(is.character),factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_hp1_2007_ranganathan_1 ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "Predicts_hp1_2007_ranganathan_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%
  dplyr::select(Body_mass_g # size
                ,Trophic_level # diet (troph 4 and 7 only zero)
                ,Habitat_breadth_IUCN # habitat
                ,Hand.Wing.Index) %>% #mobility
  mutate(across(where(is.character),factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_hp1_2007_shahabuddin_1 ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "Predicts_hp1_2007_shahabuddin_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%
  dplyr::select(Body_mass_g # size
                ,Trophic_level # diet (troph 4 and 7 only zero)
                ,Habitat_breadth_IUCN # habitat
                ,Hand.Wing.Index) %>% #mobility
  mutate(across(where(is.character),factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_hp1_2008_farwig_1 ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "Predicts_hp1_2008_farwig_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%
  dplyr::select(Body_mass_g # size
                ,Trophic_level # diet (troph 4 and 7 only zero)
                ,Habitat_breadth_IUCN # habitat
                ,Hand.Wing.Index) %>% #mobility
  mutate(across(where(is.character),factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_hp1_2008_ranganathan_1 ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "Predicts_hp1_2008_ranganathan_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%
  dplyr::select(Body_mass_g # size
                ,Trophic_level # diet (troph 4 and 7 only zero)
                ,Habitat_breadth_IUCN # habitat
                ,Hand.Wing.Index) %>% #mobility
  mutate(across(where(is.character),factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_hp1_2009_kessler_5 ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "Predicts_hp1_2009_kessler_5"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%
  dplyr::select(Body_mass_g # size
                ,Trophic_level # diet (troph 4 and 7 only zero)
                ,Habitat_breadth_IUCN # habitat
                ,Hand.Wing.Index) %>% #mobility
  mutate(across(where(is.character),factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_hp1_2010_bicknell_1 ####
#'-----------------------------------------------------------------------------------------------------------------

#setwd("/Volumes/Graco-Roza SS/Work/Ongoing manuscripts/HIATES")
focal_dataset <- "Predicts_hp1_2010_bicknell_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%  
  dplyr::select(Body_mass_g # size
                ,Trophic_level # diet (troph 4 and 7 only zero)
                ,Habitat_breadth_IUCN # habitat
                ,Hand.Wing.Index) %>% #mobility
  mutate(across(where(is.character),factor)) %>% 
  dplyr::select(-Trophic_level) %>% 
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:3))
plot(trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_hp1_2010_bicknell_2 ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "Predicts_hp1_2010_bicknell_2"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%  
  dplyr::select(Body_mass_g # size
                ,Trophic_level # diet (troph 4 and 7 only zero)
                ,Habitat_breadth_IUCN # habitat
                ,Hand.Wing.Index) %>% #mobility
  mutate(across(where(is.character),factor)) %>% 
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_hw1_2005_baldi_1 ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "Predicts_hw1_2005_baldi_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%  
  dplyr::select(Body_mass_g # size
                ,Trophic_level # diet (troph 4 and 7 only zero)
                ,Habitat_breadth_IUCN # habitat
                ,Hand.Wing.Index) %>% #mobility
  mutate(across(where(is.character),factor)) %>% 
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_hw1_2007_chapman_1 ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "Predicts_hw1_2007_chapman_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%  
  dplyr::select(Body_mass_g # size
                ,Trophic_level # diet (troph 4 and 7 only zero)
                ,Habitat_breadth_IUCN # habitat
                ,Hand.Wing.Index) %>% #mobility
  mutate(across(where(is.character),factor)) %>% 
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_hw1_2008_lantschner_1 ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "Predicts_hw1_2008_lantschner_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%  
  dplyr::select(Body_mass_g # size
                ,Trophic_level # diet (troph 4 and 7 only zero)
                ,Habitat_breadth_IUCN # habitat
                ,Hand.Wing.Index) %>% #mobility
  mutate(across(where(is.character),factor)) %>% 
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_hw1_2011_cerezo_1 ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "Predicts_hw1_2011_cerezo_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%  
  dplyr::select(Body_mass_g # size
                ,Trophic_level # diet (troph 4 and 7 only zero)
                ,Habitat_breadth_IUCN # habitat
                ,Hand.Wing.Index) %>% #mobility
  mutate(across(where(is.character),factor)) %>% 
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_hw1_2012_naoe_1 ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "Predicts_hw1_2012_naoe_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%  
  dplyr::select(Body_mass_g # size
                ,Trophic_level # diet (troph 4 and 7 only zero)
                ,Habitat_breadth_IUCN # habitat
                ,Hand.Wing.Index) %>% #mobility
  mutate(across(where(is.character),factor)) %>% 
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_hw1_2012_naoe_2 ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "Predicts_hw1_2012_naoe_2"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%  
  dplyr::select(Body_mass_g # size
                ,Trophic_level # diet (troph 4 and 7 only zero)
                ,Habitat_breadth_IUCN # habitat
                ,Hand.Wing.Index) %>% #mobility
  mutate(across(where(is.character),factor)) %>% 
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_hz1_2012_kutt_1 ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "Predicts_hz1_2012_kutt_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%  
  dplyr::select(Body_mass_g # size
                ,Trophic_level # diet (troph 4 and 7 only zero)
                ,Habitat_breadth_IUCN # habitat
                ,Hand.Wing.Index) %>% #mobility
  mutate(across(where(is.character),factor)) %>% 
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_jd1_2002_pearman_1 ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "Predicts_jd1_2002_pearman_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%  
  dplyr::select(Body_mass_g # size
                ,Trophic_level # diet (troph 4 and 7 only zero)
                ,Habitat_breadth_IUCN # habitat
                ,Hand.Wing.Index) %>% #mobility
  mutate(across(where(is.character),factor)) %>% 
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_jd1_2010_wang_1 ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "Predicts_jd1_2010_wang_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%  
  dplyr::select(Body_mass_g # size
                ,Trophic_level # diet (troph 4 and 7 only zero)
                ,Habitat_breadth_IUCN # habitat
                ,Hand.Wing.Index) %>% #mobility
  mutate(across(where(is.character),factor)) %>% 
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_ks1_2005_pons_1 ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "Predicts_ks1_2005_pons_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%  
  dplyr::select(Body_mass_g # size
                ,Trophic_level # diet (troph 4 and 7 only zero)
                ,Habitat_breadth_IUCN # habitat
                ,Hand.Wing.Index) %>% #mobility
  mutate(across(where(is.character),factor)) %>% 
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_ks1_2009_suarez_rubio_1 ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "Predicts_ks1_2009_suarez_rubio_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%  
  dplyr::select(Body_mass_g # size
                ,Trophic_level # diet (troph 4 and 7 only zero)
                ,Habitat_breadth_IUCN # habitat
                ,Hand.Wing.Index) %>% #mobility
  mutate(across(where(is.character),factor)) %>% 
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_lk1_2009_hayward_1 ####
#'-----------------------------------------------------------------------------------------------------------------
# focal_dataset <- "Predicts_lk1_2009_hayward_1"
# data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
# data_trait <- data_clean %>%  pluck("trait")
# trait_subset<- data_trait %>%  
#   dplyr::select(Body_mass_g # size
#                 ,Trophic_level # diet (troph 4 and 7 only zero)
#                 ,Habitat_breadth_IUCN # habitat
#                 ,Hand.Wing.Index) %>% #mobility
#   mutate(across(where(is.character),factor)) %>% 
#   dplyr::select(-Trophic_level) %>% 
#   slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 
# 
# trait_dissimilarity <-  gawdis::gawdis(trait_subset
#                                        , w.type="optimized")
# 
# trait_syndrome <- trait_dissimilarity %>% cmdscale(k=3) %>%  set_colnames(paste0("trait_",1:4))
# 
# output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)
# 
# saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_mh1_2010_sheldon_1 ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "Predicts_mh1_2010_sheldon_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%  
  dplyr::select(Body_mass_g # size
                ,Trophic_level # diet (troph 4 and 7 only zero)
                ,Habitat_breadth_IUCN # habitat
                ,Hand.Wing.Index) %>% #mobility
  mutate(across(where(is.character),factor)) %>% 
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_mh1_2011_phalan_1 ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "Predicts_mh1_2011_phalan_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%  
  dplyr::select(Body_mass_g # size
                ,Trophic_level # diet (troph 4 and 7 only zero)
                ,Habitat_breadth_IUCN # habitat
                ,Hand.Wing.Index) %>% #mobility
  mutate(across(where(is.character),factor)) %>% 
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_mj1_2008_munyekenye_1 ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "Predicts_mj1_2008_munyekenye_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%  
  dplyr::select(Body_mass_g # size
                ,Trophic_level # diet (troph 4 and 7 only zero)
                ,Habitat_breadth_IUCN # habitat
                ,Hand.Wing.Index) %>% #mobility
  mutate(across(where(is.character),factor)) %>% 
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_mj1_2009_lehouck_1 ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "Predicts_mj1_2009_lehouck_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%  
  dplyr::select(Body_mass_g # size
                ,Trophic_level # diet (troph 4 and 7 only zero)
                ,Habitat_breadth_IUCN # habitat
                ,Hand.Wing.Index) %>% #mobility
  mutate(across(where(is.character),factor)) %>% 
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_mj1_2009_lehouck_2 ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "Predicts_mj1_2009_lehouck_2"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%  
  dplyr::select(Body_mass_g # size
                ,Trophic_level # diet (troph 4 and 7 only zero)
                ,Habitat_breadth_IUCN # habitat
                ,Hand.Wing.Index) %>% #mobility
  mutate(across(where(is.character),factor)) %>% 
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_mj1_2009_lehouck_3 ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "Predicts_mj1_2009_lehouck_3"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%  
  dplyr::select(Body_mass_g # size
                ,Trophic_level # diet (troph 4 and 7 only zero)
                ,Habitat_breadth_IUCN # habitat
                ,Hand.Wing.Index) %>% #mobility
  mutate(across(where(is.character),factor)) %>% 
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_mj1_2009_lehouck_4 ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "Predicts_mj1_2009_lehouck_4"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%  
  dplyr::select(Body_mass_g # size
                ,Trophic_level # diet (troph 4 and 7 only zero)
                ,Habitat_breadth_IUCN # habitat
                ,Hand.Wing.Index) %>% #mobility
  mutate(across(where(is.character),factor)) %>% 
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_mj1_2009_lehouck_5 ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "Predicts_mj1_2009_lehouck_5"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%  
  dplyr::select(Body_mass_g # size
                ,Trophic_level # diet (troph 4 and 7 only zero)
                ,Habitat_breadth_IUCN # habitat
                ,Hand.Wing.Index) %>% #mobility
  mutate(across(where(is.character),factor)) %>% 
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_mj1_2013_ndanganga_1 ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "Predicts_mj1_2013_ndanganga_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%  
  dplyr::select(Body_mass_g # size
                ,Trophic_level # diet (troph 4 and 7 only zero)
                ,Habitat_breadth_IUCN # habitat
                ,Hand.Wing.Index) %>% #mobility
  mutate(across(where(is.character),factor)) %>% 
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_mj1_2013_ndanganga_2 ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "Predicts_mj1_2013_ndanganga_2"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%  
  dplyr::select(Body_mass_g # size
                ,Trophic_level # diet (troph 4 and 7 only zero)
                ,Habitat_breadth_IUCN # habitat
                ,Hand.Wing.Index) %>% #mobility
  mutate(across(where(is.character),factor)) %>% 
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_mj1_2013_reynolds_1 ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "Predicts_mj1_2013_reynolds_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%  
  dplyr::select(Body_mass_g # size
                ,Trophic_level # diet (troph 4 and 7 only zero)
                ,Habitat_breadth_IUCN # habitat
                ,Hand.Wing.Index) %>% #mobility
  mutate(across(where(is.character),factor)) %>% 
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_sc1_2005_marsh_1 ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "Predicts_sc1_2005_marsh_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%  
  dplyr::select(Body_mass_g # size
                ,Trophic_level # diet (troph 4 and 7 only zero)
                ,Habitat_breadth_IUCN # habitat
                ,Hand.Wing.Index) %>% #mobility
  mutate(across(where(is.character),factor)) %>% 
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_sc1_2010_marsh_2 ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "Predicts_sc1_2010_marsh_2"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%  
  dplyr::select(Body_mass_g # size
                ,Trophic_level # diet (troph 4 and 7 only zero)
                ,Habitat_breadth_IUCN # habitat
                ,Hand.Wing.Index) %>% #mobility
  mutate(across(where(is.character),factor)) %>% 
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_sc1_2010_rey_benayas_1 ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "Predicts_sc1_2010_rey_benayas_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%  
  dplyr::select(Body_mass_g # size
                ,Trophic_level # diet (troph 4 and 7 only zero)
                ,Habitat_breadth_IUCN # habitat
                ,Hand.Wing.Index) %>% #mobility
  mutate(across(where(is.character),factor)) %>% 
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_sc1_2011_stouffer_1 ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "Predicts_sc1_2011_stouffer_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%  
  dplyr::select(Body_mass_g # size
                ,Trophic_level # diet (troph 4 and 7 only zero)
                ,Habitat_breadth_IUCN # habitat
                ,Hand.Wing.Index) %>% #mobility
  mutate(across(where(is.character),factor)) %>% 
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_sc2_2012_santana_1 ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "Predicts_sc2_2012_santana_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%  
  dplyr::select(Body_mass_g # size
                ,Trophic_level # diet (troph 4 and 7 only zero)
                ,Habitat_breadth_IUCN # habitat
                ,Hand.Wing.Index) %>% #mobility
  mutate(across(where(is.character),factor)) %>% 
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_se2_2010_mc_carthy_1 ####
#'-----------------------------------------------------------------------------------------------------------------
# focal_dataset <- "Predicts_se2_2010_mc_carthy_1"
# data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
# data_trait <- data_clean %>%  pluck("trait")
# trait_subset<- data_trait %>%  
#   dplyr::select(Body_mass_g # size
#                 ,Trophic_level # diet (troph 4 and 7 only zero)
#                 ,Habitat_breadth_IUCN # habitat
#                 ,Hand.Wing.Index) %>% #mobility
#   mutate(across(where(is.character),factor)) %>% 
#   slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 
# 
# trait_dissimilarity <-  gawdis::gawdis(trait_subset
#                                        , w.type="optimized")
# 
# trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
# 
# output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)
# 
# saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
# rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_se2_2013_brandt_1 ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "Predicts_se2_2013_brandt_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%  
  dplyr::select(Body_mass_g # size
                ,Trophic_level # diet (troph 4 and 7 only zero)
                ,Habitat_breadth_IUCN # habitat
                ,Hand.Wing.Index) %>% #mobility
  mutate(across(where(is.character),factor)) %>% 
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_se2_2013_hassan_1 ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "Predicts_se2_2013_hassan_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%  
  dplyr::select(Body_mass_g # size
                ,Trophic_level # diet (troph 4 and 7 only zero)
                ,Habitat_breadth_IUCN # habitat
                ,Hand.Wing.Index) %>% #mobility
  mutate(across(where(is.character),factor)) %>% 
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_sh1_2012_ims_1 ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "Predicts_sh1_2012_ims_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%  
  dplyr::select(Body_mass_g # size
                ,Trophic_level # diet (troph 4 and 7 only zero)
                ,Habitat_breadth_IUCN # habitat
                ,Hand.Wing.Index) %>% #mobility
  mutate(across(where(is.character),factor)) %>% 
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_tn1_2007_o_dea_1 ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "Predicts_tn1_2007_o_dea_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%  
  dplyr::select(Body_mass_g # size
                ,Trophic_level # diet (troph 4 and 7 only zero)
                ,Habitat_breadth_IUCN # habitat
                ,Hand.Wing.Index) %>% #mobility
  mutate(across(where(is.character),factor)) %>% 
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_vk1_2007_st_laurent_1 ####
#'-----------------------------------------------------------------------------------------------------------------
# focal_dataset <- "Predicts_vk1_2007_st_laurent_1"
# data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
# data_trait <- data_clean %>%  pluck("trait")
# trait_subset<- data_trait %>%  
#   dplyr::select(Body_mass_g # size
#                 ,Trophic_level # diet (troph 4 and 7 only zero)
#                 ,Habitat_breadth_IUCN # habitat
#                 ,Hand.Wing.Index) %>% #mobility
#   mutate(across(where(is.character),factor)) %>% 
#   slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 
# 
# trait_dissimilarity <-  gawdis::gawdis(trait_subset
#                                        , w.type="optimized")
# 
# trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
# 
# output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)
# 
# saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
# rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_vk1_2007_st_laurent_3 ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "Predicts_vk1_2007_st_laurent_3"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%  
  dplyr::select(Body_mass_g # size
                ,Trophic_level # diet (troph 4 and 7 only zero)
                ,Habitat_breadth_IUCN # habitat
                ,Hand.Wing.Index) %>% #mobility
  mutate(across(where(is.character),factor)) %>% 
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_vk1_2011_edenius_1 ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "Predicts_vk1_2011_edenius_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%  
  dplyr::select(Body_mass_g # size
                ,Trophic_level # diet (troph 4 and 7 only zero)
                ,Habitat_breadth_IUCN # habitat
                ,Hand.Wing.Index) %>% #mobility
  mutate(across(where(is.character),factor)) %>% 
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_vk1_2011_zimmerman_1 ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "Predicts_vk1_2011_zimmerman_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%  
  dplyr::select(Body_mass_g # size
                ,Trophic_level # diet (troph 4 and 7 only zero)
                ,Habitat_breadth_IUCN # habitat
                ,Hand.Wing.Index) %>% #mobility
  mutate(across(where(is.character),factor)) %>% 
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_vk1_2012_otto_1 ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "Predicts_vk1_2012_otto_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%  
  dplyr::select(Body_mass_g # size
                ,Trophic_level # diet (troph 4 and 7 only zero)
                ,Habitat_breadth_IUCN # habitat
                ,Hand.Wing.Index) %>% #mobility
  mutate(across(where(is.character),factor)) %>% 
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ Predicts_vk1_2012_otto_2 ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "Predicts_vk1_2012_otto_2"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%  
  dplyr::select(Body_mass_g # size
                ,Trophic_level # diet (troph 4 and 7 only zero)
                ,Habitat_breadth_IUCN # habitat
                ,Hand.Wing.Index) %>% #mobility
  mutate(across(where(is.character),factor)) %>% 
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ raine_beetles ####
#'-----------------------------------------------------------------------------------------------------------------
# focal_dataset <- "raine_beetles"
# data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
# data_trait <- data_clean %>%  pluck("trait")
# trait_subset<- data_trait %>%  
#   dplyr::select(Body_mass_g # size
#                 ,Trophic_level # diet (troph 4 and 7 only zero)
#                 ,Habitat_breadth_IUCN # habitat
#                 ,Hand.Wing.Index) %>% #mobility
#   mutate(across(where(is.character),factor)) %>% 
#   slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 
# 
# trait_dissimilarity <-  gawdis::gawdis(trait_subset
#                                        , w.type="optimized")
# 
# trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))
# 
# output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)
# 
# saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
# rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ ribera_beetles ####
#'-----------------------------------------------------------------------------------------------------------------
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "ribera_beetles"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%  
  dplyr::select(contains("SIZE") # size
                ,contains("DIET") # diet (troph 4 and 7 only zero)
                ,contains("HABITAT") # habitat
                ,contains("MOBILITY")) %>% #mobility
  mutate(across(where(is.character),factor)) %>% 
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 

foo()
trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ romero_waterinsects ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "romero_waterinsects"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%  
  dplyr::select(contains("SIZE") # size
                ,contains("DIET") # diet (troph 4 and 7 only zero)
                ,contains("HABITAT") # habitat
                ,contains("MOBILITY")) %>% #mobility
  mutate(across(where(is.character),factor)) %>% 
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 


foo()
trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       ,groups=c(rep(1,1),rep(2,1),rep(3,1),rep(4,7))
                                       ,fuzzy=4
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ S29TTP ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "S29TTP"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%  
  dplyr::select(contains("SIZE") # size
                ,contains("DIET") # diet (troph 4 and 7 only zero)
                ,StemDens.mean # habitat
                ,contains("MOBILITY")) %>% #mobility
  mutate(across(where(is.character),factor)) %>% 
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 


foo()
trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       ,groups=c(rep(1,1),rep(2,1),rep(3,1),rep(4,1))
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ S47TTP ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "S47TTP"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%  
  dplyr::select(Height # size
                ,SLA # diet (troph 4 and 7 only zero)
                ,T:N # habitat
                ,Seed_mass) %>% #mobility
  mutate(across(where(is.character),factor)) %>% 
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       ,groups=c(rep(1,1),rep(2,1),rep(3,5),rep(4,1))
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ shieh_waterinsects ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "shieh_waterinsects"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%  
  dplyr::select(Size5:Size50 # size
                ,FP:Mai # diet (troph 4 and 7 only zero)
                ,SW:PA # habitat
                ,SL:SP) %>% #mobility
  mutate(across(where(is.character),factor)) %>% 
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 

foo()


trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       ,groups=c(rep(1,5),rep(2,6),rep(3,5),rep(4,4))
                                       ,fuzzy=c(1:4)
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())

#'-----------------------------------------------------------------------------------------------------------------
# @ stavert_bees ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "stavert_bees"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%  
  dplyr::select(contains("SIZE") # size
                ,contains("DIET") # diet (troph 4 and 7 only zero)
                ,contains("HABITAT") # habitat
                ,contains("MOBILITY")) %>% #mobility
  mutate(across(where(is.character),factor)) %>% 
  mutate(across(contains("MOBILITY"), ordered)) %>% 
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 

foo()
trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ toxywa_beetles_nonsaproxylic ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "toxywa_beetles_nonsaproxylic"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%  
  dplyr::select(contains("SIZE") # size
                ,contains("DIET") # diet 
                ,contains("HABITAT") # habitat
                #,contains("MOBILITY") #mobility -- only one class
                ) %>% 
  mutate(across(where(is.character),factor)) %>% 
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 

foo()


trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ toxywa_beetles_saproxylic ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "toxywa_beetles_saproxylic"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%  
  dplyr::select(contains("SIZE") # size
                ,contains("DIET") # diet (troph 4 and 7 only zero)
                ,contains("HABITAT") # habitat
                ,contains("MOBILITY")) %>% #mobility
  mutate(across(where(is.character),factor)) %>% 
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 

foo()

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ vanklink_bees ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "vanklink_bees"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%  
  dplyr::select(contains("SIZE") # size
                ,contains("feeding") # diet (troph 4 and 7 only zero)
                ,contains("Nesting") # habitat
                ,contains("flight")) %>% #mobility
  mutate(across(where(is.character),factor)) %>% 
  mutate(across(contains("feeding"), ordered)) %>% 
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 

foo()

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       ,groups=c(rep(1,4),2,3,4)
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ vanklink_carabids ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "vanklink_carabids"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%  
  dplyr::select(contains("SIZE") # size
                ,contains("trophic") # diet (troph 4 and 7 only zero)
                ,contains("HABITAT") # habitat
                ,contains("wing")) %>% #mobility
  mutate(across(where(is.character),factor)) %>% 
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 

foo()

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       ,groups=c(rep(1,2),2,3,4)
                                       , w.type="optimized")



trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ vanklink_hoverflies ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "vanklink_hoverflies"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%  
  dplyr::select(contains("SIZE") # size
                ,contains("feeding") # diet (troph 4 and 7 only zero)
                ,contains("substrate") # habitat
                ,contains("wing")) %>% #mobility
  mutate(across(where(is.character),factor)) %>% 
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 

foo()

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       ,groups=c(rep(1,2),2,3,4)
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ vanklink_staphylinids ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "vanklink_staphylinids"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%  
  dplyr::select(contains("length") # size
                ,contains("feeding") # diet (troph 4 and 7 only zero)
                ,contains("humidity") # habitat
                ,contains("dispersal")) %>% #mobility
  mutate(across(where(is.character),factor)) %>% 
  mutate(across(contains("dispersal"), ordered)) %>% 
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 

foo()

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
#'-----------------------------------------------------------------------------------------------------------------
# @ yates_ants ####
#'-----------------------------------------------------------------------------------------------------------------
focal_dataset <- "yates_ants"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset<- data_trait %>%  
  dplyr::select(contains("thorax") # size
                ,contains("ml") # diet (troph 4 and 7 only zero)
                ,contains("sl") # habitat
                ,contains("femur")) %>% #mobility
  slice(which(apply(., 1, function(x) !all(is.na(x)))))  #remove species where all traits are empty 


foo() 

trait_dissimilarity <-  gawdis::gawdis(trait_subset
                                       , w.type="optimized")

trait_syndrome <- trait_dissimilarity %>% cmdscale(k=4) %>%  set_colnames(paste0("trait_",1:4))

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome)

saveRDS(output, glue::glue("betadiv_input/{focal_dataset}_beta_Input.rds"))
rm(list=ls())
