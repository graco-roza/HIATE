#'##############################################################################
#  Dataset 1: barbaro_birds_1 ####
#'##############################################################################
focal_dataset <- "barbaro_birds_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(mass, diet, forag, range) %>%
  mutate(across(everything(), ordered)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis::gawdis(trait_subset, ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(comm = data_clean$comm,
               trait_syndrome = trait_syndrome,
               quality = quality,
               trait_dissimilarity = trait_dissimilarity)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 2: barbaro_birds_2 ####
#'##############################################################################
focal_dataset <- "barbaro_birds_2"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(mass, diet, forag, move) %>%
  mutate(across(everything(), ordered)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis::gawdis(trait_subset, ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(comm = data_clean$comm,
               trait_syndrome = trait_syndrome,
               quality = quality,
               trait_dissimilarity = trait_dissimilarity)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 3: barber_beetles ####
#'##############################################################################
focal_dataset <- "barber_beetles"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(SIZE_mean.body.length,
                DIET_diet,
                HABITAT_grasslands, HABITAT_cultivated.fields, HABITAT_open.forests.and.other.open.habitats,
                MOBILITY_wing.morphology) %>%
  mutate(across(c(DIET_diet, MOBILITY_wing.morphology), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis::gawdis(trait_subset,
                                      groups = c(1, 2, 3, 3, 3, 4),
                                      w.type = "optimized",
                                      ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(comm = data_clean$comm,
               trait_syndrome = trait_syndrome,
               quality = quality,
               trait_dissimilarity = trait_dissimilarity)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 4: bdm_birds ####
#'##############################################################################
focal_dataset <- "bdm_birds"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis::gawdis(trait_subset,
                                      groups = c(1, rep(2, 10), rep(3, 5), 4),
                                      fuzzy = c(2, 3),
                                      w.type = "optimized",
                                      ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(comm = data_clean$comm %>% dplyr::select(any_of(rownames(trait_syndrome))),
               trait_syndrome = trait_syndrome,
               quality = quality,
               trait_dissimilarity = trait_dissimilarity)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 5: bdm_butterflies ####
#'##############################################################################
focal_dataset <- "bdm_butterflies"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis::gawdis(trait_subset,
                                      groups = c(1, 2, rep(3, 24), rep(4, 4)),
                                      w.type = "optimized",
                                      ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(comm = data_clean$comm %>% dplyr::select(any_of(rownames(trait_syndrome))),
               trait_syndrome = trait_syndrome,
               quality = quality,
               trait_dissimilarity = trait_dissimilarity)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))
#'##############################################################################
#  Dataset 6: bettergardens_bowltrap_apidae ####
#'##############################################################################
focal_dataset <- "bettergardens_bowltrap_apidae"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(ITD_mean_ffww, feeding_spec, habitat_generalism, feeddist_km_ffww) %>%
  mutate(across(where(is.character), factor)) %>%
  mutate(across(c(feeding_spec, habitat_generalism), ordered)) %>%
  slice(which(apply(., 1, function(x) all(!is.na(x)))))

trait_dissimilarity <- gawdis(trait_subset, w.type = "analytic", ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(comm = data_clean$comm %>% dplyr::select(any_of(rownames(trait_syndrome))),
               trait_syndrome = trait_syndrome,
               quality = quality,
               trait_dissimilarity = trait_dissimilarity)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 7: bettergardens_bowltrap_carabidae ####
#'##############################################################################
focal_dataset <- "bettergardens_bowltrap_carabidae"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(bodylength, feeding_guild, habitatgeneralism, wing_form) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis(trait_subset, w.type = "optimized", ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(comm = data_clean$comm %>% dplyr::select(any_of(rownames(trait_syndrome))),
               trait_syndrome = trait_syndrome,
               quality = quality,
               trait_dissimilarity = trait_dissimilarity)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 8: bettergardens_bowltrap_heteroptera ####
#'##############################################################################
focal_dataset <- "bettergardens_bowltrap_heteroptera"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(bodylength, Feeding_guild, Stratum, Wing_form) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis(trait_subset, w.type = "optimized", ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(comm = data_clean$comm %>% dplyr::select(any_of(rownames(trait_syndrome))),
               trait_syndrome = trait_syndrome,
               quality = quality,
               trait_dissimilarity = trait_dissimilarity)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 9: bettergardens_bowltrap_homoptera ####
#'##############################################################################
focal_dataset <- "bettergardens_bowltrap_homoptera"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Size_mean, Diet_width, Stratum, Wing_form) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis(trait_subset, ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(comm = data_clean$comm %>% dplyr::select(any_of(rownames(trait_syndrome))),
               trait_syndrome = trait_syndrome,
               quality = quality,
               trait_dissimilarity = trait_dissimilarity)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 10: bettergardens_bowltrap_syrphidae ####
#'##############################################################################
focal_dataset <- "bettergardens_bowltrap_syrphidae"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(ITD, food_specialization_adult, saproxylic_larvae, WingLoading) %>%
  mutate(across(food_specialization_adult, ordered)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis(trait_subset, ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(comm = data_clean$comm %>% dplyr::select(any_of(rownames(trait_syndrome))),
               trait_syndrome = trait_syndrome,
               quality = quality,
               trait_dissimilarity = trait_dissimilarity)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 11: bettergardens_gastropoda ####
#'##############################################################################
focal_dataset <- "bettergardens_gastropoda"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(shell_size, feeding_guild, humidity, snail_mobility_index) %>%
  mutate(across(.fns = as.numeric)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis(trait_subset, w.type = "optimized", ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(
  comm = data_clean$comm %>% dplyr::select(any_of(rownames(trait_syndrome))),
  trait_syndrome = trait_syndrome,
  trait_dissimilarity = trait_dissimilarity,
  quality = quality
)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))


#'##############################################################################
#  Dataset 12: bettergardens_pitfall_araneae ####
#'##############################################################################
focal_dataset <- "bettergardens_pitfall_araneae"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(bodylength, food_preference, habitat_width, mobility_mode) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis(trait_subset, ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(
  comm = data_clean$comm %>% dplyr::select(any_of(rownames(trait_syndrome))),
  trait_syndrome = trait_syndrome,
  trait_dissimilarity = trait_dissimilarity,
  quality = quality
)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))


#'##############################################################################
#  Dataset 13: bettergardens_pitfall_carabidae ####
#'##############################################################################
focal_dataset <- "bettergardens_pitfall_carabidae"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(bodylength, feeding_guild, habitatgeneralism, wing_form) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis(trait_subset, ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(
  comm = data_clean$comm %>% dplyr::select(any_of(rownames(trait_syndrome))),
  trait_syndrome = trait_syndrome,
  trait_dissimilarity = trait_dissimilarity,
  quality = quality
)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))


#'##############################################################################
#  Dataset 14: bettergardens_pitfall_heteroptera ####
#'##############################################################################
focal_dataset <- "bettergardens_pitfall_heteroptera"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(bodylength, Feeding_guild, Stratum, Wing_form) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis(trait_subset, ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(
  comm = data_clean$comm %>% dplyr::select(any_of(rownames(trait_syndrome))),
  trait_syndrome = trait_syndrome,
  trait_dissimilarity = trait_dissimilarity,
  quality = quality
)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))


#'##############################################################################
#  Dataset 15: bettergardens_pitfall_homoptera ####
#'##############################################################################
focal_dataset <- "bettergardens_pitfall_homoptera"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Size_mean, Diet_width, Stratum, Wing_form) %>%
  mutate(across(where(is.character), factor)) %>%
  mutate(across(Diet_width, ordered)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis(trait_subset, ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(
  comm = data_clean$comm %>% dplyr::select(any_of(rownames(trait_syndrome))),
  trait_syndrome = trait_syndrome,
  trait_dissimilarity = trait_dissimilarity,
  quality = quality
)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 16: bexis_forest_birds ####
#'##############################################################################
focal_dataset <- "bexis_forest_birds"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(size_bodymass, diet, contains("habitat_"), contains("mobility")) %>%
  mutate(across(diet, factor)) %>%
  mutate(across(!diet, as.numeric)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis(trait_subset,
                              groups = c(1, 2, rep(3, 7), 4),
                              fuzzy = c(3),
                              w.type = "optimized",
                              ord = "metric")

trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(
  comm = data_clean$comm %>% dplyr::select(any_of(rownames(trait_syndrome))),
  trait_syndrome = trait_syndrome,
  trait_dissimilarity = trait_dissimilarity,
  quality = quality
)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))


#'##############################################################################
#  Dataset 17: bexis_forest_flightinterceptioncanopy_beetles ####
#'##############################################################################
focal_dataset <- "bexis_forest_flightinterceptioncanopy_beetles"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Mean_BodySize, Feeding_guild_short, Stratum_use_short, Dispersal_ability) %>%
  mutate(across(Dispersal_ability, as.numeric)) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis(trait_subset, w.type = "optimized", ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(
  comm = data_clean$comm %>% dplyr::select(any_of(rownames(trait_syndrome))),
  trait_syndrome = trait_syndrome,
  trait_dissimilarity = trait_dissimilarity,
  quality = quality
)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))


#'##############################################################################
#  Dataset 18: bexis_forest_flightinterceptionunderstory_beetles ####
#'##############################################################################
focal_dataset <- "bexis_forest_flightinterceptionunderstory_beetles"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Mean_BodySize, Feeding_guild_short, Stratum_use_short, Dispersal_ability) %>%
  mutate(across(Dispersal_ability, as.numeric)) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis(trait_subset, w.type = "optimized", ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(
  comm = data_clean$comm %>% dplyr::select(any_of(rownames(trait_syndrome))),
  trait_syndrome = trait_syndrome,
  trait_dissimilarity = trait_dissimilarity,
  quality = quality
)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))


#'##############################################################################
#  Dataset 19: bexis_forest_hemiptera ####
#'##############################################################################
focal_dataset <- "bexis_forest_hemiptera"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Mean_BodySize, Feeding_guild_short, Stratum_use_short, Dispersal_ability) %>%
  mutate(across(Dispersal_ability, ordered)) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis(trait_subset, w.type = "optimized", ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(
  comm = data_clean$comm %>% dplyr::select(any_of(rownames(trait_syndrome))),
  trait_syndrome = trait_syndrome,
  trait_dissimilarity = trait_dissimilarity,
  quality = quality
)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))


#'##############################################################################
#  Dataset 20: bexis_forest_pitfall_beetles ####
#'##############################################################################
focal_dataset <- "bexis_forest_pitfall_beetles"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Mean_BodySize, Feeding_guild_short, Stratum_use_short, Dispersal_ability) %>%
  mutate(across(Dispersal_ability, ordered)) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis(trait_subset, w.type = "optimized", ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(
  comm = data_clean$comm %>% dplyr::select(any_of(rownames(trait_syndrome))),
  trait_syndrome = trait_syndrome,
  trait_dissimilarity = trait_dissimilarity,
  quality = quality
)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 21: bexis_forest_plants ####
#'##############################################################################
focal_dataset <- "bexis_forest_plants"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Height, SLA_all, T:N, Seed_mass) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x))))) 


names_remove <- FD::gowdis(trait_subset) %>%
  as.matrix %>%
  as.data.frame.table() %>%
  filter(is.na(Freq)) %>% distinct(Var2) %>% pull(Var2) %>% as.character()

trait_subset2 <- trait_subset %>%
  rownames_to_column("species") %>%
  filter(!species %in% names_remove) %>%
  column_to_rownames("species")


trait_dissimilarity <- gawdis(trait_subset2,
                              groups = c(1, 2, rep(3, 5), 4),
                              w.type = "optimized",
                              ord = "metric")

trait_syndrome <- ape::pcoa(trait_dissimilarity, correction="lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(comm = data_clean$comm %>% dplyr::select(any_of(rownames(trait_syndrome))),
               trait_syndrome = trait_syndrome,
               trait_dissimilarity = trait_dissimilarity,
               quality = quality)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))


#'##############################################################################
#  Dataset 22: bexis_forest_spiders ####
#'##############################################################################
focal_dataset <- "bexis_forest_spiders"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Mean_BodySize, contains("diet"), Stratum_use_short, Dispersal_ability) %>%
  mutate(across(where(is.character), as.factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis(trait_subset, w.type = "optimized", ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(comm = data_clean$comm %>% dplyr::select(any_of(rownames(trait_syndrome))),
               trait_syndrome = trait_syndrome,
               trait_dissimilarity = trait_dissimilarity,
               quality = quality)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))


#'##############################################################################
#  Dataset 23: bexis_forest_waterinsects ####
#'##############################################################################
focal_dataset <- "bexis_forest_waterinsects"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  mutate(across(where(is.character), as.factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis(trait_subset, w.type = "optimized", ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(comm = data_clean$comm %>% dplyr::select(any_of(rownames(trait_syndrome))),
               trait_syndrome = trait_syndrome,
               trait_dissimilarity = trait_dissimilarity,
               quality = quality)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))


#'##############################################################################
#  Dataset 24: bexis_grassland_birds ####
#'##############################################################################
focal_dataset <- "bexis_grassland_birds"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(size_bodymass, diet, contains("habitat"), contains("mobility")) %>%
  mutate(across(where(is.character), as.factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis(trait_subset,
                              groups = c(1, 2, rep(3, 7), 4),
                              fuzzy = c(3),
                              w.type = "optimized",
                              ord = "metric")

trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(comm = data_clean$comm %>% dplyr::select(any_of(rownames(trait_syndrome))),
               trait_syndrome = trait_syndrome,
               trait_dissimilarity = trait_dissimilarity,
               quality = quality)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))


#'##############################################################################
#  Dataset 25: bexis_grassland_pitfall_beetles ####
#'##############################################################################
focal_dataset <- "bexis_grassland_pitfall_beetles"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  mutate(across(where(is.character), as.factor)) %>%
  mutate(Wingclass = ordered(Wingclass)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis(trait_subset, w.type = "optimized", ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(comm = data_clean$comm %>% dplyr::select(any_of(rownames(trait_syndrome))),
               trait_syndrome = trait_syndrome,
               trait_dissimilarity = trait_dissimilarity,
               quality = quality)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 26: bexis_grassland_pitfall_spiders ####
#'##############################################################################
focal_dataset <- "bexis_grassland_pitfall_spiders"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset <- data_trait %>%
  select(Size_mean, contains("diet"), Stratum, Mobility) %>%
  mutate(across(where(is.character), as.factor)) %>%
  mutate(Mobility = ordered(Mobility)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized", ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(comm = data_clean$comm %>% dplyr::select(any_of(rownames(trait_syndrome))),
               trait_syndrome = trait_syndrome,
               quality = quality,
               trait_dissimilarity = trait_dissimilarity)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))


#'##############################################################################
#  Dataset 27: bexis_grassland_plants ####
#'##############################################################################
focal_dataset <- "bexis_grassland_plants"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")

trait_subset <- data_trait %>%
  mutate(across(where(is.character), as.factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x))))) 

names_remove <- FD::gowdis(trait_subset) %>%
  as.matrix %>%
  as.data.frame.table() %>%
  filter(is.na(Freq)) %>% distinct(Var2) %>% pull(Var2) %>% as.character()

trait_subset2 <- trait_subset %>%
  rownames_to_column("species") %>%
  filter(!species %in% names_remove) %>%
  column_to_rownames("species")

trait_dissimilarity <- gawdis::gawdis(trait_subset2,
                                      groups = c(1, 2, rep(3, 5), 4, 3),
                                      groups.weight = TRUE,
                                      w.type = "optimized",
                                      ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(comm = data_clean$comm %>% dplyr::select(any_of(rownames(trait_syndrome))),
               trait_syndrome = trait_syndrome,
               quality = quality,
               trait_dissimilarity = trait_dissimilarity)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))


#'##############################################################################
#  Dataset 28: bexis_grassland_sweepnet_beetles ####
#'##############################################################################
focal_dataset <- "bexis_grassland_sweepnet_beetles"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset <- data_trait %>%
  mutate(across(where(is.character), as.factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized", ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(comm = data_clean$comm %>% dplyr::select(any_of(rownames(trait_syndrome))),
               trait_syndrome = trait_syndrome,
               quality = quality,
               trait_dissimilarity = trait_dissimilarity)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))


#'##############################################################################
#  Dataset 29: bexis_grassland_sweepnet_heteroptera ####
#'##############################################################################
focal_dataset <- "bexis_grassland_sweepnet_heteroptera"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Mean_BodySize, Feeding_guild_short, Stratum_use_short, Dispersal_ability) %>%
  mutate(across(where(is.character), as.factor)) %>%
  mutate(Dispersal_ability = ordered(Dispersal_ability)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized", ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(comm = data_clean$comm %>% dplyr::select(any_of(rownames(trait_syndrome))),
               trait_syndrome = trait_syndrome,
               quality = quality,
               trait_dissimilarity = trait_dissimilarity)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))


#'##############################################################################
#  Dataset 30: bexis_grassland_sweepnet_homoptera ####
#'##############################################################################
focal_dataset <- "bexis_grassland_sweepnet_homoptera"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Mean_BodySize, diet_feedingtyp, Stratum_use_short, Dispersal_ability) %>%
  mutate(across(where(is.character), as.factor)) %>%
  mutate(Dispersal_ability = ordered(Dispersal_ability)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized", ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(comm = data_clean$comm %>% dplyr::select(any_of(rownames(trait_syndrome))),
               trait_syndrome = trait_syndrome,
               quality = quality,
               trait_dissimilarity = trait_dissimilarity)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))
#'##############################################################################
#  Dataset 31: bexis_grassland_sweepnet_orthoptera ####
#'##############################################################################
focal_dataset <- "bexis_grassland_sweepnet_orthoptera"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Mean_BodySize,Feeding_guild_short, Dispersal_ability) %>%
  mutate(across(where(is.character),as.factor)) %>%
  mutate(Dispersal_ability = ordered(Dispersal_ability)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis(trait_subset, w.type="optimized", ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(
  comm = data_clean$comm %>% dplyr::select(any_of(rownames(trait_syndrome))),
  trait_syndrome = trait_syndrome,
  trait_dissimilarity = trait_dissimilarity,
  quality = quality
)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 32: bexis_grassland_sweepnet_spiders ####
#'##############################################################################
focal_dataset <- "bexis_grassland_sweepnet_spiders"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Mean_BodySize,contains("diet"), Stratum_use_short, Dispersal_ability) %>%
  mutate(across(where(is.character),as.factor)) %>%
  mutate(Dispersal_ability = ordered(Dispersal_ability)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis(trait_subset, w.type="optimized", ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(
  comm = data_clean$comm %>% dplyr::select(any_of(rownames(trait_syndrome))),
  trait_syndrome = trait_syndrome,
  trait_dissimilarity = trait_dissimilarity,
  quality = quality
)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 33: bonada_waterinsects ####
#'##############################################################################
focal_dataset <- "bonada_waterinsects"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(a1:a7, h1:h9, u1:u8, f1:f4) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis(trait_subset,
                              groups = c(rep(1, 7), rep(2, 9), rep(3, 8), rep(4, 4)),
                              fuzzy = 1:4,
                              w.type = "optimized",
                              ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(comm = data_clean$comm %>% dplyr::select(any_of(rownames(trait_syndrome))),
               trait_syndrome = trait_syndrome,
               quality = quality,
               trait_dissimilarity = trait_dissimilarity)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))


#'##############################################################################
#  Dataset 34: brindamour_fish ####
#'##############################################################################
focal_dataset <- "brindamour_fish"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(max_length_fishbase, contains("Diet"), contains("Feed"), contains("migr")) %>%
  mutate(across(where(is.character), as.factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis(trait_subset,
                              groups = c(1, rep(2, 4), rep(3, 4), rep(4, 2)),
                              fuzzy = 2:4,
                              w.type = "optimized",
                              ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(comm = data_clean$comm %>% dplyr::select(any_of(rownames(trait_syndrome))),
               trait_syndrome = trait_syndrome,
               quality = quality,
               trait_dissimilarity = trait_dissimilarity)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))


#'##############################################################################
#  Dataset 35: brodersen_fish ####
#'##############################################################################
focal_dataset <- "brodersen_fish"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Weight_g, contains("diet"), freshwater, brackish, marine, demersal,
                benthopelagic, pelagic, anadromous, potamodromous, catadromous,
                aspect.ratio.of.caudal.fin) %>%
  mutate(across(.fns = as.numeric)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis(trait_subset,
                              groups = c(1, rep(2, 7), rep(3, 9), 4),
                              fuzzy = 2:3,
                              w.type = "optimized",
                              ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(comm = data_clean$comm,
               trait_syndrome = trait_syndrome,
               quality = quality,
               trait_dissimilarity = trait_dissimilarity)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))


#'##############################################################################
#  Dataset 36: carvalho_fish ####
#'##############################################################################
focal_dataset <- "carvalho_fish"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(BM, paste0("troph", c(1:3, 5:6)), contains("watpos"), contains("migra")) %>%
  mutate(across(.fns = as.numeric)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis(trait_subset,
                              groups = c(1, rep(2, 5), rep(3, 3), rep(4, 2)),
                              fuzzy = 2:4,
                              w.type = "optimized",
                              ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(comm = data_clean$comm,
               trait_syndrome = trait_syndrome,
               quality = quality,
               trait_dissimilarity = trait_dissimilarity)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))


#'##############################################################################
#  Dataset 37: chapman_birds ####
#'##############################################################################
focal_dataset <- "chapman_birds"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(pca.bodysize, pca.trophic.shape, habitat, pca.dispersal.shape) %>%
  mutate(across(habitat, .fns = ordered)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis(trait_subset,
                              w.type = "optimized",
                              ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(comm = data_clean$comm,
               trait_syndrome = trait_syndrome,
               quality = quality,
               trait_dissimilarity = trait_dissimilarity)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 38: decastro_waterinsects ####
#'##############################################################################
focal_dataset <- "decastro_waterinsects"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(contains("size"), contains("diet"), contains("habitat"), contains("mobility")) %>%
  mutate(across(.fns = as.numeric)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis::gawdis(trait_subset,
                                      groups = c(rep(1, 6), rep(2, 5), rep(3, 3), rep(4, 5)),
                                      fuzzy = 2:4,
                                      w.type = "optimized",
                                      ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(comm = data_clean$comm,
               trait_syndrome = trait_syndrome,
               trait_dissimilarity = trait_dissimilarity,
               quality = quality)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))


#'##############################################################################
#  Dataset 39: desouzaqueiroz_amphibians ####
#'##############################################################################
focal_dataset <- "desouzaqueiroz_amphibians"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(size_max, Feeding_behaviour, Neustonic:Nektonic, WTM.TL) %>%
  mutate(across(Feeding_behaviour, factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis::gawdis(trait_subset,
                                      groups = c(1, 2, 3, 3, 3, 4),
                                      fuzzy = 3,
                                      w.type = "equal",
                                      ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(comm = data_clean$comm,
               trait_syndrome = trait_syndrome,
               trait_dissimilarity = trait_dissimilarity,
               quality = quality)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))


#'##############################################################################
#  Dataset 40: drose_ants ####
#'##############################################################################
focal_dataset <- "drose_ants"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(contains("SIZE"), contains("DIET"), contains("HABITAT"), contains("MOBILITY")) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis::gawdis(trait_subset,
                                      w.type = "optimized",
                                      ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(comm = data_clean$comm,
               trait_syndrome = trait_syndrome,
               trait_dissimilarity = trait_dissimilarity,
               quality = quality)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))
#'##############################################################################
#  Dataset 41: gallardo_waterinsects ####
#'##############################################################################
focal_dataset <- "gallardo_waterinsects"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(a1:a7, k1:k9, j1:j9, h1:h8) %>%
  mutate(across(.fns = as.numeric)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis::gawdis(trait_subset,
                                      groups = c(rep(1, 7),
                                                 rep(2, 9),
                                                 rep(3, 9),
                                                 rep(4, 8)),
                                      fuzzy = 1:4,
                                      w.type = "analytic",
                                      ord = "metric")

trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, quality = quality, trait_dissimilarity = trait_dissimilarity)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))


#'##############################################################################
#  Dataset 42: gracoroza_phytoplankton_1 ####
#'##############################################################################
focal_dataset <- "gracoroza_phytoplankton_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(MLD, Mucilage:Aerotopes, Heterocyte:Silica, Flagella) %>%
  mutate(across(.fns = as.numeric)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis::gawdis(trait_subset,
                                      groups = c(1, 2, 2, 3, 3, 4),
                                      w.type = "optimized",
                                      ord = "metric")

trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, quality = quality, trait_dissimilarity = trait_dissimilarity)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))


#'##############################################################################
#  Dataset 43: gracoroza_phytoplankton_2 ####
#'##############################################################################
focal_dataset <- "gracoroza_phytoplankton_2"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(MLD, Mucilage:Aerotopes, Heterocyte:Silica, Flagella) %>%
  mutate(across(.fns = as.numeric)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis::gawdis(trait_subset,
                                      groups = c(1, 2, 2, 3, 3, 4),
                                      w.type = "optimized",
                                      ord = "metric")

trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, quality = quality, trait_dissimilarity = trait_dissimilarity)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))


#'##############################################################################
#  Dataset 44: gracoroza_phytoplankton_3 ####
#'##############################################################################
focal_dataset <- "gracoroza_phytoplankton_3"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(MLD, sv, Mucilage:Aerotopes, Silica, Flagella) %>%
  mutate(across(.fns = as.numeric)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis::gawdis(trait_subset,
                                      groups = c(rep(1, 2), rep(2, 2), 3, 4),
                                      w.type = "optimized",
                                      ord = "metric")

trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, quality = quality, trait_dissimilarity = trait_dissimilarity)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))


#'##############################################################################
#  Dataset 45: jeliazkov_waterinsects ####
#'##############################################################################
focal_dataset <- "jeliazkov_waterinsects"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(contains("Size"), contains("Food"), contains("Sub"), contains("Loco")) %>%
  mutate(across(.fns = as.numeric)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis::gawdis(trait_subset,
                                      groups = c(rep(1, 7), rep(2, 9), rep(3, 9), rep(4, 8)),
                                      w.type = "analytic",
                                      ord = "metric")

trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, quality = quality, trait_dissimilarity = trait_dissimilarity)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 46: lee_ants ####
#'##############################################################################
focal_dataset <- "lee_ants"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(WL, ML, HW, LL) %>%
  mutate(across(.fns = as.numeric)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "analytic", ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome,
               quality = quality, trait_dissimilarity = trait_dissimilarity)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))


#'##############################################################################
#  Dataset 47: lehikoinen_birds ####
#'##############################################################################
focal_dataset <- "lehikoinen_birds"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Mass, Diet.Inv:Diet.PlantO, contains("ForStrat"), contains("Mig")) %>%
  mutate(across(Mig, factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis::gawdis(trait_subset,
                                      groups = c(1, rep(2, 10), rep(3, 7), 4),
                                      fuzzy = c(2, 3),
                                      w.type = "optimized",
                                      ord = "metric")

trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome,
               quality = quality, trait_dissimilarity = trait_dissimilarity)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))


#'##############################################################################
#  Dataset 48: lowe_spiders ####
#'##############################################################################
focal_dataset <- "lowe_spiders"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(contains("SIZE"), contains("DIET"), contains("HABITAT"), contains("MOBILITY")) %>%
  mutate(across(contains("DIET"), factor),
         SIZE_Size = ordered(SIZE_Size, levels = c("small", "medium", "large"))) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis::gawdis(trait_subset,
                                      w.type = "optimized",
                                      ord = "metric")

trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome,
               quality = quality, trait_dissimilarity = trait_dissimilarity)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))


#'##############################################################################
#  Dataset 49: luoto_butterflies ####
#'##############################################################################
focal_dataset <- "luoto_butterflies"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset <- data_trait %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis::gawdis(trait_subset,
                                      groups = c(2,2,3,4,1),
                                      w.type = "optimized",
                                      ord = "metric")

trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome,
               quality = quality, trait_dissimilarity = trait_dissimilarity)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))


#'##############################################################################
#  Dataset 50: meffert_birds ####
#'##############################################################################
focal_dataset <- "meffert_birds"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>%  read_rds()
data_trait <- data_clean %>%  pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(body.mass, food, foraging.technique, mig1.strategy) %>%
  mutate(across(where(is.character), factor),
         mig1.strategy = ordered(mig1.strategy, c("sedentary", "short", "long"))) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis::gawdis(trait_subset,
                                      w.type = "optimized",
                                      ord = "metric")

trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome,
               quality = quality, trait_dissimilarity = trait_dissimilarity)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 51: N21FMI ####
#'##############################################################################
focal_dataset <- "N21FMI"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(contains("SIZE"), contains("DIET"), contains("HABITAT"), contains("MOBILITY")) %>%
  mutate(across(.fns = ordered)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized", ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome,
               quality = quality, trait_dissimilarity = trait_dissimilarity)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))


#'##############################################################################
#  Dataset 52: N25FBD ####
#'##############################################################################
focal_dataset <- "N25FBD"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(contains("SIZE"), High.profile.guild:Planktonic, contains("HABITAT"), contains("MOBILITY")) %>%
  mutate(SIZE_size.class = ordered(SIZE_size.class)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis::gawdis(trait_subset,
                                      groups = c(1, rep(2, 4), rep(3, 4), 4),
                                      w.type = "optimized",
                                      ord = "metric")

trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome,
               quality = quality, trait_dissimilarity = trait_dissimilarity)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))


#'##############################################################################
#  Dataset 53: N25FMI ####
#'##############################################################################
focal_dataset <- "N25FMI"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(contains("SIZE"), contains("DIET"), contains("HABITAT"), contains("MOBILITY")) %>%
  mutate(across(.fns = ordered)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized", ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome,
               quality = quality, trait_dissimilarity = trait_dissimilarity)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))


#'##############################################################################
#  Dataset 54: N29FBD ####
#'##############################################################################
focal_dataset <- "N29FBD"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(contains("SIZE"), contains("DIET"), contains("HABITAT"), contains("MOBILITY")) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis::gawdis(trait_subset,
                                      groups = c(1, 2, rep(3, 4), 4),
                                      w.type = "optimized",
                                      ord = "metric")

trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome,
               quality = quality, trait_dissimilarity = trait_dissimilarity)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))


#'##############################################################################
#  Dataset 55: N29FMI ####
#'##############################################################################
focal_dataset <- "N29FMI"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(contains("SIZE"), contains("DIET"), contains("HABITAT"), contains("MOBILITY")) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis::gawdis(trait_subset,
                                      w.type = "optimized",
                                      ord = "metric")

trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome,
               quality = quality, trait_dissimilarity = trait_dissimilarity)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))


#'##############################################################################
#  Dataset 56: N34FBD ####
#'##############################################################################
focal_dataset <- "N34FBD"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(contains("SIZE"), contains("DIET"), contains("HABITAT"), contains("MOBILITY")) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis::gawdis(trait_subset,
                                      groups = c(1, 2, rep(3, 4), 4),
                                      w.type = "optimized",
                                      ord = "metric")

trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(comm = data_clean$comm,
               trait_syndrome = trait_syndrome,
               quality = quality,
               trait_dissimilarity = trait_dissimilarity)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))


#'##############################################################################
#  Dataset 57: N34FMI ####
#'##############################################################################
focal_dataset <- "N34FMI"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(contains("SIZE"), contains("DIET"), contains("HABITAT"), contains("MOBILITY")) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis::gawdis(trait_subset,
                                      w.type = "optimized",
                                      ord = "metric")

trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(comm = data_clean$comm,
               trait_syndrome = trait_syndrome,
               quality = quality,
               trait_dissimilarity = trait_dissimilarity)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))


#'##############################################################################
#  Dataset 58: N38TTP ####
#'##############################################################################
focal_dataset <- "N38TTP"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(contains("SIZE"), contains("DIET"), T:N, contains("MOBILITY")) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

names_remove <- FD::gowdis(trait_subset) %>%
  as.matrix() %>%
  as.data.frame.table() %>%
  slice(which(is.na(.$Freq))) %>%
  select(Var2) %>%
  t() %>%
  c() %>%
  unique()

trait_subset2 <- trait_subset %>%
  rownames_to_column("species") %>%
  filter(!species %in% names_remove) %>%
  column_to_rownames("species")

trait_dissimilarity <- gawdis::gawdis(trait_subset2,
                                      groups = c(1, 2, rep(3, 5), 4),
                                      w.type = "optimized",
                                      ord = "metric")

trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(comm = data_clean$comm,
               trait_syndrome = trait_syndrome,
               quality = quality,
               trait_dissimilarity = trait_dissimilarity)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))


#'##############################################################################
#  Dataset 59: N39FMI_2 ####
#'##############################################################################
focal_dataset <- "N39FMI_2"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(contains("SIZE"), contains("food"), contains("feedhabit"), contains("loc")) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis::gawdis(trait_subset,
                                      groups = c(rep(1, 2), rep(2, 2), rep(3, 2), rep(4, 2)),
                                      w.type = "optimized",
                                      ord = "metric")

trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(comm = data_clean$comm,
               trait_syndrome = trait_syndrome,
               quality = quality,
               trait_dissimilarity = trait_dissimilarity)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))


#'##############################################################################
#  Dataset 60: N39FMI ####
#'##############################################################################
focal_dataset <- "N39FMI"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(contains("SIZE"), contains("food"), contains("feedhabit"), contains("loc")) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis::gawdis(trait_subset,
                                      groups = c(rep(1, 2), rep(2, 2), rep(3, 2), rep(4, 2)),
                                      w.type = "optimized",
                                      ord = "metric")

trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(comm = data_clean$comm,
               trait_syndrome = trait_syndrome,
               quality = quality,
               trait_dissimilarity = trait_dissimilarity)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))
#'##############################################################################
#  Dataset 61: N39TTP ####
#'##############################################################################
focal_dataset <- "N39TTP"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")

trait_subset <- data_trait %>%
  dplyr::select(contains("SIZE"), contains("DIET"), T:N, contains("MOBILITY")) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

names_remove <- FD::gowdis(trait_subset) %>%
  as.matrix() %>%
  as.data.frame.table() %>%
  slice(which(is.na(.$Freq))) %>%
  select(Var2) %>%
  t() %>%
  c() %>%
  unique()

trait_subset2 <- trait_subset %>%
  rownames_to_column("species") %>%
  filter(!species %in% names_remove) %>%
  column_to_rownames("species")

trait_dissimilarity <- gawdis::gawdis(trait_subset2,
                                      groups = c(1, 2, rep(3, 5), 4),
                                      w.type = "optimized",
                                      ord = "metric")

trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(comm = data_clean$comm,
               trait_syndrome = trait_syndrome,
               quality = quality,
               trait_dissimilarity = trait_dissimilarity)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))


#'##############################################################################
#  Dataset 62: N40FMI ####
#'##############################################################################
focal_dataset <- "N40FMI"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")

trait_subset <- data_trait %>%
  rename(Loc_sec = Log_sec) %>%
  dplyr::select(contains("SIZE"), contains("Food"), contains("Feed"), contains("Loc")) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis::gawdis(trait_subset,
                                      groups = c(rep(1, 2), rep(2, 2), rep(3, 2), rep(4, 2)),
                                      w.type = "optimized",
                                      ord = "metric")

trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(comm = data_clean$comm,
               trait_syndrome = trait_syndrome,
               quality = quality,
               trait_dissimilarity = trait_dissimilarity)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))


#'##############################################################################
#  Dataset 63: N41TTP ####
#'##############################################################################
focal_dataset <- "N41TTP"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")

trait_subset <- data_trait %>%
  dplyr::select(contains("SIZE"), contains("DIET"), T:N, contains("MOBILITY")) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

names_remove <- FD::gowdis(trait_subset) %>%
  as.matrix() %>%
  as.data.frame.table() %>%
  slice(which(is.na(.$Freq))) %>%
  select(Var2) %>%
  t() %>%
  c() %>%
  unique()

trait_subset2 <- trait_subset %>%
  rownames_to_column("species") %>%
  filter(!species %in% names_remove) %>%
  column_to_rownames("species")

trait_dissimilarity <- gawdis::gawdis(trait_subset2,
                                      groups = c(1, 2, rep(3, 5), 4),
                                      w.type = "optimized",
                                      ord = "metric")

trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(comm = data_clean$comm,
               trait_syndrome = trait_syndrome,
               quality = quality,
               trait_dissimilarity = trait_dissimilarity)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))


#'##############################################################################
#  Dataset 64: N42TTP_2 ####
#'##############################################################################
focal_dataset <- "N42TTP_2"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")

trait_subset <- data_trait %>%
  dplyr::select(contains("SIZE"), contains("DIET"), T:N, contains("MOBILITY")) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

names_remove <- FD::gowdis(trait_subset) %>%
  as.matrix() %>%
  as.data.frame.table() %>%
  slice(which(is.na(.$Freq))) %>%
  select(Var2) %>%
  t() %>%
  c() %>%
  unique()

trait_subset2 <- trait_subset %>%
  rownames_to_column("species") %>%
  filter(!species %in% names_remove) %>%
  column_to_rownames("species")

trait_dissimilarity <- gawdis::gawdis(trait_subset2,
                                      groups = c(1, 2, rep(3, 5), 4),
                                      w.type = "optimized",
                                      ord = "metric")

trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(comm = data_clean$comm,
               trait_syndrome = trait_syndrome,
               quality = quality,
               trait_dissimilarity = trait_dissimilarity)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))


#'##############################################################################
#  Dataset 65: N42TTP_3 ####
#'##############################################################################
focal_dataset <- "N42TTP_3"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")

trait_subset <- data_trait %>%
  dplyr::select(contains("SIZE"), contains("DIET"), T:N, contains("MOBILITY")) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

names_remove <- FD::gowdis(trait_subset) %>%
  as.matrix() %>%
  as.data.frame.table() %>%
  slice(which(is.na(.$Freq))) %>%
  select(Var2) %>%
  t() %>%
  c() %>%
  unique()

trait_subset2 <- trait_subset %>%
  rownames_to_column("species") %>%
  filter(!species %in% names_remove) %>%
  column_to_rownames("species")

trait_dissimilarity <- gawdis::gawdis(trait_subset2,
                                      groups = c(1, 2, rep(3, 5), 4),
                                      w.type = "optimized",
                                      ord = "metric")

trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(comm = data_clean$comm,
               trait_syndrome = trait_syndrome,
               quality = quality,
               trait_dissimilarity = trait_dissimilarity)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))
#'##############################################################################
#  Dataset 66: N42TTP ####
#'##############################################################################
focal_dataset <- "N42TTP"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(contains("SIZE"), contains("DIET"), T:N, contains("MOBILITY")) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

names_remove <- FD::gowdis(trait_subset) %>%
  as.matrix() %>%
  as.data.frame.table() %>%
  slice(which(is.na(.$Freq))) %>%
  select(Var2) %>%
  t() %>%
  c() %>%
  unique()

trait_subset2 <- trait_subset %>%
  rownames_to_column("species") %>%
  filter(!species %in% names_remove) %>%
  column_to_rownames("species")

trait_dissimilarity <- gawdis::gawdis(trait_subset2,
                                      groups = c(1, 2, rep(3, 5), 4),
                                      ord = "metric",
                                      w.type = "optimized"
)

trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(
  comm = data_clean$comm,
  trait_syndrome = trait_syndrome,
  quality = quality,
  trait_dissimilarity = trait_dissimilarity
)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))
#'##############################################################################
#  Dataset 67: N43TTP_1 ####
#'##############################################################################
focal_dataset <- "N43TTP_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(contains("SIZE"), contains("DIET"), T:N, contains("MOBILITY")) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

names_remove <- FD::gowdis(trait_subset) %>%
  as.matrix() %>%
  as.data.frame.table() %>%
  slice(which(is.na(.$Freq))) %>%
  select(Var2) %>%
  t() %>%
  c() %>%
  unique()

trait_subset2 <- trait_subset %>%
  rownames_to_column("species") %>%
  filter(!species %in% names_remove) %>%
  column_to_rownames("species")

trait_dissimilarity <- gawdis::gawdis(trait_subset2,
                                      groups = c(1, 2, rep(3, 5), 4),
                                      ord = "metric",
                                      w.type = "optimized"
)

trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(
  comm = data_clean$comm,
  trait_syndrome = trait_syndrome,
  quality = quality,
  trait_dissimilarity = trait_dissimilarity
)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))
#'##############################################################################
#  Dataset 68: N43TTP_2 ####
#'##############################################################################
focal_dataset <- "N43TTP_2"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(contains("SIZE"), contains("DIET"), T:N, contains("MOBILITY")) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

names_remove <- FD::gowdis(trait_subset) %>%
  as.matrix() %>%
  as.data.frame.table() %>%
  slice(which(is.na(.$Freq))) %>%
  select(Var2) %>%
  t() %>%
  c() %>%
  unique()

trait_subset2 <- trait_subset %>%
  rownames_to_column("species") %>%
  filter(!species %in% names_remove) %>%
  column_to_rownames("species")

trait_dissimilarity <- gawdis::gawdis(trait_subset2,
                                      groups = c(1, 2, rep(3, 5), 4),
                                      ord = "metric",
                                      w.type = "optimized"
)

trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(
  comm = data_clean$comm,
  trait_syndrome = trait_syndrome,
  quality = quality,
  trait_dissimilarity = trait_dissimilarity
)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 69: N43TTP_3 ####
#'##############################################################################
focal_dataset <- "N43TTP_3"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(contains("SIZE"), contains("DIET"), T:N, contains("MOBILITY")) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

names_remove <- FD::gowdis(trait_subset) %>%
  as.matrix() %>%
  as.data.frame.table() %>%
  slice(which(is.na(.$Freq))) %>%
  select(Var2) %>%
  t() %>%
  c() %>%
  unique()

trait_subset2 <- trait_subset %>%
  rownames_to_column("species") %>%
  filter(!species %in% names_remove) %>%
  column_to_rownames("species")

trait_dissimilarity <- gawdis::gawdis(trait_subset2,
                                      groups = c(1, 2, rep(3, 5), 4),
                                      ord = "metric",
                                      w.type = "optimized"
)

trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(
  comm = data_clean$comm,
  trait_syndrome = trait_syndrome,
  quality = quality,
  trait_dissimilarity = trait_dissimilarity
)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 70: N44TBT ####
#'##############################################################################
focal_dataset <- "N44TBT"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(contains("SIZE"), contains("DIET"), contains("HABITAT"), contains("MOBILITY")) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis::gawdis(trait_subset,
                                      ord = "metric",
                                      w.type = "optimized"
)

trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)

output <- list(
  comm = data_clean$comm,
  trait_syndrome = trait_syndrome,
  quality = quality,
  trait_dissimilarity = trait_dissimilarity
)

saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))
#'##############################################################################
#  Dataset 71: N44TTP ####
#'##############################################################################
focal_dataset <- "N44TTP"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(PlantHeight_mean, SLA_mean, e_light:e_temp, SeedMass_mean) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

names_remove <- FD::gowdis(trait_subset) %>%
  as.matrix() %>%
  as.data.frame.table() %>%
  slice(which(is.na(.$Freq))) %>%
  select(Var2) %>%
  t() %>%
  c() %>%
  unique()

trait_subset2 <- trait_subset %>%
  rownames_to_column("species") %>%
  filter(!species %in% names_remove) %>%
  column_to_rownames("species")

trait_dissimilarity <- gawdis::gawdis(trait_subset2,
                                      groups = c(1, 2, rep(3, 6), 4),
                                      ord = "metric",
                                      w.type = "optimized"
)
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, quality = quality,
               trait_dissimilarity = trait_dissimilarity)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 72: N47FFI ####
#'##############################################################################
focal_dataset <- "N47FFI"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(contains("SIZE_"), contains("DIET"), contains("HABITAT"), contains("MOBILITY")) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis::gawdis(trait_subset,
                                      groups = c(1, rep(2, 6), rep(3, 3), rep(4, 7)),
                                      fuzzy = c(2:4),
                                      ord = "metric",
                                      w.type = "optimized"
)
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, quality = quality,
               trait_dissimilarity = trait_dissimilarity)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 73: N47TTP_2 ####
#'##############################################################################
focal_dataset <- "N47TTP_2"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(contains("SIZE"), contains("DIET"), T:N, contains("MOBILITY")) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

names_remove <- FD::gowdis(trait_subset) %>%
  as.matrix() %>%
  as.data.frame.table() %>%
  slice(which(is.na(.$Freq))) %>%
  select(Var2) %>%
  t() %>%
  c() %>%
  unique()

trait_subset2 <- trait_subset %>%
  rownames_to_column("species") %>%
  filter(!species %in% names_remove) %>%
  column_to_rownames("species")

trait_dissimilarity <- gawdis::gawdis(trait_subset2,
                                      groups = c(1, 2, rep(3, 5), 4),
                                      ord = "metric",
                                      w.type = "optimized"
)
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, quality = quality,
               trait_dissimilarity = trait_dissimilarity)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 74: N47TTP ####
#'##############################################################################
focal_dataset <- "N47TTP"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(contains("SIZE"), contains("DIET"), T:N, contains("MOBILITY")) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

names_remove <- FD::gowdis(trait_subset) %>%
  as.matrix() %>%
  as.data.frame.table() %>%
  slice(which(is.na(.$Freq))) %>%
  select(Var2) %>%
  t() %>%
  c() %>%
  unique()

trait_subset2 <- trait_subset %>%
  rownames_to_column("species") %>%
  filter(!species %in% names_remove) %>%
  column_to_rownames("species")

trait_dissimilarity <- gawdis::gawdis(trait_subset2,
                                      groups = c(1, 2, rep(3, 5), 4),
                                      ord = "metric",
                                      w.type = "optimized"
)
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, quality = quality,
               trait_dissimilarity = trait_dissimilarity)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 75: N48FBD ####
#'##############################################################################
focal_dataset <- "N48FBD"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(contains("SIZE"), contains("DIET"), contains("HABITAT"), contains("MOBILITY")) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis::gawdis(trait_subset,
                                      groups = c(1, 2, rep(3, 4), 4),
                                      ord = "metric",
                                      w.type = "optimized"
)
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, quality = quality,
               trait_dissimilarity = trait_dissimilarity)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))
#'##############################################################################
#  Dataset 76: N48FMI ####
#'##############################################################################
focal_dataset <- "N48FMI"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(contains("SIZE"), contains("DIET"), contains("HABITAT"), contains("MOBILITY")) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis::gawdis(trait_subset, ord = "metric", w.type = "optimized")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, quality = quality)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 77: N49TBI ####
#'##############################################################################
focal_dataset <- "N49TBI"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(contains("SIZE"), contains("DIET"), contains("HABITAT"), contains("MOBILITY")) %>%
  mutate(across(where(is.character), factor)) %>%
  mutate(MOBILITY_home = ordered(MOBILITY_home, levels = c("small", "mid", "large"))) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis::gawdis(trait_subset, ord = "metric", w.type = "optimized")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, quality = quality,
               trait_dissimilarity = trait_dissimilarity)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 78: N49TBT ####
#'##############################################################################
focal_dataset <- "N49TBT"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(contains("SIZE"), contains("DIET"), contains("HABITAT"), contains("MOBILITY")) %>%
  mutate(across(where(is.character), factor)) %>%
  mutate(MOBILITY_home = ordered(MOBILITY_home, levels = c("small", "mid", "large"))) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis::gawdis(trait_subset, ord = "metric", w.type = "optimized")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, quality = quality,
               trait_dissimilarity = trait_dissimilarity)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 79: N50TFG ####
#'##############################################################################
focal_dataset <- "N50TFG"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(contains("SIZE"), contains("DIET"), contains("HABITAT"), contains("MOBILITY")) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis::gawdis(trait_subset, ord = "metric", w.type = "optimized")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, quality = quality,
               trait_dissimilarity = trait_dissimilarity)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 80: N51TTP ####
#'##############################################################################
focal_dataset <- "N51TTP"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(contains("SIZE"), contains("DIET"), T:N, contains("MOBILITY")) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

names_remove <- FD::gowdis(trait_subset) %>%
  as.matrix() %>%
  as.data.frame.table() %>%
  slice(which(is.na(.$Freq))) %>%
  select(Var2) %>%
  t() %>%
  c() %>%
  unique()

trait_subset2 <- trait_subset %>%
  rownames_to_column("species") %>%
  filter(!species %in% names_remove) %>%
  column_to_rownames("species")

trait_dissimilarity <- gawdis::gawdis(trait_subset2,
                                      groups = c(1, 2, rep(3, 5), 4),
                                      ord = "metric",
                                      w.type = "optimized"
)
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, quality = quality,
               trait_dissimilarity = trait_dissimilarity)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))
#'##############################################################################
#  Dataset 81: N55TAT ####
#'##############################################################################
focal_dataset <- "N55TAT"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(contains("SIZE"), contains("DIET"), contains("HABITAT"), contains("MOBILITY")) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis::gawdis(trait_subset, ord = "metric", w.type = "optimized")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, quality = quality,
               trait_dissimilarity = trait_dissimilarity)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 82: N55TTP ####
#'##############################################################################
focal_dataset <- "N55TTP"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(contains("SIZE"), contains("DIET"), T:N, contains("MOBILITY")) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

names_remove <- FD::gowdis(trait_subset) %>%
  as.matrix() %>%
  as.data.frame.table() %>%
  slice(which(is.na(.$Freq))) %>%
  select(Var2) %>%
  t() %>%
  c() %>%
  unique()

trait_subset2 <- trait_subset %>%
  rownames_to_column("species") %>%
  filter(!species %in% names_remove) %>%
  column_to_rownames("species")

trait_dissimilarity <- gawdis::gawdis(trait_subset2,
                                      groups = c(1, 2, rep(3, 5), 4),
                                      ord = "metric",
                                      w.type = "optimized"
)
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, quality = quality,
               trait_dissimilarity = trait_dissimilarity)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 83: N62FMI ####
#'##############################################################################
focal_dataset <- "N62FMI"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(contains("SIZE"), contains("DIET"), contains("HABITAT"), contains("MOBILITY")) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis::gawdis(trait_subset,
                                      groups = c(1, 2, 3, rep(4, 8)),
                                      ord = "metric",
                                      w.type = "optimized"
)
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, quality = quality,
               trait_dissimilarity = trait_dissimilarity)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 84: N67TTP ####
#'##############################################################################
focal_dataset <- "N67TTP"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Height, SLA, T:N, Seed_mass) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

names_remove <- FD::gowdis(trait_subset) %>%
  as.matrix() %>%
  as.data.frame.table() %>%
  slice(which(is.na(.$Freq))) %>%
  select(Var2) %>%
  t() %>%
  c() %>%
  unique()

trait_subset2 <- trait_subset %>%
  rownames_to_column("species") %>%
  filter(!species %in% names_remove) %>%
  column_to_rownames("species")

trait_dissimilarity <- gawdis::gawdis(trait_subset2,
                                      groups = c(1, 2, rep(3, 5), 4),
                                      ord = "metric",
                                      w.type = "optimized"
)
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, quality = quality,
               trait_dissimilarity = trait_dissimilarity)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 85: N69FMI ####
#'##############################################################################
focal_dataset <- "N69FMI"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(contains("SIZE"), contains("DIET"), contains("HABITAT"), contains("MOBILITY")) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis::gawdis(trait_subset,
                                      groups = c(1, 2, 3, rep(4, 4)),
                                      ord = "metric",
                                      w.type = "optimized"
)
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, quality = quality,
               trait_dissimilarity = trait_dissimilarity)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))
#'##############################################################################
#  Dataset 86: N78TTP ####
#'##############################################################################
focal_dataset <- "N78TTP"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Height, SLA, T:N, Seed_mass) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

names_remove <- FD::gowdis(trait_subset) %>%
  as.matrix() %>%
  as.data.frame.table() %>%
  slice(which(is.na(.$Freq))) %>%
  select(Var2) %>%
  t() %>%
  c() %>%
  unique()

trait_subset2 <- trait_subset %>%
  rownames_to_column("species") %>%
  filter(!species %in% names_remove) %>%
  column_to_rownames("species")

trait_dissimilarity <- gawdis::gawdis(trait_subset2, groups = c(1, 2, rep(3, 5), 4), ord = "metric", w.type = "optimized")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, quality = quality,
               trait_dissimilarity = trait_dissimilarity)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 87: ossola_ants ####
#'##############################################################################
focal_dataset <- "ossola_ants"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(contains("SIZE"), contains("DIET_"), contains("HABITAT"), contains("MOBILITY")) %>%
  dplyr::select(-DIET_Diet) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis::gawdis(trait_subset, groups = c(1, rep(2, 6), 3, 4), ord = "metric", w.type = "optimized")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, quality = quality,
               trait_dissimilarity = trait_dissimilarity)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 88: penone_orthoptera ####
#'##############################################################################
focal_dataset <- "penone_orthoptera"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(contains("SIZE"), contains("DIET"), contains("HABITAT"), contains("dispersal")) %>%
  select(-diet.from.pyrgus.de, -habitat.from.pyrgus.de) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis::gawdis(trait_subset, groups = c(rep(1, 3), rep(2, 5), 3, 4), ord = "metric", w.type = "optimized")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, quality = quality,
               trait_dissimilarity = trait_dissimilarity)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 89: piano_beetles ####
#'##############################################################################
focal_dataset <- "piano_beetles"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(contains("SIZE"), contains("DIET"), contains("HABITAT"), contains("MOBILITY")) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis::gawdis(trait_subset, ord = "metric", w.type = "optimized")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, quality = quality,
               trait_dissimilarity = trait_dissimilarity)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 90: piano_spiders ####
#'##############################################################################
focal_dataset <- "piano_spiders"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(contains("SIZE"), contains("DIET"), contains("HABITAT"), contains("MOBILITY")) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis::gawdis(trait_subset, ord = "metric", w.type = "optimized")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, quality = quality,
               trait_dissimilarity = trait_dissimilarity)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 91: Predicts_cc1_2013_waite_1 ####
#'##############################################################################
focal_dataset <- "Predicts_cc1_2013_waite_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Body_mass_g, Trophic_level, Habitat_breadth_IUCN, Hand.Wing.Index) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized", ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, quality = quality,
               trait_dissimilarity = trait_dissimilarity)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))
 

#'##############################################################################
#  Dataset 92: Predicts_db1_2010_dures_1 ####
#'##############################################################################
focal_dataset <- "Predicts_db1_2010_dures_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Body_mass_g, Trophic_level, Habitat_breadth_IUCN, Hand.Wing.Index) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized", ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, quality = quality,
               trait_dissimilarity = trait_dissimilarity)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))
 

#'##############################################################################
#  Dataset 93: Predicts_di1_2004_naidoo_1 ####
#'##############################################################################
focal_dataset <- "Predicts_di1_2004_naidoo_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Body_mass_g, Trophic_level, Habitat_breadth_IUCN, Hand.Wing.Index) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized", ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, quality = quality,
               trait_dissimilarity = trait_dissimilarity)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))
 

#'##############################################################################
#  Dataset 94: Predicts_di1_2010_milder_2 ####
#'##############################################################################
focal_dataset <- "Predicts_di1_2010_milder_2"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Body_mass_g, Trophic_level, Habitat_breadth_IUCN, Hand.Wing.Index) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized", ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, quality = quality,
               trait_dissimilarity = trait_dissimilarity)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))
 

#'##############################################################################
#  Dataset 95: Predicts_di1_2011_dawson_1 ####
#'##############################################################################
focal_dataset <- "Predicts_di1_2011_dawson_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Body_mass_g, Trophic_level, Habitat_breadth_IUCN, Hand.Wing.Index) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized", ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, quality = quality,
               trait_dissimilarity = trait_dissimilarity)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))
#'##############################################################################
#  Dataset 96: Predicts_di1_2011_neuschulz_1 ####
#'##############################################################################
focal_dataset <- "Predicts_di1_2011_neuschulz_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Body_mass_g, Trophic_level, Habitat_breadth_IUCN, Hand.Wing.Index) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized", ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, quality = quality,
               trait_dissimilarity = trait_dissimilarity)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))
 

#'##############################################################################
#  Dataset 97: Predicts_di1_2012_reid_1 ####
#'##############################################################################
focal_dataset <- "Predicts_di1_2012_reid_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Body_mass_g, Trophic_level, Habitat_breadth_IUCN, Hand.Wing.Index) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized", ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, quality = quality,
               trait_dissimilarity = trait_dissimilarity)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))
 

#'##############################################################################
#  Dataset 98: Predicts_di1_2013_azhar_1 ####
#'##############################################################################
focal_dataset <- "Predicts_di1_2013_azhar_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Body_mass_g, Trophic_level, Habitat_breadth_IUCN, Hand.Wing.Index) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized", ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, quality = quality,
               trait_dissimilarity = trait_dissimilarity)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))
 

#'##############################################################################
#  Dataset 99: Predicts_di1_2013_de_lima_1 ####
#'##############################################################################
focal_dataset <- "Predicts_di1_2013_de_lima_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Body_mass_g, Trophic_level, Habitat_breadth_IUCN, Hand.Wing.Index) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized", ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, quality = quality,
               trait_dissimilarity = trait_dissimilarity)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))
 

#'##############################################################################
#  Dataset 100: Predicts_dl1_2009_woinarski_1 ####
#'##############################################################################
focal_dataset <- "Predicts_dl1_2009_woinarski_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Body_mass_g, Trophic_level, Habitat_breadth_IUCN, Hand.Wing.Index) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))

trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized", ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, quality = quality,
               trait_dissimilarity = trait_dissimilarity)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 101: Predicts_dl1_2010_proenca_2 ####
#'##############################################################################
focal_dataset <- "Predicts_dl1_2010_proenca_2"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Body_mass_g, Trophic_level, Habitat_breadth_IUCN, Hand.Wing.Index) %>%
  mutate(across(where(is.character), factor)) %>%
  filter(if_any(everything(), ~ !is.na(.)))

trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized", ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, quality = quality, trait_dissimilarity = trait_dissimilarity)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 102: Predicts_dl1_2011_mallari_1 ####
#'##############################################################################
focal_dataset <- "Predicts_dl1_2011_mallari_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Body_mass_g, Trophic_level, Habitat_breadth_IUCN, Hand.Wing.Index) %>%
  mutate(across(where(is.character), factor)) %>%
  filter(if_any(everything(), ~ !is.na(.)))

trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized", ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, quality = quality, trait_dissimilarity = trait_dissimilarity)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 103: Predicts_dl1_2011_moreno_mateos_1 ####
#'##############################################################################
focal_dataset <- "Predicts_dl1_2011_moreno_mateos_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Body_mass_g, Trophic_level, Habitat_breadth_IUCN, Hand.Wing.Index) %>%
  mutate(across(where(is.character), factor)) %>%
  dplyr::select(-Trophic_level) %>%
  filter(if_any(everything(), ~ !is.na(.)))

trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized", ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, quality = quality, trait_dissimilarity = trait_dissimilarity)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 104: Predicts_dl1_2012_dallimer_1 ####
#'##############################################################################
focal_dataset <- "Predicts_dl1_2012_dallimer_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Body_mass_g, Trophic_level, Habitat_breadth_IUCN, Hand.Wing.Index) %>%
  mutate(across(where(is.character), factor)) %>%
  filter(if_any(everything(), ~ !is.na(.)))

trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized", ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, quality = quality, trait_dissimilarity = trait_dissimilarity)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))
#'##############################################################################
#  Dataset 105: Predicts_dl1_2013_bartolommei_1 ####
#'##############################################################################
focal_dataset <- "Predicts_dl1_2013_bartolommei_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Body_mass_g, Trophic_level, Habitat_breadth_IUCN, Hand.Wing.Index) %>%
  mutate(across(where(is.character), factor)) %>%
  dplyr::select(-Trophic_level) %>%
  filter(if_any(everything(), ~ !is.na(.)))

trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized", ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, quality = quality, trait_dissimilarity = trait_dissimilarity)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 106: Predicts_gp1_2007_kutt_2 ####
#'##############################################################################
focal_dataset <- "Predicts_gp1_2007_kutt_2"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Body_mass_g, Trophic_level, Habitat_breadth_IUCN, Hand.Wing.Index) %>%
  mutate(across(where(is.character), factor)) %>%
  filter(if_any(everything(), ~ !is.na(.)))

trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized", ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, quality = quality, trait_dissimilarity = trait_dissimilarity)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 107: Predicts_hb1_2001_aumann_1 ####
#'##############################################################################
focal_dataset <- "Predicts_hb1_2001_aumann_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Body_mass_g, Trophic_level, Habitat_breadth_IUCN, Hand.Wing.Index) %>%
  mutate(across(where(is.character), factor)) %>%
  dplyr::select(-Trophic_level) %>%
  filter(if_any(everything(), ~ !is.na(.)))

trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized", ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, quality = quality, trait_dissimilarity = trait_dissimilarity)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 108: Predicts_hp1_2006_wunderle_1 ####
#'##############################################################################
focal_dataset <- "Predicts_hp1_2006_wunderle_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Body_mass_g, Trophic_level, Habitat_breadth_IUCN, Hand.Wing.Index) %>%
  mutate(across(where(is.character), factor)) %>%
  filter(if_any(everything(), ~ !is.na(.)))

trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized", ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, quality = quality, trait_dissimilarity = trait_dissimilarity)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 109: Predicts_hp1_2007_borges_1 ####
#'##############################################################################
focal_dataset <- "Predicts_hp1_2007_borges_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Body_mass_g, Trophic_level, Habitat_breadth_IUCN, Hand.Wing.Index) %>%
  mutate(across(where(is.character), factor)) %>%
  filter(if_any(everything(), ~ !is.na(.)))

trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized", ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, quality = quality, trait_dissimilarity = trait_dissimilarity)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 110: Predicts_hp1_2007_ranganathan_1 ####
#'##############################################################################
focal_dataset <- "Predicts_hp1_2007_ranganathan_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Body_mass_g, Trophic_level, Habitat_breadth_IUCN, Hand.Wing.Index) %>%
  mutate(across(where(is.character), factor)) %>%
  filter(if_any(everything(), ~ !is.na(.)))

trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized", ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, quality = quality, trait_dissimilarity = trait_dissimilarity)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 111: Predicts_hp1_2007_shahabuddin_1 ####
#'##############################################################################
focal_dataset <- "Predicts_hp1_2007_shahabuddin_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Body_mass_g, Trophic_level, Habitat_breadth_IUCN, Hand.Wing.Index) %>%
  mutate(across(where(is.character), factor)) %>%
  filter(if_any(everything(), ~ !is.na(.)))

trait_dissimilarity <- gawdis(trait_subset, w.type = "optimized", ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, quality = quality, trait_dissimilarity = trait_dissimilarity)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 112: Predicts_hp1_2008_farwig_1 ####
#'##############################################################################
focal_dataset <- "Predicts_hp1_2008_farwig_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Body_mass_g, Trophic_level, Habitat_breadth_IUCN, Hand.Wing.Index) %>%
  mutate(across(where(is.character), factor)) %>%
  filter(if_any(everything(), ~ !is.na(.)))

trait_dissimilarity <- gawdis(trait_subset, w.type = "optimized", ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, quality = quality, trait_dissimilarity = trait_dissimilarity)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 113: Predicts_hp1_2008_ranganathan_1 ####
#'##############################################################################
focal_dataset <- "Predicts_hp1_2008_ranganathan_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Body_mass_g, Trophic_level, Habitat_breadth_IUCN, Hand.Wing.Index) %>%
  mutate(across(where(is.character), factor)) %>%
  filter(if_any(everything(), ~ !is.na(.)))

trait_dissimilarity <- gawdis(trait_subset, w.type = "optimized", ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, quality = quality, trait_dissimilarity = trait_dissimilarity)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 114: Predicts_hp1_2009_kessler_5 ####
#'##############################################################################
focal_dataset <- "Predicts_hp1_2009_kessler_5"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Body_mass_g, Trophic_level, Habitat_breadth_IUCN, Hand.Wing.Index) %>%
  mutate(across(where(is.character), factor)) %>%
  filter(if_any(everything(), ~ !is.na(.)))

trait_dissimilarity <- gawdis(trait_subset, w.type = "optimized", ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, quality = quality, trait_dissimilarity = trait_dissimilarity)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 115: Predicts_hp1_2010_bicknell_1 ####
#'##############################################################################
focal_dataset <- "Predicts_hp1_2010_bicknell_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Body_mass_g, Trophic_level, Habitat_breadth_IUCN, Hand.Wing.Index) %>%
  dplyr::select(-Trophic_level) %>%
  mutate(across(where(is.character), factor)) %>%
  filter(if_any(everything(), ~ !is.na(.)))

trait_dissimilarity <- gawdis(trait_subset, w.type = "optimized", ord = "metric")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, quality = quality, trait_dissimilarity = trait_dissimilarity)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 116: Predicts_hp1_2010_bicknell_2 ####
#'##############################################################################
focal_dataset <- "Predicts_hp1_2010_bicknell_2"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Body_mass_g, Trophic_level, Habitat_breadth_IUCN, Hand.Wing.Index) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))
trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(
  comm = data_clean$comm,
  trait_syndrome = trait_syndrome,
  trait_dissimilarity = trait_dissimilarity,
  quality = quality
)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 117: Predicts_hw1_2005_baldi_1 ####
#'##############################################################################
focal_dataset <- "Predicts_hw1_2005_baldi_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Body_mass_g, Trophic_level, Habitat_breadth_IUCN, Hand.Wing.Index) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))
trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(
  comm = data_clean$comm,
  trait_syndrome = trait_syndrome,
  trait_dissimilarity = trait_dissimilarity,
  quality = quality
)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 118: Predicts_hw1_2007_chapman_1 ####
#'##############################################################################
focal_dataset <- "Predicts_hw1_2007_chapman_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Body_mass_g, Trophic_level, Habitat_breadth_IUCN, Hand.Wing.Index) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))
trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(
  comm = data_clean$comm,
  trait_syndrome = trait_syndrome,
  trait_dissimilarity = trait_dissimilarity,
  quality = quality
)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 119: Predicts_hw1_2008_lantschner_1 ####
#'##############################################################################
focal_dataset <- "Predicts_hw1_2008_lantschner_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Body_mass_g, Trophic_level, Habitat_breadth_IUCN, Hand.Wing.Index) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))
trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(
  comm = data_clean$comm,
  trait_syndrome = trait_syndrome,
  trait_dissimilarity = trait_dissimilarity,
  quality = quality
)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 120: Predicts_hw1_2011_cerezo_1 ####
#'##############################################################################
focal_dataset <- "Predicts_hw1_2011_cerezo_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Body_mass_g, Trophic_level, Habitat_breadth_IUCN, Hand.Wing.Index) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))
trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(
  comm = data_clean$comm,
  trait_syndrome = trait_syndrome,
  trait_dissimilarity = trait_dissimilarity,
  quality = quality
)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 121: Predicts_hw1_2012_naoe_1 ####
#'##############################################################################
focal_dataset <- "Predicts_hw1_2012_naoe_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Body_mass_g, Trophic_level, Habitat_breadth_IUCN, Hand.Wing.Index) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))
trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(
  comm = data_clean$comm,
  trait_syndrome = trait_syndrome,
  trait_dissimilarity = trait_dissimilarity,
  quality = quality
)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 122: Predicts_hw1_2012_naoe_2 ####
#'##############################################################################
focal_dataset <- "Predicts_hw1_2012_naoe_2"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Body_mass_g, Trophic_level, Habitat_breadth_IUCN, Hand.Wing.Index) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))
trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(
  comm = data_clean$comm,
  trait_syndrome = trait_syndrome,
  trait_dissimilarity = trait_dissimilarity,
  quality = quality
)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 123: Predicts_hz1_2012_kutt_1 ####
#'##############################################################################
focal_dataset <- "Predicts_hz1_2012_kutt_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Body_mass_g, Trophic_level, Habitat_breadth_IUCN, Hand.Wing.Index) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))
trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(
  comm = data_clean$comm,
  trait_syndrome = trait_syndrome,
  trait_dissimilarity = trait_dissimilarity,
  quality = quality
)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 124: Predicts_jd1_2002_pearman_1 ####
#'##############################################################################
focal_dataset <- "Predicts_jd1_2002_pearman_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Body_mass_g, Trophic_level, Habitat_breadth_IUCN, Hand.Wing.Index) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))
trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(
  comm = data_clean$comm,
  trait_syndrome = trait_syndrome,
  trait_dissimilarity = trait_dissimilarity,
  quality = quality
)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 125: Predicts_jd1_2010_wang_1 ####
#'##############################################################################
focal_dataset <- "Predicts_jd1_2010_wang_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Body_mass_g, Trophic_level, Habitat_breadth_IUCN, Hand.Wing.Index) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))
trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(
  comm = data_clean$comm,
  trait_syndrome = trait_syndrome,
  trait_dissimilarity = trait_dissimilarity,
  quality = quality
)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 126: Predicts_ks1_2005_pons_1 ####
#'##############################################################################
focal_dataset <- "Predicts_ks1_2005_pons_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Body_mass_g, Trophic_level, Habitat_breadth_IUCN, Hand.Wing.Index) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))
trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, trait_dissimilarity = trait_dissimilarity, quality = quality)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 127: Predicts_ks1_2009_suarez_rubio_1 ####
#'##############################################################################
focal_dataset <- "Predicts_ks1_2009_suarez_rubio_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Body_mass_g, Trophic_level, Habitat_breadth_IUCN, Hand.Wing.Index) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))
trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, trait_dissimilarity = trait_dissimilarity, quality = quality)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 128: Predicts_mh1_2010_sheldon_1 ####
#'##############################################################################
focal_dataset <- "Predicts_mh1_2010_sheldon_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Body_mass_g, Trophic_level, Habitat_breadth_IUCN, Hand.Wing.Index) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))
trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, trait_dissimilarity = trait_dissimilarity, quality = quality)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 129: Predicts_mh1_2011_phalan_1 ####
#'##############################################################################
focal_dataset <- "Predicts_mh1_2011_phalan_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Body_mass_g, Trophic_level, Habitat_breadth_IUCN, Hand.Wing.Index) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))
trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, trait_dissimilarity = trait_dissimilarity, quality = quality)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 130: Predicts_mj1_2008_munyekenye_1 ####
#'##############################################################################
focal_dataset <- "Predicts_mj1_2008_munyekenye_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Body_mass_g, Trophic_level, Habitat_breadth_IUCN, Hand.Wing.Index) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))
trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, trait_dissimilarity = trait_dissimilarity, quality = quality)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 131: Predicts_mj1_2009_lehouck_1 ####
#'##############################################################################
focal_dataset <- "Predicts_mj1_2009_lehouck_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Body_mass_g, Trophic_level, Habitat_breadth_IUCN, Hand.Wing.Index) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))
trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, trait_dissimilarity = trait_dissimilarity, quality = quality)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 132: Predicts_mj1_2009_lehouck_2 ####
#'##############################################################################
focal_dataset <- "Predicts_mj1_2009_lehouck_2"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Body_mass_g, Trophic_level, Habitat_breadth_IUCN, Hand.Wing.Index) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))
trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, trait_dissimilarity = trait_dissimilarity, quality = quality)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 133: Predicts_mj1_2009_lehouck_3 ####
#'##############################################################################
focal_dataset <- "Predicts_mj1_2009_lehouck_3"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Body_mass_g, Trophic_level, Habitat_breadth_IUCN, Hand.Wing.Index) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))
trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, trait_dissimilarity = trait_dissimilarity, quality = quality)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 134: Predicts_mj1_2009_lehouck_4 ####
#'##############################################################################
focal_dataset <- "Predicts_mj1_2009_lehouck_4"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Body_mass_g, Trophic_level, Habitat_breadth_IUCN, Hand.Wing.Index) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))
trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, trait_dissimilarity = trait_dissimilarity, quality = quality)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 135: Predicts_mj1_2009_lehouck_5 ####
#'##############################################################################
focal_dataset <- "Predicts_mj1_2009_lehouck_5"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Body_mass_g, Trophic_level, Habitat_breadth_IUCN, Hand.Wing.Index) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))
trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, trait_dissimilarity = trait_dissimilarity, quality = quality)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 136: Predicts_mj1_2013_ndanganga_1 ####
#'##############################################################################
focal_dataset <- "Predicts_mj1_2013_ndanganga_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Body_mass_g, Trophic_level, Habitat_breadth_IUCN, Hand.Wing.Index) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))
trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, trait_dissimilarity = trait_dissimilarity, quality = quality)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 137: Predicts_mj1_2013_ndanganga_2 ####
#'##############################################################################
focal_dataset <- "Predicts_mj1_2013_ndanganga_2"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Body_mass_g, Trophic_level, Habitat_breadth_IUCN, Hand.Wing.Index) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))
trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, trait_dissimilarity = trait_dissimilarity, quality = quality)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 138: Predicts_mj1_2013_reynolds_1 ####
#'##############################################################################
focal_dataset <- "Predicts_mj1_2013_reynolds_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Body_mass_g, Trophic_level, Habitat_breadth_IUCN, Hand.Wing.Index) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))
trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, trait_dissimilarity = trait_dissimilarity, quality = quality)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 139: Predicts_sc1_2005_marsh_1 ####
#'##############################################################################
focal_dataset <- "Predicts_sc1_2005_marsh_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Body_mass_g, Trophic_level, Habitat_breadth_IUCN, Hand.Wing.Index) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))
trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, trait_dissimilarity = trait_dissimilarity, quality = quality)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 140: Predicts_sc1_2010_marsh_2 ####
#'##############################################################################
focal_dataset <- "Predicts_sc1_2010_marsh_2"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Body_mass_g, Trophic_level, Habitat_breadth_IUCN, Hand.Wing.Index) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))
trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, trait_dissimilarity = trait_dissimilarity, quality = quality)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 141: Predicts_sc1_2010_rey_benayas_1 ####
#'##############################################################################
focal_dataset <- "Predicts_sc1_2010_rey_benayas_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Body_mass_g, Trophic_level, Habitat_breadth_IUCN, Hand.Wing.Index) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))
trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, trait_dissimilarity = trait_dissimilarity, quality = quality)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 142: Predicts_sc1_2011_stouffer_1 ####
#'##############################################################################
focal_dataset <- "Predicts_sc1_2011_stouffer_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Body_mass_g, Trophic_level, Habitat_breadth_IUCN, Hand.Wing.Index) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))
trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, trait_dissimilarity = trait_dissimilarity, quality = quality)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 143: Predicts_sc2_2012_santana_1 ####
#'##############################################################################
focal_dataset <- "Predicts_sc2_2012_santana_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Body_mass_g, Trophic_level, Habitat_breadth_IUCN, Hand.Wing.Index) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))
trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, trait_dissimilarity = trait_dissimilarity, quality = quality)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 144: Predicts_se2_2013_brandt_1 ####
#'##############################################################################
focal_dataset <- "Predicts_se2_2013_brandt_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Body_mass_g, Trophic_level, Habitat_breadth_IUCN, Hand.Wing.Index) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))
trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, trait_dissimilarity = trait_dissimilarity, quality = quality)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 145: Predicts_se2_2013_hassan_1 ####
#'##############################################################################
focal_dataset <- "Predicts_se2_2013_hassan_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Body_mass_g, Trophic_level, Habitat_breadth_IUCN, Hand.Wing.Index) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))
trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, trait_dissimilarity = trait_dissimilarity, quality = quality)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))
#'##############################################################################
#  Dataset 146: Predicts_sh1_2012_ims_1 ####
#'##############################################################################
focal_dataset <- "Predicts_sh1_2012_ims_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Body_mass_g, Trophic_level, Habitat_breadth_IUCN, Hand.Wing.Index) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))
trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, trait_dissimilarity = trait_dissimilarity, quality = quality)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 147: Predicts_tn1_2007_o_dea_1 ####
#'##############################################################################
focal_dataset <- "Predicts_tn1_2007_o_dea_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Body_mass_g, Trophic_level, Habitat_breadth_IUCN, Hand.Wing.Index) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))
trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, trait_dissimilarity = trait_dissimilarity, quality = quality)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 148: Predicts_vk1_2007_st_laurent_3 ####
#'##############################################################################
focal_dataset <- "Predicts_vk1_2007_st_laurent_3"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Body_mass_g, Trophic_level, Habitat_breadth_IUCN, Hand.Wing.Index) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))
trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, trait_dissimilarity = trait_dissimilarity, quality = quality)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 149: Predicts_vk1_2011_edenius_1 ####
#'##############################################################################
focal_dataset <- "Predicts_vk1_2011_edenius_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Body_mass_g, Trophic_level, Habitat_breadth_IUCN, Hand.Wing.Index) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))
trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, trait_dissimilarity = trait_dissimilarity, quality = quality)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 150: Predicts_vk1_2011_zimmerman_1 ####
#'##############################################################################
focal_dataset <- "Predicts_vk1_2011_zimmerman_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Body_mass_g, Trophic_level, Habitat_breadth_IUCN, Hand.Wing.Index) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))
trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, trait_dissimilarity = trait_dissimilarity, quality = quality)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 151: Predicts_vk1_2012_otto_1 ####
#'##############################################################################
focal_dataset <- "Predicts_vk1_2012_otto_1"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Body_mass_g, Trophic_level, Habitat_breadth_IUCN, Hand.Wing.Index) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))
trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, trait_dissimilarity = trait_dissimilarity, quality = quality)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 152: Predicts_vk1_2012_otto_2 ####
#'##############################################################################
focal_dataset <- "Predicts_vk1_2012_otto_2"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Body_mass_g, Trophic_level, Habitat_breadth_IUCN, Hand.Wing.Index) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))
trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, trait_dissimilarity = trait_dissimilarity, quality = quality)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 153: ribera_beetles ####
#'##############################################################################
focal_dataset <- "ribera_beetles"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(contains("SIZE"), contains("DIET"), contains("HABITAT"), contains("MOBILITY")) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))
trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, trait_dissimilarity = trait_dissimilarity, quality = quality)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 154: romero_waterinsects ####
#'##############################################################################
focal_dataset <- "romero_waterinsects"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(contains("SIZE"), contains("DIET"), contains("HABITAT"), contains("MOBILITY")) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))
trait_dissimilarity <- gawdis::gawdis(trait_subset,
                                      groups = c(rep(1, 1), rep(2, 1), rep(3, 1), rep(4, 7)),
                                      fuzzy = 4,
                                      w.type = "optimized")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, trait_dissimilarity = trait_dissimilarity, quality = quality)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 155: S29TTP ####
#'##############################################################################
focal_dataset <- "S29TTP"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(contains("SIZE"), contains("DIET"), StemDens.mean, contains("MOBILITY")) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))
trait_dissimilarity <- gawdis::gawdis(trait_subset,
                                      groups = c(rep(1, 1), rep(2, 1), rep(3, 1), rep(4, 1)),
                                      w.type = "optimized")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, trait_dissimilarity = trait_dissimilarity, quality = quality)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))
#'##############################################################################
#  Dataset 156: S47TTP ####
#'##############################################################################
focal_dataset <- "S47TTP"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Height, SLA, T:N, Seed_mass) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))
trait_dissimilarity <- gawdis::gawdis(trait_subset, groups = c(rep(1,1), rep(2,1), rep(3,5), rep(4,1)), w.type = "optimized")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, trait_dissimilarity = trait_dissimilarity, quality = quality)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 157: shieh_waterinsects ####
#'##############################################################################
focal_dataset <- "shieh_waterinsects"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(Size5:Size50, FP:Mai, SW:PA, SL:SP) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))
trait_dissimilarity <- gawdis::gawdis(trait_subset, groups = c(rep(1,5), rep(2,6), rep(3,5), rep(4,4)), fuzzy = 1:4, w.type = "optimized")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, trait_dissimilarity = trait_dissimilarity, quality = quality)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 158: stavert_bees ####
#'##############################################################################
focal_dataset <- "stavert_bees"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(contains("SIZE"), contains("DIET"), contains("HABITAT"), contains("MOBILITY")) %>%
  mutate(across(where(is.character), factor)) %>%
  mutate(across(contains("MOBILITY"), ordered)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))
trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, trait_dissimilarity = trait_dissimilarity, quality = quality)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 159: toxywa_beetles_nonsaproxylic ####
#'##############################################################################
focal_dataset <- "toxywa_beetles_nonsaproxylic"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(contains("SIZE"), contains("DIET"), contains("HABITAT")) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))
trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, trait_dissimilarity = trait_dissimilarity, quality = quality)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 160: toxywa_beetles_saproxylic ####
#'##############################################################################
focal_dataset <- "toxywa_beetles_saproxylic"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(contains("SIZE"), contains("DIET"), contains("HABITAT"), contains("MOBILITY")) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))
trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, trait_dissimilarity = trait_dissimilarity, quality = quality)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 161: vanklink_bees ####
#'##############################################################################
focal_dataset <- "vanklink_bees"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(contains("SIZE"), contains("feeding"), contains("Nesting"), contains("flight")) %>%
  mutate(across(where(is.character), factor)) %>%
  mutate(across(contains("feeding"), ordered)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))
trait_dissimilarity <- gawdis::gawdis(trait_subset, groups = c(rep(1,4), 2, 3, 4), w.type = "optimized")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, trait_dissimilarity = trait_dissimilarity, quality = quality)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 162: vanklink_carabids ####
#'##############################################################################
focal_dataset <- "vanklink_carabids"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(contains("SIZE"), contains("trophic"), contains("HABITAT"), contains("wing")) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))
trait_dissimilarity <- gawdis::gawdis(trait_subset, groups = c(rep(1,2), 2, 3, 4), w.type = "optimized")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, trait_dissimilarity = trait_dissimilarity, quality = quality)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 163: vanklink_hoverflies ####
#'##############################################################################
focal_dataset <- "vanklink_hoverflies"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(contains("SIZE"), contains("feeding"), contains("substrate"), contains("wing")) %>%
  mutate(across(where(is.character), factor)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))
trait_dissimilarity <- gawdis::gawdis(trait_subset, groups = c(rep(1,2), 2, 3, 4), w.type = "optimized")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, trait_dissimilarity = trait_dissimilarity, quality = quality)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 164: vanklink_staphylinids ####
#'##############################################################################
focal_dataset <- "vanklink_staphylinids"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(contains("length"), contains("feeding"), contains("humidity"), contains("dispersal")) %>%
  mutate(across(where(is.character), factor)) %>%
  mutate(across(contains("dispersal"), ordered)) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))
trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, trait_dissimilarity = trait_dissimilarity, quality = quality)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 165: yates_ants ####
#'##############################################################################
focal_dataset <- "yates_ants"
data_clean <- glue::glue("S1_Preprocessing/Processed/{focal_dataset}_clean.rds") %>% read_rds()
data_trait <- data_clean %>% pluck("trait")
trait_subset <- data_trait %>%
  dplyr::select(contains("thorax"), contains("ml"), contains("sl"), contains("femur")) %>%
  slice(which(apply(., 1, function(x) !all(is.na(x)))))
trait_dissimilarity <- gawdis::gawdis(trait_subset, w.type = "optimized")
trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, trait_dissimilarity = trait_dissimilarity, quality = quality)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))

#'##############################################################################
#  Dataset 166: N21FBD_2 ####
#'##############################################################################
focal_dataset <- "N21FBD_2"
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
                                       , w.type="optimized"
                                       ,ord = "metric")

trait_syndrome <- ape::pcoa(trait_dissimilarity, correction = "lingoes")$vectors.cor[, 1:3]
quality <- calculate_trait_space_quality(trait_dissimilarity, trait_syndrome)
output <- list(comm = data_clean$comm, trait_syndrome = trait_syndrome, trait_dissimilarity = trait_dissimilarity, quality = quality)
saveRDS(output, glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds"))