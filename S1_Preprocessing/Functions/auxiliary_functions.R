###############################################################################
# SCRIPT NAME: auxiliary_functions.R
#
# DESCRIPTION:
#   This script contains auxiliary functions used throughout the project for
#   data cleaning, processing, and extraction of predictors such as human
#   footprint, MODIS land cover, and climate data.
#
# USAGE:
#   This script is sourced by other scripts in the workflow.
#
# INPUTS:
#   - Raw data files: species data, trait data, coordinates, and environment data
#   - Spatial raster files: Human Footprint (2009), WorldClim bioclimatic variables
#   - MODIS product: MCD12Q1 (land cover data)
#
# OUTPUTS:
#   - Cleaned datasets: community, traits, coordinates, and environment
#   - Predictors extracted for each site: Human Footprint, MODIS, and climate data
#
# AUTHOR: [Caio Graco-Roza]
# LAST UPDATED: [20.11.2024]
#
# NOTES:
#   - This script uses libraries such as `dplyr`, `raster`, `sf`, and `MODISTools`
#   - Ensure all input files are properly formatted and projections are consistent.
###############################################################################


# Function to clean all raw datasets
clean_data <- function(focal_dataset) { 
  # Define the path to the dataset
  dataset_path <- glue::glue("S1_Preprocessing/raw_data/{focal_dataset}.xlsx")
  
  # Load data from the dataset
  comm_original <- dataset_path %>% read_excel(sheet = "species", na = c("NA", "", "#N/A")) %>% column_to_rownames("site") # Load community data
  trait_original <- dataset_path %>% read_excel(sheet = "traits", na = c("NA", "", "#N/A")) # Load trait data
  coord_original <- dataset_path %>% read_excel(sheet = "coordinates", na = c("NA", "", "#N/A")) # Load coordinates
  env_original <- dataset_path %>% read_excel(sheet = "environment", na = c("NA", "", "#N/A")) # Load environmental data
  
  # Filter and clean community data
  comm <- comm_original %>%  
    slice(which(rownames(comm_original) %in% coord_original$site)) %>% # Keep only sites with coordinates
    dplyr::select(any_of(trait_original$species)) %>%  # Keep only species with trait data
    slice(which(specnumber(.) > 2)) %>% # Remove sites with fewer than 2 species because functional diversity metrics require at least 3 species to compute a hypervolume
    dplyr::select(which(!colSums(., na.rm = TRUE) %in% 0)) # Remove species with zero abundance across all sites because they are irrelevant to the study
  
  # Filter and clean trait data
  trait <- trait_original %>%  
    filter(species %in% colnames(comm)) %>%  # Keep only species present in the community data
    column_to_rownames("species") %>%  # Set species as row names
    data.frame()
  
  # Filter and clean coordinate data
  coord <- coord_original %>%  
    filter(site %in% rownames(comm))  # Keep only coordinates for the retained sites
  
  # Filter and clean environmental data
  env <- env_original %>%  
    filter(site %in% rownames(comm))  # Keep only environmental data for the retained sites
  
  # Combine cleaned data into a list
  data_parsed <- list(comm = comm,
                      trait = trait, 
                      coord = coord, 
                      env = env)
  
  return(data_parsed) # Return cleaned data
}

# Function to extract Human Footprint data
hfp_get <- function(coordinates, buffers = c(1000, 1500, 2000)) {
  require(raster)
  require(sf)
  
  # Load the Human Footprint raster
  humanfootprint <- raster::raster("S1_Preprocessing/Miscellaneous/Human_footprint/wildareas-v3-2009-human-footprint_geo.tif")
  
  # Convert coordinates to spatial object and set projection
  coordinates(coordinates) <- ~x + y
  proj4string(coordinates) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") # WGS84 ensures consistent projection across all data
  coordinates <- spTransform(coordinates, raster::projection(humanfootprint))
  
  # Extract mean human footprint for each buffer
  # Buffers allow testing how different spatial extents of human footprint affect predictions
  hfp <- lapply(buffers, function(i) {
    raster::extract(humanfootprint, coordinates, buffer = i, fun = mean, na.rm = TRUE)
  }) %>%
    set_names(glue::glue("hfp_{buffers}")) %>%
    bind_cols() %>% 
    add_column(site = coordinates$site) %>% 
    relocate(site)
  
  return(hfp)
}

# Function to extract MODIS land cover data
modis_get <- function(coordinates, dataset, buffers = seq(500, 2000, by = 500)) { 
  # Retrieve start and end years from dataset info
  start_year <- data_info %>% filter(dataset_name == dataset) %>% pull(start_year)
  final_year <- data_info %>% filter(dataset_name == dataset) %>% pull(end_year)
  
  # Download MODIS raster data
  write("download modis raster", stderr())
  modis_df <- mt_batch_subset(
    df = coordinates %>% rename(site_name = "site", lat = "y", lon = "x"),
    product = "MCD12Q1",
    band = "LC_Type1",
    start = as.character(glue::glue("{start_year}-01-01")),
    end = as.character(glue::glue("{final_year}-12-30")),
    km_lr = 3, km_ab = 3, internal = TRUE
  )
  
  # Convert MODIS data into rasters grouped by year
  write("make raster modis", stderr())
  raster_modis <- modis_df %>% 
    group_split(calendar_date) %>% 
    purrr::map(~ .x %>% group_split(site) %>% purrr::map(~MODISTools::mt_to_terra(df = .x))) %>% 
    c() %>% 
    map(~terra::merge(terra::sprc(.x)))
  
  # Convert coordinates to spatial object and reproject
  write("convert coordinates into spatial object", stderr())
  coordinates(coordinates) <- ~x + y
  proj4string(coordinates) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  coordinates <- spTransform(coordinates, crs(raster_modis[[1]], proj = TRUE))
  
  # Extract land cover proportions for human-altered classes (e.g., buildings, croplands)
  # Multiple buffers allow sensitivity analysis to spatial scale
  write("final calculation of MODIS", stderr())
  modis <- list()
  for (i in seq_along(raster_modis)) {
    modis[[i]] <- lapply(buffers, function(x) {
      raster::extract(raster(raster_modis[[i]]), coordinates, buffer = x) %>%
        purrr::map(~data.frame(class = .x)) %>%
        bind_rows(.id = "site") %>%
        group_by(site) %>%
        summarise(modis = sum(freq, na.rm = TRUE), heterogeneity = mean(heterogeneity, na.rm = TRUE)) %>%
        rename(!!glue("modis_{x}") := modis, !!glue("het_{x}") := heterogeneity)
    }) %>% Reduce(f = full_join, x = .)
  }
  
  modis_mean <- modis %>% bind_rows() %>% group_by(site) %>% 
    summarise(across(starts_with(c("modis", "het")), ~mean(.x, na.rm = TRUE)))
  
  return(modis_mean)
}

# Function to extract climate data
climate_get <- function(coord) { 
  require(raster)
  require(sp)
  
  # Load climate rasters
  temp <- raster::raster("S1_Preprocessing/Miscellaneous/wc10/wc2.1_10m/wc2.1_10m_bio_1.tif")
  prec <- raster::raster("S1_Preprocessing/Miscellaneous/wc10/wc2.1_10m/wc2.1_10m_bio_12.tif")
  
  # Convert coordinates to spatial object and reproject
  sp::coordinates(coord) <- ~x + y
  sp::proj4string(coord) <- sp::CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  coord <- sp::spTransform(coord, raster::projection(temp))
  
  # Extract temperature and precipitation data
  temp_extracted <- raster::extract(temp, coord) %>% data.frame()
  prec_extracted <- raster::extract(prec, coord) %>% data.frame()
  
  # Combine into a single data frame
  climate <- data.frame(Temp = temp_extracted[,1], Prec = prec_extracted[,1])
  return(climate)
}

# Function to retrieve all predictors (HFP, MODIS, Climate)
get_predictors <- function(coordinates, dataset_name) {
  write("estimate HFP", stderr())
  hfp <- hfp_get(coordinates) %>% mutate(site = as.character(site))
  
  write("estimate MODIS", stderr())
  modis <- modis_get(coordinates, dataset_name) %>% mutate(site = as.character(site))
  
  write("estimate climate", stderr()) 
  climate <- coordinates %>%
    dplyr::select(x, y) %>% 
    climate_get() %>% 
    add_column(site = coordinates$site) %>% 
    mutate(site = as.character(site))
  
  # Combine all predictors
  predictors <- Reduce(full_join, list(coordinates %>% mutate(site = as.character(site)), hfp, modis, climate))
  return(predictors)
}






