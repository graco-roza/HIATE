#Handle data from predicts 
library(tidyverse)
library(readxl)
library(xlsx)
library(lubridate)

predicts_data <- readRDS("PREDICTs/Predicts_merged_sites.rds")
predicts_traits <- readRDS("Ongoing manuscripts/HIATES/PREDICTs/VertebrateTraits_Caio.rds")
AVONET_traits <- read_excel("PREDICTs/AVONET_traits.xlsx", sheet =4)
imputed_traits <- readRDS("PREDICTs/ImputedBirds.rds")

birds_traits <-
  imputed_traits %>%
  pluck(5) %>% 
    left_join(
    AVONET_traits %>%
      dplyr::select(Species3, `Hand-Wing.Index`)
    ,by = c("Best_guess_binomial" = "Species3")
  ) %>%
  mutate(species = janitor::make_clean_names(Best_guess_binomial))

birds_traits %>% dplyr::select("Hand-Wing.Index") %>%  summary()

birds_datasets <- predicts_data %>%
  filter(Class == "Aves") %>%
  dplyr::select(
    SS,
    Site_name,
    Predominant_land_use,
    Sample_start_earliest,
    Sample_end_latest,
    Longitude,
    Latitude,
    Site_name,
    Best_guess_binomial,
    Measurement
  ) %>%
  mutate(SS = as.character(SS)
         #,Measurement = ifelse(Measurement > 0, 1, 0) #to convert all to presence-absence (not needed!)
         ) %>%
  split(.$SS) %>%
  map(
    ~ .x %>%
      pivot_wider(
       # id_cols = c(Site_name, Latitude, Longitude, Predominant_land_use,Sample_start_earliest,Sample_end_latest), #(not needed!)
        names_from = Best_guess_binomial,
        values_from = Measurement,
        values_fn = max
      ) %>%
      janitor::clean_names()
  )

birds_datasets_final <- birds_datasets %>% 
  keep(function(x) nrow(x) > 10) #only more than 10 sites 


for (i in 1:length(birds_datasets_final)){
print(i)
comm <-   birds_datasets_final %>% 
  pluck(i) %>% 
  dplyr::select(-c(ss,predominant_land_use,sample_start_earliest,sample_end_latest, longitude, latitude)) %>% 
  rename(site = site_name) %>% 
  data.frame()

trait <- birds_traits %>% 
  filter(species %in% names(comm)[-1]) %>% 
  relocate(species, .before = "Trophic_level")

coord <- birds_datasets_final %>% 
  pluck(i) %>% 
  dplyr::select(c(site_name, longitude, latitude)) %>% 
  rename(site = site_name, x= longitude, y= latitude) %>% 
  data.frame()


env <-  birds_datasets_final %>% 
  pluck(i) %>% 
  dplyr::select(c(site_name, predominant_land_use,sample_start_earliest,sample_end_latest)) %>% 
  rename(site = site_name, sample_start = sample_start_earliest, sample_end = sample_end_latest) %>% 
  data.frame() %>% 
  mutate(across(c(sample_start,sample_end),year))

wb = createWorkbook()

sheet = createSheet(wb, "species")
addDataFrame(comm, sheet=sheet, startColumn=1, row.names=FALSE)

sheet = createSheet(wb, "traits")
addDataFrame(trait, sheet=sheet, startColumn=1, row.names=FALSE)

sheet = createSheet(wb, "coordinates")
addDataFrame(coord, sheet=sheet, startColumn=1, row.names=FALSE)

sheet = createSheet(wb, "environment")
addDataFrame(env, sheet=sheet, startColumn=1, row.names=FALSE)


saveWorkbook(wb, glue::glue("raw_predicts/Predicts_{janitor::make_clean_names(names(birds_datasets_final)[i])}.xlsx"))

}


data_info <- read_csv("S1_Preprocessing/dataset_info_all.csv")

start <- birds_datasets_final %>% 
  map( ~ .x %>%  pull(sample_start_earliest) %>% year %>%   min) %>% 
  bind_rows(.id = "study") %>% 
  t %>%  c

end <- birds_datasets_final %>% 
  map( ~ .x %>%  pull(sample_end_latest) %>% year %>%   max) %>% 
  bind_rows(.id = "study") %>% 
  t %>%  c

new_info <- data.frame(dataset = paste("pb",1:66),
                       dataset.name = c(glue::glue("Predicts_{janitor::make_clean_names(names(birds_datasets_final))}")),
                       taxa = "birds",
                       group1 = "vertebrates",
                       group2 = "terrestrial vertebrates",
                       Start_year = start ,
                       End_year = end
                       )
data_info %>% 
  full_join(new_info) %>% 
  write_csv("dataset_info_all.csv")



new_varying <- read_csv("S1_Preprocessing/dataset_info_all.csv") %>% 
  filter(grepl("Predicts",dataset_name)) %>% 
  pull(dataset_name) %>% 
  data.frame(dataset= .)
i=3
read.xlsx(file = glue::glue("S1_Preprocessing/raw_data/{new_varying$dataset[i]}.xlsx"), sheetName="environment") %>% 
  pull(predominant_land_use) %>%  table

new_varying$disturbance[i] <- "forestry"



files <- tools::file_path_sans_ext(list.files("S1_Preprocessing/raw_data"))

test <- files[grepl("Predicts",files)] |> 
  map(~ glue::glue("S1_Preprocessing/raw_data/{.x}.xlsx") |> 
        read_excel(sheet = "environment",na=c("NA","","#N/A")) |>  
        mutate(dataset = .x) |> 
        group_by(dataset) |> 
        count(predominant_land_use)) |> 
  bind_rows() |> 
  pivot_wider(names_from = "predominant_land_use", values_from="n", values_fill = 0) |> 
  mutate(Forestry = 
           `Secondary vegetation (indeterminate age)` +
           `Intermediate secondary vegetation` +
           `Young secondary vegetation` + 
           `Plantation forest`+
           `Mature secondary vegetation`,
         Agriculture = Cropland + Pasture ,
         Varying = `Cannot decide`,
         neutral =  `Primary vegetation`,
         .keep = "unused"
           ) |> 
  mutate(Total = Agriculture + Forestry + Urban + Varying) |> 
  mutate(across(!Total, ~ (.x/Total)*100)) 

test |>  print(n=66)

lu_pred<-test |> 
  pivot_longer(!dataset, names_to = "lu", values_to = "pct") |> 
  filter(lu != "Total") |> 
  group_by(dataset) |> 
  slice_max(pct) |> 
  mutate(highest = case_when(pct >= 100 ~ lu, TRUE ~"Varying")) |> 
  select(dataset,highest) |> 
  set_names(c("dataset_name","highest"))

lu_pred |>  ungroup() |> count(highest)

dataset_info <- read_csv("S1_Preprocessing/dataset_info_all.csv")

new_datainfo <- dataset_info |> 
  full_join(lu_pred) |> 
  mutate(disturbance_type = case_when(!is.na(highest) ~ tolower(highest),TRUE ~  disturbance_type)) 

write_csv(new_datainfo,"S1_Preprocessing/dataset_info_all.csv")



files <- tools::file_path_sans_ext(list.files("S1_Preprocessing/raw_data"))

test <- files |> 
  map(~ glue::glue("S1_Preprocessing/raw_data/{.x}.xlsx") |> 
        read_excel(sheet = "coordinates",na=c("NA","","#N/A")) |>  
        mutate(dataset = .x) |> 
        group_by(dataset))

test[-158] |> 
  map(~ write_csv(x= .x, file= glue::glue("csv_coord/{unique(.x$dataset)}.csv")))
