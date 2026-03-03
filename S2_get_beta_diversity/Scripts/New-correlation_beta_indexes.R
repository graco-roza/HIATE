pacman::p_load(
  BAT #make taxonomic beta diversity
  ,hypervolume #make functional beta diversity
  ,tidyverse #manage data
  ,magrittr #use set names
  ,glue #create strings
  ,readxl #read rds files
  ,future #set up cores available
  ,doSNOW # make parallel procedure
)

files <- tools::file_path_sans_ext(list.files("S2_get_beta_diversity/betadiv_input")) #get a vector with dataset names

final_res <- data.frame(
  dataset = character(),
  PA_corr = numeric(),
  ABUND_corr = numeric(),
  Morisita_corr = numeric(),
  pa_mean =numeric(),
  pa_25 = numeric(),
  pa_75 = numeric(),
  abund_mean = numeric(),
  abund_25 = numeric(),
  abund_75 = numeric()
)

for (ii in 1:length(files)){
print(ii)
focal_dataset <- gsub("_beta_Input","",files[ii]) #chose one dataset 
species_diff <- glue::glue("S2_get_beta_diversity/betadiv_input/{focal_dataset}_beta_Input.rds") %>%  read_rds() #read the dataset pre processed file

comm <- species_diff %>%  pluck("comm") #extract the community data
traits <- species_diff %>%  pluck("trait_syndrome") #extract the synthetic traits (PCoA Axes)

#make sure the community only has species with trait values. 
comm <- comm %>%  select(any_of(rownames(traits)))

comm<-comm[rowSums(comm) > 0,]
#Estimate Taxonomic Beta diversity ---------------------------------------------------------------------------------
podani_pa <- BAT::beta(vegan::decostand(comm,"pa"), abund = FALSE, func = "sorensen")$Brepl
podani_abund <- BAT::beta(comm, abund = TRUE, func = "sorensen")$Brepl 

baselga_pa <- betapart::beta.pair(vegan::decostand(comm,"pa"), index.family="sorensen")$beta.sim
baselga_abun <- betapart::beta.pair.abund(comm, index.family="bray")$beta.bray.bal

morisita_horn <- vegan::vegdist(comm, "horn")

comparison_pa <-vegan::mantel(podani_pa,baselga_pa)
comparison_abun <- vegan::mantel(podani_abund,baselga_abun)
comparison_index <- vegan::mantel(podani_abund, morisita_horn)

df<- data.frame(
 dataset =  focal_dataset,
 PA_corr = ifelse(comparison_pa$signif <= 0.05, comparison_pa$statistic, NA),
 ABUND_corr = ifelse(comparison_abun$signif <= 0.05, comparison_abun$statistic, NA),
 Morisita_corr = ifelse(comparison_index$signif <= 0.05, comparison_index$statistic, NA),
 pa_mean = mean(podani_pa),
 pa_25 = quantile(podani_pa,.025),
 pa_75 = quantile(podani_pa,.975),
 abund_mean = mean(podani_abund),
 abund_25 = quantile(podani_abund,.025),
 abund_75 = quantile(podani_abund,.975)
)

final_res[ii,] <-  df

}


final_res |> 
  pivot_longer(cols = c(PA_corr,ABUND_corr,Morisita_corr), names_to = "feature", values_to="corr") |> 
  ggplot(aes(x=corr, fill=feature))+
  geom_density(alpha=.5, colour="white")+
  theme_minimal()+
  labs(y="Density", x = "Mantel correlation")+


  scale_fill_viridis_d(labels=c("Abundance","Morisita","Presence-absence"), direction=-1)
