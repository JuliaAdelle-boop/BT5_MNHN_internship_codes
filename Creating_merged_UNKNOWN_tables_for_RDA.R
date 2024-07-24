library(dplyr)
library(VIM)
library(stringr)
#Opening the Size Fraction
SSUU <- read.table("SSUU_all_mags_indexed_GeneID_corr_final_UNKNOWNS", header = T, quote = "", sep = "\t", row.names = NULL)

#Make a subset of MMQQ to make it lighter and faster to handle
comp_SSUU <- subset(SSUU, select = c(Genes, Water, MAG_Cov_Y, Gene_Expression, Gene_Expression_Index, CC))

#Decomposing the 'Water' column
comp_SSUU_scission <- comp_SSUU %>%
  mutate(Water = as.character(Water),
         Station = str_extract(Water, "^[0-9]+"),
         Rest = str_extract(Water, "[A-Za-z]+$"),
         Depth = str_sub(Rest, 1, 3),
         Size = str_sub(Rest, 4, 7),
         Water = NULL,
         Rest = NULL)
Final_comp_SSUU_scission <- comp_SSUU_scission %>%
  mutate(Genes = NULL,
         Depth = NULL,
         Size = NULL)

#Opening that contains the environmental data
Tara_insitu <- read.table("SMAGS_insitu_WOA_Env.tsv", sep="\t", header=T)

#removal of the last row, because it doesn't contain any data
Tara_insitu <- Tara_insitu[1:244,]

#removal of the environmental parameters that we don't want (these are not yet the ones decided by the PCA)
drop.cols <- c('Ocean.region','Marine.biome', 'Moon.phase.prop','Season','DepthBathy','Depth.top.min','Depth.bottom.max','PO4_WOA', 'Si_WOA', 'NO3_WOA','Part.beam.att.coef', 'Latitude_WOA', 'Longitude_WOA', 'Temperature_WOA', 'Salinity_WOA', 'Part.backscat.coef')
Tara_insitu <- Tara_insitu %>%
  select(-one_of(drop.cols))

#Preparing Tara_insitu for kNN Imputation: making sure the values of each column are correctly formatted and of the right class
sapply(Tara_insitu, class)
Tara_insitu <- Tara_insitu %>%
  mutate(Station = str_sub(Station, 6),
         Station = str_extract(Station, "^[0-9]+"),
         Station = as.numeric(Station),
         Station = as.character(Station),
         Depth.nominal = as.character(Depth.nominal),
         Sunshine.duration = as.character(Sunshine.duration),
         Depth.Mixed.Layer = as.character(Depth.Mixed.Layer),
         Depth.chloro.max = as.character(Depth.chloro.max),
         Depth.max.Brunt.V..is..l.. = as.character(Depth.max.Brunt.V..is..l..),
         Depth.Max.O2 = as.character(Depth.Max.O2),
         Depth.Min.O2 = as.character(Depth.Min.O2))

#kNN Imputation of Tara_insitu's missing values (no z-centering done, at this stage; done later, just before the RDA)
Imputed_Tara_insitu <- kNN(Tara_insitu, imp_var = FALSE)

#Creation of a table containing both Abundance data and Environmental data
Full_SSUU_trimmed <- merge(Final_comp_SSUU_scission, Imputed_Tara_insitu, by="Station", all.x=F)
write.table(Full_SSUU_trimmed, "Merged_UNKNOWNS_SSUU_for_RDA", quote = FALSE, sep = "\t")
