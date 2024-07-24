#!/usr/bin/env Rscript

#install.packages("dplyr")
#install.packages("tidyr")
library("dplyr")
library("tidyr")

getwd()
setwd("/Users/julia/Documents/Documents/SupBiotech/BT5/Stage/Codestuff/Faure_coding_material/Homogeneity_Factor")
getwd()

CC_Abund_80=read.table('SMAGs_CC_Abund_80_all_noOLD_SAGs', header=F)
colnames(CC_Abund_80)<-c("CC_ID","Genes")
annots_egg = read.delim("All_eggNOG_Functional_description_CC80_NA", header = T, sep="\t", quote="",  na.strings = "")

CC_Annot = merge(x=CC_Abund_80, y=annots_egg, by="Genes", all.x=T)
CC_Annot = CC_Annot[,c(2,3)]

#####SIZE
cc_size = CC_Annot %>%
  group_by(CC_ID) %>%
  do(data.frame(size = nrow(.)))
summary(cc_size)

#creating a new, singleton-free version of CC_Abund_80
list_of_singletons = cc_size[cc_size$size == 1,] #is a data frame
unwanted_CCs = list_of_singletons$CC_ID

CC_Abund_80_filtered = CC_Abund_80[!(CC_Abund_80$CC_ID %in% unwanted_CCs), ]
CC_Abund_80_filtered_second_attempt <- CC_Abund_80 %>%
  filter(!CC_ID %in% unwanted_CCs)

write.table(CC_Abund_80_filtered_second_attempt, "SMAGs_CC_Abund_80_all_noOLD_SAGs_singleton_free", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


#### Taxonomical indices ####
cc_taxo=read.table('CC_80_taxo', header=T, sep=" ")
cc_taxo_filtered = cc_taxo %>% 
  filter(!CC_ID %in% unwanted_CCs)
write.table(cc_taxo_filtered, "CC_80_taxo_singleton_free", sep = " ", quote = FALSE, row.names = FALSE, col.names = TRUE)
#with " " as the sep, because that was the separator in CC_80_taxo, not a tab

