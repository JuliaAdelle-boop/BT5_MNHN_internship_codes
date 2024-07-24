library(dplyr)
library(tidyr)

CC_Abund_80=read.table('SMAGs_CC_Abund_80_all_noOLD_SAGs', header=F)
colnames(CC_Abund_80)<-c("CC_ID","Genes")
annots_egg = read.delim("All_eggNOG_Functional_description_CC80_NA", header = T, sep="\t", quote="",  na.strings = "")
KOprot = read.delim("CC_KEGGKO_8080_NA", header = T, na.strings = "", sep = "\t")
PFAM = read.delim("CC_PFAM_8080_NA", header = T, na.strings = "", sep = "\t")
colnames(PFAM)<-c("Genes","Annot")


CC_Annot = merge(x=CC_Abund_80, y=annots_egg, by="Genes", all.x=T)
CC_KEGG = merge(x=CC_Abund_80, y=KOprot, by="Genes", all.x=T)
sum(is.na(CC_KEGG$Annot))
CC_PFAM = merge(x=CC_Abund_80, y=PFAM, by="Genes", all.x=T)
sum(is.na(CC_PFAM$Annot))

CC_Annot = CC_Annot[,c(2,3)]
CC_PFAM = CC_PFAM[,c(2,3)]
CC_KEGG = CC_KEGG[,c(2,3)]

#####SIZE
cc_size = CC_Annot %>%
  group_by(CC_ID) %>%
  do(data.frame(size = nrow(.)))
summary(cc_size)

# print(sum(CC_Annot$CC_ID %in% "CC_10")) (just to check on a few ccs that cc_size's result was indeed the number of occurences of each cc in CC_Annot)

#### Functional indices 80 ####
func_homscore = function(x) {
  annot = as.character(x$Annot)
  nprot = nrow(x)
  annot = unique(annot)
  annot = annot[!is.na(annot)]
  if(length(annot) == 0) {
    homscore = NA
  } else if (length(annot) == 1) {
    homscore = 1  #because your ratio of "interesting stuff" is 100%
  } else {
    homscore = 1 - length(annot)/nprot #your ratio of "UNinteresting stuff"     
  }
  return(homscore)
}

homscore_kegg = CC_KEGG %>%
  group_by(CC_ID) %>%
  do(data.frame(Homogeneity.kegg =func_homscore(.)))
homscore_kegg = as.data.frame(homscore_kegg)
summary(homscore_kegg)

homscore_eggnog = CC_Annot %>%
  group_by(CC_ID) %>%
  do(data.frame(Homogeneity.eggnog=func_homscore(.)))
homscore_eggnog = as.data.frame(homscore_eggnog)
summary(homscore_eggnog)

homscore_pfam = CC_PFAM %>%
  group_by(CC_ID) %>%
  do(data.frame(Homogeneity.pfam=func_homscore(.)))
homscore_pfam = as.data.frame(homscore_pfam)
summary(homscore_pfam)

#### Functional indices UNKNOWN ####

func_unkscore = function(x) {
  annot = as.character(x$Annot)
  nprot = nrow(x)
  annot = annot[!is.na(annot)]
  unkscore = (nprot-length(annot))/nprot
  return(unkscore)
}

unkscore_eggnog = CC_Annot %>%
  group_by(CC_ID) %>%
  do(data.frame(Unknowns.eggnog=func_unkscore(.)))
unkscore_eggnog = as.data.frame(unkscore_eggnog)
summary(unkscore_eggnog)

unkscore_kegg = CC_KEGG %>%
  group_by(CC_ID) %>%
  do(data.frame(Unknowns.kegg=func_unkscore(.)))
unkscore_kegg = as.data.frame(unkscore_kegg)
summary(unkscore_kegg)

unkscore_pfam = CC_PFAM %>%
  group_by(CC_ID) %>%
  do(data.frame(Unknowns.pfam=func_unkscore(.)))
unkscore_pfam = as.data.frame(unkscore_pfam)
summary(unkscore_pfam)

#### Taxonomical indices ####
cc_taxo=read.table('CC_80_taxo', header=T, sep=" ")


##### Domain level #####
func_homscore = function(x) {
  annot = as.character(x$Domain)
  nprot = nrow(x)
  annot = unique(annot)
  annot = annot[!is.na(annot)]
  if(length(annot) == 0) {
    homscore = NA
  } else if (length(annot) == 1) {
    homscore = 1
  } else {
    homscore = 1 - length(annot)/nprot    
  }
  return(homscore)
}

homscore_taxo_Domain = cc_taxo %>%
  group_by(CC_ID) %>%
  do(data.frame(Homogeneity.taxo.Domain =func_homscore(.)))
homscore_taxo_Domain = as.data.frame(homscore_taxo_Domain)

##### Phylum level #####
func_homscore = function(x) {
  annot = as.character(x$Phylum)
  nprot = nrow(x)
  annot = unique(annot)
  annot = annot[!is.na(annot)]
  if(length(annot) == 0) {
    homscore = NA
  } else if (length(annot) == 1) {
    homscore = 1
  } else {
    homscore = 1 - length(annot)/nprot    
  }
  return(homscore)
}

homscore_taxo_Phylum = cc_taxo %>%
  group_by(CC_ID) %>%
  do(data.frame(Homogeneity.taxo.Phylum =func_homscore(.)))
homscore_taxo_Phylum = as.data.frame(homscore_taxo_Phylum)

func_unkscore = function(x) {
  annot = as.character(x$Phylum)
  nprot = nrow(x)
  annot = annot[!is.na(annot)]
  unkscore = (nprot-length(annot))/nprot
  return(unkscore)
}

unkscore_taxo_Phylum = cc_taxo %>%
  group_by(CC_ID) %>%
  do(data.frame(Unknowns.taxo.Phylum=func_unkscore(.)))
unkscore_taxo_Phylum = as.data.frame(unkscore_taxo_Phylum)

##### Class level #####
func_homscore = function(x) {
  annot = as.character(x$Class)
  nprot = nrow(x)
  annot = unique(annot)
  annot = annot[!is.na(annot)]
  if(length(annot) == 0) {
    homscore = NA
  } else if (length(annot) == 1) {
    homscore = 1
  } else {
    homscore = 1 - length(annot)/nprot    
  }
  return(homscore)
}

homscore_taxo_Class = cc_taxo %>%
  group_by(CC_ID) %>%
  do(data.frame(Homogeneity.taxo.Class=func_homscore(.)))
homscore_taxo_Class = as.data.frame(homscore_taxo_Class)

func_unkscore = function(x) {
  annot = as.character(x$Class)
  nprot = nrow(x)
  annot = annot[!is.na(annot)]
  unkscore = (nprot-length(annot))/nprot
  return(unkscore)
}

unkscore_taxo_Class = cc_taxo %>%
  group_by(CC_ID) %>%
  do(data.frame(Unknowns.taxo.Class=func_unkscore(.)))
unkscore_taxo_Class = as.data.frame(unkscore_taxo_Class)

##### Order level #####
func_homscore = function(x) {
  annot = as.character(x$Order)
  nprot = nrow(x)
  annot = unique(annot)
  annot = annot[!is.na(annot)]
  if(length(annot) == 0) {
    homscore = NA
  } else if (length(annot) == 1) {
    homscore = 1
  } else {
    homscore = 1 - length(annot)/nprot    
  }
  return(homscore)
}

homscore_taxo_Order = cc_taxo %>%
  group_by(CC_ID) %>%
  do(data.frame(Homogeneity.taxo.Order =func_homscore(.)))
homscore_taxo_Order = as.data.frame(homscore_taxo_Order)

func_unkscore = function(x) {
  annot = as.character(x$Order)
  nprot = nrow(x)
  annot = annot[!is.na(annot)]
  unkscore = (nprot-length(annot))/nprot
  return(unkscore)
}

unkscore_taxo_Order = cc_taxo %>%
  group_by(CC_ID) %>%
  do(data.frame(Unknowns.taxo.Order=func_unkscore(.)))
unkscore_taxo_Order = as.data.frame(unkscore_taxo_Order)

##### Family level #####
func_homscore = function(x) {
  annot = as.character(x$Family)
  nprot = nrow(x)
  annot = unique(annot)
  annot = annot[!is.na(annot)]
  if(length(annot) == 0) {
    homscore = NA
  } else if (length(annot) == 1) {
    homscore = 1
  } else {
    homscore = 1 - length(annot)/nprot    
  }
  return(homscore)
}

homscore_taxo_Family = cc_taxo %>%
  group_by(CC_ID) %>%
  do(data.frame(Homogeneity.taxo.Family =func_homscore(.)))
homscore_taxo_Family = as.data.frame(homscore_taxo_Family)

func_unkscore = function(x) {
  annot = as.character(x$Family)
  nprot = nrow(x)
  annot = annot[!is.na(annot)]
  unkscore = (nprot-length(annot))/nprot
  return(unkscore)
}

unkscore_taxo_Family = cc_taxo %>%
  group_by(CC_ID) %>%
  do(data.frame(Unknowns.taxo.Family=func_unkscore(.)))
unkscore_taxo_Family = as.data.frame(unkscore_taxo_Family)

##### Genus level #####
func_homscore = function(x) {
  annot = as.character(x$Genus)
  nprot = nrow(x)
  annot = unique(annot)
  annot = annot[!is.na(annot)]
  if(length(annot) == 0) {
    homscore = NA
  } else if (length(annot) == 1) {
    homscore = 1
  } else {
    homscore = 1 - length(annot)/nprot    
  }
  return(homscore)
}

homscore_taxo_Genus = cc_taxo %>%
  group_by(CC_ID) %>%
  do(data.frame(Homogeneity.taxo.Genus =func_homscore(.)))
homscore_taxo_Genus = as.data.frame(homscore_taxo_Genus)

func_unkscore = function(x) {
  annot = as.character(x$Genus)
  nprot = nrow(x)
  annot = annot[!is.na(annot)]
  unkscore = (nprot-length(annot))/nprot
  return(unkscore)
}

unkscore_taxo_Genus = cc_taxo %>%
  group_by(CC_ID) %>%
  do(data.frame(Unknowns.taxo.Genus=func_unkscore(.)))
unkscore_taxo_Genus = as.data.frame(unkscore_taxo_Genus)

##### MAG level #####
func_homscore = function(x) {
  annot = as.character(x$Genome_Id)
  nprot = nrow(x)
  annot = unique(annot)
  annot = annot[!is.na(annot)]
  if(length(annot) == 0) {
    homscore = NA
  } else if (length(annot) == 1) {
    homscore = 1
  } else {
    homscore = 1 - length(annot)/nprot    
  }
  return(homscore)
}

homscore_taxo_MAG = cc_taxo %>%
  group_by(CC_ID) %>%
  do(data.frame(Homogeneity.taxo.MAG =func_homscore(.)))
homscore_taxo_MAG = as.data.frame(homscore_taxo_MAG)

func_unkscore = function(x) {
  annot = as.character(x$Genome_Id)
  nprot = nrow(x)
  annot = annot[!is.na(annot)]
  unkscore = (nprot-length(annot))/nprot
  return(unkscore)
}

unkscore_taxo_MAG = cc_taxo %>%
  group_by(CC_ID) %>%
  do(data.frame(Unknowns.taxo.MAG=func_unkscore(.)))
unkscore_taxo_MAG = as.data.frame(unkscore_taxo_MAG)
#####

# Regroup all levels
cc_stats = Reduce(function(x,y) merge(x,y),list(cc_size, homscore_kegg, homscore_eggnog, homscore_pfam, homscore_taxo_Phylum, 
                                                homscore_taxo_Class, homscore_taxo_Order, homscore_taxo_Family,
                                                homscore_taxo_Genus, homscore_taxo_MAG, unkscore_eggnog, unkscore_kegg, unkscore_pfam, homscore_pfam,
                                                unkscore_taxo_Phylum, unkscore_taxo_Class, unkscore_taxo_Order,
                                                unkscore_taxo_Family, unkscore_taxo_Genus, unkscore_taxo_MAG))
write.table(cc_stats, "80CC_stats_with_taxo", quote = FALSE)
