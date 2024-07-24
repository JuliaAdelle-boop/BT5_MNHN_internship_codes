small <- read.table("GGMM_all_mags_indexed_GeneID_corr_final_corrected", header = T, quote = "", sep = "\t", row.names = NULL)
med_small <- read.table("MMQQ_all_mags_indexed_GeneID_corr_final_corrected", header = T, quote = "", sep = "\t", row.names = NULL)
med_big <- read.table("QQSS_all_mags_indexed_GeneID_corr_final_corrected", header = T, quote = "", sep = "\t", row.names = NULL)
big <- read.table("SSUU_all_mags_indexed_GeneID_corr_final_corrected", header = T, quote = "", sep = "\t", row.names = NULL)



#SF_grouped_CCs <- split(small, small$CC)
#SF_cc_size <- lapply(SF_grouped_CCs, function(x) data.frame(size = nrow(x)))
#SF_cc_size <- do.call(rbind, SF_cc_size)
#summary(SF_cc_size)

#Minimum_SF = min(SF_cc_size)
#Maximum_SF = max(SF_cc_size)
#Mean_SF = round(mean(SF_cc_size), digits = 2)

#SMF_grouped_CCs <- split(med_small, med_small$CC)
#SMF_cc_size <- lapply(SMF_grouped_CCs, function(x) data.frame(size = nrow(x))) 
#SMF_cc_size <- do.call(rbind, SMF_cc_size)

#Minimum_SMF = min(SMF_cc_size)
#Maximum_SMF = max(SMF_cc_size)
#Mean_SMF = round(mean(SMF_cc_size), digits = 2)


#BMF_grouped_CCs <- split(med_big, med_big$CC) 
#BMF_cc_size <- lapply(BMF_grouped_CCs, function(x) data.frame(size = nrow(x)))  
#BMF_cc_size <- do.call(rbind, BMF_cc_size)

#Minimum_BMF = min(BMF_cc_size)
#Maximum_BMF = max(BMF_cc_size)
#Mean_BMF = round(mean(BMF_cc_size), digits = 2)

#BF_grouped_CCs <- split(big, big$CC)  # Split the data frame into groups based on 'CC'
#BF_cc_size <- lapply(BF_grouped_CCs, function(x) data.frame(size = nrow(x)))
#BF_cc_size <- do.call(rbind, BF_cc_size)

#Minimum_BF = min(BF_cc_size)
#Maximum_BF = max(BF_cc_size)
#Mean_BF = round(mean(BF_cc_size), digits = 2)

#homscore functions, for each database
func_homscore_EggNOG = function(x) {
  annot = as.character(x$Description)
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

func_homscore_KEGG = function(x) {
  annot = as.character(x$KEGG_ko)
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

func_homscore_PFAMs = function(x) {
  annot = as.character(x$PFAMs)
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

#condensation of the size fractions by CC ID
grouped_CCs_small <- split(small, small$CC)
grouped_CCs_med_small <- split(med_small, med_small$CC)
grouped_CCs_med_big <- split(med_big, med_big$CC)
grouped_CCs_big <- split(big, big$CC)

#EggNOG homscores, for each size fraction
homscore_EggNOG_small <- lapply(grouped_CCs_small, func_homscore_EggNOG)
homscore_EggNOG_small <- data.frame(CC_ID = names(homscore_EggNOG_small), Homogeneity.EggNOG = unlist(homscore_EggNOG_small))
rownames(homscore_EggNOG_small) <- 1:nrow(homscore_EggNOG_small)

homscore_EggNOG_med_small <- lapply(grouped_CCs_med_small, func_homscore_EggNOG)
homscore_EggNOG_med_small <- data.frame(CC_ID = names(homscore_EggNOG_med_small), Homogeneity.EggNOG = unlist(homscore_EggNOG_med_small))
rownames(homscore_EggNOG_med_small) <- 1:nrow(homscore_EggNOG_med_small)

homscore_EggNOG_med_big <- lapply(grouped_CCs_med_big, func_homscore_EggNOG)
homscore_EggNOG_med_big <- data.frame(CC_ID = names(homscore_EggNOG_med_big), Homogeneity.EggNOG = unlist(homscore_EggNOG_med_big))
rownames(homscore_EggNOG_med_big) <- 1:nrow(homscore_EggNOG_med_big)

homscore_EggNOG_big <- lapply(grouped_CCs_big, func_homscore_EggNOG)
homscore_EggNOG_big <- data.frame(CC_ID = names(homscore_EggNOG_big), Homogeneity.EggNOG = unlist(homscore_EggNOG_big))
rownames(homscore_EggNOG_big) <- 1:nrow(homscore_EggNOG_big)

#KEGG homscores, for each size fraction
homscore_KEGG_small <- lapply(grouped_CCs_small, func_homscore_KEGG)
homscore_KEGG_small <- data.frame(CC_ID = names(homscore_KEGG_small), Homogeneity.KEGG = unlist(homscore_KEGG_small))
rownames(homscore_KEGG_small) <- 1:nrow(homscore_KEGG_small)

homscore_KEGG_med_small <- lapply(grouped_CCs_med_small, func_homscore_KEGG)
homscore_KEGG_med_small <- data.frame(CC_ID = names(homscore_KEGG_med_small), Homogeneity.KEGG = unlist(homscore_KEGG_med_small))
rownames(homscore_KEGG_med_small) <- 1:nrow(homscore_KEGG_med_small)

homscore_KEGG_med_big <- lapply(grouped_CCs_med_big, func_homscore_KEGG)
homscore_KEGG_med_big <- data.frame(CC_ID = names(homscore_KEGG_med_big), Homogeneity.KEGG = unlist(homscore_KEGG_med_big))
rownames(homscore_KEGG_med_big) <- 1:nrow(homscore_KEGG_med_big)

homscore_KEGG_big <- lapply(grouped_CCs_big, func_homscore_KEGG)
homscore_KEGG_big <- data.frame(CC_ID = names(homscore_KEGG_big), Homogeneity.KEGG = unlist(homscore_KEGG_big))
rownames(homscore_KEGG_big) <- 1:nrow(homscore_KEGG_big)

#PFAM homscores, for each size fraction
homscore_PFAM_small <- lapply(grouped_CCs_small, func_homscore_PFAMs)
homscore_PFAM_small <- data.frame(CC_ID = names(homscore_PFAM_small), Homogeneity.PFAM = unlist(homscore_PFAM_small))
rownames(homscore_PFAM_small) <- 1:nrow(homscore_PFAM_small)

homscore_PFAM_med_small <- lapply(grouped_CCs_med_small, func_homscore_PFAMs)
homscore_PFAM_med_small <- data.frame(CC_ID = names(homscore_PFAM_med_small), Homogeneity.PFAM = unlist(homscore_PFAM_med_small))
rownames(homscore_PFAM_med_small) <- 1:nrow(homscore_PFAM_med_small)

homscore_PFAM_med_big <- lapply(grouped_CCs_med_big, func_homscore_PFAMs)
homscore_PFAM_med_big <- data.frame(CC_ID = names(homscore_PFAM_med_big), Homogeneity.PFAM = unlist(homscore_PFAM_med_big))
rownames(homscore_PFAM_med_big) <- 1:nrow(homscore_PFAM_med_big)

homscore_PFAM_big <- lapply(grouped_CCs_big, func_homscore_PFAMs)
homscore_PFAM_big <- data.frame(CC_ID = names(homscore_PFAM_big), Homogeneity.PFAM = unlist(homscore_PFAM_big))
rownames(homscore_PFAM_big) <- 1:nrow(homscore_PFAM_big)

#Unknown score functions, for each database
func_unkscore_EggNOG = function(x) {
  annot = as.character(x$Description)
  nprot = nrow(x)
  annot = annot[!is.na(annot)]
  unkscore = (nprot-length(annot))/nprot
  return(unkscore)
}

func_unkscore_KEGG = function(x) {
  annot = as.character(x$KEGG_ko)
  nprot = nrow(x)
  annot = annot[!is.na(annot)]
  unkscore = (nprot-length(annot))/nprot
  return(unkscore)
}

func_unkscore_PFAM = function(x) {
  annot = as.character(x$PFAMs)
  nprot = nrow(x)
  annot = annot[!is.na(annot)]
  unkscore = (nprot-length(annot))/nprot
  return(unkscore)
}

#EggNOG Unkscores, for each size fraction
unkscore_EggNOG_small <- lapply(grouped_CCs_small, func_unkscore_EggNOG)
unkscore_EggNOG_small <- data.frame(CC_ID = names(unkscore_EggNOG_small), Unknowns.EggNOG = unlist(unkscore_EggNOG_small))
rownames(unkscore_EggNOG_small) <- 1:nrow(unkscore_EggNOG_small)

unkscore_EggNOG_med_small <- lapply(grouped_CCs_med_small, func_unkscore_EggNOG)
unkscore_EggNOG_med_small <- data.frame(CC_ID = names(unkscore_EggNOG_med_small), Unknowns.EggNOG = unlist(unkscore_EggNOG_med_small))
rownames(unkscore_EggNOG_med_small) <- 1:nrow(unkscore_EggNOG_med_small)

unkscore_EggNOG_med_big <- lapply(grouped_CCs_med_big, func_unkscore_EggNOG)
unkscore_EggNOG_med_big <- data.frame(CC_ID = names(unkscore_EggNOG_med_big), Unknowns.EggNOG = unlist(unkscore_EggNOG_med_big))
rownames(unkscore_EggNOG_med_big) <- 1:nrow(unkscore_EggNOG_med_big)

unkscore_EggNOG_big <- lapply(grouped_CCs_big, func_unkscore_EggNOG)
unkscore_EggNOG_big <- data.frame(CC_ID = names(unkscore_EggNOG_big), Unknowns.EggNOG = unlist(unkscore_EggNOG_big))
rownames(unkscore_EggNOG_big) <- 1:nrow(unkscore_EggNOG_big)

#KEGG Unkscores, for each size fraction
unkscore_KEGG_small <- lapply(grouped_CCs_small, func_unkscore_KEGG)
unkscore_KEGG_small <- data.frame(CC_ID = names(unkscore_KEGG_small), Unknowns.KEGG = unlist(unkscore_KEGG_small))
rownames(unkscore_KEGG_small) <- 1:nrow(unkscore_KEGG_small)

unkscore_KEGG_med_small <- lapply(grouped_CCs_med_small, func_unkscore_KEGG)
unkscore_KEGG_med_small <- data.frame(CC_ID = names(unkscore_KEGG_med_small), Unknowns.KEGG = unlist(unkscore_KEGG_med_small))
rownames(unkscore_KEGG_med_small) <- 1:nrow(unkscore_KEGG_med_small)

unkscore_KEGG_med_big <- lapply(grouped_CCs_med_big, func_unkscore_KEGG)
unkscore_KEGG_med_big <- data.frame(CC_ID = names(unkscore_KEGG_med_big), Unknowns.KEGG = unlist(unkscore_KEGG_med_big))
rownames(unkscore_KEGG_med_big) <- 1:nrow(unkscore_KEGG_med_big)

unkscore_KEGG_big <- lapply(grouped_CCs_big, func_unkscore_KEGG)
unkscore_KEGG_big <- data.frame(CC_ID = names(unkscore_KEGG_big), Unknowns.KEGG = unlist(unkscore_KEGG_big))
rownames(unkscore_KEGG_big) <- 1:nrow(unkscore_KEGG_big)

#PFAM Unkscores, for each size fraction
unkscore_PFAM_small <- lapply(grouped_CCs_small, func_unkscore_PFAM)
unkscore_PFAM_small <- data.frame(CC_ID = names(unkscore_PFAM_small), Unknowns.PFAM = unlist(unkscore_PFAM_small))
rownames(unkscore_PFAM_small) <- 1:nrow(unkscore_PFAM_small)

unkscore_PFAM_med_small <- lapply(grouped_CCs_med_small, func_unkscore_PFAM)
unkscore_PFAM_med_small <- data.frame(CC_ID = names(unkscore_PFAM_med_small), Unknowns.PFAM = unlist(unkscore_PFAM_med_small))
rownames(unkscore_PFAM_med_small) <- 1:nrow(unkscore_PFAM_med_small)

unkscore_PFAM_med_big <- lapply(grouped_CCs_big, func_unkscore_PFAM)
unkscore_PFAM_med_big <- data.frame(CC_ID = names(unkscore_PFAM_med_big), Unknowns.PFAM = unlist(unkscore_PFAM_med_big))
rownames(unkscore_PFAM_med_big) <- 1:nrow(unkscore_PFAM_med_big)

unkscore_PFAM_big <- lapply(grouped_CCs_big, func_unkscore_PFAM)
unkscore_PFAM_big <- data.frame(CC_ID = names(unkscore_PFAM_big), Unknowns.PFAM = unlist(unkscore_PFAM_big))
rownames(unkscore_PFAM_big) <- 1:nrow(unkscore_PFAM_big)

#put all of the homscores and unkscores into four tables, one for each size fraction, like Funk_R2, which will allow me to calculate min, max and mean homscore and unkscore
GGMM_fractions_scores = Reduce(function(x,y) merge(x,y),list(homscore_EggNOG_small, 
                                                homscore_KEGG_small,
                                                homscore_PFAM_small,
						unkscore_EggNOG_small,
                                                unkscore_KEGG_small,
                                                unkscore_PFAM_small))
write.table(GGMM_fractions_scores, "Funk_R2_GGMM_Fraction_corrected", quote = FALSE, sep = "\t", row.names = F)

MMQQ_fractions_scores = Reduce(function(x,y) merge(x,y),list(homscore_EggNOG_med_small,
						homscore_KEGG_med_small,
                                                homscore_PFAM_med_small,
                                                unkscore_EggNOG_med_small,
                                                unkscore_KEGG_med_small,
                                                unkscore_PFAM_med_small))
write.table(MMQQ_fractions_scores, "Funk_R2_MMQQ_Fraction_corrected", quote = FALSE, sep = "\t", row.names = F)

QQSS_fractions_scores = Reduce(function(x,y) merge(x,y),list(homscore_EggNOG_med_big, 
						homscore_KEGG_big,
                                                homscore_PFAM_med_big,
                                                unkscore_EggNOG_med_big,
                                                unkscore_KEGG_med_big,
                                                unkscore_PFAM_med_big))
write.table(QQSS_fractions_scores, "Funk_R2_QQSS_Fraction_corrected", quote = FALSE, sep = "\t", row.names = F)

SSUU_fractions_scores = Reduce(function(x,y) merge(x,y),list(homscore_EggNOG_big,
						homscore_KEGG_big,
                                                homscore_PFAM_big,
                                                unkscore_EggNOG_big,
                                                unkscore_KEGG_big,
                                                unkscore_PFAM_big))
write.table(SSUU_fractions_scores, "Funk_R2_SSUU_Fraction_corrected", quote = FALSE, sep = "\t", row.names = F)
