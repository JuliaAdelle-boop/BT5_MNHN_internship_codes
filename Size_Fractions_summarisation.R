small <- read.table("Funk_R2_GGMM_Fraction_corrected", header = T, sep = "\t", row.names = NULL)
med_small <- read.table("Funk_R2_MMQQ_Fraction_corrected", header = T, sep = "\t", row.names = NULL)
med_big <- read.table("Funk_R2_QQSS_Fraction_corrected", header = T, sep = "\t", row.names = NULL)
big <- read.table("Funk_R2_SSUU_Fraction_corrected", header = T, sep = "\t", row.names = NULL)


SF <- read.table("GGMM_all_mags_indexed_GeneID_corr_final_corrected", header = T, quote = "", sep = "\t", row.names = NULL)
SMF <- read.table("MMQQ_all_mags_indexed_GeneID_corr_final_corrected", header = T, quote = "", sep = "\t", row.names = NULL)
BMF <- read.table("QQSS_all_mags_indexed_GeneID_corr_final_corrected", header = T, quote = "", sep = "\t", row.names = NULL)
BF <- read.table("SSUU_all_mags_indexed_GeneID_corr_final_corrected", header = T, quote = "", sep = "\t", row.names = NULL)

SF_grouped_CCs <- split(SF, SF$CC)
SF_cc_size <- lapply(SF_grouped_CCs, function(x) data.frame(size = nrow(x)))
SF_cc_size <- do.call(rbind, SF_cc_size)
summary(SF_cc_size)

Minimum_SF = min(SF_cc_size$size)
Maximum_SF = max(SF_cc_size$size)
Mean_SF = round(mean(SF_cc_size$size), digits = 2)

SMF_grouped_CCs <- split(SMF, SMF$CC)
SMF_cc_size <- lapply(SMF_grouped_CCs, function(x) data.frame(size = nrow(x)))
SMF_cc_size <- do.call(rbind, SMF_cc_size)

Minimum_SMF = min(SMF_cc_size$size)
Maximum_SMF = max(SMF_cc_size$size)
Mean_SMF = round(mean(SMF_cc_size$size), digits = 2)


BMF_grouped_CCs <- split(BMF, BMF$CC)
BMF_cc_size <- lapply(BMF_grouped_CCs, function(x) data.frame(size = nrow(x)))
BMF_cc_size <- do.call(rbind, BMF_cc_size)

Minimum_BMF = min(BMF_cc_size$size)
Maximum_BMF = max(BMF_cc_size$size)
Mean_BMF = round(mean(BMF_cc_size$size), digits = 2)

BF_grouped_CCs <- split(BF, BF$CC)
BF_cc_size <- lapply(BF_grouped_CCs, function(x) data.frame(size = nrow(x)))
BF_cc_size <- do.call(rbind, BF_cc_size)

Minimum_BF = min(BF_cc_size$size)
Maximum_BF = max(BF_cc_size$size)
Mean_BF = round(mean(BF_cc_size$size), digits = 2)

#### SF's FUNCTIONAL SCORES

#### Mean homscores (for each annotation database)
SF_pre_eggnog_homsc = small$Homogeneity.EggNOG[!is.na(small$Homogeneity.EggNOG)]
SF_nb_NA_homscores_eggnog = sum(is.na(small$Homogeneity.EggNOG))
SF_mean_eggnog_homsc = paste(round(mean(SF_pre_eggnog_homsc), digits = 2), " (", SF_nb_NA_homscores_eggnog, ")", sep = "")

SF_pre_kegg_homsc = small$Homogeneity.KEGG[!is.na(small$Homogeneity.KEGG)]
SF_nb_NA_homscores_kegg = sum(is.na(small$Homogeneity.KEGG))
SF_mean_kegg_homsc = paste(round(mean(SF_pre_kegg_homsc), digits = 2), " (", SF_nb_NA_homscores_kegg, ")", sep="")

SF_pre_pfam_homsc = small$Homogeneity.PFAM[!is.na(small$Homogeneity.PFAM)]
SF_nb_NA_homscores_pfam = sum(is.na(small$Homogeneity.PFAM))
SF_mean_pfam_homsc = paste(round(mean(SF_pre_pfam_homsc), digits = 2), " (", SF_nb_NA_homscores_pfam, ")", sep="")


#### Mean unkscores (for each annotation database)
SF_mean_eggnog_unksc = paste("Eggnog:", round(mean(small$Unknowns.EggNOG), digits = 2))
SF_mean_kegg_unksc = paste("Kegg:", round(mean(small$Unknowns.KEGG), digits = 2))
SF_mean_pfam_unksc = paste("PFAM:", round(mean(small$Unknowns.PFAM), digits = 2))

#### Unknown quantifs
SF_eggnogs_noNAs = nrow(small[small$Unknowns.EggNOG == 0,])
SF_eggnogs_onlyNAs = nrow(small[small$Unknowns.EggNOG == 1,])
SF_eggnogs_atleast_1_anno = nrow(small[small$Unknowns.EggNOG !=1,])

SF_kegg_noNAs = nrow(small[small$Unknowns.KEGG == 0,])
SF_kegg_onlyNAs = nrow(small[small$Unknowns.KEGG == 1,])
SF_kegg_atleast_1_anno = nrow(small[small$Unknowns.KEGG !=1,])

SF_pfam_noNAs = nrow(small[small$Unknowns.PFAM == 0,])
SF_pfam_onlyNAs = nrow(small[small$Unknowns.PFAM == 1,])
SF_pfam_atleast_1_anno = nrow(small[small$Unknowns.PFAM !=1,])

# Number (and percentage) of CCs which don't contain any annotated protein with regards to any database
SF_rowindexes_with_only_completely_unannotated_prots <- which(small$Unknowns.EggNOG == 1 & small$Unknowns.KEGG == 1 & small$Unknowns.PFAM == 1)
SF_number_of_CCs_with_only_completely_unannotated_prots <- length(SF_rowindexes_with_only_completely_unannotated_prots)
SF_percentage_of_CCs_with_only_completely_unannotated_prots_over_total = round(((SF_number_of_CCs_with_only_completely_unannotated_prots / nrow(small))*100), digits = 2)

#### Percentages of unknown quantifs over the total number of CCs
SF_eggnogs_noNAs_percentage = round(((SF_eggnogs_noNAs) / nrow(small) * 100), digits = 2)
SF_eggnogs_onlyNAs_percentage = round(((SF_eggnogs_onlyNAs) / nrow(small) * 100), digits = 2)
SF_eggnogs_atleast_1_anno_percentage = round(((SF_eggnogs_atleast_1_anno) / nrow(small) * 100), digits = 2)

SF_kegg_noNAs_percentage = round(((SF_kegg_noNAs) / nrow(small) * 100), digits = 2)
SF_kegg_onlyNAs_percentage = round(((SF_kegg_onlyNAs) / nrow(small) * 100), digits = 2)
SF_kegg_atleast_1_anno_percentage = round(((SF_kegg_atleast_1_anno) / nrow(small) * 100), digits = 2)

SF_pfam_noNAs_percentage = round(((SF_pfam_noNAs) / nrow(small) * 100), digits = 2)
SF_pfam_onlyNAs_percentage = round(((SF_pfam_onlyNAs) / nrow(small) * 100), digits = 2)
SF_pfam_atleast_1_anno_percentage = round(((SF_pfam_atleast_1_anno) / nrow(small) * 100), digits = 2)

####Combining the unknown quantifs and percentages of unknown quantifs into unique variables, so that they both fit into a single 'Unknown quantifications' column
SF_eggnogs_noNAs_both = paste(SF_eggnogs_noNAs, " (", SF_eggnogs_noNAs_percentage, "%)", sep = "")
SF_eggnogs_onlyNAs_both = paste(SF_eggnogs_onlyNAs, " (", SF_eggnogs_onlyNAs_percentage, "%)", sep = "")
SF_eggnogs_atleast_1_anno_both = paste(SF_eggnogs_atleast_1_anno, " (", SF_eggnogs_atleast_1_anno_percentage, "%)", sep = "")
SF_kegg_noNAs_both = paste(SF_kegg_noNAs, " (", SF_kegg_noNAs_percentage, "%)", sep = "")
SF_kegg_onlyNAs_both = paste(SF_kegg_onlyNAs, " (", SF_kegg_onlyNAs_percentage, "%)", sep = "")
SF_kegg_atleast_1_anno_both = paste(SF_kegg_atleast_1_anno, " (", SF_kegg_atleast_1_anno_percentage, "%)", sep = "")
SF_pfam_noNAs_both = paste(SF_pfam_noNAs, " (", SF_pfam_noNAs_percentage, "%)", sep = "")
SF_pfam_onlyNAs_both = paste(SF_pfam_onlyNAs, " (", SF_pfam_onlyNAs_percentage, "%)", sep = "")
SF_pfam_atleast_1_anno_both = paste(SF_pfam_atleast_1_anno, " (", SF_pfam_atleast_1_anno_percentage, "%)", sep = "")
SF_only_NAs_throughout_all_databases = paste(SF_number_of_CCs_with_only_completely_unannotated_prots, " (", SF_percentage_of_CCs_with_only_completely_unannotated_prots_over_total, "%)", sep = "")

#### Final table
SF_summary <- data.frame(
  Sizes_description = c("Mean", "Min", "Max", rep("",7)),
  Sizes = c(Mean_SF,Minimum_SF,Maximum_SF,rep("", 7)),
  Homogeneity = c("Mean homogeneity with EggNOG annotations", "Mean homogeneity with PFAM annotations","Mean homogeneity with KEGG annotations", rep("",7)),
  Mean_Homscores = c(SF_mean_eggnog_homsc,SF_mean_pfam_homsc,SF_mean_kegg_homsc, rep("", 7)),
  Mean_Unkscores = c(SF_mean_eggnog_unksc,SF_mean_pfam_unksc,SF_mean_kegg_unksc, rep("", 7)),
  Unknown_quantif_descriptions = c("eggnogs_noNAs","eggnogs_atleast_1_anno","eggnogs_onlyNAs", "pfam_noNAs","pfam_atleast_1_anno","pfam_onlyNAs", "kegg_noNAs","kegg_atleast_1_anno","kegg_onlyNAs","Only_NAs_throughout_all_databases"),  
  Unknown_quantifications = c(SF_eggnogs_noNAs_both, SF_eggnogs_atleast_1_anno_both, SF_eggnogs_onlyNAs_both, SF_pfam_noNAs_both, SF_pfam_atleast_1_anno_both, SF_pfam_onlyNAs_both, SF_kegg_noNAs_both, SF_kegg_atleast_1_anno_both, SF_kegg_onlyNAs_both, SF_only_NAs_throughout_all_databases)
)

write.table(SF_summary, "GGMM_summary",quote = FALSE, sep = "\t", row.names = F)

#### SMF's FUNCTIONAL SCORES

#### Mean homscores (for each annotation database)
SMF_pre_eggnog_homsc = med_small$Homogeneity.EggNOG[!is.na(med_small$Homogeneity.EggNOG)]
SMF_nb_NA_homscores_eggnog = sum(is.na(med_small$Homogeneity.EggNOG))
SMF_mean_eggnog_homsc = paste(round(mean(SMF_pre_eggnog_homsc), digits = 2), " (", SMF_nb_NA_homscores_eggnog, ")", sep = "")

SMF_pre_kegg_homsc = med_small$Homogeneity.KEGG[!is.na(med_small$Homogeneity.KEGG)]
SMF_nb_NA_homscores_kegg = sum(is.na(med_small$Homogeneity.KEGG))
SMF_mean_kegg_homsc = paste(round(mean(SMF_pre_kegg_homsc), digits = 2), " (", SMF_nb_NA_homscores_kegg, ")", sep="")

SMF_pre_pfam_homsc = med_small$Homogeneity.PFAM[!is.na(med_small$Homogeneity.PFAM)]
SMF_nb_NA_homscores_pfam = sum(is.na(med_small$Homogeneity.PFAM))
SMF_mean_pfam_homsc = paste(round(mean(SMF_pre_pfam_homsc), digits = 2), " (", SMF_nb_NA_homscores_pfam, ")", sep="")


#### Mean unkscores (for each annotation database)
SMF_mean_eggnog_unksc = paste("Eggnog:", round(mean(med_small$Unknowns.EggNOG), digits = 2))
SMF_mean_kegg_unksc = paste("Kegg:", round(mean(med_small$Unknowns.KEGG), digits = 2))
SMF_mean_pfam_unksc = paste("PFAM:", round(mean(med_small$Unknowns.PFAM), digits = 2))

#### Unknown quantifs
SMF_eggnogs_noNAs = nrow(med_small[med_small$Unknowns.EggNOG == 0,])
SMF_eggnogs_onlyNAs = nrow(med_small[med_small$Unknowns.EggNOG == 1,])
SMF_eggnogs_atleast_1_anno = nrow(med_small[med_small$Unknowns.EggNOG !=1,])

SMF_kegg_noNAs = nrow(med_small[med_small$Unknowns.KEGG == 0,])
SMF_kegg_onlyNAs = nrow(med_small[med_small$Unknowns.KEGG == 1,])
SMF_kegg_atleast_1_anno = nrow(med_small[med_small$Unknowns.KEGG !=1,])

SMF_pfam_noNAs = nrow(med_small[med_small$Unknowns.PFAM == 0,]) 
SMF_pfam_onlyNAs = nrow(med_small[med_small$Unknowns.PFAM == 1,])
SMF_pfam_atleast_1_anno = nrow(med_small[med_small$Unknowns.PFAM !=1,])

# Number (and percentage) of CCs which don't contain any annotated protein with regards to any database
SMF_rowindexes_with_only_completely_unannotated_prots <- which(med_small$Unknowns.EggNOG == 1 & med_small$Unknowns.KEGG == 1 & med_small$Unknowns.PFAM == 1)
SMF_number_of_CCs_with_only_completely_unannotated_prots <- length(SMF_rowindexes_with_only_completely_unannotated_prots)
SMF_percentage_of_CCs_with_only_completely_unannotated_prots_over_total = round(((SMF_number_of_CCs_with_only_completely_unannotated_prots / nrow(med_small))*100), digits = 2)

#### Percentages of unknown quantifs over the total number of CCs
SMF_eggnogs_noNAs_percentage = round(((SMF_eggnogs_noNAs) / nrow(med_small) * 100), digits = 2)
SMF_eggnogs_onlyNAs_percentage = round(((SMF_eggnogs_onlyNAs) / nrow(med_small) * 100), digits = 2)
SMF_eggnogs_atleast_1_anno_percentage = round(((SMF_eggnogs_atleast_1_anno) / nrow(med_small) * 100), digits = 2)

SMF_kegg_noNAs_percentage = round(((SMF_kegg_noNAs) / nrow(med_small) * 100), digits = 2)
SMF_kegg_onlyNAs_percentage = round(((SMF_kegg_onlyNAs) / nrow(med_small) * 100), digits = 2)
SMF_kegg_atleast_1_anno_percentage = round(((SMF_kegg_atleast_1_anno) / nrow(med_small) * 100), digits = 2)

SMF_pfam_noNAs_percentage = round(((SMF_pfam_noNAs) / nrow(med_small) * 100), digits = 2)
SMF_pfam_onlyNAs_percentage = round(((SMF_pfam_onlyNAs) / nrow(med_small) * 100), digits = 2)
SMF_pfam_atleast_1_anno_percentage = round(((SMF_pfam_atleast_1_anno) / nrow(med_small) * 100), digits = 2)

####Combining the unknown quantifs and percentages of unknown quantifs into unique variables, so that they both fit into a single 'Unknown quantifications' column
SMF_eggnogs_noNAs_both = paste(SMF_eggnogs_noNAs, " (", SMF_eggnogs_noNAs_percentage, "%)", sep = "")
SMF_eggnogs_onlyNAs_both = paste(SMF_eggnogs_onlyNAs, " (", SMF_eggnogs_onlyNAs_percentage, "%)", sep = "")
SMF_eggnogs_atleast_1_anno_both = paste(SMF_eggnogs_atleast_1_anno, " (", SMF_eggnogs_atleast_1_anno_percentage, "%)", sep = "")
SMF_kegg_noNAs_both = paste(SMF_kegg_noNAs, " (", SMF_kegg_noNAs_percentage, "%)", sep = "")
SMF_kegg_onlyNAs_both = paste(SMF_kegg_onlyNAs, " (", SMF_kegg_onlyNAs_percentage, "%)", sep = "")
SMF_kegg_atleast_1_anno_both = paste(SMF_kegg_atleast_1_anno, " (", SMF_kegg_atleast_1_anno_percentage, "%)", sep = "")
SMF_pfam_noNAs_both = paste(SMF_pfam_noNAs, " (", SMF_pfam_noNAs_percentage, "%)", sep = "")
SMF_pfam_onlyNAs_both = paste(SMF_pfam_onlyNAs, " (", SMF_pfam_onlyNAs_percentage, "%)", sep = "")
SMF_pfam_atleast_1_anno_both = paste(SMF_pfam_atleast_1_anno, " (", SMF_pfam_atleast_1_anno_percentage, "%)", sep = "")
SMF_only_NAs_throughout_all_databases = paste(SMF_number_of_CCs_with_only_completely_unannotated_prots, " (", SMF_percentage_of_CCs_with_only_completely_unannotated_prots_over_total, "%)", sep = "")

#### Final table
SMF_summary <- data.frame(
  Sizes_description = c("Mean", "Min", "Max", rep("", 7)),
  Sizes = c(Mean_SMF, Minimum_SMF, Maximum_SMF, rep("", 7)),
  Homogeneity = c("Mean homogeneity with EggNOG annotations", "Mean homogeneity with PFAM annotations","Mean homogeneity with KEGG annotations", rep("",7)),
  Mean_Homscores = c(SMF_mean_eggnog_homsc, SMF_mean_pfam_homsc, SMF_mean_kegg_homsc, rep("", 7)),
  Mean_Unkscores = c(SMF_mean_eggnog_unksc, SMF_mean_pfam_unksc, SMF_mean_kegg_unksc, rep("", 7)),
  Unknown_quantif_descriptions = c("eggnogs_noNAs", "eggnogs_atleast_1_anno", "eggnogs_onlyNAs", "pfam_noNAs", "pfam_atleast_1_anno", "pfam_onlyNAs", "kegg_noNAs","kegg_atleast_1_anno","kegg_onlyNAs","Only_NAs_throughout_all_databases"),  
  Unknown_quantifications = c(SMF_eggnogs_noNAs_both, SMF_eggnogs_atleast_1_anno_both, SMF_eggnogs_onlyNAs_both, SMF_pfam_noNAs_both, SMF_pfam_atleast_1_anno_both, SMF_pfam_onlyNAs_both, SMF_kegg_noNAs_both, SMF_kegg_atleast_1_anno_both, SMF_kegg_onlyNAs_both, SMF_only_NAs_throughout_all_databases)
)

write.table(SMF_summary, "MMQQ_summary",quote = FALSE, sep = "\t", row.names = F)

#### BMF's FUNCTIONAL SCORES

#### Mean homscores (for each annotation database)
BMF_pre_eggnog_homsc = med_big$Homogeneity.EggNOG[!is.na(med_big$Homogeneity.EggNOG)]
BMF_nb_NA_homscores_eggnog = sum(is.na(med_big$Homogeneity.EggNOG))
BMF_mean_eggnog_homsc = paste(round(mean(BMF_pre_eggnog_homsc), digits = 2), " (", BMF_nb_NA_homscores_eggnog, ")", sep = "")

BMF_pre_kegg_homsc = med_big$Homogeneity.KEGG[!is.na(med_big$Homogeneity.KEGG)]
BMF_nb_NA_homscores_kegg = sum(is.na(med_big$Homogeneity.KEGG))
BMF_mean_kegg_homsc = paste(round(mean(BMF_pre_kegg_homsc), digits = 2), " (", BMF_nb_NA_homscores_kegg, ")", sep="")

BMF_pre_pfam_homsc = med_big$Homogeneity.PFAM[!is.na(med_big$Homogeneity.PFAM)]
BMF_nb_NA_homscores_pfam = sum(is.na(med_big$Homogeneity.PFAM))
BMF_mean_pfam_homsc = paste(round(mean(BMF_pre_pfam_homsc), digits = 2), " (", BMF_nb_NA_homscores_pfam, ")", sep="")


#### Mean unkscores (for each annotation database)
BMF_mean_eggnog_unksc = paste("Eggnog:", round(mean(med_big$Unknowns.EggNOG), digits = 2))
BMF_mean_kegg_unksc = paste("Kegg:", round(mean(med_big$Unknowns.KEGG), digits = 2))
BMF_mean_pfam_unksc = paste("PFAM:", round(mean(med_big$Unknowns.PFAM), digits = 2))

#### Unknown quantifs
BMF_eggnogs_noNAs = nrow(med_big[med_big$Unknowns.EggNOG == 0,])
BMF_eggnogs_onlyNAs = nrow(med_big[med_big$Unknowns.EggNOG == 1,])
BMF_eggnogs_atleast_1_anno = nrow(med_big[med_big$Unknowns.EggNOG !=1,])

BMF_kegg_noNAs = nrow(med_big[med_big$Unknowns.KEGG == 0,])
BMF_kegg_onlyNAs = nrow(med_big[med_big$Unknowns.KEGG == 1,])
BMF_kegg_atleast_1_anno = nrow(med_big[med_big$Unknowns.KEGG !=1,])

BMF_pfam_noNAs = nrow(med_big[med_big$Unknowns.PFAM == 0,])
BMF_pfam_onlyNAs = nrow(med_big[med_big$Unknowns.PFAM == 1,])
BMF_pfam_atleast_1_anno = nrow(med_big[med_big$Unknowns.PFAM !=1,])

# Number (and percentage) of CCs which don't contain any annotated protein with regards to any database
BMF_rowindexes_with_only_completely_unannotated_prots <- which(med_big$Unknowns.EggNOG == 1 & med_big$Unknowns.KEGG == 1 & med_big$Unknowns.PFAM == 1)
BMF_number_of_CCs_with_only_completely_unannotated_prots <- length(BMF_rowindexes_with_only_completely_unannotated_prots)
BMF_percentage_of_CCs_with_only_completely_unannotated_prots_over_total = round(((BMF_number_of_CCs_with_only_completely_unannotated_prots / nrow(med_big))*100), digits = 2)

#### Percentages of unknown quantifs over the total number of CCs
BMF_eggnogs_noNAs_percentage = round(((BMF_eggnogs_noNAs) / nrow(med_big) * 100), digits = 2)
BMF_eggnogs_onlyNAs_percentage = round(((BMF_eggnogs_onlyNAs) / nrow(med_big) * 100), digits = 2)
BMF_eggnogs_atleast_1_anno_percentage = round(((BMF_eggnogs_atleast_1_anno) / nrow(med_big) * 100), digits = 2)

BMF_kegg_noNAs_percentage = round(((BMF_kegg_noNAs) / nrow(med_big) * 100), digits = 2)
BMF_kegg_onlyNAs_percentage = round(((BMF_kegg_onlyNAs) / nrow(med_big) * 100), digits = 2)
BMF_kegg_atleast_1_anno_percentage = round(((BMF_kegg_atleast_1_anno) / nrow(med_big) * 100), digits = 2)

BMF_pfam_noNAs_percentage = round(((BMF_pfam_noNAs) / nrow(med_big) * 100), digits = 2)
BMF_pfam_onlyNAs_percentage = round(((BMF_pfam_onlyNAs) / nrow(med_big) * 100), digits = 2)
BMF_pfam_atleast_1_anno_percentage = round(((BMF_pfam_atleast_1_anno) / nrow(med_big) * 100), digits = 2)

####Combining the unknown quantifs and percentages of unknown quantifs into unique variables, so that they both fit into a single 'Unknown quantifications' column
BMF_eggnogs_noNAs_both = paste(BMF_eggnogs_noNAs, " (", BMF_eggnogs_noNAs_percentage, "%)", sep = "")
BMF_eggnogs_onlyNAs_both = paste(BMF_eggnogs_onlyNAs, " (", BMF_eggnogs_onlyNAs_percentage, "%)", sep = "")
BMF_eggnogs_atleast_1_anno_both = paste(BMF_eggnogs_atleast_1_anno, " (", BMF_eggnogs_atleast_1_anno_percentage, "%)", sep = "")
BMF_kegg_noNAs_both = paste(BMF_kegg_noNAs, " (", BMF_kegg_noNAs_percentage, "%)", sep = "")
BMF_kegg_onlyNAs_both = paste(BMF_kegg_onlyNAs, " (", BMF_kegg_onlyNAs_percentage, "%)", sep = "")
BMF_kegg_atleast_1_anno_both = paste(BMF_kegg_atleast_1_anno, " (", BMF_kegg_atleast_1_anno_percentage, "%)", sep = "")
BMF_pfam_noNAs_both = paste(BMF_pfam_noNAs, " (", BMF_pfam_noNAs_percentage, "%)", sep = "")
BMF_pfam_onlyNAs_both = paste(BMF_pfam_onlyNAs, " (", BMF_pfam_onlyNAs_percentage, "%)", sep = "")
BMF_pfam_atleast_1_anno_both = paste(BMF_pfam_atleast_1_anno, " (", BMF_pfam_atleast_1_anno_percentage, "%)", sep = "")
BMF_only_NAs_throughout_all_databases = paste(BMF_number_of_CCs_with_only_completely_unannotated_prots, " (", BMF_percentage_of_CCs_with_only_completely_unannotated_prots_over_total, "%)", sep = "")

#### Final table
BMF_summary <- data.frame(
  Sizes_description = c("Mean", "Min", "Max", rep("",7)),
  Sizes = c(Mean_BMF,Minimum_BMF,Maximum_BMF,rep("", 7)),
  Homogeneity = c("Mean homogeneity with EggNOG annotations", "Mean homogeneity with PFAM annotations","Mean homogeneity with KEGG annotations", rep("",7)),
  Mean_Homscores = c(BMF_mean_eggnog_homsc,BMF_mean_pfam_homsc,BMF_mean_kegg_homsc, rep("", 7)),
  Mean_Unkscores = c(BMF_mean_eggnog_unksc,BMF_mean_pfam_unksc,BMF_mean_kegg_unksc, rep("", 7)),
  Unknown_quantif_descriptions = c("eggnogs_noNAs","eggnogs_atleast_1_anno","eggnogs_onlyNAs", "pfam_noNAs","pfam_atleast_1_anno","pfam_onlyNAs", "kegg_noNAs","kegg_atleast_1_anno","kegg_onlyNAs","Only_NAs_throughout_all_databases"),  
  Unknown_quantifications = c(BMF_eggnogs_noNAs_both, BMF_eggnogs_atleast_1_anno_both, BMF_eggnogs_onlyNAs_both, BMF_pfam_noNAs_both, BMF_pfam_atleast_1_anno_both, BMF_pfam_onlyNAs_both, BMF_kegg_noNAs_both, BMF_kegg_atleast_1_anno_both, BMF_kegg_onlyNAs_both, BMF_only_NAs_throughout_all_databases)
)

write.table(BMF_summary, "QQSS_summary",quote = FALSE, sep = "\t", row.names = F)

#### BF's FUNCTIONAL SCORES

#### Mean homscores (for each annotation database)
BF_pre_eggnog_homsc = big$Homogeneity.EggNOG[!is.na(big$Homogeneity.EggNOG)]
BF_nb_NA_homscores_eggnog = sum(is.na(big$Homogeneity.EggNOG))
BF_mean_eggnog_homsc = paste(round(mean(BF_pre_eggnog_homsc), digits = 2), " (", BF_nb_NA_homscores_eggnog, ")", sep = "")

BF_pre_kegg_homsc = big$Homogeneity.KEGG[!is.na(big$Homogeneity.KEGG)]
BF_nb_NA_homscores_kegg = sum(is.na(big$Homogeneity.KEGG))
BF_mean_kegg_homsc = paste(round(mean(BF_pre_kegg_homsc), digits = 2), " (", BF_nb_NA_homscores_kegg, ")", sep="")

BF_pre_pfam_homsc = big$Homogeneity.PFAM[!is.na(big$Homogeneity.PFAM)]
BF_nb_NA_homscores_pfam = sum(is.na(big$Homogeneity.PFAM))
BF_mean_pfam_homsc = paste(round(mean(BF_pre_pfam_homsc), digits = 2), " (", BF_nb_NA_homscores_pfam, ")", sep="")


#### Mean unkscores (for each annotation database)
BF_mean_eggnog_unksc = paste("Eggnog:", round(mean(big$Unknowns.EggNOG), digits = 2))
BF_mean_kegg_unksc = paste("Kegg:", round(mean(big$Unknowns.KEGG), digits = 2))
BF_mean_pfam_unksc = paste("PFAM:", round(mean(big$Unknowns.PFAM), digits = 2))

#### Unknown quantifs
BF_eggnogs_noNAs = nrow(big[big$Unknowns.EggNOG == 0,])
BF_eggnogs_onlyNAs = nrow(big[big$Unknowns.EggNOG == 1,])
BF_eggnogs_atleast_1_anno = nrow(big[big$Unknowns.EggNOG !=1,])

BF_kegg_noNAs = nrow(big[big$Unknowns.KEGG == 0,])
BF_kegg_onlyNAs = nrow(big[big$Unknowns.KEGG == 1,])
BF_kegg_atleast_1_anno = nrow(big[big$Unknowns.KEGG !=1,])

BF_pfam_noNAs = nrow(big[big$Unknowns.PFAM == 0,])
BF_pfam_onlyNAs = nrow(big[big$Unknowns.PFAM == 1,])
BF_pfam_atleast_1_anno = nrow(big[big$Unknowns.PFAM !=1,])

# Number (and percentage) of CCs which don't contain any annotated protein with regards to any database
BF_rowindexes_with_only_completely_unannotated_prots <- which(big$Unknowns.EggNOG == 1 & big$Unknowns.KEGG == 1 & big$Unknowns.PFAM == 1)
BF_number_of_CCs_with_only_completely_unannotated_prots <- length(BF_rowindexes_with_only_completely_unannotated_prots)
BF_percentage_of_CCs_with_only_completely_unannotated_prots_over_total = round(((BF_number_of_CCs_with_only_completely_unannotated_prots / nrow(big))*100), digits = 2)

#### Percentages of unknown quantifs over the total number of CCs
BF_eggnogs_noNAs_percentage = round(((BF_eggnogs_noNAs) / nrow(big) * 100), digits = 2)
BF_eggnogs_onlyNAs_percentage = round(((BF_eggnogs_onlyNAs) / nrow(big) * 100), digits = 2)
BF_eggnogs_atleast_1_anno_percentage = round(((BF_eggnogs_atleast_1_anno) / nrow(big) * 100), digits = 2)

BF_kegg_noNAs_percentage = round(((BF_kegg_noNAs) / nrow(big) * 100), digits = 2)
BF_kegg_onlyNAs_percentage = round(((BF_kegg_onlyNAs) / nrow(big) * 100), digits = 2)
BF_kegg_atleast_1_anno_percentage = round(((BF_kegg_atleast_1_anno) / nrow(big) * 100), digits = 2)

BF_pfam_noNAs_percentage = round(((BF_pfam_noNAs) / nrow(big) * 100), digits = 2)
BF_pfam_onlyNAs_percentage = round(((BF_pfam_onlyNAs) / nrow(big) * 100), digits = 2)
BF_pfam_atleast_1_anno_percentage = round(((BF_pfam_atleast_1_anno) / nrow(big) * 100), digits = 2)

####Combining the unknown quantifs and percentages of unknown quantifs into unique variables, so that they both fit into a single 'Unknown quantifications' column
BF_eggnogs_noNAs_both = paste(BF_eggnogs_noNAs, " (", BF_eggnogs_noNAs_percentage, "%)", sep = "")
BF_eggnogs_onlyNAs_both = paste(BF_eggnogs_onlyNAs, " (", BF_eggnogs_onlyNAs_percentage, "%)", sep = "")
BF_eggnogs_atleast_1_anno_both = paste(BF_eggnogs_atleast_1_anno, " (", BF_eggnogs_atleast_1_anno_percentage, "%)", sep = "")
BF_kegg_noNAs_both = paste(BF_kegg_noNAs, " (", BF_kegg_noNAs_percentage, "%)", sep = "")
BF_kegg_onlyNAs_both = paste(BF_kegg_onlyNAs, " (", BF_kegg_onlyNAs_percentage, "%)", sep = "")
BF_kegg_atleast_1_anno_both = paste(BF_kegg_atleast_1_anno, " (", BF_kegg_atleast_1_anno_percentage, "%)", sep = "")
BF_pfam_noNAs_both = paste(BF_pfam_noNAs, " (", BF_pfam_noNAs_percentage, "%)", sep = "")
BF_pfam_onlyNAs_both = paste(BF_pfam_onlyNAs, " (", BF_pfam_onlyNAs_percentage, "%)", sep = "")
BF_pfam_atleast_1_anno_both = paste(BF_pfam_atleast_1_anno, " (", BF_pfam_atleast_1_anno_percentage, "%)", sep = "")
BF_only_NAs_throughout_all_databases = paste(BF_number_of_CCs_with_only_completely_unannotated_prots, " (", BF_percentage_of_CCs_with_only_completely_unannotated_prots_over_total, "%)", sep = "")

#### Final table
BF_summary <- data.frame(
  Sizes_description = c("Mean", "Min", "Max", rep("",7)),
  Sizes = c(Mean_BF,Minimum_BF,Maximum_BF,rep("", 7)),
  Homogeneity = c("Mean homogeneity with EggNOG annotations", "Mean homogeneity with PFAM annotations","Mean homogeneity with KEGG annotations", rep("",7)),
  Mean_Homscores = c(BF_mean_eggnog_homsc,BF_mean_pfam_homsc,BF_mean_kegg_homsc, rep("", 7)),
  Mean_Unkscores = c(BF_mean_eggnog_unksc,BF_mean_pfam_unksc,BF_mean_kegg_unksc, rep("", 7)),
  Unknown_quantif_descriptions = c("eggnogs_noNAs","eggnogs_atleast_1_anno","eggnogs_onlyNAs", "pfam_noNAs","pfam_atleast_1_anno","pfam_onlyNAs", "kegg_noNAs","kegg_atleast_1_anno","kegg_onlyNAs","Only_NAs_throughout_all_databases"),  
  Unknown_quantifications = c(BF_eggnogs_noNAs_both, BF_eggnogs_atleast_1_anno_both, BF_eggnogs_onlyNAs_both, BF_pfam_noNAs_both, BF_pfam_atleast_1_anno_both, BF_pfam_onlyNAs_both, BF_kegg_noNAs_both, BF_kegg_atleast_1_anno_both, BF_kegg_onlyNAs_both, BF_only_NAs_throughout_all_databases)
)

write.table(BF_summary, "SSUU_summary", quote = FALSE, sep = "\t", row.names = F)