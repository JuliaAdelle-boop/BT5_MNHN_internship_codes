library(dplyr)
library(tidyr)

#Opening the Size Fraction
GGMM <- read.table("GGMM_all_mags_indexed_GeneID_corr_final", header = T, quote = "", sep = "\t", row.names = NULL)
MMQQ <- read.table("MMQQ_all_mags_indexed_GeneID_corr_final", header = T, quote = "", sep = "\t", row.names = NULL)
QQSS <- read.table("QQSS_all_mags_indexed_GeneID_corr_final", header = T, quote = "", sep = "\t", row.names = NULL)
SSUU <- read.table("SSUU_all_mags_indexed_GeneID_corr_final", header = T, quote = "", sep = "\t", row.names = NULL)

#Filter out the rows where Description, KEGG_ko and PFAMs are NA
subset_GGMM <- GGMM %>% filter(is.na(Description) & is.na(KEGG_ko) & is.na(PFAMs))
subset_MMQQ <- MMQQ %>% filter(is.na(Description) & is.na(KEGG_ko) & is.na(PFAMs))
subset_QQSS <- QQSS %>% filter(is.na(Description) & is.na(KEGG_ko) & is.na(PFAMs))
subset_SSUU <- SSUU %>% filter(is.na(Description) & is.na(KEGG_ko) & is.na(PFAMs))

#Create files which contain the unknown-unknown size fractions
write.table(subset_GGMM, "GGMM_all_mags_indexed_GeneID_corr_final_UNKNOWNS", quote = FALSE, sep = "\t", row.names=F)
write.table(subset_MMQQ, "MMQQ_all_mags_indexed_GeneID_corr_final_UNKNOWNS", quote = FALSE, sep = "\t", row.names=F)
write.table(subset_QQSS, "QQSS_all_mags_indexed_GeneID_corr_final_UNKNOWNS", quote = FALSE, sep = "\t", row.names=F)
write.table(subset_SSUU, "SSUU_all_mags_indexed_GeneID_corr_final_UNKNOWNS", quote = FALSE, sep = "\t", row.names=F)
