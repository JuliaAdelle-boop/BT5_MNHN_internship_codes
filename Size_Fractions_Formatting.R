CC_ALL <- read.csv('quoted_sanstiretsCC_80_withallannot.csv', header = T, row.names = NULL, na.strings=c("NA"))

small_file <- read.table("GGMM_all_mags_indexed_GeneID_corr", header = T, sep = " ", row.names = NULL)
small_med_file <- read.table("MMQQ_all_mags_indexed_GeneID_corr", header = T, sep = " ", row.names = NULL)
big_med_file <- read.table("QQSS_all_mags_indexed_GeneID_corr", header = T, sep = " ", row.names = NULL)
big_file <- read.table("SSUU_all_mags_indexed_GeneID_corr", header = T, sep = " ", row.names = NULL)

full_small <- merge(x=small_file, y=CC_ALL, by="Genes", all.x=T)
full_med_small <- merge(x=small_med_file, y=CC_ALL, by="Genes", all.x=T)
full_med_big <- merge(x=big_med_file, y=CC_ALL, by="Genes", all.x=T)
full_big <- merge(x=big_file, y=CC_ALL, by="Genes", all.x=T)

#removal of the GO column from the 'full' tables
full_small <- subset(full_small, select = - GOs)
full_med_small <- subset(full_med_small, select = - GOs)
full_med_big <- subset(full_med_big, select = - GOs)
full_big <- subset(full_big, select = - GOs)

#borttagning av alla dem raderna som inte är sammankopplade med CCer
full_small <- full_small[!is.na(full_small$CC),]
full_med_small <- full_med_small[!is.na(full_med_small$CC),]
full_med_big <- full_med_big[!is.na(full_med_big$CC),]
full_big <- full_big[!is.na(full_big$CC),]

write.table(full_small, "GGMM_all_mags_indexed_GeneID_corr_final_corrected", quote = FALSE, sep = "\t", row.names = F)
write.table(full_med_small, "MMQQ_all_mags_indexed_GeneID_corr_final_corrected", quote = FALSE, sep = "\t", row.names = F)
write.table(full_med_big, "QQSS_all_mags_indexed_GeneID_corr_final_corrected", quote = FALSE, sep = "\t", row.names = F)
write.table(full_big, "SSUU_all_mags_indexed_GeneID_corr_final_corrected", quote = FALSE, sep = "\t", row.names = F)
