all_mags_indexed <- read.table("all_mags_indexed_GeneID_corr", sep = " ")
header_row <- c("Genes","Gene_ID","MAG_ID","Water","MAG_Cov_X","MAG_Cov_Y","Gene_Expression","Gene_Expression_Index")
colnames(all_mags_indexed) <- header_row

small_orgs_i_tabellen <- grep("GGMM$", all_mags_indexed$Water)
small_orgs <- all_mags_indexed[small_orgs_i_tabellen,]
#small_orgs <- filter(all_mags_indexed, str_ends(Water, "GGMM"))

medium_small_orgs_i_tabellen <- grep("MMQQ$", all_mags_indexed$Water)
medium_small_orgs <- all_mags_indexed[medium_small_orgs_i_tabellen,]

medium_big_orgs_i_tabellen <- grep("QQSS$", all_mags_indexed$Water)
medium_big_orgs <- all_mags_indexed[medium_big_orgs_i_tabellen,]

big_orgs_i_tabellen <- grep("SSUU$", all_mags_indexed$Water)
big_orgs <- all_mags_indexed[big_orgs_i_tabellen,]

#Create four table files: one for small_orgs, one for medium_small_orgs, etc. (write.table())
write.table(small_orgs, "GGMM_all_mags_indexed_GeneID_corr", quote = FALSE, sep = " ", row.names=F)
write.table(medium_small_orgs, "MMQQ_all_mags_indexed_GeneID_corr", quote = FALSE, sep = " ", row.names=F)
write.table(medium_big_orgs, "QQSS_all_mags_indexed_GeneID_corr", quote = FALSE, sep = " ", row.names=F)
write.table(big_orgs, "SSUU_all_mags_indexed_GeneID_corr", quote = FALSE, sep = " ", row.names=F)
