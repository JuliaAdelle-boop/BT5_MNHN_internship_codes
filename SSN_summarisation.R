library("dplyr")
library("tidyr")

taxostats = read.table('80CC_stats_with_taxo', header=T)

#### Size of the CCs
Minimum = min(taxostats$size)
Maximum = max(taxostats$size)
Mean = round(mean(taxostats$size), digits = 2)

#### FUNCTIONAL SCORES

#### Mean homscores (for each annotation database)
pre_eggnog_homsc = taxostats$Homogeneity.eggnog[!is.na(taxostats$Homogeneity.eggnog)]
nb_NA_homscores_eggnog = sum(is.na(taxostats$Homogeneity.eggnog))
Mean_eggnog_homsc = paste(round(mean(pre_eggnog_homsc), digits = 2), " (", nb_NA_homscores_eggnog, ")", sep = "")

pre_pfam_homsc = taxostats$Homogeneity.pfam[!is.na(taxostats$Homogeneity.pfam)]
nb_NA_homscores_pfam = sum(is.na(taxostats$Homogeneity.pfam))
Mean_pfam_homsc = paste(round(mean(pre_pfam_homsc), digits = 2), " (", nb_NA_homscores_pfam, ")", sep="")

pre_kegg_homsc = taxostats$Homogeneity.kegg[!is.na(taxostats$Homogeneity.kegg)]
nb_NA_homscores_kegg = sum(is.na(taxostats$Homogeneity.kegg))
Mean_kegg_homsc = paste(round(mean(pre_kegg_homsc), digits = 2), " (", nb_NA_homscores_kegg, ")", sep="")


#### Mean unkscores (for each annotation database)
Mean_eggnog_unksc = paste("Eggnog:", round(mean(taxostats$Unknowns.eggnog), digits = 2))
Mean_pfam_unksc = paste("PFAM:", round(mean(taxostats$Unknowns.pfam), digits = 2))
Mean_kegg_unksc = paste("Kegg:", round(mean(taxostats$Unknowns.kegg), digits = 2))


#### Unknown quantifs
eggnogs_noNAs = nrow(taxostats[taxostats$Unknowns.eggnog == 0,])
eggnogs_onlyNAs = nrow(taxostats[taxostats$Unknowns.eggnog == 1,])
eggnogs_atleast_1_anno = nrow(taxostats[taxostats$Unknowns.eggnog !=1,])

pfam_noNAs = nrow(taxostats[taxostats$Unknowns.pfam == 0,]) 
pfam_onlyNAs = nrow(taxostats[taxostats$Unknowns.pfam == 1,])
pfam_atleast_1_anno = nrow(taxostats[taxostats$Unknowns.pfam !=1,])

kegg_noNAs = nrow(taxostats[taxostats$Unknowns.kegg == 0,])
kegg_onlyNAs = nrow(taxostats[taxostats$Unknowns.kegg == 1,])
kegg_atleast_1_anno = nrow(taxostats[taxostats$Unknowns.kegg !=1,])


# Number (and percentage) of CCs which don't contain any annotated protein with regards to any database
rowindexes_with_only_completely_unannotated_prots <- which(taxostats$Unknowns.eggnog == 1 & taxostats$Unknowns.pfam == 1 & taxostats$Unknowns.kegg == 1)
test <- which(is.na(taxostats$Homogeneity.eggnog) & is.na(taxostats$Homogeneity.kegg) & is.na(taxostats$Homogeneity.pfam))
number_of_CCs_with_only_completely_unannotated_prots <- length(rowindexes_with_only_completely_unannotated_prots)
test_svit <- length(test)
print(taxostats$CC_ID[rowindexes_with_only_completely_unannotated_prots])
print(taxostats$CC_ID[test])
percentage_of_CCs_with_only_completely_unannotated_prots_over_total = round(((number_of_CCs_with_only_completely_unannotated_prots / nrow(taxostats))*100), digits = 2)
test_svit_svit <- round(((test_svit / nrow(taxostats))*100), digits = 2)

#### Percentages of unknown quantifs over the total number of CCs
eggnogs_noNAs_percentage = round(((eggnogs_noNAs) / nrow(taxostats) * 100), digits = 2)
eggnogs_onlyNAs_percentage = round(((eggnogs_onlyNAs) / nrow(taxostats) * 100), digits = 2)
eggnogs_atleast_1_anno_percentage = round(((eggnogs_atleast_1_anno) / nrow(taxostats) * 100), digits = 2)

pfam_noNAs_percentage = round(((pfam_noNAs) / nrow(taxostats) * 100), digits = 2)
pfam_onlyNAs_percentage = round(((pfam_onlyNAs) / nrow(taxostats) * 100), digits = 2)
pfam_atleast_1_anno_percentage = round(((pfam_atleast_1_anno) / nrow(taxostats) * 100), digits = 2)

kegg_noNAs_percentage = round(((kegg_noNAs) / nrow(taxostats) * 100), digits = 2)
kegg_onlyNAs_percentage = round(((kegg_onlyNAs) / nrow(taxostats) * 100), digits = 2)
kegg_atleast_1_anno_percentage = round(((kegg_atleast_1_anno) / nrow(taxostats) * 100), digits = 2)


####Combining the unknown quantifs and percentages of unknown quantifs into unique variables, so that they both fit into a single 'Unknown quantifications' column
eggnogs_noNAs_both = paste(eggnogs_noNAs, " (", eggnogs_noNAs_percentage, "%)", sep = "")
eggnogs_onlyNAs_both = paste(eggnogs_onlyNAs, " (", eggnogs_onlyNAs_percentage, "%)", sep = "")
eggnogs_atleast_1_anno_both = paste(eggnogs_atleast_1_anno, " (", eggnogs_atleast_1_anno_percentage, "%)", sep = "")
pfam_noNAs_both = paste(pfam_noNAs, " (", pfam_noNAs_percentage, "%)", sep = "")
pfam_onlyNAs_both = paste(pfam_onlyNAs, " (", pfam_onlyNAs_percentage, "%)", sep = "")
pfam_atleast_1_anno_both = paste(pfam_atleast_1_anno, " (", pfam_atleast_1_anno_percentage, "%)", sep = "")
kegg_noNAs_both = paste(kegg_noNAs, " (", kegg_noNAs_percentage, "%)", sep = "")
kegg_onlyNAs_both = paste(kegg_onlyNAs, " (", kegg_onlyNAs_percentage, "%)", sep = "")
kegg_atleast_1_anno_both = paste(kegg_atleast_1_anno, " (", kegg_atleast_1_anno_percentage, "%)", sep = "")
only_NAs_throughout_all_databases = paste(number_of_CCs_with_only_completely_unannotated_prots, " (", percentage_of_CCs_with_only_completely_unannotated_prots_over_total, "%)", sep = "")

#### TAXONOMY SCORES

#### number of CCs associated to only 1 "species level"
CCs_associated_to_only_phylum = nrow(taxostats[taxostats$Homogeneity.taxo.Phylum == 1,])
CCs_associated_to_only_class = nrow(taxostats[taxostats$Homogeneity.taxo.Class == 1,])
CCs_associated_to_only_order = nrow(taxostats[taxostats$Homogeneity.taxo.Order == 1,])
CCs_associated_to_only_family = nrow(taxostats[taxostats$Homogeneity.taxo.Family == 1,])
CCs_associated_to_only_genus = nrow(taxostats[taxostats$Homogeneity.taxo.Genus == 1,])

#### percentage of CCs associated to only 1 "species level" over the total number of CCs
percentage_of_CCs_associated_to_only_phylum_over_total = round(((CCs_associated_to_only_phylum / nrow(taxostats)) * 100), digits = 2)
percentage_of_CCs_associated_to_only_class_over_total = round(((CCs_associated_to_only_class / nrow(taxostats)) * 100), digits = 2)
percentage_of_CCs_associated_to_only_order_over_total = round(((CCs_associated_to_only_order / nrow(taxostats)) * 100), digits = 2)
percentage_of_CCs_associated_to_only_family_over_total = round(((CCs_associated_to_only_family / nrow(taxostats)) * 100), digits = 2)
percentage_of_CCs_associated_to_only_genus_over_total = round(((CCs_associated_to_only_genus / nrow(taxostats)) * 100), digits = 2)

#### percentage of CCs associated to only 1 "species level" over the number of CCs with at least 1 "species level" annotation
percentage_of_CCs_associated_to_only_phylum_over_atleast_1_annot = round(((CCs_associated_to_only_phylum / nrow(taxostats[!is.na(taxostats$Homogeneity.taxo.Phylum),])) * 100), digits = 2)
percentage_of_CCs_associated_to_only_class_over_atleast_1_annot = round(((CCs_associated_to_only_class / nrow(taxostats[!is.na(taxostats$Homogeneity.taxo.Class),])) * 100), digits = 2)
percentage_of_CCs_associated_to_only_order_over_atleast_1_annot = round(((CCs_associated_to_only_order / nrow(taxostats[!is.na(taxostats$Homogeneity.taxo.Order),])) * 100), digits = 2)
percentage_of_CCs_associated_to_only_family_over_atleast_1_annot = round(((CCs_associated_to_only_family / nrow(taxostats[!is.na(taxostats$Homogeneity.taxo.Family),])) * 100), digits = 2)
percentage_of_CCs_associated_to_only_genus_over_atleast_1_annot = round(((CCs_associated_to_only_genus / nrow(taxostats[!is.na(taxostats$Homogeneity.taxo.Genus),])) * 100), digits = 2)


#### Number of CCs which exclusively contain proteins from annotated MAGS (so, not functional annotations, but taxonomic annotations)
only_prots_from_annotated_MAGs_phylum = nrow(taxostats[taxostats$Unknowns.taxo.Phylum == 0,])
only_prots_from_annotated_MAGs_class = nrow(taxostats[taxostats$Unknowns.taxo.Class == 0,])
only_prots_from_annotated_MAGs_order = nrow(taxostats[taxostats$Unknowns.taxo.Order == 0,])
only_prots_from_annotated_MAGs_family = nrow(taxostats[taxostats$Unknowns.taxo.Family == 0,])
only_prots_from_annotated_MAGs_genus = nrow(taxostats[taxostats$Unknowns.taxo.Genus == 0,])

#### Number of CCs which exclusively contain proteins from UNannotated MAGs (so, not functional annotations, but taxonomic annotations)
only_prots_from_unannotated_MAGs_phylum = nrow(taxostats[taxostats$Unknowns.taxo.Phylum == 1,])
only_prots_from_unannotated_MAGs_class = nrow(taxostats[taxostats$Unknowns.taxo.Class == 1,])
only_prots_from_unannotated_MAGs_order = nrow(taxostats[taxostats$Unknowns.taxo.Order == 1,])
only_prots_from_unannotated_MAGs_family = nrow(taxostats[taxostats$Unknowns.taxo.Family == 1,])
only_prots_from_unannotated_MAGs_genus = nrow(taxostats[taxostats$Unknowns.taxo.Genus == 1,])

####Percentages of the number of the 10 previous objects
percentage_only_prots_from_annotated_MAGs_phylum = round(((only_prots_from_annotated_MAGs_phylum / nrow(taxostats)) * 100), digits = 2)
percentage_only_prots_from_annotated_MAGs_class = round(((only_prots_from_annotated_MAGs_class / nrow(taxostats)) * 100), digits = 2)
percentage_only_prots_from_annotated_MAGs_order = round(((only_prots_from_annotated_MAGs_order / nrow(taxostats)) * 100), digits = 2)
percentage_only_prots_from_annotated_MAGs_family = round(((only_prots_from_annotated_MAGs_family / nrow(taxostats)) * 100), digits = 2)
percentage_only_prots_from_annotated_MAGs_genus = round(((only_prots_from_annotated_MAGs_genus / nrow(taxostats))* 100), digits = 2)
percentage_only_prots_from_unannotated_MAGs_phylum = round(((only_prots_from_unannotated_MAGs_phylum / nrow(taxostats)) * 100), digits = 2)
percentage_only_prots_from_unannotated_MAGs_class = round(((only_prots_from_unannotated_MAGs_class / nrow(taxostats))* 100), digits = 2)
percentage_only_prots_from_unannotated_MAGs_order = round(((only_prots_from_unannotated_MAGs_order / nrow(taxostats)) * 100), digits = 2)
percentage_only_prots_from_unannotated_MAGs_family = round(((only_prots_from_unannotated_MAGs_family / nrow(taxostats)) * 100), digits = 2)
percentage_only_prots_from_unannotated_MAGs_genus = round(((only_prots_from_unannotated_MAGs_genus / nrow(taxostats)) * 100), digits = 2)

#### Final table
Mini_Faure <- data.frame(
  Sizes_description = c("Mean", "Min", "Max", rep("",7)),
  Sizes = c(Mean,Minimum,Maximum,rep("", 7)),
  Homogeneity = c("Mean homogeneity with EggNOG annotations", "Mean homogeneity with PFAM annotations","Mean homogeneity with KEGG annotations", rep("",7)),
  Mean_Homscores = c(Mean_eggnog_homsc,Mean_pfam_homsc,Mean_kegg_homsc, rep("", 7)),
  Mean_Unkscores = c(Mean_eggnog_unksc,Mean_pfam_unksc,Mean_kegg_unksc, rep("", 7)),
  Unknown_quantif_descriptions = c("eggnogs_noNAs","eggnogs_atleast_1_anno","eggnogs_onlyNAs", "pfam_noNAs","pfam_atleast_1_anno","pfam_onlyNAs", "kegg_noNAs","kegg_atleast_1_anno","kegg_onlyNAs","Only_NAs_throughout_all_databases"),
  Unknown_quantifications = c(eggnogs_noNAs_both, eggnogs_atleast_1_anno_both, eggnogs_onlyNAs_both, pfam_noNAs_both, pfam_atleast_1_anno_both, pfam_onlyNAs_both, kegg_noNAs_both, kegg_atleast_1_anno_both, kegg_onlyNAs_both, only_NAs_throughout_all_databases),
  Species_levels_for_the_next_columns = c("Phylum","Class","Order","Family","Genus",rep("", 5)),
  Number_of_CCs_associated_to_only_1_species_level = c(CCs_associated_to_only_phylum,CCs_associated_to_only_class,CCs_associated_to_only_order,CCs_associated_to_only_family,CCs_associated_to_only_genus, rep("", 5)),
  Percentage_CCs_associated_to_only_1_specieslevel_over_total_number_CCs = c(percentage_of_CCs_associated_to_only_phylum_over_total,percentage_of_CCs_associated_to_only_class_over_total,percentage_of_CCs_associated_to_only_order_over_total,percentage_of_CCs_associated_to_only_family_over_total,percentage_of_CCs_associated_to_only_genus_over_total,rep("", 5)),
  percentage_of_CCs_associated_to_only_1_specieslevel_over_atleast_1_annot = c(percentage_of_CCs_associated_to_only_phylum_over_atleast_1_annot,percentage_of_CCs_associated_to_only_class_over_atleast_1_annot,percentage_of_CCs_associated_to_only_order_over_atleast_1_annot,percentage_of_CCs_associated_to_only_family_over_atleast_1_annot,percentage_of_CCs_associated_to_only_genus_over_atleast_1_annot,rep("",5)),
  Taxonomic_unknown_quantif_descriptions = c("Phylum: Number of CCs which exclusively contain proteins from annotated MAGS","Phylum: Number of CCs which exclusively contain proteins from UNannotated MAGs","Class: Number of CCs which exclusively contain proteins from annotated MAGs","Class: Number of CCs which exclusively contain proteins from UNannotated MAGs","Order: Number of CCs which exclusively contain proteins from annotated MAGs","Order: Number of CCs which exclusively contain proteins from UNannotated MAGs","Family: Number of CCs which exclusively contain proteins from annotated MAGs","Family: Number of CCs which exclusively contain proteins from UNannotated MAGs","Genus: Number of CCs which exclusively contain proteins from annotated MAGs","Genus: Number of CCs which exclusively contain proteins from UNannotated MAGs"),
  Taxonomic_Unknown_quantifs = c(only_prots_from_annotated_MAGs_phylum,only_prots_from_unannotated_MAGs_phylum, only_prots_from_annotated_MAGs_class,only_prots_from_unannotated_MAGs_class, only_prots_from_annotated_MAGs_order,only_prots_from_unannotated_MAGs_order, only_prots_from_annotated_MAGs_family, only_prots_from_unannotated_MAGs_family, only_prots_from_annotated_MAGs_genus, only_prots_from_unannotated_MAGs_genus),                                           
  Percentages_of_taxonomic_unknown_quantifs = c(percentage_only_prots_from_annotated_MAGs_phylum, percentage_only_prots_from_unannotated_MAGs_phylum, percentage_only_prots_from_annotated_MAGs_class, percentage_only_prots_from_unannotated_MAGs_class, percentage_only_prots_from_annotated_MAGs_order, percentage_only_prots_from_unannotated_MAGs_order, percentage_only_prots_from_annotated_MAGs_family, percentage_only_prots_from_unannotated_MAGs_family, percentage_only_prots_from_annotated_MAGs_genus, percentage_only_prots_from_unannotated_MAGs_genus)
)

library(openxlsx)
write.xlsx(Mini_Faure, file = "SSN_summary")

