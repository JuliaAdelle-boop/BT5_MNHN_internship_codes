library(vegan)
library(data.table)
library(dplyr)
library(caret)
library(tidyr)

#### Opening the merged table that contains both the Abundance and the Envi data, necessary to create the response and explanatory matrices

See <- read.table("Merged_UNKNOWNS_GGMM_for_RDA", sep = "\t", header = T)


#### Creation of the response matrix
## Resp_Matrix is a version of See with the mean value of MAG_Cov_Y abundance of a CC at every station
Resp_Matrix <- See %>%
  group_by(CC) %>%
  group_by(Station, .add = TRUE) %>%
  summarise(
    MAG_Cov_Y = mean(MAG_Cov_Y),
    .groups = 'keep')


Resp_Matrix <- as.data.table(Resp_Matrix)
pivot_Resp_Matrix <- dcast(Resp_Matrix, Station ~ CC, value.var = "MAG_Cov_Y", fill = 0)
pivot_Resp_Matrix <- as.data.frame(pivot_Resp_Matrix)
##response_matrix is the final response matrix, which will be fed to rda()
response_matrix <- pivot_Resp_Matrix[,-1]

#### Creation of the explanatory matrix

### One-hot encoding of character Environmental values

## One-hot encoding of Depth values
#To check how many possible values of Depth there are
Depth_values <- unique(See$Depth)
print(Depth_values) 
#only 2 possible values (as expected)
#One-hot encoding of Depth values within a new table called See_recode
Depth_values_mod <- c("SRF" = 0, "DCM" = 1)
See_recode <- See %>%
  mutate(Depth = recode(Depth, !!!Depth_values_mod),
  )


### Expl_Matrix is a modified version of See_recode which only contains the environmental parameters. Like with Resp_Matrix, the environmental parameters have been "compressed" so that there is only one value of each per Station.
#For Depth, the mean is computed because there are only two possible values of Depth, which means that the mean gives us a score representing "in general, this CC is mostly associated to the SRF or the DCM" 

#pre_Expl is a pre-version of Expl_Matrix, where the envi parameters we don't want in the final Expl-Matrix are preemptively removed
pre_Expl <- See_recode %>% select(-CC, -MAG_Cov_Y, -Gene_Expression, -Gene_Expression_Index, -Biogeographical.province, -Season.moment)

Expl_Matrix <- pre_Expl %>%
  group_by(Station) %>%
  summarise(across(everything(), mean))

###explanatory_matrix is the final explanatory matrix, which will be fed to rda(): it's Expl_Matrix without the Abundance columns, without the 'CC' column, and without 'Station'
explanatory_matrix <- Expl_Matrix[,-1]

#### The RDA (on a Hellinger-transformed resp. matrix and z-centred expl. matrix)

hellinger_resp <- decostand(response_matrix, method = "hellinger")
#To identify columns and remove with 0 variance, as they can cause error during the z-standardization
zero_var_columns <- sapply(explanatory_matrix, function(x) sd(x) == 0)
print(names(explanatory_matrix)[zero_var_columns])
explanatory_matrix_cleaned <- explanatory_matrix[, !zero_var_columns]

zcentred_expl<- decostand(explanatory_matrix_cleaned, method = "standardize")

hellinger_rda <- rda(hellinger_resp ~., data = zcentred_expl)

summary(hellinger_rda)
print("Distinction between summary(rda) and print(hellinger_rda)")
print(hellinger_rda)

png("RDA_Hellinger_UNKNOWNS_GGMM_MAGCovY.png", width = 800, height = 1000)

plot(hellinger_rda)
row_names <- pivot_Resp_Matrix$Station
#To add the real names of the stations
text(hellinger_rda, display = "sites", labels = row_names, cex = 0.7, col = "black")

dev.off()

#permutation test to test the significance of the RDA model (using default value of 999 permutation executed)
print("Permutation test on the RDA")
anova.cca(hellinger_rda)