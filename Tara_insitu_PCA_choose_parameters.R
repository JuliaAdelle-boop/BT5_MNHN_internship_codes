library(caret)
library(MASS)
library(dplyr)
library(VIM)

Tara_insitu <- read.table("SMAGS_insitu_WOA_Env.tsv", sep = "\t", header = T)
#of the next three lines, choose the one appropriate to the current round of PCA performed
cols_to_remove <- c("Ocean.region", "Marine.biome", "Moon.phase.prop", "Season", "DepthBathy", "Depth.top.min", "Depth.bottom.max")
cols_to_remove <- c("Ocean.region", "Marine.biome", "Moon.phase.prop", "Season", "DepthBathy", "Depth.top.min", "Depth.bottom.max", "Temperature_WOA", "Salinity_WOA", "Longitude_WOA", "Latitude_WOA", "Si_WOA", "NO3_WOA")
cols_to_remove <- c("Ocean.region", "Marine.biome", "Moon.phase.prop", "Season", "DepthBathy", "Depth.top.min", "Depth.bottom.max", "Temperature_WOA", "Salinity_WOA", "Longitude_WOA", "Latitude_WOA", "Si_WOA", "NO3_WOA", "Part.beam.att.coef")
Tara_insitu <- Tara_insitu[, -which(names(Tara_insitu) %in% cols_to_remove)]

#Numerical encoding of character values
Tara_insitu$Biogeographical.province <- as.numeric(factor(Tara_insitu$Biogeographical.province))
Tara_insitu$Depth <- as.numeric(factor(Tara_insitu$Depth))
Tara_insitu$Season.moment <- as.numeric(factor(Tara_insitu$Season.moment))


#To see if some columns still have non-numeric values
column_classes <- sapply(Tara_insitu, class)
column_classes <- names(column_classes[column_classes != "numeric"])
print(column_classes)
#To see what columns contain integer-class values
column_classes <- sapply(Tara_insitu, class)
integer_columns <- names(column_classes[column_classes == "integer"])
print(integer_columns)

#To change integer columns into numeric columns
Tara_insitu$Depth.nominal <- as.numeric(Tara_insitu$Depth.nominal)
Tara_insitu$Depth.chloro.max <- as.numeric(Tara_insitu$Depth.chloro.max)
Tara_insitu$Depth.Min.O2 <- as.numeric(Tara_insitu$Depth.Min.O2)
Tara_insitu$Sunshine.duration <- as.numeric(Tara_insitu$Sunshine.duration)
Tara_insitu$Depth.max.Brunt.V..is..l.. <- as.numeric(Tara_insitu$Depth.max.Brunt.V..is..l..)
Tara_insitu$Depth.Mixed.Layer <- as.numeric(Tara_insitu$Depth.Mixed.Layer)
Tara_insitu$Depth.Max.O2 <- as.numeric(Tara_insitu$Depth.Max.O2)


##KNN-Imputation (without preProcess)
KNNImputed <- kNN(Tara_insitu, imp_var = FALSE)

####Imputing the missing values (with different methods)
#KNNImputation <- preProcess(Tara_insitu, method = "knnImpute") ##knnImpute fuckar upp värden som inte fattas, här
#KNNImputed <- predict(KNNImputation, Tara_insitu)
#SVDImputation <- preProcess(Tara_insitu, method = "svdImpute")
#SVDImputed <- predict(SVDImputation, Tara_insitu)

#### Finding outliers in the data (KNNImputed)
##With Interquartile Range (IQR) method 
# Calculate Q1 and Q3 (=first and third quartiles)
numeric_KNNImputed <- KNNImputed %>% select_if(is.numeric)
Q1_KNN <- apply(numeric_KNNImputed, 2, quantile, 0.25, na.rm = TRUE)
Q3_KNN <- apply(numeric_KNNImputed, 2, quantile, 0.75, na.rm = TRUE)
# Compute the IQR
IQR_KNN <- Q3_KNN - Q1_KNN
# Detect outliers (values which deviate "too much" (1.5 is quite arbitrarily chosen) from the IQR)
outliers_KNN <- t(apply(numeric_KNNImputed, 1, function(x) (x < (Q1_KNN - 1.5 * IQR_KNN)) | (x > (Q3_KNN + 1.5 * IQR_KNN))))
# View outliers
outliers_detected_KNN <- numeric_KNNImputed[apply(outliers_KNN, 1, any), ]
head(outliers_detected_KNN)
##==> 217 outliers! This method is too stringent and/or I shouldn't have applied it
#to the one-hot encoded categorical values

##With Isolation Forest
# Install and load required package
#install.packages("IsolationForest")
#library(IsolationForest)
# Fit the model
#iso_forest <- isolation.forest(KNNImputed, ntrees = 100)
# Get anomaly scores
#anomaly_scores <- predict(iso_forest, KNNImputed)
# Detect outliers (using a threshold of 0.6, adjust as necessary)
#outliers <- anomaly_scores > 0.6
# View outliers
#outliers_detected <- KNNImputed[outliers, ]
#head(outliers_detected)

#### z-standardization and PCA
pca<-prcomp(numeric_KNNImputed,retx=T,center=T,scale.=T) 

screeplot(pca,npcs=length(pca$sdev),type = "lines") 
#have merely copied Pavla's code here, so I don't know if the parameters of the screeplot are optimal, here.
#Not that the screeplot is super important here (at least I think ^^)

# Plot proportion of variance (to better see what PCs are worth looking at)
pca.pct<-100*round(summary(pca)$importance[2,],3)          
barplot(pca.pct)

####Get the loadings of PC1 and PC2, by order of weight
pca.scores<-pca$x
#x is pca's rotation matrix
head(pca.scores)
(pca.scores12<-cbind(pca$x[,1],pca$x[,2]))
(pca.scores34<-cbind(pca$x[,3],pca$x[,4]))

eigenvec<-pca$rotation		 	  
eigenvec123<-pca$rotation[,1:3]
eigenvec123<-as.data.frame(eigenvec123)
PC1_loadings <- eigenvec123$PC1
PC2_loadings <- eigenvec123$PC2
PC3_loadings <- eigenvec123$PC3
ordered_PC1 <- order(PC1_loadings, decreasing = TRUE)
ordered_PC2 <- order(PC2_loadings, decreasing = TRUE)
ordered_PC3 <- order(PC3_loadings, decreasing = TRUE)
PC1_df_reordered <- eigenvec123[ordered_PC1, , drop = FALSE]
PC1_df_reordered <- subset(PC1_df_reordered, select = "PC1")
PC2_df_reordered <- eigenvec123[ordered_PC2, , drop = FALSE]
PC2_df_reordered <- subset(PC2_df_reordered, select = "PC2")
PC3_df_reordered <- eigenvec123[ordered_PC3, , drop = FALSE]
PC3_df_reordered <- subset(PC3_df_reordered, select = "PC3")

write.table(PC1_df_reordered, "round2_PC1_corr", sep = "\t")
write.table(PC2_df_reordered, "round2_PC2_corr", sep = "\t")
write.table(PC3_df_reordered, "round2_PC3_corr", sep = "\t")
