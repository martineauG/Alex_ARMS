############################################
########  A. Load working directory ########
############################################


##A. 1 Set working directory
setwd("~/U_of_H/SPONGES/Alex/R_clean_tables/Metadata_organization/")

##A.2 Install libraires
#Install necessary packages and load libraries
library ("data.table")
library("CHNOSZ")
library ("data.table")
library("RFLPtools")
library ("DECIPHER")
library ("magrittr")
library ("dplyr")
library ("tidyverse")
library ("plyr")
library("taxize")
library("dplyr")
library("purrr")
library("reshape2")
library("seqinr")
library("Biostrings")


## A.3 Read files
ASV_table <- as.data.frame(fread ("Alex_abundance_P2_no_singletons_040824.csv", header = TRUE))
haplosclerida_metadata <- as.data.frame(fread ("Martineau_Haliclona_Vicente_metadata_110123_addedVicente.csv", header = TRUE))
metadataMolly <- as.data.frame(fread ("ARMS_Locations_Gabby - ARMS_Locations_Gabby.csv", header = TRUE))
metadatafreezer <- as.data.frame(fread ("ARMS_freezer_samples_ID_MT_GM.csv", header = TRUE))



#####################################
########  B. Data wrangling  ########
#####################################

## B.1 Change the samples name in haplosclerida metadata
real_samples_names <- colnames(ASV_table)[2:ncol(ASV_table)]
new_samples_names <- sub("_", "", real_samples_names)
new_samples_names <-  sub("_.*", "", new_samples_names)

# verify if there are duplicates
length(unique(new_samples_names)) == length(real_samples_names)

# Change samples name for Haplosclerida data
real_samples_haplosclerida <- sub("_.*", "", haplosclerida_metadata$Sample)
haplosclerida_metadata$Sample <- real_samples_haplosclerida

#Delete duplicates Haplosclerida metadata
haplosclerida_clean <- haplosclerida_metadata %>%
  dplyr::distinct()

## B.2 Create meta data file new
Alexdf <- data.frame(new_samples_names, real_samples_names)
colnames(Alexdf)[1] = "Sample"


# Merge with existing df
Alexdf_2 <- left_join(Alexdf, haplosclerida_clean)
Alex_df_3 <- left_join(Alexdf_2, metadataMolly)

# Change sample name back
Alex_df_3$Sample <- NULL
Alex_df_3$V8 <- NULL
colnames(Alex_df_3)[1] = "Sample"
rownames(Alex_df_3) = Alex_df_3$Sample


# Save Alexs data
write.csv(Alex_df_3, "Alex_metadata_040824.csv",  row.names = F)
