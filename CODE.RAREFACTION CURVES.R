rm(list = ls())

## ----------------------1. PREPARE DATA----------------------##

#########################################################
### 1.A. Set working directory and load packages  #######
#########################################################

setwd("C:/Users/Usuario/Desktop/2025 FINAL PAPER MASTER/0.CSV_DATA_EXCELLS")

library(dplyr)
library(vegan)
library(tibble)
library(phyloseq)

########################
### 1.B. Load data #####
########################

# Taxonomy table
taxonomy <- read.csv("TAXONOMY_TABLE_ORIGINAL.csv", header = TRUE) %>%
  column_to_rownames(var = "ID")

# Abundance table
abundance <- read.csv("ABUNDANCE_TABLE_ORIGINAL.csv", header = TRUE) %>%
  column_to_rownames(var = "ID") %>%
  mutate(across(everything(), as.numeric))

# Metadata table
metadata <- read.csv("METADATA_TABLE_ORIGINAL.csv", header = TRUE) %>%
  column_to_rownames(var = "Sample")

#####################################
### 1.C. Create phyloseq object  ####
#####################################

phylo <- phyloseq(
  otu_table(as.matrix(abundance), taxa_are_rows = TRUE),
  tax_table(as.matrix(taxonomy)),
  sample_data(metadata)
)

######################################
## 1.D. Filter only Porifera taxa ###
######################################

# Keep only OTUs belonging to phylum Porifera
phylo_porifera <- subset_taxa(phylo, phylum == "Porifera")

###########################################
## 1.E. Extract abundance table ###########
###########################################

# Convert to data frame
# Rows = ARMS samples
# Columns = Porifera OTUs

abundance_porifera <- as.data.frame(t(otu_table(phylo_porifera)))

# Check structure
str(abundance_porifera)

# Optional view
View(abundance_porifera)

## ----------------------2. PREPARE METADATA----------------------##

########################
## 2.A. Keep Region ####
########################

metadata_region <- metadata %>%
  select(Region)

##########################################################
## 2.B. Ensure samples match between metadata and OTUs ###
##########################################################

# Confirm sample IDs match
all(rownames(abundance_porifera) %in% rownames(metadata_region))

# Reorder metadata to match abundance table
metadata_region <- metadata_region[rownames(abundance_porifera), , drop = FALSE]

###########################################################
## 2.C. Split abundance matrix by region ##################
###########################################################

MHI_samples <- abundance_porifera[
  metadata_region$Region == "Main Hawaiian Islands", ]

NWHI_samples <- abundance_porifera[
  metadata_region$Region == "Northwestern Hawaiian Islands", ]

############################################################
## 2.D. OPTIONAL: Convert to presence/absence ##############
############################################################

# Recommended for species accumulation curves

MHI_samples[MHI_samples > 0] <- 1
NWHI_samples[NWHI_samples > 0] <- 1

## ----------------------3. SPECIES ACCUMULATION----------------------##

###########################################################
## 3.A. Generate species accumulation curves ##############
###########################################################

spec_MHI <- specaccum(
  MHI_samples,
  method = "random"
)

spec_NWHI <- specaccum(
  NWHI_samples,
  method = "random"
)

###########################################################
## 3.B. Plot curves #######################################
###########################################################

plot(
  spec_MHI,
  col = "blue",
  lwd = 3,
  ci.type = "bar",
  ci.col = rgb(0, 0, 1, 0.2),
  ylim = c(0, 120),
  xlab = "Number of ARMS units",
  ylab = "Cumulative Porifera OTUs",
  main = "Species Accumulation Curves"
)

plot(
  spec_NWHI,
  add = TRUE,
  col = "red",
  lwd = 3,
  ci.type = "bar",
  ci.col = rgb(1, 0, 0, 0.2)
)

###########################################################
## 3.C. Add legend ########################################
###########################################################

legend(
  "bottomright",
  legend = c(
    "Main Hawaiian Islands",
    "Northwestern Hawaiian Islands"
  ),
  col = c("blue", "red"),
  lwd = 3,
  bty = "n"
)
