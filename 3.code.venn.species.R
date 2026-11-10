rm(list = ls())




#_________________SHARED AND UNIQUE OTUs______________________#

###################################################
########### A. IMPORT AND LOAD PACKAGES ###########
###################################################


# Package import 
PkgTest <- function(x)
{
  if(!require(x,character.only = TRUE))
  {
    install.packages(x, dep=TRUE, repos = "http://cran.us.r-project.org")
    
    if (!require(x, character.only = TRUE))
    {
      if (!require("BiocManager", quietly = TRUE))
      {
        install.packages("BiocManager")
        BiocManager::install(version = "3.18")
      }
      if(!require(x, character.only=TRUE))
      {stop(x, ":Package not found")
      }
    }
  }
}


# Create list of required packages
list.of.pack <- c("vegan", "plyr", "dplyr", "reshape2", "ComplexHeatmap","ggplot2","ggVennDiagram",
                  "circlize", "EcolUtils", "data.table", "RColorBrewer", "ggpubr","pairwiseAdonis", "car","funrar","emmeans", "VennDiagram")

# Create list of installed packages
pkges <- installed.packages() [,"Package"]
# Install and load all required packages
for (pk in list.of.pack) {
  PkgTest(pk)
}

#download missing packages (los que no he conseguido)
if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/ggpubr")
install.packages("emmeans",repos = "http://cran.us.r-project.org")
devtools::install_github("GuillemSalazar/EcolUtils")
devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
BiocManager::install("systemPipeR")


library(systemPipeR)
library(EcolUtils)
library(pairwiseAdonis)
library(emmeans)
library(ggpubr)
library(devtools)
library(vegan)
library(plyr)
library(dplyr)
library(reshape2)
library(ComplexHeatmap) 
library (ggplot2)
library(ggVennDiagram)
library(circlize)
library(data.table)
library(RColorBrewer)
library(ggpubr) 
library(car) 
library(funrar)


## set directory
setwd("C:/Users/Usuario/Desktop/2025 FINAL PAPER MASTER/0. CSV_DATA_EXCELLS")


################################################
######## B. VENN DIAGRAM ORDERS/REGION##########
################################################

## B.1 Load data taxonomy and abundance for only porifera
library(systemPipeR)
porifera<-read.csv("ABUND_TAX.csv")


# B.2 Keep only abundance and group by order and island
# exclude taxonomic columns that are not "order"
new_df <- porifera[, -c(1, 2, 3, 4, 6, 7, 8)]
# Extract column groups (first three letters of column names)
column_groups <- substr(names(new_df)[-1], 1, 3)
# Aggregate the columns and sum their values for each row
aggregated_df <- aggregate(. ~ order, data = new_df, FUN = sum)
# Print the aggregated dataframe
print(aggregated_df)
# Create a new dataframe to store aggregated values by group
grouped_df <- data.frame(order = aggregated_df$order)
# Iterate over unique groups
for (group_name in unique(substr(names(aggregated_df)[-1], 1, 3))) {
  # Subset columns that start with the same three letters
  group_columns <- grep(paste0("^", group_name), names(aggregated_df)[-1], value = TRUE)
  # Sum values for each row in the group and store in new dataframe
  grouped_df[[group_name]] <- rowSums(aggregated_df[group_columns])
}
# Print the new dataframe
print(grouped_df)


## B.3 Create a presence/absence dataframe 
# Create a copy of the original dataframe
presence_absence_df <- grouped_df
# Replace non-zero values with 1 in columns 2 to 9
presence_absence_df[, 2:9][presence_absence_df[, 2:9] != 0] <- 1
# Set the first column as row names
rownames(presence_absence_df) <- presence_absence_df$order
# Remove the first column
presence_absence_df <- presence_absence_df[, -1]


## B.4 Extract presence-absence data for each category
MHI <- presence_absence_df[, c("HAW", "KAU", "OAH", "MAI")]
NWHI <- presence_absence_df[, c("FFS", "PHR", "KUR", "LIS")]


## B.5 Create df with only presence for each category
# B.5.1  MHI
# Find rows (orders) where at least one island has a presence
orders_with_presence_MHI <- rownames(MHI)[apply(MHI, 1, any)]
# Subset the original dataframe to get only orders with presence
MHI_orders_with_presence <- MHI[orders_with_presence_MHI, ]
# Extract the row names (orders) from the dataframe
orders_list_MHI <- rownames(MHI_orders_with_presence)
# Print the list of orders
print(orders_list_MHI)
# Convert the orders_list to a character vector
MHI <- as.character(orders_list_MHI)

# B.5.2  NWHI
# Find rows (orders) where at least one island has a presence
orders_with_presence_NWHI <- rownames(NWHI)[apply(NWHI, 1, any)]
# Subset the original dataframe to get only orders with presence
NWHI_orders_with_presence <- NWHI[orders_with_presence_NWHI, ]
# Extract the row names (orders) from the dataframe
orders_list_NWHI <- rownames(NWHI_orders_with_presence)
# Print the list of orders
print(orders_list_NWHI)
# Convert the orders_list to a character vector
NWHI <- as.character(orders_list_NWHI)




install.packages("VennDiagram")
library(VennDiagram)
library(grid)

# Calculate values
area1 <- length(MHI)
area2 <- length(NWHI)
cross.area <- length(intersect(MHI, NWHI))

# Draw Venn diagram
venn.plot <- draw.pairwise.venn(
  area1 = length(MHI),
  area2 = length(NWHI),
  cross.area = length(intersect(MHI, NWHI)),
  
  category = c("MHI", "NWHI"),
  
  lwd = 4,
  col = c("red", "blue"),
  fill = c(NA, NA),
  
  cex = 2,
  cat.cex = 2,
  cat.col = c("red", "blue"),
  
  scaled = FALSE   
)

grid.draw(venn.plot)






# Unique orders in each region
unique_MHI <- setdiff(MHI, NWHI)
unique_NWHI <- setdiff(NWHI, MHI)

# Print results
cat("Unique orders in Main Hawaiian Islands:\n")
print(unique_MHI)

cat("\nUnique orders in Northwestern Hawaiian Islands:\n")
print(unique_NWHI)


################################################
######## C. VENN DIAGRAM SPECIES/REGION#########
################################################

## C.1 Load data and packages
rm(list=ls())
porifera<-read.csv("ABUND_TAX.csv")


## C.2 Keep only abundance and group by species and island
# Create a new dataframe excluding columns that are not "species"
df_species <- porifera[, -c(1, 2, 3, 4, 5, 6, 7)]
# Extract column groups (first three letters of column names)
column_groups <- substr(names(df_species)[-1], 1, 3)
# Aggregate the columns and sum their values for each row
aggregated_df_sp <- aggregate(. ~ species, data = df_species, FUN = sum)
# Create a new dataframe to store aggregated values by group
grouped_df_sp <- data.frame(species = aggregated_df_sp$species)
# Iterate over unique groups
for (group_name in unique(substr(names(aggregated_df_sp)[-1], 1, 3))) {
  # Subset columns that start with the same three letters
  group_columns <- grep(paste0("^", group_name), names(aggregated_df_sp)[-1], value = TRUE)
  # Sum values for each row in the group and store in new dataframe
  grouped_df_sp[[group_name]] <- rowSums(aggregated_df_sp[group_columns])
}


## C.3 Create a presence/absence dataframe 
# Create a copy of the original dataframe
df <- grouped_df_sp
# Replace non-zero values with 1 in columns 2 to 9
df[, 2:9][df[, 2:9] != 0] <- 1
# Print the presence/absence dataframe
print(df)
# Set the first column as row names
rownames(df) <- df$species
# Remove the first column
df <- df[, -1]
# Print the updated dataframe
print(df)


## C.4 Extract presence-absence data for each region
MHI_sp <- df[, c("HAW", "KAU", "OAH", "MAI")]
NWHI_sp <- df[, c("FFS", "PHR", "KUR", "LIS")]


## C.5 Create df with only presence for each category
# C.5.1  MHI
# Find rows (species) where at least one island has a presence 
sp_with_presence_MHI <- rownames(MHI_sp)[apply(MHI_sp, 1, any)]
# Subset the original dataframe to get only orders with presence
MHI_sp_with_presence <- MHI_sp[sp_with_presence_MHI, ]
# Extract the row names (orders) from the dataframe
sp_list_MHI <- rownames(MHI_sp_with_presence)
# Print the list of sp in MHI
print(sp_list_MHI)
# Convert the orders_list to a character vector
MHI <- as.character(sp_list_MHI)

# C.5.2  NWHI
# Find rows (species) where at least one island has a presence NWHI
sp_with_presence_NWHI <- rownames(NWHI_sp)[apply(NWHI_sp, 1, any)]
# Subset the original dataframe to get only orders with presence
NWHI_sp_with_presence <- NWHI_sp[sp_with_presence_NWHI, ]
# Extract the row names (orders) from the dataframe
sp_list_NWHI <- rownames(NWHI_sp_with_presence)
# Print the list of sp in NWHI
print(sp_list_NWHI)
# Convert the orders_list to a character vector
NWHI <- as.character(sp_list_NWHI)



library(VennDiagram)
library(grid)

# Calculate values
area1 <- length(MHI)
area2 <- length(NWHI)
cross.area <- length(intersect(MHI, NWHI))

# Draw Venn diagram
venn.plot <- draw.pairwise.venn(
  area1 = area1,
  area2 = area2,
  cross.area = cross.area,
  
  category = c("MHI", "NWHI"),
  
  lwd = 4,
  col = c("red", "blue"),
  fill = c(NA, NA),
  
  cex = 2,
  cat.cex = 2,
  
  cat.col = c("red", "blue"),
  scaled = FALSE 
  
)

grid.draw(venn.plot)


# Unique orders in each region
unique_MHI <- setdiff(MHI, NWHI)
unique_NWHI <- setdiff(NWHI, MHI)

# Print results
cat("Unique species in Main Hawaiian Islands:\n")
print(unique_MHI)

cat("Unique species in Northwestern Hawaiian Islands:\n")
print(unique_NWHI)



#############################################################
######## D. HEATMAP SHARED AND UNIQUE ORDERS/ISLANDS#########
#############################################################

## D.1 Load data
rm(list=ls())
porifera<-read.csv("tax_abund_porifera.csv")
metaporifera<-read.csv("metaD.csv")


## D.2 Keep only abundance and group by orders and island
# Create a new dataframe excluding specified columns
new_df <- porifera[, -c(1, 2, 3, 4, 6, 7, 8)]
new_df
length(new_df)
# Extract column groups (first three letters of column names)
column_groups <- substr(names(new_df)[-1], 1, 3)
length(column_groups)
# Aggregate the columns and sum their values for each row
aggregated_df <- aggregate(. ~ order, data = new_df, FUN = sum)
# Print the aggregated dataframe
print(aggregated_df)
# Create a new dataframe to store aggregated values by group
grouped_df <- data.frame(order = aggregated_df$order)
# Iterate over unique groups
for (group_name in unique(substr(names(aggregated_df)[-1], 1, 3))) {
  # Subset columns that start with the same three letters
  group_columns <- grep(paste0("^", group_name), names(aggregated_df)[-1], value = TRUE)
  # Sum values for each row in the group and store in new dataframe
  grouped_df[[group_name]] <- rowSums(aggregated_df[group_columns])
}


## D.3 Create a presence/absence dataframe
# Create a copy of the original dataframe
presence_absence_df <- grouped_df
# Replace non-zero values with 1 in columns 2 to 9
presence_absence_df[, 2:9][presence_absence_df[, 2:9] != 0] <- 1
# Print the presence/absence dataframe
print(presence_absence_df)
# Set the first column as row names
rownames(presence_absence_df) <- presence_absence_df$order
# Remove the first column
presence_absence_df <- presence_absence_df[, -1]


## D.4 Do Heatmap
library(ComplexHeatmap)
# Define colors for heatmap
colors <- c("red3", "palegreen3")
# Create heatmap
Heatmap(as.matrix(presence_absence_df), 
        col = colors,
        name = "Presence/Absence",
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 8),
        show_row_names = TRUE,
        show_column_names = TRUE)


