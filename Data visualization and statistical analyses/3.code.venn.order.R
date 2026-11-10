
####   VENN ORDERS  ####

rm(list = ls())
## set directory
setwd("C:/Users/Usuario/Desktop/2025 FINAL PAPER MASTER/0. CSV_DATA_EXCELLS")


#_________________SHARED AND UNIQUE OTUs______________________#

###################################################
########### A. IMPORT AND LOAD PACKAGES ###########
###################################################
install.packages("VennDiagram", "dplyr", "grid")

library(dplyr)
library(VennDiagram)
library(grid)


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
