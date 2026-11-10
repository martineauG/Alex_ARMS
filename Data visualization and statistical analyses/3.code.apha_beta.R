rm(list = ls())
#########################################################
####  A. Set working directory and download packages ####
#########################################################

## A.1 Set WD 
setwd("C:/Users/Usuario/Desktop/2025 FINAL PAPER MASTER/0. CSV_DATA_EXCELLS")

## A.2 Load packages
library(tidyverse)
library(patchwork)
library(phyloseq)
library(ggpubr)


#########################################################
####  B. Import all files (abund, tax and metadata) #####
#########################################################

## B.1 Taxonomy
taxonomy <- read.csv("TAXONOMY_TABLE_ORIGINAL.csv", header = TRUE) %>%
  column_to_rownames(var = "ID")

## B.2 Taxonomy and abundance with only Porifera
abundance <- read.csv("ABUNDANCE_TABLE_ORIGINAL.csv", header = TRUE) %>%
  column_to_rownames(var = "ID")%>%
  mutate_if(is.integer, as.numeric)

combined_table <- merge(abundance, taxonomy, by = "row.names") %>%
  filter(phylum == "Porifera") %>%
  column_to_rownames(var = "Row.names")

#extract numeric columns (abundance)
abundance.only.porifera <- combined_table %>%
  dplyr::select_if(is.numeric)

#extract character columns (taxonomy)
taxonomy.only.porfiera <- combined_table %>%
  dplyr::select_if(is.character)

## B.3 Meta data
metadata <- read.csv("./METADATA_TABLE_ORIGINAL.csv", header = TRUE) 


#######################################
####  C. Phyloseq: data wrangling  ####
#######################################

## C.1 Create phyloseq object
# Format tables so that rownames are Samples or IDs
metadata1 <- metadata %>%
  column_to_rownames(var = "Sample")

phylo = phyloseq(otu_table(abundance, taxa_are_rows = TRUE), tax_table (as.matrix(taxonomy)), sample_data(metadata1))
phylo

# Select only sponges
spongeshawaii <- subset_taxa(phylo, phylum == "Porifera")
spongeshawaii


## C.2 Compute chao diversity
richness_phyloseq <- estimate_richness(spongeshawaii, measures = c("Chao1")) %>%
  rename(chao_diversity = Chao1) %>%
  rownames_to_column(var = "Sample") %>%
  # merge meta data 
  left_join(metadata)

print(richness_phyloseq)


## C.3 Plot everything
phyloseq_alpha <- ggplot(richness_phyloseq, 
                         aes(x = Region, y = chao_diversity, colour = Region)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  coord_cartesian(ylim = c(2.5, 50)) +
  theme_bw() +
  scale_colour_discrete(labels = c(
    "Main Hawaiian Islands" = "MHI",
    "Northwestern Hawaiian Islands" = "NWHI"
  )) +
  annotate(
    "text",
    x = 1.5,
    y = 45,
    label = "Kruskal–Wallis test\np = 0.005516",
    size = 6
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18)
  )

phyloseq_alpha

###########################################
####  D. Statistical Testing of Groups ####
###########################################


#---------------------first check normality-----------#
boxplot(chao_diversity ~ Region, data = richness_phyloseq, main = "Boxplot of Chao Diversity by Region")

#Check for normality
by(richness_phyloseq$chao_diversity, richness_phyloseq$Region, shapiro.test)

#Plot histograms or Q-Q plots:
library(ggpubr)
ggqqplot(richness_phyloseq, x = "chao_diversity", facet.by = "Region")

#Check for homogeneity of variances (needed for ANOVA)
library(car)
leveneTest(chao_diversity ~ Region, data = richness_phyloseq)


# Perform one-way ANOVA (because all assumptions met)
anova_result <- aov(chao_diversity ~ Region, data = richness_phyloseq)
summary(anova_result)


###########################################
####  E. beta diversity / ordination ######
########## COMMUNIT COMPOSITION ###########
###########################################

## E.1 Perform NMDS ordination using Bray-Curtis dissimilarity
spongeshawaii_ord <- ordinate(spongeshawaii, method = "NMDS", distance = "bray")

# E.1.1 Plot the ordination colored by region

ordination_plot <- plot_ordination(spongeshawaii, spongeshawaii_ord, type = "samples", color = "Region", title = "NMDS Ordination (Bray-Curtis)") +
  geom_point(size = 3) +
  stat_ellipse(type = "t", level = 0.95, aes(group = Region)) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )

print(ordination_plot)



## E.1.2 ADD ELLIPSE AND LABEL OUTLIERS (mahalanobis distances)

library(tidyverse)
library(phyloseq)
library(vegan)
library(ggrepel)

# Perform NMDS ordination using Bray-Curtis dissimilarity
spongeshawaii_ord <- ordinate(spongeshawaii, method = "NMDS", distance = "bray")

# Extract ordination scores
ord_scores <- as.data.frame(scores(spongeshawaii_ord, display = "sites"))
ord_scores$Sample <- rownames(ord_scores)
ord_scores <- left_join(ord_scores, metadata, by = c("Sample" = "Sample"))

# Function to determine points outside the ellipses
find_outside_points <- function(data, region_col, coord_cols) {
  outside_points <- data %>%
    group_by(!!sym(region_col)) %>%
    do({
      region_data <- .
      cov_matrix <- cov(region_data[, coord_cols])
      center <- colMeans(region_data[, coord_cols])
      mahalanobis_distances <- mahalanobis(region_data[, coord_cols], center, cov_matrix)
      critical_value <- qchisq(0.95, df = ncol(region_data[, coord_cols]))
      region_data[mahalanobis_distances > critical_value, ]
    })
  return(outside_points)
}

# Find points outside the ellipses
outside_points <- find_outside_points(ord_scores, "Region", c("NMDS1", "NMDS2"))
outside_points

# Manually inspect and include additional points that are outside the ellipses
additional_outside_points <- ord_scores %>%
  filter((Region == "MHI" & NMDS1 < -0.67) | (Region == "NWHI" & NMDS1 > 0.96))
additional_outside_points

# Combine the automatically identified and manually inspected points
all_outside_points <- bind_rows(outside_points, additional_outside_points) %>%
  distinct()

# Plot the ordination with ellipses and labels for points outside ellipses
ordination_plot <- plot_ordination(
  spongeshawaii, 
  spongeshawaii_ord, 
  type = "samples", 
  color = "Region", 
  title = "NMDS Ordination (Bray-Curtis)"
) +
  geom_point(size = 3) +
  stat_ellipse(type = "t", level = 0.95, aes(group = Region)) +
  geom_text_repel(
    data = all_outside_points, 
    aes(x = NMDS1, y = NMDS2, label = Sample), 
    size = 4
  ) +
  scale_color_discrete(labels = c(
    "Main Hawaiian Islands" = "MHI",
    "Northwestern Hawaiian Islands" = "NWHI"
  )) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18)
  )

print(ordination_plot)

library(vegan)



## E.1.3. check for statistical differences

# Extract the OTU table and metadata
otu_table_matrix <- as(otu_table(spongeshawaii), "matrix")
metadata_df <- as(sample_data(spongeshawaii), "data.frame")

# Transpose the OTU table to have samples as rows
otu_table_matrix <- t(otu_table_matrix)

# Make sure the sample names match between the OTU table and metadata
matching_samples <- intersect(rownames(otu_table_matrix), rownames(metadata_df))
otu_table_matrix <- otu_table_matrix[matching_samples, ]
metadata_df <- metadata_df[matching_samples, ]

# Perform PERMANOVA
permanova_results <- adonis2(otu_table_matrix ~ Region, data = metadata_df, method = "bray")

# Display PERMANOVA results
print(permanova_results)

