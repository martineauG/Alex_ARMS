rm(list = ls())
#########################################################
####  A. Set working directory and download packages ####
#########################################################

## A.1 Working directory
setwd("C:/Users/Usuario/Desktop/2025 FINAL PAPER MASTER/0. CSV_DATA_EXCELLS")

## A.2 Load packages
library(tidyverse)
library(patchwork)
library(phyloseq)
library(ggpubr)
library(ggplot2)

## B.1 Taxonomy
taxonomy <- read.csv("TAXONOMY_TABLE_ORIGINAL.csv", header = TRUE) %>%
  column_to_rownames(var = "ID")

## B.2 Abundance
abundance <- read.csv("ABUNDANCE_TABLE_ORIGINAL.csv", header = TRUE) %>%
  column_to_rownames(var = "ID") %>%
  mutate_if(is.integer, as.numeric)

combined_table <- merge(abundance, taxonomy, by = "row.names") %>%
  filter(phylum == "Porifera") %>%
  column_to_rownames(var = "Row.names")

abundance.only.porifera <- combined_table %>%
  dplyr::select_if(is.numeric)

taxonomy.only.porfiera <- combined_table %>%
  dplyr::select_if(is.character)

## B.3 Metadata
metadata <- read.csv("METADATA_TABLE_ORIGINAL.csv", header = TRUE)

metadata1 <- metadata %>%
  column_to_rownames(var = "Sample")

## C.1 Create phyloseq object
phylo <- phyloseq(
  otu_table(abundance, taxa_are_rows = TRUE),
  tax_table(as.matrix(taxonomy)),
  sample_data(metadata1)
)

## Select only Porifera
spongeshawaii <- subset_taxa(phylo, phylum == "Porifera")

#######################################
####  D. Relabel REGION / CLASS  #######
#######################################

my_colors_class <- c("#d62728", "#98df8a", "#1f77b4")

## Merge samples by Region
spongeshawaii_fraction_Region <- merge_samples(spongeshawaii, "Region")

## Glom by class
spongesbyclass <- tax_glom(
  spongeshawaii_fraction_Region,
  "class",
  NArm = TRUE
)

## Convert to relative abundance (%)
spongesbyclass_relabun_region <- transform_sample_counts(
  spongesbyclass,
  function(OTU) OTU / sum(OTU) * 100
)

## =========================
## Panel A: Relative abundance by class
## =========================
pA <- plot_bar(spongesbyclass_relabun_region, fill = "class") +
  geom_bar(stat = "identity", position = "stack", color = NA) +
  scale_fill_manual(values = my_colors_class) +
  coord_cartesian(ylim = c(0, 100)) +
  scale_x_discrete(labels = c(
    "Northwestern Hawaiian Islands" = "NWHI",
    "Main Hawaiian Islands" = "MHI"
  )) +
  ylab("Relative Abundance (%)") +
  theme(
    panel.border = element_blank(),
    axis.text.x  = element_text(size = 16, hjust = 0.5, angle = 0),  # ← FIX
    axis.text.y  = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.position = "none"
  )


## =========================
## Panel B: Relative abundance of OTUs by class
## =========================

## Melt phyloseq object
df <- psmelt(spongeshawaii) %>%
  filter(Abundance > 0)

## Count unique OTUs per class per region
otu_counts_class <- df %>%
  group_by(Region, class, OTU) %>%
  summarise(Present = 1, .groups = "drop") %>%
  group_by(Region, class) %>%
  summarise(NumUniqueOTUs = sum(Present), .groups = "drop")

## Convert to relative abundance (%)
otu_rel_abundance_class <- otu_counts_class %>%
  group_by(Region) %>%
  mutate(RelAbundance = NumUniqueOTUs / sum(NumUniqueOTUs) * 100) %>%
  ungroup()

pB <- ggplot(
  otu_rel_abundance_class,
  aes(x = Region, y = RelAbundance, fill = class)
) +
  geom_bar(stat = "identity", position = "stack", color = NA) +
  scale_fill_manual(values = my_colors_class) +
  coord_cartesian(ylim = c(0, 100)) +
  scale_x_discrete(labels = c(
    "Northwestern Hawaiian Islands" = "NWHI",
    "Main Hawaiian Islands" = "MHI"
  )) +
  ylab("Relative Abundance of OTUs (%)") +
  theme(
    axis.text.x  = element_text(size = 16, hjust = 0.5),
    axis.text.y  = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  )

## =========================
## Combine panels
## =========================
combined_plot <- (pA | pB) +
  plot_layout(guides = "collect") +
  plot_annotation(
    tag_levels = "A",
    tag_suffix = ")",
    theme = theme(
      plot.tag = element_text(size = 18, face = "bold"),
      legend.position = "right",
      legend.title = element_text(size = 14),
      legend.text  = element_text(size = 12)
    )
  )

## Display figure
combined_plot
