rm(list = ls())

#########################################################
####  A. Set working directory and load packages      ####
#########################################################

setwd("C:/Users/Usuario/Desktop/2025 FINAL PAPER MASTER/0.CSV_DATA_EXCELLS")

library(tidyverse)
library(phyloseq)
library(ggplot2)
library(patchwork)
library(RColorBrewer)

#########################################################
####  B. Load data                                   ####
#########################################################

# Taxonomy
taxonomy <- read.csv("TAXONOMY_TABLE_ORIGINAL.csv", header = TRUE) %>%
  column_to_rownames(var = "ID")

# Abundance
abundance <- read.csv("ABUNDANCE_TABLE_ORIGINAL.csv", header = TRUE) %>%
  column_to_rownames(var = "ID") %>%
  mutate_if(is.integer, as.numeric)

# Metadata
metadata <- read.csv("METADATA_TABLE_ORIGINAL.csv", header = TRUE) %>%
  column_to_rownames(var = "Sample")

#########################################################
####  C. Create phyloseq object                     ####
#########################################################

phylo <- phyloseq(
  otu_table(abundance, taxa_are_rows = TRUE),
  tax_table(as.matrix(taxonomy)),
  sample_data(metadata)
)

#########################################################
####  D. Filter unwanted taxa                       ####
#########################################################

phylo_filtered <- subset_taxa(
  phylo,
  !(kingdom %in% c("Bacillati", "Pseudomonadati", "Fungi")) &
    !(phylum %in% c("Candidatus Binatota", "unclassified"))
)

#########################################################
####  E. Plot A – Relative sequence abundance        ####
#########################################################

taxonomy_region <- merge_samples(phylo_filtered, "Region")
taxonomy_phylum <- tax_glom(taxonomy_region, "phylum", NArm = TRUE)

phylum_relabun_region <- transform_sample_counts(
  taxonomy_phylum,
  function(OTU) OTU / sum(OTU) * 100
)

# Shared color palette
my_colors <- colorRampPalette(brewer.pal(12, "Paired"))(21)

pA <- plot_bar(phylum_relabun_region, fill = "phylum") +
  geom_bar(stat = "identity", position = "stack", color = NA) +
  scale_fill_manual(values = my_colors) +
  scale_x_discrete(labels = c(
    "Northwestern Hawaiian Islands" = "NWHI",
    "Main Hawaiian Islands" = "MHI"
  )) +
  ylab("Relative Abundance (%)") +
  theme(
    axis.text.x  = element_text(size = 14, angle = 0, hjust = 0.5),
    axis.text.y  = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    panel.border = element_blank(),   # removes black rectangle
    legend.position = "none"
  )

#########################################################
####  F. Plot B – Relative abundance of OTUs         ####
#########################################################

# Melt phyloseq
df <- psmelt(phylo_filtered) %>%
  filter(Abundance > 0)

# Count unique OTUs per phylum per region
otu_counts <- df %>%
  group_by(Region, phylum, OTU) %>%
  summarise(Present = 1, .groups = "drop") %>%
  group_by(Region, phylum) %>%
  summarise(NumUniqueOTUs = sum(Present), .groups = "drop")

# Convert to relative abundance
otu_rel_abundance <- otu_counts %>%
  group_by(Region) %>%
  mutate(RelAbundance = NumUniqueOTUs / sum(NumUniqueOTUs) * 100) %>%
  ungroup()

pB <- ggplot(
  otu_rel_abundance,
  aes(x = Region, y = RelAbundance, fill = phylum)
) +
  geom_bar(stat = "identity", position = "stack", color = NA) +
  scale_fill_manual(values = my_colors) +
  scale_x_discrete(labels = c(
    "Northwestern Hawaiian Islands" = "NWHI",
    "Main Hawaiian Islands" = "MHI"
  )) +
  scale_y_continuous(limits = c(0, 100)) +
  ylab("Relative Abundance of OTUs (%)") +
  theme(
    axis.text.x  = element_text(size = 14, angle = 0, hjust = 0.5),
    axis.text.y  = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16)
    # legend kept here
  )

#########################################################
####  G. Combine plots – shared legend & labels     ####
#########################################################

(pA | pB) +
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

