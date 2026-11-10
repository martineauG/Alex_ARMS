#######################################
####  D. Relabel REGION / ORDER  #######
#######################################

library(dplyr)
library(ggplot2)
library(patchwork)
library(phyloseq)

## -------------------------------
## Glom by ORDER and rel. abundance
## -------------------------------
spongesbyorder <- tax_glom(spongeshawaii_fraction_Region, "order", NArm = TRUE)

spongesbyorder_relabun_region <- transform_sample_counts(
  spongesbyorder,
  function(OTU) OTU / sum(OTU) * 100
)

## Melt to dataframe
df <- psmelt(spongesbyorder_relabun_region)

## -------------------------------
## Function to add "Others"
## -------------------------------
add_others_category <- function(df_region, top_orders) {
  df_region %>%
    mutate(order = ifelse(order %in% top_orders, order, "Others")) %>%
    group_by(order) %>%
    summarise(Abundance = sum(Abundance), .groups = "drop") %>%
    mutate(Region = unique(df_region$Sample))
}

## ===============================
## ==========  MHI  ==============
## ===============================

top6_mhi <- c(
  "Poecilosclerida", "Haplosclerida", "Tetractinellida",
  "Dictyoceratida", "Suberitida", "Tethyida"
)

df_mhi <- df %>% filter(Sample == "Main Hawaiian Islands")
df_mhi_combined <- add_others_category(df_mhi, top6_mhi)

my_colors_mhi <- c(
  "Poecilosclerida" = "#ffbb78",
  "Haplosclerida"   = "#dbdb8d",
  "Tetractinellida" = "#f7b6d2",
  "Dictyoceratida"  = "#666666",
  "Suberitida"      = "#d2b4de",
  "Tethyida"        = "#d8a29d",
  "Others"          = "#c0c0c0"
)

stack_order_mhi <- c(
  "Others", "Tethyida", "Suberitida",
  "Dictyoceratida", "Haplosclerida",
  "Tetractinellida", "Poecilosclerida"
)

df_mhi_combined$order <- factor(df_mhi_combined$order, levels = stack_order_mhi)

p_mhi <- ggplot(df_mhi_combined, aes(x = "MHI", y = Abundance, fill = order)) +
  geom_bar(stat = "identity", color = NA) +
  scale_fill_manual(values = my_colors_mhi, name = "Order") +
  coord_cartesian(ylim = c(0, 100)) +
  ylab("Relative Abundance (%)") +
  theme_grey(base_size = 16) +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_line(color = "white", linewidth = 0.6),
    panel.grid.minor = element_blank(),
    axis.title.y = element_text(size = 20),
    axis.text.y  = element_text(size = 16),
    axis.text.x  = element_text(size = 20),
    axis.title.x = element_blank(),
    legend.position = "none"
  )

## ===============================
## ==========  NWHI  =============
## ===============================

top6_nwhi <- c(
  "Poecilosclerida", "Haplosclerida", "Tetractinellida",
  "Scopalinida", "Clathrinida", "Tethyida"
)

df_nwhi <- df %>% filter(Sample == "Northwestern Hawaiian Islands")
df_nwhi_combined <- add_others_category(df_nwhi, top6_nwhi)

my_colors_nwhi <- c(
  "Poecilosclerida" = "#ffbb78",
  "Haplosclerida"   = "#dbdb8d",
  "Tetractinellida" = "#f7b6d2",
  "Scopalinida"     = "#ea8983",
  "Clathrinida"     = "#9e77b7",
  "Tethyida"        = "#d8a29d",
  "Others"          = "#c0c0c0"
)

stack_order_nwhi <- c(
  "Others", "Tethyida", "Clathrinida",
  "Scopalinida", "Haplosclerida",
  "Tetractinellida", "Poecilosclerida"
)

df_nwhi_combined$order <- factor(df_nwhi_combined$order, levels = stack_order_nwhi)

p_nwhi <- ggplot(df_nwhi_combined, aes(x = "NWHI", y = Abundance, fill = order)) +
  geom_bar(stat = "identity", color = NA) +
  scale_fill_manual(values = my_colors_nwhi, name = "Order") +
  coord_cartesian(ylim = c(0, 100)) +
  theme_grey(base_size = 16) +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_line(color = "white", linewidth = 0.6),
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_blank(),
    legend.position = "none"
  )

## ===============================
## ========  COMBINE  ============
## ===============================

order_plot <- (p_mhi | p_nwhi) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "right",
    legend.title = element_text(size = 20),
    legend.text  = element_text(size = 16)
  )

## Display final figure
order_plot



