# ============================================================================
# Title: UMAP + Cell Fraction Bar Plot
# ============================================================================

# ============================================================================
# 1. INSTALL AND LOAD LIBRARIES
# ============================================================================

# remotes::install_version("SeuratObject", "4.1.4", 
#                          repos = c("https://satijalab.r-universe.dev", 
#                                    getOption("repos")))
# remotes::install_version("Seurat", "4.4.0", 
#                          repos = c("https://satijalab.r-universe.dev", 
#                                    getOption("repos")))

library(Seurat)      # Note: Need Seurat V4
library(tidyverse)
library(patchwork)
library(RColorBrewer)

# ============================================================================
# 2. LOAD DATA
# ============================================================================

combined <- readRDS('final_combined.RDS')

# ============================================================================
# 3. DEFINE CELL TYPE ORDER AND COLORS
# ============================================================================

# ---- Define cell type order (35 cell types) ----
cell_order <- c(
  "celltype1",   # celltype1
  "celltype2",   # celltype2
  "celltype3",   # celltype3
  "celltype4",   # celltype4
  "celltype5",   # celltype5
  "celltype6",   # celltype6
  "celltype7",   # celltype7
  "celltype8",   # celltype8
  "celltype9",   # celltype9
  "celltype10",  # celltype10
  "celltype11",  # celltype11
  "celltype12",  # celltype12
  "celltype13",  # celltype13
  "celltype14",  # celltype14
  "celltype15",  # celltype15
  "celltype16",  # celltype16
  "celltype17",  # celltype17
  "celltype18",  # celltype18
  "celltype19",  # celltype19
  "celltype20",  # celltype20
  "celltype21",  # celltype21
  "celltype22",  # celltype22
  "celltype23",  # celltype23
  "celltype24",  # celltype24
  "celltype25",  # celltype25
  "celltype26",  # celltype26
  "celltype27",  # celltype27
  "celltype28",  # celltype28
  "celltype29",  # celltype29
  "celltype30",  # celltype30
  "celltype31",  # celltype31
  "celltype32",  # celltype32
  "celltype33",  # celltype33
  "celltype34",  # celltype34
  "celltype35"   # celltype35
)

# ---- Define color palette (35 colors) ----
cell_colors <- c(
  "#CF9FFF", "#E7C7DC", "#CAA7DD", "#953553", "#a15891", 
  "#9C58A1", "#79127F", "#BF40BF", "#A8A2D2", "#B6D0E2", 
  "#2874A6", "#5599C8", "#6495ED", "#96C5D7", "#B0DFE6", 
  "#4682B4", "#40B5AD", "#8FCACA", "#C6DBDA", "#CCE2CB", 
  "#97C1A9", "#63BA97", "#019477", "#7DB954", "#64864A", 
  "#FFBF00", "#ff9d5c", "#DD3F4E", "#FFC8A2", "#FFD580", 
  "#F89880", "#FA8072", "#F3B0C3", "#FED7C3", "#FFC5BF"
)
names(cell_colors) <- cell_order

# ---- Define lineage level 1 (7 broad lineages) ----
lineage_order <- c(
  "lineage1",   # lineage1
  "lineage2",   # lineage2
  "lineage3",   # lineage3
  "lineage4",   # lineage4
  "lineage5",   # lineage5
  "lineage6",   # lineage6
  "lineage7"    # lineage7
)

lineage_colors <- c("#E0B0FF", "#A7C7E7", "#AFE1AF", "#BDB5D5", 
                    "#FFB6C1", "#F28C28", "#DD3F4E")
names(lineage_colors) <- lineage_order

# ============================================================================
# 4. FIGURE 1B - UMAP Showing All Cell Types
# ============================================================================

p1 <- DimPlot(combined,
              group.by = 'celltype',
              raster = TRUE, 
              pt.size = 2,
              raster.dpi = c(1028, 1028),
              label = TRUE, 
              repel = TRUE, 
              reduction = "UMAP_dim30", 
              cols = cell_colors) + 
  coord_fixed() + 
  ggtitle("Comprehensive scRNA-Seq Atlas\nn=XX,XXX cells") +
  SeuratExtend::theme_umap_arrows(
    x_label = "UMAP1",
    y_label = "UMAP2", 
    text_offset_x = unit(2, 'mm'),
    text_offset_y = unit(2, 'mm'),
    text_size = 15
  ) +
  NoLegend()

# Display and save
p1
ggsave(p1, file = "Fig6SA-scRNA_UMAP.pdf", width = 6, height = 6)

# ============================================================================
# 5. FIGURE 1C - Bar Chart Showing Lineage Frequencies
# ============================================================================

# ---- Figure 1C (left) - Barplot of cell counts per lineage ----
p2 <- combined@meta.data %>% 
  ggplot(aes(y = forcats::fct_rev(forcats::fct_infreq(cluster_anno)), 
             fill = cluster_anno)) + 
  geom_bar(stat = 'count') +
  labs(x = 'Cell count', y = NULL) +
  scale_fill_manual(
    name = "Cell lineage", 
    values = lineage_colors, 
    labels = lineage_order
  ) + 
  theme_bw(base_size = 16) + 
  theme(axis.text = element_text(size = 16, color = 'black')) 

# ---- Figure 1C (right) - Sample composition barplot ----
cell_counts <- as.data.frame(table(combined$cluster_anno, combined$orig.ident))

p3 <- ggplot(data = cell_counts, 
             aes(x = forcats::fct_rev(Var2), y = Freq, fill = Var1)) +
  geom_bar(position = "fill", stat = "identity") + 
  coord_flip() +
  labs(x = NULL, y = 'Cell Lineage Frequency') + 
  scale_x_discrete(labels = c("Sample1", "Sample2", "Sample3", "Sample4", 
                              "Sample5", "Sample6", "Sample7", "Sample8")) +
  scale_fill_manual(
    name = "Cell lineage", 
    values = lineage_colors
  ) + 
  theme_bw(base_size = 16) + 
  theme(axis.text = element_text(size = 16, color = 'black')) 

# ---- Combine plots ----
fig1c <- p2 + p3 + 
  plot_layout(guides = "collect") &
  plot_annotation(
    title = "Atlas Composition",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'))
  )

# Display and save
fig1c
ggsave("Fig6SB_Barplots_Lineage_Composition.pdf", fig1c, width = 9, height = 4)

# ============================================================================
# 6. SESSION INFORMATION
# ============================================================================

sessionInfo()

# ============================================================================
# END OF SCRIPT
# ============================================================================