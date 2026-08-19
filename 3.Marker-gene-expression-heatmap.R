# ============================================================================
# Title: Marker gene expression heatmap
# ============================================================================

# ============================================================================
# 1. LOAD LIBRARIES
# ============================================================================

library(ComplexHeatmap)
library(circlize)
library(Seurat)      # Note: Need Seurat V4
library(tidyverse)

# ============================================================================
# 2. LOAD DATA
# ============================================================================

combined <- readRDS('final_combined.RDS')
combined_markers_anno <- read.csv('combined_markers_annotated.csv')

# ============================================================================
# 3. DEFINE CELL TYPES AND LINEAGES
# ============================================================================

# ---- Define 35 cell types ----
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

# ---- Define 7 broad lineages ----
cell_lineages <- c(
  "lineage1",   # lineage1 (e.g., celltype1)
  "lineage1",   # lineage1 (e.g., celltype2)
  "lineage1",   # lineage1
  "lineage2",   # lineage2 (e.g., celltype4)
  "lineage2",   # lineage2
  "lineage2",   # lineage2
  "lineage2",   # lineage2
  "lineage2",   # lineage2
  "lineage3",   # lineage3 (e.g., celltype9)
  "lineage3",   # lineage3
  "lineage3",   # lineage3
  "lineage3",   # lineage3
  "lineage3",   # lineage3
  "lineage3",   # lineage3
  "lineage3",   # lineage3
  "lineage3",   # lineage3
  "lineage3",   # lineage3
  "lineage4",   # lineage4 (e.g., celltype?)
  "lineage4",   # lineage4
  "lineage4",   # lineage4
  "lineage4",   # lineage4
  "lineage4",   # lineage4
  "lineage4",   # lineage4
  "lineage4",   # lineage4
  "lineage4",   # lineage4
  "lineage5",   # lineage5 (e.g., celltype?)
  "lineage5",   # lineage5
  "lineage6",   # lineage6 (e.g., celltype?)
  "lineage7",   # lineage7 (e.g., celltype?)
  "lineage7",   # lineage7
  "lineage7",   # lineage7
  "lineage7",   # lineage7
  "lineage7",   # lineage7
  "lineage7",   # lineage7
  "lineage7"    # lineage7
)

# ============================================================================
# 4. GET TOP MARKER GENES
# ============================================================================

# Get top 10 most significant genes for each cell type
top10 <- combined_markers_anno %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = 1/p_val_adj)

# ============================================================================
# 5. SET CELL TYPE ORDER AND CALCULATE AVERAGE EXPRESSION
# ============================================================================

combined$cluster_anno<- factor(combined$cluster_anno, levels = cell_order)
combined <- SetIdent(combined, value = "cluster_anno")

# Calculate average expression for top genes
avg <- AverageExpression(
  object = combined, 
  group.by = 'cluster_anno', 
  slot = 'data',
  features = top10$gene
)

# ============================================================================
# 6. DEFINE COLORS
# ============================================================================

# ---- Cell type colors (35 colors) ----
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

# ---- Lineage colors (7 lineages) ----
lineage_color_map <- c(
  "lineage1" = "#E0B0FF",   # e.g., Mono/Macro
  "lineage2" = "#BDB5D5",   # e.g., CAF
  "lineage3" = "#A7C7E7",   # e.g., cDC
  "lineage4" = "#AFE1AF",   # e.g., Neutrophil
  "lineage5" = "#F28C28",   # e.g., NK
  "lineage6" = "#DD3F4E",   # e.g., T_cell
  "lineage7" = "#FFB6C1"    # e.g., Mast
)

# Assign lineage colors to each cell type
cell_lineage_colors <- lineage_color_map[cell_lineages]

# ---- Heatmap color function ----
col_fun <- colorRamp2(c(-4, 0, 4), c("blue", "white", "red"))

# ============================================================================
# 7. DEFINE GENE SETS AND COLORS
# ============================================================================

# ---- Define marker genes to highlight (27 genes) ----
genes_show <- c(
  "gene1",   # gene1
  "gene2",   # gene2
  "gene3",   # gene3
  "gene4",   # gene4
  "gene5",   # gene5
  "gene6",   # gene6
  "gene7",   # gene7
  "gene8",   # gene8
  "gene9",   # gene9
  "gene10",  # gene10
  "gene11",  # gene11
  "gene12",  # gene12
  "gene13",  # gene13
  "gene14",  # gene14
  "gene15",  # gene15
  "gene16",  # gene16
  "gene17",  # gene17
  "gene18",  # gene18
  "gene19",  # gene19
  "gene20",  # gene20
  "gene21",  # gene21
  "gene22",  # gene22
  "gene23",  # gene23
  "gene24",  # gene24
  "gene25",  # gene25
  "gene26",  # gene26
  "gene27"   # gene27
)

# ---- Assign colors to genes (matching lineage colors) ----
gene_colors <- c(
  "gene1"  = "#E0B0FF",   # lineage1
  "gene2"  = "#E0B0FF",   # lineage1
  "gene3"  = "#AFE1AF",   # lineage4
  "gene4"  = "#AFE1AF",   # lineage4
  "gene5"  = "#AFE1AF",   # lineage4
  "gene6"  = "#AFE1AF",   # lineage4
  "gene7"  = "#AFE1AF",   # lineage4
  "gene8"  = "#AFE1AF",   # lineage4
  "gene9"  = "#A7C7E7",   # lineage3
  "gene10" = "#A7C7E7",   # lineage3
  "gene11" = "#A7C7E7",   # lineage3
  "gene12" = "#A7C7E7",   # lineage3
  "gene13" = "#A7C7E7",   # lineage3
  "gene14" = "#A7C7E7",   # lineage3
  "gene15" = "#BDB5D5",   # lineage2
  "gene16" = "#BDB5D5",   # lineage2
  "gene17" = "#BDB5D5",   # lineage2
  "gene18" = "#BDB5D5",   # lineage2
  "gene19" = "#FFB6C1",   # lineage7
  "gene20" = "#FFB6C1",   # lineage7
  "gene21" = "#F28C28",   # lineage5
  "gene22" = "#F28C28",   # lineage5
  "gene23" = "#F28C28",   # lineage5
  "gene24" = "#DD3F4E",   # lineage6
  "gene25" = "#DD3F4E",   # lineage6
  "gene26" = "#A7C7E7",   # lineage3
  "gene27" = "#BDB5D5"    # lineage2
)

# ============================================================================
# 8. PREPARE HEATMAP ANNOTATIONS
# ============================================================================

# ---- Get cell counts for each cell type ----
cell_counts <- as.numeric(table(factor(combined$cluster_anno_l2, levels = cell_order)))

# ---- Define alignment for lineage blocks ----
# Each lineage maps to specific cell type indices
align_to <- list(
  "lineage1" = 1:3,     # celltype1-3
  "lineage2" = 4:8,     # celltype4-8
  "lineage3" = 9:17,    # celltype9-17
  "lineage4" = 18:25,   # celltype18-25
  "lineage5" = 26:27,   # celltype26-27
  "lineage6" = 28,      # celltype28
  "lineage7" = 29:35    # celltype29-35
)

# ---- Panel function for lineage annotation ----
panel_fun <- function(index, nm) {
  grid.rect(gp = gpar(fill = cell_lineage_colors[index], col = 'black'))
  txt <- cell_lineages[index]
  grid.text(txt, 0.5, 0.5, gp = gpar(fontface = "bold"))
}

# ---- Create column annotation ----
ha <- HeatmapAnnotation(
  `Cell Lineage` = anno_block(
    align_to = align_to,
    panel_fun = panel_fun,
    labels_gp = gpar(col = "black", fontsize = 10),
    show_name = TRUE
  ),
  `Cell Type` = anno_simple(
    cell_order, 
    col = cell_colors,
    border = TRUE
  ),
  annotation_name_gp = gpar(fontsize = 14), 
  which = 'col', 
  show_legend = FALSE,
  border = TRUE
)

# ============================================================================
# 9. GENERATE HEATMAP
# ============================================================================

# Convert average expression to data frame
avg <- as.data.frame(avg$RNA)

# Create PDF output
pdf("Heatmap_Marker_Genes.pdf", height = 8, width = 10)

# Main heatmap
Heatmap(
  t(scale(t(avg))), 
  name = "mat", 
  rect_gp = gpar(col = "transparent", lwd = 0), 
  col = col_fun, 
  column_title = "scRNA Average Scaled Cell Type Marker Expression", 
  clustering_method_rows = "single", 
  row_dend_gp = gpar(lwd = 0.2),
  cluster_columns = FALSE,
  show_row_names = FALSE, 
  show_heatmap_legend = FALSE,
  column_names_gp = gpar(fontsize = 15), 
  column_title_gp = gpar(fontsize = 18, fontface = "bold"),
  top_annotation = ha
) +
  rowAnnotation(
    link = anno_mark(
      at = which(rownames(avg) %in% genes_show),
      labels = rownames(avg)[rownames(avg) %in% genes_show], 
      labels_gp = gpar(col = gene_colors, fontsize = 15), 
      which = "rows"
    )
  )

# Add legend
lgd <- Legend(
  col_fun = col_fun,
  title = "Normalized\nExpression", 
  direction = "horizontal", 
  title_position = "topcenter", 
  border = FALSE, 
  grid_width = unit(1, "cm")
)
draw(lgd, x = unit(0.94, "npc"), y = unit(0.15, "npc"))

dev.off()

# ============================================================================
# 10. SESSION INFORMATION
# ============================================================================

sessionInfo()

# ============================================================================
# END OF SCRIPT
# ============================================================================