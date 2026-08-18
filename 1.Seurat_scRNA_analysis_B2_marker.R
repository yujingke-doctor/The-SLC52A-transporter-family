
###方法1：分面小提琴图（最直观）####
# 创建标记基因列表
markers_list <- list(
  Tcell = c("Cd3e", "Cd3d"),
  NKcell = c("Klrb1c", "Klrd1", "Nkg7", "Xcl1"),
  Bcell = c("Cd79a", "Iglc2", "Mzb1", "Jchain"),
  MonoMacro = c("Fcer1g", "Tyrobp", "S100a11", "Fn1"),
  Neutrophil = c("Clec4d", "Mmp9", "S100a8", "Camp"),
  Basophil = c("Cpa3", "Gata2", "Ms4a2", "Mcpt8"),
  cDC = c("Clec9a", "Clec10a", "Cd209a"),
  MigratorycDC = c("Ccr7", "Ccl22", "Fscn1"),
  pDC = c("Siglech", "Ccr9", "Upb1"),
  Fibroblast = c("Ly6c1", "Col3a1", "Dcn", "Col1a1")
)

# 转换为向量方便绘图
all_markers <- unlist(markers_list)

# 分面小提琴图（按细胞类型分组）
VlnPlot(combined, features = all_markers, 
        group.by = "seurat_clusters",
        stack = TRUE, flip = TRUE,
        ncol = 4) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 分面小提琴图（按细胞类型分组）
library(RColorBrewer)
set3_colors <- brewer.pal(min(12, length(unique(combined$seurat_clusters))), "Set3")

zm_colors <- colorRampPalette(brewer.pal(12, "Set3"))(35)

# 应用并保存
VlnPlot(combined, features = all_markers, 
        group.by = "seurat_clusters",
        stack = TRUE, #小提琴堆叠
        flip = TRUE,#小提琴翻转
        cols = zm_colors) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



markers_list_mouse <- list(
  # T细胞系列
  Tcell = c("Cd3e", "Cd3d", "Cd3g", "Cd8a"),
  
  CD8Tcell = c("Cd8a", "Cd8b1", "Gzma", "Gzmb", "Gzmk", "Prf1", 
               "Ifng"),  # CXC13 → Cxcl13
  
  CD4Tcell = c("Cd4", "Il7r", "Ccr7", "Sell"),
  
  Treg = c("Foxp3", "Il2ra", "Ctla4", "Tigit"),
  
  Trm = c("Cd69", "Itgae", "Cxcr6"),
  
  Tem <- c("Ccl5","Ccr5","Sell", "Il7r", "Ccr7", "Tcf7"),
  # NK细胞
  NKcell = c("Klrb1c", "Klrd1", "Nkg7", "Xcl1"),
  
  # B细胞系列
  Bcell = c("Cd79a", "Cd79b", "Ms4a1"),  
  
  plasma = c("Mzb1", "Jchain", "Igha", "Ighg1", "Sdc1"),  # IGHA1 → Igha, IGHG1 → Ighg1
  
  # 髓系细胞
  MonoMacro = c("Fcer1g", "Tyrobp", "S100a11", "Fn1"),
  
  Neutrophil = c("Clec4d", "Mmp9", "S100a8", "Camp"),
  
  Basophil = c("Cpa3", "Gata2", "Ms4a2", "Mcpt8"),
  
  # 树突状细胞
  cDC = c("Clec9a", "Clec10a", "Cd209a","Cd74","Cd83"),
  
  MigratorycDC = c( "Ccl22", "Fscn1"),
  
  pDC = c("Siglech", "Ccr9", "Upp1"),  # Upb1 → Upp1
  
  # 基质细胞
  Fibroblast = c("Ly6c1", "Col3a1", "Dcn", "Col1a1")
)


###方法2：点图（DotPlot）- 推荐####
DotPlot(combined, features = markers_list_mouse,
        group.by = "RNA_snn_res.0.5",
        cols = c("lightblue", "red"),
        dot.scale = 6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("细胞类型标记基因表达模式")


###方法3：热图（Heatmap）####
DoHeatmap(combined, features = all_markers,
          group.by = "seurat_clusters",
          group.colors = rainbow(length(unique(combined$seurat_clusters)))) +
  theme(axis.text.y = element_text(size = 8))


###方法4：分cluster的小提琴图####
# 按cluster查看标记基因表达
for (cluster in sort(unique(combined$seurat_clusters))) {
  print(
    VlnPlot(combined, features = all_markers,
            idents = cluster) +
      ggtitle(paste("Cluster", cluster))
  )
}

###方法5：特征图（FeaturePlot）分面显示####
# 在UMAP/tSNE上查看标记基因表达
FeaturePlot(combined, features = all_markers,
            reduction = "umap_harmony",
            ncol = 4, order = TRUE,
            min.cutoff = 'q10', max.cutoff = 'q90')


###方法6：自动评分系统####
# 为每个细胞类型计算得分
library(Seurat)

# 计算每个细胞类型的模块得分
for(cell_type in names(markers_list)) {
  combined <- AddModuleScore(combined,
                             features = list(markers_list[[cell_type]]),
                             name = cell_type,
                             ctrl = 100)  # 使用100个控制基因
}

# 查看每个cluster的得分
score_plot_data <- FetchData(combined, 
                             vars = c(paste0(names(markers_list), "1"), "seurat_clusters"))

# 汇总每个cluster的平均得分
cluster_scores <- aggregate(. ~ seurat_clusters, data = score_plot_data, mean)
rownames(cluster_scores) <- cluster_scores$seurat_clusters
cluster_scores$seurat_clusters <- NULL

# 热图显示得分
pheatmap::pheatmap(t(cluster_scores),
                   color = viridis::viridis(100),
                   main = "细胞类型得分热图",
                   cluster_rows = TRUE,
                   cluster_cols = TRUE)

###方法7：快速分类函数####
# 自动识别每个cluster的主要细胞类型
identify_cluster_celltypes <- function(seurat_obj, markers_list, min_frac = 0.2) {
  # 获取表达矩阵
  exp_mat <- GetAssayData(seurat_obj, slot = "data")
  
  # 为每个cluster计算每个标记集的表达比例
  cluster_ids <- as.character(seurat_obj$seurat_clusters)
  clusters <- unique(cluster_ids)
  
  results <- data.frame()
  
  for(cluster in clusters) {
    # 获取该cluster的细胞
    cluster_cells <- colnames(seurat_obj)[cluster_ids == cluster]
    
    cluster_result <- data.frame(cluster = cluster)
    
    for(cell_type in names(markers_list)) {
      markers <- markers_list[[cell_type]]
      markers <- markers[markers %in% rownames(exp_mat)]
      
      if(length(markers) > 0) {
        # 计算该cluster中表达标记的细胞比例
        exp_subset <- exp_mat[markers, cluster_cells, drop = FALSE]
        cell_frac <- mean(colSums(exp_subset > 0) > 0)
        
        # 计算平均表达水平
        avg_exp <- mean(exp_subset[exp_subset > 0])
        if(is.na(avg_exp)) avg_exp <- 0
        
        cluster_result[[paste0(cell_type, "_frac")]] <- cell_frac
        cluster_result[[paste0(cell_type, "_exp")]] <- avg_exp
      }
    }
    
    results <- rbind(results, cluster_result)
  }
  
  # 找出每个cluster的主要细胞类型
  for(cluster in clusters) {
    row_idx <- which(results$cluster == cluster)
    
    # 获取所有frac列
    frac_cols <- grep("_frac$", colnames(results), value = TRUE)
    frac_values <- as.numeric(results[row_idx, frac_cols])
    
    # 如果最高比例大于阈值，则赋值
    if(max(frac_values) > min_frac) {
      main_type <- gsub("_frac", "", frac_cols[which.max(frac_values)])
      results$predicted_type[row_idx] <- main_type
    } else {
      results$predicted_type[row_idx] <- "Unknown"
    }
  }
  
  return(results)
}

# 使用函数
cluster_predictions <- identify_cluster_celltypes(combined, markers_list, min_frac = 0.3)
print(cluster_predictions[, c("cluster", "predicted_type")])

###方法8：交互式可视化####
# 安装并运行CellBrowser
library(cellbrowser)

# 导出数据到CellBrowser
ExportToCellbrowser(combined, 
                    dataset.name = "MyDataset",
                    out.dir = "cellbrowser",
                    meta.fields = c("seurat_clusters"),
                    markers = all_markers)

# 在浏览器中打开
cbBuild("cellbrowser")

####9.快速筛查流程建议：####
# 9.1. 先用点图看整体模式
p1 <- DotPlot(combined, features = unlist(markers_list),
              group.by = "seurat_clusters")
print(p1)

# 9.2. 查看疑似淋巴细胞的cluster（T/NK/B）
immune_clusters <- c()  # 根据点图确定
if(length(immune_clusters) > 0) {
  VlnPlot(combined, features = c("Cd3e", "Cd3d", "Nkg7", "Cd79a"),
          idents = immune_clusters)
}

# 9.3. 查看疑似髓系细胞的cluster
myeloid_clusters <- c()  # 根据点图确定
if(length(myeloid_clusters) > 0) {
  VlnPlot(combined, features = c("Fcer1g", "Tyrobp", "S100a8", "Clec9a"),
          idents = myeloid_clusters)
}

# 9.4. 查看成纤维细胞
VlnPlot(combined, features = c("Col1a1", "Dcn"),
        group.by = "seurat_clusters")

