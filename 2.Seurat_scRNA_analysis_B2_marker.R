
###方法1：分面小提琴图（最直观）####
# 创建标记基因列表
markers_list_mouse<- list(
  # T细胞
  Tcell = c("Cd3e", "Cd3d", "Cd3g", "Cd2", "Cd5", "Cd28", "Cd4", "Cd8a", "Cd8b1", 
            "Lat", "Lck", "Zap70", "Il2rb", "Il7r"),
  
  # NK细胞
  NKcell = c("Klrb1c", "Klrd1", "Nkg7", "Xcl1", "Klra8", "Klre1", "Gzma", "Gzmb", 
             "Prf1", "Ncr1", "Cd49b", "Cd69", "Il2rb", "Fcer1g"),
  
  # B细胞
  Bcell = c("Cd79a", "Iglc2", "Mzb1", "Jchain", "Cd19", "Ms4a1", "Cd22", "Cd24a",
            "Ighm", "Ighd", "Pax5", "Bank1", "Cd37"),
  
  # 单核/巨噬细胞
  MonoMacro = c("Fcer1g", "Tyrobp", "S100a11", "Fn1", "Cd68", "Ly6c2", "Adgre1", 
                "C1qa", "C1qb", "C1qc", "Cd14", "Itgam", "Csf1r", "Aif1"),
  
  # 中性粒细胞
  Neutrophil = c("Clec4d", "Mmp9", "S100a8", "Camp", "S100a9", "S100a12", "Cxcr2",
                 "Csf3r", "Fpr1", "Fpr2", "Mpo", "Elane", "Lcn2", "Mmp8", "Pglyrp1"),
  
  # 嗜碱性粒细胞
  Basophil = c("Cpa3", "Gata2", "Ms4a2", "Mcpt8", "Fcer1a", "Fcer1g", "Hdc",
               "Cd200r3", "Il1rl1"),
  
  # 常规树突状细胞 (cDC)
  cDC = c("Clec9a", "Clec10a", "Cd209a", "Flt3", "Itgax", "Cd24a", "Xcr1", 
          "Zbtb46", "Cd8a", "Sirpa"),
  
  # 迁移性树突状细胞 (Migratory cDC)
  MigratorycDC = c("Ccr7", "Ccl22", "Fscn1", "Cd40", "Cd80", "Cd86", "Il12b",
                   "Tnf", "Il6", "Ccl5", "Ebi3"),
  
  # 浆细胞样树突状细胞 (pDC)
  pDC = c("Siglech", "Ccr9", "Upb1", "Tcf4", "Bst2", "Il3ra", "Clec4c", 
          "Ptprc", "Pkn1"),
  
  # 肿瘤相关成纤维细胞 (CAF)
  CAF = c("Ly6c1", "Col3a1", "Dcn", "Col1a1", "Col1a2", "Col5a2", "Fap", 
          "Pdpn", "Acta2", "Tagln", "Vim", "Pdgfra", "Pdgfrb", "Cxcl12"),
  
  # 浆细胞 (Plasma Cell) 
  Plasma = c("Mzb1", "Jchain", "Prdm1", "Xbp1", "Ighg1", "Igha", "Sdc1", 
             "Cd38", "Tnfrsf17", "Slpi"),
  
  # 经典单核细胞 (Classical Monocyte) 
  MonoClassical = c("Ly6c2", "Ccr2", "Cd14", "Fcgr1", "S100a9", "S100a8",
                    "Csf1r", "Il1b", "Tnf", "Nfil3"),
  
  # 非经典单核细胞 (Non-classical Monocyte) 
  MonoNonClassical = c("Ace", "Cx3cr1", "Itga4", "Ms4a7", "Fcgr4",
                       "Lyz2", "Treml4", "Cd43"),
  
  # 巨噬细胞 (Macrophage) - 新增亚群
  Macrophage = c("Adgre1", "C1qa", "C1qb", "Trem2", "Cd68", "Mrc1", "Il10",
                 "Tgfbi", "Thbs1", "Lgals3", "Arg1", "Csf1r")
)

# 转换为向量方便绘图
all_markers <- unlist(markers_list_mouse)

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



markers_list_Tcell <- list(
  # ===== 泛T细胞标记 =====
  # 核心：T细胞受体复合物、共受体与信号分子
  Tcell = c("Cd3e", "Cd3d", "Cd3g", "Cd2", "Cd5", "Lat", "Lck", "Zap70"),
  
  # ===== 主要功能性T细胞亚群 =====
  # CD8+ T细胞：细胞毒性T细胞，需同时表达Cd8a与Cd8b1才能形成功能性共受体
  CD8Tcell = c("Cd8a", "Cd8b1", "Gzma", "Gzmb", "Gzmk", "Prf1", "Ifng", "Tnf", "Fasl"),
  
  # CD4+ T细胞：辅助性T细胞
  CD4Tcell = c("Cd4", "Il7r"),
  
  # ===== 初始与记忆T细胞亚群（Nature Reviews Immunology 2025定义） =====
  # 初始T细胞（Naive T）
  Tnaive = c("Sell", "Ccr7", "Il7r", "Tcf7", "Lef1", "S1pr1", "Klf2"),
  
  # 中央记忆T细胞 (TCM) —— CD44hi, CD62L+, CCR7+, TCF1+ 
  Tcm = c("Sell", "Ccr7", "Cd44", "Tcf7", "Bcl2"),
  
  # 效应记忆T细胞 (TEM) —— CD44hi, CD62L-, CCR7- 
  # 注：Sell、Ccr7、Tcf7在TEM中均为低表达或不表达，不应作为TEM标记
  Tem = c("Cd44", "Ccl5", "Cxcr3", "Gzma", "Gzmb"),
  
  # 终末分化效应记忆T细胞 (TEMRA) —— 人TEMRA特征为CCR7-CD45RA+，小鼠中此群体较少且具异质性 
  Temra = c("Cx3cr1", "Klrg1", "Gzmb", "Prf1"),
  
  # 干细胞样记忆T细胞 (TSCM) —— 具自我更新能力，Tcf7、Lef1高表达 
  Tscm = c("Sell", "Ccr7", "Cd44", "Tcf7", "Lef1", "Bcl2", "Il7r"),
  
  # 组织驻留记忆T细胞 (TRM) —— CD69+, CD103(Itgae)+, CXCR6+ 
  Trm = c("Cd69", "Itgae", "Cxcr6", "Cd49a", "Zfp683", "Bhlhe40", "Klf2"),
  
  # ===== 调节性T细胞 (Treg) =====
  # 核心：Foxp3+ 为金标准；表面标记使用CD25(Il2ra)+CD127(Il7r)low为优选 
  Treg = c("Foxp3", "Il2ra", "Il7r", "Ctla4", "Tigit", "Tnfrsf18", "Entpd1"),
  
  # ===== CD4+ T辅助细胞亚群（效应性） =====
  # 滤泡辅助性T细胞 (TFH)
  Tfh = c("Cxcr5", "Pdcd1", "Bcl6", "Icos", "Il21", "Slamf6"),
  
  # Th1细胞 —— 特征性转录因子Tbx21（T-bet），效应分子Ifng
  Th1 = c("Tbx21", "Ifng", "Stat4", "Cxcr3", "Il12rb2"),
  
  # Th2细胞 —— 特征性转录因子Gata3，效应分子Il4、Il13
  Th2 = c("Gata3", "Il4", "Il13", "Stat6"),
  
  # Th17细胞 —— 特征性转录因子Rorc，效应分子Il17a
  Th17 = c("Rorc", "Il17a", "Il17f", "Il22", "Ccr6"),
  
  # Th9细胞
  Th9 = c("Irf4", "Il9", "Il10"),
  
  # Th22细胞
  Th22 = c("Ahr", "Il22"),
  
  # ===== 耗竭T细胞 (TEX)与耗竭前体 (TPEX) =====
  # 耗竭T细胞 (TEX) —— 高表达多个抑制性受体，缺乏Tcf7 
  Tex = c("Pdcd1", "Tigit", "Ctla4", "Lag3", "Havcr2", "Tox", "Eomes", "Cd244", "Entpd1"),
  
  # 耗竭前体T细胞 (TPEX) —— 表达Tcf7，具自我更新能力，是对免疫检查点抑制剂应答的关键群体 
  Tpex = c("Tcf7", "Sell", "Pdcd1", "Tox"),
  
  # ===== 非常规T细胞亚群 =====
  # γδ T细胞 —— 使用Trdc和Trgc标识TCRγδ 
  Tgd = c("Trdc", "Trgc", "Cd27", "Cd44"),
  
  # NKT细胞（自然杀伤T细胞）
  NKT = c("Zbtb16", "Cd1d1", "Klrb1c", "Ifng", "Il4"),
  
  # MAIT细胞（粘膜相关恒定T细胞）
  MAIT = c("Rorc", "Cxcr6", "Il7r", "Klrd1", "Ifng", "Il17a")
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

