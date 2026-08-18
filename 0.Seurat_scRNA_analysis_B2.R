# 1.加载必要的库 ----
# Seurat V5已经整合了许多功能，不再需要某些单独的包
library(dplyr)
library(Seurat) # 使用V5版本
library(patchwork)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(forcats)
library(ggrastr)
library(cowplot)
library(SeuratDisk) # 用于h5ad文件转换
library(ComplexHeatmap)
library(harmony)

# 注意：DoubletFinder在Seurat V5中可能需要特别处理
# 如果需要，请先安装：install.packages("remotes"); remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")
library(DoubletFinder)
# 检查Seurat版本
packageVersion("Seurat")

# 2.修改后的FindDoublets函数，适配Seurat V5 ----
FindDoublets <- function(seurat.rna, PCs = 1:50, exp_rate = 0.02, sct = FALSE){
  # sct--是否进行SCTransform
  
  ## pK identification
  sweep.res.list <- paramSweep(seurat.rna, PCs = PCs, sct = sct)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  # 选择最优pK值
  optimal_pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  
  ## Homotypic Doublet proportion Estimate
  # Seurat V5中，聚类结果可能在active.id中
  if (!"seurat_clusters" %in% colnames(seurat.rna@meta.data)) {
    stop("Please run FindClusters first")
  }
  annotations <- seurat.rna@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(exp_rate * length(seurat.rna$seurat_clusters))
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  
  ## Run DoubletFinder
  # 使用最优pK值而不是硬编码的0.09
  seurat.rna <- doubletFinder(
    seurat.rna, 
    PCs = PCs, 
    pN = 0.25,
    pK = optimal_pK, 
    nExp = nExp_poi, 
    reuse.pANN = FALSE, 
    sct = sct
  )
  
  # 获取第一次运行创建的pANN列名
  pANN_name <- colnames(seurat.rna@meta.data)[
    grep("pANN", colnames(seurat.rna@meta.data))
  ][1]
  
  seurat.rna <- doubletFinder(
    seurat.rna, 
    PCs = PCs, 
    pN = 0.25, 
    pK = optimal_pK, 
    nExp = nExp_poi.adj,
    reuse.pANN = pANN_name, 
    sct = sct
  )
  
  # 获取双检测结果列名
  doublet_col <- colnames(seurat.rna@meta.data)[
    grep(paste0("DF.classifications_0.25_", optimal_pK, "_", nExp_poi.adj), 
         colnames(seurat.rna@meta.data))
  ]
  
  if (length(doublet_col) > 0) {
    seurat.rna[['Doublet_Singlet']] <- seurat.rna[[doublet_col]]
  } else {
    warning("Doublet column not found. Using default naming.")
    doublet_col <- paste0('DF.classifications_0.25_', optimal_pK, '_', nExp_poi.adj)
    if (doublet_col %in% colnames(seurat.rna@meta.data)) {
      seurat.rna[['Doublet_Singlet']] <- seurat.rna[[doublet_col]]
    }
  }
  
  # 清理中间列
  pann_cols <- grep("^pANN_|^DF.classifications_", colnames(seurat.rna@meta.data), value = TRUE)
  cols_to_remove <- setdiff(pann_cols, "Doublet_Singlet")
  seurat.rna@meta.data[, cols_to_remove] <- NULL
  return(seurat.rna)
}

# 定义样本名称和组别
samples <- c("C1", "C2", "C3", "C4", "T1", "T2", "T3", "T4")
sample_dirs <- c("C1", "C2", "C3", "C4", "T1", "T2", "T3", "T4")
groups <- c(rep("Scramble", 4), rep("shSLC52A2", 4))

# 检查所有样本目录是否存在
missing_dirs <- sample_dirs[!dir.exists(sample_dirs)]
if (length(missing_dirs) > 0) {
  stop(paste("以下目录不存在：", paste(missing_dirs, collapse = ", ")))
}

# 创建样本列表
sample_list <- list()

# 读取所有样本数据
for (i in seq_along(samples)) {
  sample_name <- samples[i]
  sample_dir <- sample_dirs[i]
  
  cat("正在读取样本:", sample_name, "从目录:", sample_dir, "\n")
  
  # 读取10X数据
  data <- Read10X(data.dir = sample_dir, gene.column = 1)
  
  # 创建Seurat对象 - Seurat V5方式
  seurat_obj <- CreateSeuratObject(
    counts = data, 
    project = sample_name,
    min.cells = 3,
    min.features = 100,
    assay = "RNA"  # 明确指定assay名称
  )
  
  # 添加样本信息和组别
  seurat_obj$sample <- sample_name
  seurat_obj$group <- groups[i]
  
  sample_list[[sample_name]] <- seurat_obj
  cat("样本", sample_name, "读取完成，细胞数:", ncol(seurat_obj), "\n")
}

# 3.预处理每个样本以去除双细胞 ----
cat("\n开始预处理每个样本...\n")

for (i in seq_along(sample_list)) {
  sample_name <- names(sample_list)[i]
  cat("\n处理样本:", sample_name, "\n")
  
  obj <- sample_list[[sample_name]]
  
  # 计算线粒体基因百分比
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-|^mt-")
  
  # 质量过滤
  obj <- subset(obj, subset = nFeature_RNA > 100 & nCount_RNA > 1000 & nFeature_RNA < 10000 & percent.mt < 10)
  cat("  质量过滤后细胞数:", ncol(obj), "\n")
  
  # 标准化
  obj <- NormalizeData(obj)
  
  # 寻找高变基因
  obj <- FindVariableFeatures(obj, nfeatures = 2000)
  
  # 缩放数据
  obj <- ScaleData(obj)
  
  # PCA降维
  obj <- RunPCA(obj, features = VariableFeatures(object = obj), verbose = FALSE)
  
  # 寻找邻居和聚类
  obj <- FindNeighbors(obj, dims = 1:30)
  obj <- FindClusters(obj, resolution = 1)
  
  # UMAP可视化
  obj <- RunUMAP(obj, dims = 1:30, reduction.name = "umap_30")
  
  # 双细胞检测
  # 计算预期双细胞率（基于细胞数量）
  n_cells <- ncol(obj)
  exp_rate <- min(0.075, max(0.01, n_cells / 100000))
  cat("  预期双细胞率:", round(exp_rate * 100, 2), "%\n")
  
  obj <- FindDoublets(obj, PCs = 1:30, sct = FALSE, exp_rate = exp_rate)
  
  # 更新样本列表
  sample_list[[sample_name]] <- obj
}

# 保存处理后的单个样本（可选）
for (sample_name in names(sample_list)) {
  saveRDS(sample_list[[sample_name]], file = paste0(sample_name, "_processed.rds"))
}

# 合并所有样本 - Seurat V5方式
cat("\n合并所有样本...\n")

# 方法1：使用merge函数
combined <- merge(
  x = sample_list[[1]],
  y = sample_list[-1],
  add.cell.ids = names(sample_list),
  project = "MultiSample"
)

# 或者方法2：使用Seurat V5的集成功能（推荐用于批次校正）
# 如果样本间批次效应明显，建议使用整合方法

# 简单查看合并后数据
cat("合并后总细胞数:", ncol(combined), "\n")
cat("样本分布:\n")
print(table(combined$sample))


# 4.对整个数据集进行RNA污染检测 ----
library(Matrix)
library(glue)
library(celda)
library(Seurat)
library(ggplot2)
# 4.1. 确保 seu 中有 n_genes 信息（如果已有可以跳过这一步）
seu$is_HQ <- seu$nFeature_RNA >= 500   # 注意：Seurat 中基因数叫 nFeature_RNA

# 4.2. 提取高质量细胞的计数矩阵
counts_hq <- GetAssayData(seu, assay = "RNA", layer = "counts")[, seu$is_HQ]

# 4.3. 创建 SingleCellExperiment 对象（只包含高质量细胞）
sce <- SingleCellExperiment(
  assays = list(counts = counts_hq),
  colData = seu[[]][seu$is_HQ, ]   # 注意：seu[[]] 获取所有 meta.data
)

# 4.4. 运行 decontX
sce <- decontX(sce)

# 4.5. 将污染分数回填到 Seurat 对象的 meta.data 中
seu$ambient_frac_decontx <- NA_real_
seu[[]][colnames(sce), "ambient_frac_decontx"] <- sce$decontX_contamination

# 4.6. 查看结果
head(seu[[]])
summary(seu$ambient_frac_decontx)

FeaturePlot(seu, features = "ambient_frac_decontx", 
            reduction = "umap_harmony", label = TRUE) +
  scale_colour_gradientn(colours = c("blue", "grey90", "red"), 
                         limits = c(0, 0.25))

FeaturePlot(seu, features = "ambient_frac_decontx", 
            reduction = "Tcell_UMAP_dim50", label = TRUE) +
  scale_colour_gradientn(colours = c("blue", "grey90", "red"), 
                         limits = c(0, 0.5))


# 5.对整个数据集进行重新预处理和聚类 ----
cat("\n对整个数据集进行预处理...\n")

# 质量过滤（再次确保）
combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^MT-|^mt-")
combined <- subset(combined, subset = nFeature_RNA > 100 & nFeature_RNA < 10000 & percent.mt < 10)
combined <- subset(combined, subset = Doublet_Singlet == "Singlet")

cat("双细胞去除后细胞数:", ncol(combined), "\n")

# 标准化和特征选择
combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined, nfeatures = 3000)
combined <- ScaleData(combined, verbose = FALSE)

# PCA降维
combined <- RunPCA(combined, features = VariableFeatures(object = combined), verbose = FALSE)

# 检查PCA结果
ElbowPlot(combined, ndims = 50)

combined <- RunHarmony(
  combined,
  group.by.vars = "sample",  # 根据样本去除批次
  assay.use = "RNA",
  plot_convergence = TRUE
)


combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:30)
combined <- FindClusters(combined, resolution = 0.5)
# 后续分析使用Harmony矫正后的降维
combined <- RunUMAP(combined, dims = 1:30, reduction = "harmony", reduction.name = "umap_harmony")
# combined <- RunUMAP(combined, dims = 1:10, reduction = "harmony", reduction.name = "GCDH_UMAP_dim10")
# combined <- RunUMAP(combined, dims = 1:20, reduction = "harmony", reduction.name = "GCDH_UMAP_dim20")
# combined<- RunUMAP(combined, dims = 1:30, reduction = "harmony", reduction.name = "GCDH_UMAP_dim30")
# combined<- RunUMAP(combined, dims = 1:40, reduction = "harmony", reduction.name = "GCDH_UMAP_dim40")
# combined<- RunUMAP(combined, dims = 1:50, reduction = "harmony", reduction.name = "GCDH_UMAP_dim50")

# 可以比较批次效应去除前后的差异
# combined <- RunUMAP(combined, reduction = "pca", dims = 1:30, reduction.name = "umap_pca")

# 设置默认分辨率

# 5.1. 计算多分辨率
resolutions <- seq(0.5, 2.0, by = 0.1)
for (res in resolutions) {
  combined <- FindClusters(combined, 
                           resolution = res, 
                           cluster.name = paste0("RNA_snn_res.", res))
}

# 5.2. 使用 clustree 可视化
library(clustree)
p=clustree(combined, prefix = "RNA_snn_res.")
ggsave("clustree.png", plot = p, width = 20, height = 18, dpi = 300)


# 更简单：只选择你计算过的分辨率列表
resolutions_used <- seq(0.5, 2.0, by = 0.1)  # 改成你实际算的

plot_list <- list()
for (res in resolutions_used) {
  col_name <- paste0("RNA_snn_res.", res)
  
  # 检查这个分辨率是否存在
  if (!col_name %in% colnames(combined@meta.data)) {
    warning(paste(col_name, "not found, skipping"))
    next
  }
  
  p <- DimPlot(combined, 
               reduction = "umap_harmony", 
               group.by = col_name,
               label = TRUE, 
               repel = TRUE) +
    ggtitle(paste0("Resolution = ", res)) +
    theme(plot.title = element_text(hjust = 0.5, size = 10))
  
  plot_list[[length(plot_list) + 1]] <- p
}

# 拼接
P2=wrap_plots(plot_list, ncol = 4)

ggsave("umaptree_umap_harmony.png", plot =P2, width = 24, height = 18, dpi = 300)


# 5.3. 选一个最合理的作为最终聚类
# 假设你决定用 res = 0.5
# combined <- FindClusters(combined, resolution = 0.5)
Idents(combined) <- "RNA_snn_res.0.5"

# 4.4. 检查分群质量
DimPlot(combined, reduction = "umap_harmony", label = TRUE)


# 可视化
p1 <- DimPlot(combined, reduction = "umap_harmony", group.by = "sample", label = FALSE) +
  ggtitle("By Sample")
p2 <- DimPlot(combined, reduction = "umap_harmony", group.by = "group", label = FALSE) +
  ggtitle("By Group")
p3 <- DimPlot(combined, reduction = "umap_harmony", label = TRUE, repel = TRUE) +
  ggtitle("By Cluster")

# 保存可视化结果
pdf("combined_visualization.pdf", width = 15, height = 5)
print(p1 + p2 + p3)
dev.off()


# 保存质量控制的统计图
VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        group.by = "sample", pt.size = 0.1, ncol = 3) +
  NoLegend()
ggsave("QC_metrics_by_sample.pdf", width = 12, height = 6)

# 6.寻找所有marker基因 ----
cat("\n寻找差异表达基因...\n")
# 合并数据层

DefaultAssay(combined) <- "RNA"
combined <- JoinLayers(combined)
combined$RNA_snn_res.0.5=factor(combined$RNA_snn_res.0.5)
# 使用默认聚类分辨率
Idents(combined) <- "RNA_snn_res.0.5"

# 寻找marker基因，设置随机种子保证可重复性
set.seed(123)
combined_markers <- FindAllMarkers(
  combined,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  max.cells.per.ident = 1000,  # 下采样以提高速度
  random.seed = 123
)

# 快速测试一个cluster
# test_markers <- FindMarkers(
#   combined,
#   ident.1 = "0",
#   only.pos = TRUE,
#   min.pct = 0.1,  # 降低阈值测试
#   logfc.threshold = 0.1
# )


# 保存marker基因
write.csv(combined_markers , "combined_markers_res.csv", row.names = FALSE)
cat("Marker基因已保存到 combined_markers_res.csv\n")

# 保存完整的Seurat对象
saveRDS(combined, file = "combined_anno.rds")
cat("完整的Seurat对象已保存到 combined_seurat_object.rds\n")

# 6.去除污染细胞后的再分析（如果需要）----
# 这部分代码可以根据第一次聚类结果手动调整
cat("\n分析完成！\n")
cat("下一步建议：\n")
cat("1. 检查combined_markers_res1.csv文件，识别细胞类型\n")
cat("2. 查看UMAP图，识别可能的污染细胞群\n")
cat("3. 根据marker基因手动注释细胞类型\n")
cat("4. 去除污染细胞后重新运行分析（如果需要）\n")

# 辅助函数：去除指定聚类并重新分析
remove_clusters_and_reanalyze <- function(seurat_obj, clusters_to_remove, output_prefix = "filtered") {
  # 去除指定聚类
  seurat_obj <- subset(seurat_obj, subset = RNA_snn_res.1 %in% clusters_to_remove, invert = TRUE)
  cat("去除聚类", paste(clusters_to_remove, collapse = ", "), "后细胞数:", ncol(seurat_obj), "\n")
  
  # 重新分析
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 3000)
  seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
  seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
  seurat_obj <- FindClusters(seurat_obj, resolution = c(0.5, 1.0, 1.5))
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)
  
  # 保存结果
  saveRDS(seurat_obj, file = paste0(output_prefix, "_seurat_object.rds"))
  
  # 寻找marker基因
  Idents(seurat_obj) <- "RNA_snn_res.1"
  markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, max.cells.per.ident = 1000)
  write.csv(markers, file = paste0(output_prefix, "_markers.csv"), row.names = FALSE)
  
  return(seurat_obj)
}

# 使用方法示例（取消注释并根据需要修改）：
# combined_filtered <- remove_clusters_and_reanalyze(
#   combined, 
#   clusters_to_remove = c(10, 15, 20),  # 需要去除的聚类编号
#   output_prefix = "combined_filtered"
# )