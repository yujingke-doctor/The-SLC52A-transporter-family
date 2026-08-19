library(readxl)
library(tidyverse)
library(scales)
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)
library(eulerr)
library(grid)

# 读取数据
data = read_xlsx("XXXXXXXXXX..xlsx") 

# 对数值列进行标准化
data_normalized <- data %>%
  group_by(Group) %>%
  arrange(qval, .by_group = TRUE) %>%
  mutate(across(c(`B0-1`, `B0-2`, `B0-3`, `B50-1`,`B50-2`,`B50-3`),
                ~ rescale(., to = c(0,0.5)))) %>%  #标准化到0-0.5的范围内
  ungroup()


# 转换数据为矩阵
data_matrix <- data_normalized %>%
  select(`B0-1`, `B0-2`, `B0-3`, `B50-1`,`B50-2`,`B50-3`) %>%
  as.matrix()

rownames(data_matrix) <- data_normalized$Gene  #设置行名

colnames(data_matrix)=c("  B0-1","  B0-2", "  B0-3", "  B50-1","  B50-2","  B50-3")
##########################################################
#定义颜色
green_pink = colorRamp2(seq(0,0.5,length.out=100), 
                        colorRampPalette(rev(brewer.pal(n = 5, name = "RdBu")))(100)) #文章配色RdBu

group_colors <- c("OXPHOS" = rev(brewer.pal(n = 5, name = "RdBu"))[1], 
                  "Mitophagy"= rev(brewer.pal(n = 5, name = "RdBu"))[3], 
                  "Mitochondrial Complex" = rev(brewer.pal(n = 5, name = "RdBu"))[5])

red_black <- colorRamp2(seq(0,1,length.out=6), 
                            c("#fff1e7","#f99c70","#dd4139","#97044d","#4f2043","black"))

# 设置分组的展示顺序
data_normalized$Group = factor(data_normalized$Group, levels = c("OXPHOS","Mitophagy","Mitochondrial Complex"))


# 设置全局字体家族为 Times New Roman
par(family = "Times")

pdf("circos_heatmap2.pdf", width = 10, height = 10)


# 初始化circlize
circos.clear()
circos.par(start.degree = 60,  #起始的角度
           gap.after = c(rep(5, length(levels(data_normalized$Group)) - 1), 30), #扇区之间的间隙大小
           track.margin = c(0, 0.01),
           cell.padding = c(0, 0, 0, 0))


##########################################################
# 使用circos.heatmap进行分组热图绘制
circos.heatmap(data_matrix, 
               split = data_normalized$Group, #按照分组分割热图
               cluster = FALSE,
               bg.border = "black",
               bg.lwd = 1,
               cell.border = "white",
               cell.lwd = 0.5,
               rownames.side = "outside",
               rownames.cex = 1.5, 
               rownames.font = 3, 
               col = green_pink, 
               track.height = 0.25)

# 添加列名：在最后一个扇区的右侧
# circos.track(
#   track.index = get.current.track.index(),
#   bg.border = NA,
#   panel.fun = function(x,y){
#     #判断当前的扇区是否是环形图中的最后一个扇区
#     if(CELL_META$sector.numeric.index == length(levels(data_normalized$Group))) {
#       #计算每个标签的y坐标，确保均匀分布
#       cn <- colnames(data_matrix)
#       n <- length(cn)
#       cell_height <- (CELL_META$cell.ylim[2] - CELL_META$cell.ylim[1]) / n
#       y_coords <- seq(CELL_META$cell.ylim[1] + cell_height / 2,  
#                       CELL_META$cell.ylim[2] - cell_height / 2, 
#                       length.out = n)
#       #添加线段
#       for (i in 1:n) {
#         circos.lines(
#           c(CELL_META$cell.xlim[2], CELL_META$cell.xlim[2] + convert_x(1, "mm")), # x坐标，1mm偏移量
#           c(y_coords[i], y_coords[i]),
#           col = "black",
#           lwd = 2
#         )
#       } 
#       #添加文本
#       circos.text(
#         rep(CELL_META$cell.xlim[2], n) + convert_x(1.5, "mm"), # x坐标，1.5mm偏移量
#         y_coords, # y坐标
#         cn,
#         cex = 1.4,
#         adj = c(0, 0.5),
#         facing = "inside"
#       )
#     }
#   }
# )

circos.track(
  track.index = get.current.track.index(),
  bg.border = NA,
  panel.fun = function(x, y) {
    if(CELL_META$sector.numeric.index == length(levels(data_normalized$Group))) {
      cn <- colnames(data_matrix)
      n <- length(cn)
      cell_height <- (CELL_META$cell.ylim[2] - CELL_META$cell.ylim[1]) / n
      
      # 计算 y 坐标（从下到上）
      y_coords <- seq(CELL_META$cell.ylim[1] + cell_height / 2,  
                      CELL_META$cell.ylim[2] - cell_height / 2, 
                      length.out = n)
      
      # 反转 y_coords，使标签从上到下排列
      y_coords <- rev(y_coords)  # 关键修改！
      
      # 添加线段
      for (i in 1:n) {
        circos.lines(
          c(CELL_META$cell.xlim[2], CELL_META$cell.xlim[2] + convert_x(1, "mm")),
          c(y_coords[i], y_coords[i]),
          col = "black",
          lwd = 2
        )
      }
      
      # 添加文本
      circos.text(
        rep(CELL_META$cell.xlim[2], n) + convert_x(1.5, "mm"),
        y_coords,
        cn,  # cn 保持原顺序，但 y_coords 反转了
        cex = 1.4,
        adj = c(0, 0.5),
        facing = "inside"
      )
    }
  }
)

##########################################################
# 添加qval散点图
circos.track(
  ylim = c(0, 1), 
  track.height = 0.05, 
  bg.border = NA, 
  panel.fun = function(x, y) {
    sector_data <- data_normalized[data_normalized$Group == CELL_META$sector.index, ]
    
    #添加散点
    for (i in 1:nrow(sector_data)) {
      circos.points(
        CELL_META$xlim[1] + (CELL_META$xlim[2] - CELL_META$xlim[1]) * (i - 0.5) / nrow(sector_data), 
        0.5, 
        pch = 18,
        cex = 1.8,
        col = red_black(sector_data$qval[i]*100)
      )
    }
    
    #添加线段和文本
    if(CELL_META$sector.numeric.index == length(levels(data_normalized$Group))) {
      circos.lines(
        c(CELL_META$cell.xlim[2], CELL_META$cell.xlim[2] + convert_x(1, "mm")), 
        c(0.5, 0.5),
        col = "black",
        lwd = 2
      )
      circos.text(
        CELL_META$cell.xlim[2] + convert_x(1.5, "mm"), 
        0.5,
        "pvalue",
        cex = 1.4,
        adj = c(0, 0.5),
        facing = "inside"
      )
    }
  }
)


##########################################################
# 添加扇形分组标签
circos.track(
  ylim = c(0, 1), 
  track.height = 0.065, 
  bg.col = adjustcolor(group_colors[levels(data_normalized$Group)], alpha.f = 0.3),
  panel.fun = function(x, y) {
    circos.text(
      CELL_META$xcenter,
      CELL_META$ylim[2] - 0.75,
      CELL_META$sector.index, 
      facing = "bending.inside", 
      cex = 1.5, 
      adj = c(0.5, 0)
    )
  }
)


##########################################################
# 添加图例
heatmap_legend <- Legend(title = "Expression\nfraction", 
                         col_fun = green_pink, 
                         at = seq(0, 0.5, length.out = 6),
                         title_position = "leftcenter-rot",
                         title_gp = gpar(fontsize = 14),
                         labels_gp = gpar(fontsize = 14))

qval_legend <- Legend(title = "pvalue", 
                      col_fun = red_black, 
                      at = seq(0, 1, length.out = 6),
                      title_position = "leftcenter-rot",
                      title_gp = gpar(fontsize = 14),
                      labels_gp = gpar(fontsize = 14))

draw(heatmap_legend, x = unit(0.95, "npc") - unit(5, "mm"), y = unit(0.9, "npc") - unit(5, "mm"), just = c("right", "top"))
draw(qval_legend, x = unit(0.93, "npc") - unit(5, "mm"), y = unit(0.25, "npc") - unit(5, "mm"), just = c("right", "top"))


##########################################################
# 提取 venn 图数据：获取每个Group中的基因集合
genes_acceleration <- data %>% filter(Group =="OXPHOS") %>% pull(Gene)
genes_divergence <- data %>% filter(Group =="Mitophagy") %>% pull(Gene)
genes_curl <- data %>% filter(Group =="Mitochondrial Complex") %>% pull(Gene)

# 创建一个列表，包含每个Group的基因集合
gene_groups <- list(
  Acceleration = genes_acceleration,
  Divergence = genes_divergence,
  Curl = genes_curl
)

# 使用 eulerr 包生成按比例调整大小的 venn 图
fit <- euler(gene_groups)

# 创建一个空视口来管理Venn图的插入
pushViewport(viewport(width = 0.24, height = 0.24))

# 在最内层插入Venn图
circos.track(
  ylim = c(0, 1),
  track.height = 0.3,
  bg.border = NA,
  panel.fun = function(x, y) {
    grid.draw(plot(fit, 
                   fills = list(fill = c("#005094","#F7F7F7","#7A111F"), alpha = 0.2), 
                   edges = list(col = "black", lty = 2, lwd = 2),
                   labels = FALSE,
                   quantities = TRUE))
    grid.text("Mitochondrial\ngenes", 
              x = unit(0.5, "npc"), 
              y = unit(-0.01, "npc"), 
              gp = gpar(fontsize = 25, lineheight = 0.8))
  }
)

##########################################################
# 关闭PDF设备
dev.off()

###########重新保留3种格式的图片######
# 设置全局字体家族为 Times New Roman
par(family = "Times")

# 定义绘图函数，避免代码重复
draw_circos_heatmap <- function() {
  # 初始化circlize
  circos.clear()
  circos.par(start.degree = 60,
             gap.after = c(rep(5, length(levels(data_normalized$Group)) - 1), 30),
             track.margin = c(0, 0.01),
             cell.padding = c(0, 0, 0, 0))
  
  # 使用circos.heatmap进行分组热图绘制
  circos.heatmap(data_matrix, 
                 split = data_normalized$Group,
                 cluster = FALSE,
                 bg.border = "black",
                 bg.lwd = 1,
                 cell.border = "white",
                 cell.lwd = 0.5,
                 rownames.side = "outside",
                 rownames.cex = 1.6, 
                 rownames.font = 3,
                 col = green_pink, 
                 track.height = 0.25)
  
  # 添加列名（代码省略，保持原样）
  circos.track(
    track.index = get.current.track.index(),
    bg.border = NA,
    panel.fun = function(x, y) {
      if(CELL_META$sector.numeric.index == length(levels(data_normalized$Group))) {
        cn <- colnames(data_matrix)
        n <- length(cn)
        cell_height <- (CELL_META$cell.ylim[2] - CELL_META$cell.ylim[1]) / n
        y_coords <- seq(CELL_META$cell.ylim[1] + cell_height / 2,  
                        CELL_META$cell.ylim[2] - cell_height / 2, 
                        length.out = n)
        y_coords <- rev(y_coords)
        
        for (i in 1:n) {
          circos.lines(
            c(CELL_META$cell.xlim[2], CELL_META$cell.xlim[2] + convert_x(1, "mm")),
            c(y_coords[i], y_coords[i]),
            col = "black",
            lwd = 2
          )
        }
        
        circos.text(
          rep(CELL_META$cell.xlim[2], n) + convert_x(1.5, "mm"),
          y_coords,
          cn,
          cex = 1.4,
          adj = c(0, 0.5),
          facing = "inside"
        )
      }
    }
  )
  
  # 添加qval散点图
  circos.track(
    ylim = c(0, 1), 
    track.height = 0.05, 
    bg.border = NA, 
    panel.fun = function(x, y) {
      sector_data <- data_normalized[data_normalized$Group == CELL_META$sector.index, ]
      
      for (i in 1:nrow(sector_data)) {
        circos.points(
          CELL_META$xlim[1] + (CELL_META$xlim[2] - CELL_META$xlim[1]) * (i - 0.5) / nrow(sector_data), 
          0.5, 
          pch = 18,
          cex = 1.8,
          col = red_black(sector_data$qval[i]*100)
        )
      }
      
      if(CELL_META$sector.numeric.index == length(levels(data_normalized$Group))) {
        circos.lines(
          c(CELL_META$cell.xlim[2], CELL_META$cell.xlim[2] + convert_x(1, "mm")), 
          c(0.5, 0.5),
          col = "black",
          lwd = 2
        )
        circos.text(
          CELL_META$cell.xlim[2] + convert_x(1.5, "mm"), 
          0.5,
          "pvalue",
          cex = 1.4,
          adj = c(0, 0.5),
          facing = "inside"
        )
      }
    }
  )
  
  # 添加扇形分组标签
  circos.track(
    ylim = c(0, 1), 
    track.height = 0.065, 
    bg.col = adjustcolor(group_colors[levels(data_normalized$Group)], alpha.f = 0.3),
    panel.fun = function(x, y) {
      circos.text(
        CELL_META$xcenter,
        CELL_META$ylim[2] - 0.75,
        CELL_META$sector.index, 
        facing = "bending.inside", 
        cex = 1.5, 
        adj = c(0.5, 0)
      )
    }
  )
  
  # 添加图例
  heatmap_legend <- Legend(title = "Expression\nfraction", 
                           col_fun = green_pink, 
                           at = seq(0, 0.5, length.out = 6),
                           title_position = "leftcenter-rot",
                           title_gp = gpar(fontsize = 14),
                           labels_gp = gpar(fontsize = 14))
  
  qval_legend <- Legend(title = "pvalue", 
                        col_fun = red_black, 
                        at = seq(0, 1, length.out = 6),
                        title_position = "leftcenter-rot",
                        title_gp = gpar(fontsize = 14),
                        labels_gp = gpar(fontsize = 14))
  
  draw(heatmap_legend, x = unit(0.95, "npc") - unit(5, "mm"), y = unit(0.9, "npc") - unit(5, "mm"), just = c("right", "top"))
  draw(qval_legend, x = unit(0.93, "npc") - unit(5, "mm"), y = unit(0.25, "npc") - unit(5, "mm"), just = c("right", "top"))
  
  # 添加Venn图
  genes_acceleration <- data %>% filter(Group == "OXPHOS") %>% pull(Gene)
  genes_divergence <- data %>% filter(Group == "Mitophagy") %>% pull(Gene)
  genes_curl <- data %>% filter(Group == "Mitochondrial Complex") %>% pull(Gene)
  
  gene_groups <- list(
    Acceleration = genes_acceleration,
    Divergence = genes_divergence,
    Curl = genes_curl
  )
  
  fit <- euler(gene_groups)
  
  pushViewport(viewport(width = 0.24, height = 0.24))
  circos.track(
    ylim = c(0, 1),
    track.height = 0.3,
    bg.border = NA,
    panel.fun = function(x, y) {
      grid.draw(plot(fit, 
                     fills = list(fill = c("#005094","#F7F7F7","#7A111F"), alpha = 0.2), 
                     edges = list(col = "black", lty = 2, lwd = 2),
                     labels = FALSE,
                     quantities = TRUE))
      grid.text("Mitochondrial\ngenes", 
                x = unit(0.5, "npc"), 
                y = unit(0.1, "npc"), 
                gp = gpar(fontsize = 25, lineheight = 0.8))
    }
  )
}

# 保存为 PDF
pdf("Fig3N-circos_heatmap.pdf", width = 10, height = 10)
draw_circos_heatmap()
dev.off()

# 保存为 PNG
png("Fig3N-circos_heatmap.png", width = 10, height = 10, units = "in", res = 300)
draw_circos_heatmap()
dev.off()

# 保存为 TIFF
tiff("Fig3N-circos_heatmap.tiff", width = 10, height = 10, units = "in", res = 300, compression = "lzw")
draw_circos_heatmap()
dev.off()

