##galaxy plot showing cell density 
library(ggplot2)
library(viridis)
library(ggpubr)

galaxyTheme_black = function(base_size = 12, base_family = "") {
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(
      # Specify axis options
      axis.line = element_blank(),  
      axis.text.x = element_text(size = base_size*0.8, color = "white", lineheight = 0.9),  
      axis.text.y = element_text(size = base_size*0.8, color = "white", lineheight = 0.9),  
      axis.ticks = element_line(color = "white", size  =  0.2),  
      axis.title.x = element_text(size = base_size, color = "white", margin = margin(0, 10, 0, 0)),  
      axis.title.y = element_text(size = base_size, color = "white", angle = 90, margin = margin(0, 10, 0, 0)),  
      axis.ticks.length = unit(0.3, "lines"),   
      # Specify legend options
      legend.background = element_rect(color = NA, fill = "black"),  
      legend.key = element_rect(color = "white",  fill = "black"),  
      legend.key.size = unit(1.2, "lines"),  
      legend.key.height = NULL,  
      legend.key.width = NULL,      
      legend.text = element_text(size = base_size*0.8, color = "white"),  
      legend.title = element_text(size = base_size*0.8, face = "bold", hjust = 0, color = "white"),  
      legend.position = "none",  
      legend.text.align = NULL,  
      legend.title.align = NULL,  
      legend.direction = "vertical",  
      legend.box = NULL, 
      # Specify panel options
      panel.background = element_rect(fill = "black", color  =  NA),  
      panel.border = element_rect(fill = NA, color = "white"),  
      ##panel.grid.major = element_line(color = "grey35"),  
      panel.grid.major = element_blank(),  
      ##panel.grid.minor = element_line(color = "grey20"),  
      panel.grid.minor = element_blank(),  
      panel.spacing = unit(0.5, "lines"),   
      # Specify facetting options
      strip.background = element_rect(fill = "grey30", color = "grey10"),  
      strip.text.x = element_text(size = base_size*0.8, color = "white"),  
      strip.text.y = element_text(size = base_size*0.8, color = "white",angle = -90),  
      # Specify plot options
      plot.background = element_rect(color = "black", fill = "black"),  
      plot.title = element_text(size = base_size*1.2, color = "white"),  
      plot.margin = unit(rep(1, 4), "lines")
    )
}   


#####Load scRNA-seq

# 提取 UMAP 坐标
umap_coords <- Embeddings(Singlet, reduction ="umap_harmony")
# 或者如果是其他命名，如 "tsne"、"pca" 等
# umap_coords <- Embeddings(Singlet, reduction = "umap")  # 根据实际 reduction 名称调整

# 提取元数据（包含分组信息和细胞名称）
metadata <- Singlet@meta.data

# 添加细胞名称列（如果元数据中没有）
metadata$cell_name <- rownames(metadata)

# 合并 UMAP 坐标和元数据
umap_data <- cbind(metadata, umap_coords)

# 查看数据结构
head(umap_data)

data=umap_data
colnames(data)
data$group <- factor(data$group, levels = c("Scramble", "shSLC52A2"))

galaxy_Umap <- ggplot(data = data,
                      aes(x = umapharmony_1, y = umapharmony_2)) + 
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = F) + 
  # geom_polygon(data = hull, mapping = aes(x = V1, y = V2), 
  #              color = "white", fill = 'transparent', linetype = "dashed") +
  geom_point(color = 'white',size = .02) + 
  # annotate(geom = "text", x = 1, y = 4, label = "xxx", color = "white",size = 6) +
  # annotate(geom = "text", x = -2, y = 0, label = "xxx", color = "white",size = 6) +
  facet_wrap(~group,ncol = 2) +
  scale_fill_viridis(option="magma") +
  galaxyTheme_black(base_size = 20)+
  SeuratExtend::theme_umap_arrows(x_label = "UMAP1",y_label = "UMAP2", 
                                  text_offset_x = unit(2, 'mm'),
                                  text_offset_y = unit(2, 'mm'),
                                  text_size = 15) +NoLegend()

galaxy_Umap <- ggplot(data = data,
                      aes(x = umapharmony_1, y = umapharmony_2)) + 
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = F) + 
  geom_point(color = 'white', size = .02) + 
  facet_wrap(~group, ncol = 2) +
  scale_fill_viridis(option = "magma") +
  galaxyTheme_black(base_size = 20) +
  SeuratExtend::theme_umap_arrows(x_label = "UMAP1", y_label = "UMAP2", 
                                  text_offset_x = unit(2, 'mm'),
                                  text_offset_y = unit(2, 'mm'),
                                  text_size = 15) +
  NoLegend() +
  theme(axis.title.x = element_blank(),      # 隐藏 x 轴标题
        axis.title.y = element_blank(),      # 隐藏 y 轴标题
        axis.text.x = element_blank(),       # 隐藏 x 轴刻度数字
        axis.text.y = element_blank(),       # 隐藏 y 轴刻度数字
        axis.ticks.x = element_blank(),      # 隐藏 x 轴刻度线
        axis.ticks.y = element_blank())      # 隐藏 y 轴刻度线

galaxy_Umap 

##suggest to save in a figure with large width and height, cause the cell dots plotted in a small figure affect the visualization of the actual density
ggsave("galaxy plot.pdf", plot = galaxy_Umap,  width = 12, height = 7)
ggsave("galaxy plot.png", plot = galaxy_Umap,  width = 12, height = 7)
ggsave("galaxy plot.tiff", plot = galaxy_Umap,  width = 12, height = 7)
