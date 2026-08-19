# Load necessary libraries
library(ggplot2)
library(dplyr)
library(ggrepel)





data =read.table("B0_vs_B50.result.diff.txt",header = 1,sep = "\t")
data$logFC=data$log2FoldChange

gmt <- readLines("XXXXXXX.gmt")

# Extract gene sets
gene_set1 <- strsplit(gmt[grepl("MORF_GPX4", gmt)], "\t")[[1]][-c(1, 2)]
gene_set2 <- strsplit(gmt[grepl("MORF_GPX4", gmt)], "\t")[[1]][-c(1, 2)]
gene_sig <- c("Slc7a11")
gene_set3=as.data.frame(gene_set1)


data$external_gene_name=data$symbol
# 1. Filter data based on gene sets
data_gene_set <- data %>%
  filter(external_gene_name %in% c(gene_set1
                                   # , gene_set2
                                   , gene_sig)) %>%
  mutate(gene_set = case_when(
    external_gene_name %in% gene_set1 ~ "gene_set1",
    # external_gene_name %in% gene_set2 ~ "gene_set2",
    external_gene_name %in% gene_sig ~ "gene_sig"
  ))

data_other_genes <- data %>%
  filter(!external_gene_name %in% c(gene_set1
                                    # , gene_set2
                                    , gene_sig))
#data$logFC=data$log2FoldChange
# 2. Filter DOWN/UP DEGs
down_genes <- data %>% filter(FDR < 0.05, logFC < -0.5)
up_genes <- data %>% filter(FDR < 0.05, logFC > 0.5)

# 3. Filter DOWN/UP && top200 DEGs
top_down200_genes <- data %>% 
  filter(FDR < 0.01, logFC < 0) %>% 
  arrange(logFC) %>% 
  head(200)

top_up200_genes <- data %>% 
  filter(FDR < 0.01, logFC > 0) %>% 
  arrange(desc(logFC)) %>% 
  head(200)




# volcano
p1 <- ggplot(data, aes(x = logFC, y = -log10(FDR))) +
  # annotate("rect", xmin = min(up_genes$logFC), xmax = max(up_genes$logFC), ymin = -log10(0.05), ymax = -log10(min(up_genes$FDR)), fill = "#DCEABB") + #文章颜色#FDE7E9
  # #annotate("rect", xmin = sort(down_genes$logFC)[2], xmax = max(down_genes$logFC), ymin = -log10(0.05), ymax = -log10(min(up_genes$FDR)), fill = "#CCDFF1") + #文章颜色#E6E7FC
  # annotate("rect", xmin =-min(up_genes$logFC), xmax = -max(up_genes$logFC), ymin = -log10(0.05), ymax = -log10(min(up_genes$FDR)), fill = "#CCDFF1") +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 2) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "#A2CA55", linewidth = 0.6) +
  geom_point(data = data_other_genes, shape = 21, color = "black", alpha = 0.1, size = 1.2, stroke =1) +
  geom_point(data = data_gene_set, aes(fill = gene_set), shape = 21, color ="black", size = 4, stroke = 1) +
  # geom_label_repel(data = filter(data_gene_set, gene_set == "gene_sig"), aes(label = external_gene_name), size = 5, box.padding=unit(0.35, "lines"), segment.colour = "grey30") +
  scale_fill_manual(values = c("gene_set1" = "#5EA7D3"
                               # , "gene_set2" = "#A2CA55"
                               , "gene_sig" = "#5EA7D3"
                               )) + #文章颜色#07F1F9
  # annotate("rect", xmin = min(top_up200_genes$logFC), xmax = max(top_up200_genes$logFC), ymin = -log10(0.01), ymax = -log10(min(up_genes$FDR)), fill = "transparent", linetype = "dotted", color = "black", linewidth = 0.6) +
  # annotate("rect", xmin = sort(top_down200_genes$logFC)[2], xmax = max(top_down200_genes$logFC), ymin = -log10(0.01), ymax = -log10(min(up_genes$FDR)), fill = "transparent", linetype = "dotted", color = "black", linewidth = 0.6) +
  scale_x_continuous(limits = c(-5, 5), breaks = seq(-5, 5, 2)) +
  labs(x = "", y = "Significance (-log10 adj.p-value)") + 
  theme_classic(base_size = 24) + 
  theme (legend.position = "none")

p1



# add signature enrichment bar添加富集信号通路的位置和竖直线的粗细
p2 <- p1 + 
  scale_y_continuous(limits = c(-20, 110), breaks = seq(0, 110, 25), expand = c(0, 0)) +
  geom_linerange(data = filter(data_gene_set, gene_set == "gene_set1"), aes(x = logFC, ymin = -16, ymax = -4), color = "#5EA7D3" , size = 1, linewidth = 0.5) 

  # +geom_linerange(data = filter(data_gene_set, gene_set == "gene_set2"), aes(x = logFC, ymin = -26.25 , ymax = -16.5), color = "#A2CA55", size = 0.5, linewidth = 0.1)

p2

# add text

# p3=p2 + geom_label_repel(data = filter(data_gene_set, external_gene_name == "Slc7a11"), 
#                       aes(label = external_gene_name), 
#                       size = 5, 
#                       box.padding = unit(0.35, "lines"), 
#                       segment.colour = "grey30")
# 
# p3

# 确保加载必要的包
# 先创建带nudge参数的标签数据
# gene_labels <- data.frame(
#   gene = c("Slc7a11","Slc3a2", "Gclc", "Gclm" ,"Gpx4","Aifm2"),
#   nudge_x = c(-1,-2,-1,-1.5,-2,2),        # 统一向左偏移
#   nudge_y = c(50, 32,30, 20, 10,10)  # 垂直方向错开
# )
# 
# p3 = p2 + 
#   geom_text_repel(
#     data = dplyr::filter(data_gene_set, external_gene_name %in% gene_labels$gene) %>%
#       left_join(gene_labels, by = c("external_gene_name" = "gene")),
#     aes(label =external_gene_name,
#         nudge_x = nudge_x,
#         nudge_y = nudge_y),
#     size = 7,
#     box.padding = unit(0.35, "lines"),
#     segment.colour = "grey30",
#     point.padding = unit(0.8, "lines"),
#     min.segment.length = 0,
#     max.overlaps = Inf,
#     force = 2,
#     segment.size = 0.3
#   )
# 
# p3


gene_labels <- data.frame(
  gene = c("Slc7a11","Slc3a2", "Gclc", "Gclm" ,"Gpx4","Aifm2"),
  nudge_x = c(-1,-2,-1,-1.5,-2,2),
  nudge_y = c(55, 33,30, 20, 10,10)
)

# 准备数据
plot_data <- dplyr::filter(data_gene_set, external_gene_name %in% gene_labels$gene) %>%
  left_join(gene_labels, by = c("external_gene_name" = "gene")) %>%
  mutate(
    label_color = ifelse(external_gene_name == "Aifm2", "#F38989", "black"),
    segment_color = ifelse(external_gene_name == "Aifm2", "#F38989", "black"),
    point_color = ifelse(external_gene_name == "Aifm2", "#F38989", "#3A84B3")  # 点颜色
  )

# 方法1：在p2中直接修改
p2_modified <- p2 + 
  geom_point(data = plot_data, 
             aes(color= point_color),
              size = 3, stroke = 1,
             # 根据您的p2中的点大小调整
             show.legend = FALSE) +
  scale_color_identity()

p3 = p2_modified + 
  geom_text_repel(
    data = plot_data,
    aes(label = external_gene_name,
        nudge_x = nudge_x,
        nudge_y = nudge_y,
        color = label_color),
    size = 7,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.8, "lines"),
    min.segment.length = 0,
    max.overlaps = Inf,
    force = 2,
    segment.size = 0.3
  ) +
  scale_color_identity()

p3

p4=p3 + 
  annotate("text", x = 3, y =100, label = "Up-regulated", color ="black", size =7.5, lineheight = 0.8, vjust = 0) + #文章颜色#857CD9
  annotate("text", x = -2.5, y =100, label =  "Down-regulated", color =  "black", size = 7.5, lineheight = 0.8, vjust = 0) + #文章颜色#FF7D81
  # annotate("text", x = -5, y = 200, label = "Top 200", color = "black", size = 5) +
  # annotate("text", x = 5, y = 200, label = "Top 200", color = "black", size = 5) +
  # annotate("text", x = 8, y = 1.1, label = "α = 0.05", color =  "#FFBF00", size = 4.5) +
  annotate("text", x = -3.5, y = -10.5, label = "Anti-ferroptosis", color =  "#5EA7D3",size = 5) +
  # annotate("text", x = -7, y = -21, label = "Ferroptosis", color = "#A2CA55", size = 4.5) +
  annotate("text", x =5, y = -10.5, label = "NES:-2.29 FDR:2.01e-07", color ="#5EA7D3", size =4.5, lineheight = 1, hjust = 1) 
# +
# annotate("text", x = 10, y = -21, label = "NES:1.719 pval:2.81e-05", color = "#A2CA55", size = 4.5, lineheight = 0.8, hjust = 1)
p4

p5=p4+ 
  theme(
    # X轴和Y轴线条加粗
    axis.line = element_line(linewidth =2, color = "black"),  # 同时加粗X和Y轴线
    # 或者分别设置
    # axis.line.x = element_line(linewidth = 1.2, color = "black"),  # 只加粗X轴
    # axis.line.y = element_line(linewidth = 1.2, color = "black"),  # 只加粗Y轴
    
    # X轴和Y轴刻度数字加粗
    axis.text = element_text(face = "bold"),  # 同时加粗X和Y轴数字
    # 或者分别设置
    # axis.text.x = element_text(face = "bold"),  # 只加粗X轴数字
    # axis.text.y = element_text(face = "bold"),  # 只加粗Y轴数字
    
    # 还可以调整数字大小
    axis.text.x = element_text( family = "Arial", face = "plain",size = 28, color = "black"),  # X轴数字加粗并调大
    axis.text.y = element_text( family = "Arial", face = "plain",size = 28, color = "black"),  # Y轴数字加粗并调大
    
    # 轴标题加粗
    axis.title = element_text(face = "bold"),  # 同时加粗X和Y轴标题
    axis.title.x = element_text(family = "Arial", face = "plain",size = 20),  # X轴标题加粗
    axis.title.y = element_text(family = "Arial", face = "plain",size = 20)   # Y轴标题加粗
  )

p5
ggsave("Fig4G volcano_NES.png", width = 6.5, height = 5)
ggsave("Fig4G volcano_NES.pdf", width = 6.5, height = 5)
ggsave("Fig4G volcano_NES.tiff",width =6.5, height = 5)
