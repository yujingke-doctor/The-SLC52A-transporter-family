library("magrittr")
library("sscVis")
library("data.table")
library("R.utils")
library("ggpubr")
library("ggplot2")
library("plyr")
library("grid")

dotplot <- read.table("Tcell_Markers.txt", 
                                sep = "\t", 
                                header = TRUE
                                # , 
                                # row.names = 1
                                )


# 确保x和y是因子，并设置顺序
gene.plot.tb1=dotplot

table(gene.plot.tb1$x)


gene.plot.tb1$x <- factor(gene.plot.tb1$x, levels =c("CD8+Tcm", "CD8+Teff", "CD8+Tem", 
                                                     "CD8+Tpex", "CD8+Tex", "CD4+Tcm",
                                                     "CD4+Treg", "CD4+Tconv", "γδ T"))
table(gene.plot.tb1$x)

# gene.plot.tb1$y <- factor(gene.plot.tb1$y, levels = c(
#   # 谱系
#   "Cd3d","Cd8a", "Cd4", "Foxp3", "Tcrg-C1", "Trdc",
#   # CD4亚型
#   "Tbx21", "Gata3", "Bcl6",
#   # 干性/记忆
#   "Tcf7", "Lef1", "Ccr7", "Sell", "Il7r", "Bcl2", "Satb1", "Bach2",
#   # 增殖
#   "Mki67", "Top2a",
#   # 代谢
#   "Rpl3", "Rps3", "Npm1", "Eef1b2",
#   # 效应
#   "Gzmb", "Gzma", "Gzmk", "Prf1", "Ifng", "Nkg7", "Cxcr6", "Ccl5",
#   # 激活
#   "Cd69", "Icos", "Cd28", "Klrk1",
#   # 趋化
#   "Cxcr3", "Ccr5", "S1pr1", "Itgae",
#   # 耗竭
#   "Pdcd1", "Tox", "Lag3", "Havcr2", "Tigit", "Cd160", "Klrc1", "Klrd1",
#   # γδ特征
#   "Rorc", "Il23r", "Zbtb16"
#   
# ))

# table(gene.plot.tb1$y)


gene.plot.tb1$Group <- factor(gene.plot.tb1$Group, levels =c(
  "Lineage ",
  "CD4",
  "Stemness_Memory",
  "Proliferation",
  "Metabolic_Ribosomal",
  "Effector",
  "Activation",
  "Chemotaxis_Homing",
  "Exhaustion",
  "GammaDelta"
))

gene.plot.tb1$y <- factor(gene.plot.tb1$y, levels = c(
  # γδ特征
  "Zbtb16", "Il23r", "Rorc",
  # 耗竭
  "Klrd1", "Klrc1", "Cd160", "Tigit", "Havcr2", "Lag3", "Tox", "Pdcd1",
  # 趋化
  "Itgae", "S1pr1", "Ccr5", "Cxcr3",
  # 激活
  "Klrk1", "Cd28", "Icos", "Cd69",
  # 效应
  "Ccl5", "Cxcr6", "Nkg7", "Ifng", "Prf1", "Gzmk", "Gzma", "Gzmb",
  # 代谢
  "Eef1b2", "Npm1", "Rps3", "Rpl3",
  # 增殖
  "Top2a", "Mki67",
  # 干性/记忆
  "Bach2", "Satb1", "Bcl2", "Il7r", "Sell", "Ccr7", "Lef1", "Tcf7",
  # CD4亚型
  "Bcl6", "Gata3", "Tbx21",
  # 谱系
  "Trdc", "Tcrg-C1", "Foxp3", "Cd4", "Cd8a", "Cd3d"
))

table(gene.plot.tb1$y)

# gene.plot.tb1$Group <- factor(gene.plot.tb1$Group, levels = c(
#   "GammaDelta",
#   "Exhaustion",
#   "Chemotaxis_Homing",
#   "Activation",
#   "Effector",
#   "Metabolic_Ribosomal",
#   "Proliferation",
#   "Stemness_Memory",
#   "CD4",
#   "Lineage "
# ))
# 
# table(gene.plot.tb1$Group)

library(ggplot2)
p <- ggplot(gene.plot.tb1,aes(x,y)) +
  geom_point(aes(size=pct.exp,color=avg.exp.scaled),shape=16) +
  facet_grid(Group ~ ., scales = "free", space = "free") +
  scale_colour_distiller(palette = "RdYlBu") +
  scale_size_continuous(range = c(1, 10)) +
  labs(x="",y="") +
  #theme_pubr() +
  theme(strip.text.y = element_blank(),
        axis.line.x=element_blank(),
        axis.line.y=element_blank(),
        panel.background = element_rect(colour = "black", fill = "white"),
        # panel.grid = element_line(colour = "grey", linetype = "dashed"),
        # panel.grid.major = element_line( colour = "grey", linetype = "dashed", size = 0.2),
        axis.text.y = element_text(size=20,colour = "black", face = "italic"),
        axis.text.x = element_text(angle = 60,size=24,colour = "black", hjust = 1)
        ,panel.border = element_rect(
          colour = "black",      # 边框颜色
          fill = NA,             # 不填充
          linewidth = 2        # 边框粗细
        ))

ggsave("DotPlot_Tcell_Markers_Complete.pdf", width = 8, height = 20)
ggsave("DotPlot_Tcell_Markers_Complete.png", width = 8, height = 20)
ggsave("DotPlot_Tcell_Markers_Complete.tiff", width = 8, height =20)



