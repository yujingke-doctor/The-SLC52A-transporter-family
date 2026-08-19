
library(jstable)
library(survival)
library(grid)
library(forestploter)
library(tidyverse)

tm <- forest_theme(
  base_size = 20,                             # 设置文本的基础大小
  xaxis_gp = gpar(lwd = 2,                    # <--- 设置 X 轴主线的宽度，数值越大线越粗
                  cex = 2),                   # 可选：同时调整 X 轴刻度标签的字体大小
  arrow_gp = gpar(lwd = 2,                    # <--- 设置箭头线条的宽度，数值越大箭头越粗
                  fill = "black",
                  # cex = 1.5,                # 方法1：相对大小，1.5表示base_size的1.5倍
                  fontsize = 20),             # 可选：箭头填充颜色
  colhead=list(fg_params=list(hjust=0)),      # 设置列名对齐方式 
  core = list(padding = unit(c(6, 6), "mm")), #增加森林图每行的宽度
  
  # 设置置信区间的外观
  ci_pch = 15,                                # 置信区间点的形状
  ci_col = "#5EAE3E",                         # 置信区间的边框颜色
  ci_fill ="#71AFC4",                         # 置信区间的填充颜色
  ci_alpha = 0.8,                             # 置信区间的透明度
  ci_lty = 1,                                 # 置信区间的线型
  ci_lwd = 4,                                 # 置信区间的线宽
  ci_Theight = 0.2,                           # 设置T字在置信区间末端的高度，默认是NULL
  
  # 设置参考线的外观
  refline_gp = gpar(lwd =2,                   # 参考线的线宽
                    lty = "dashed",           # 参考线的线型
                    col = "grey20"            # 参考线的颜色
  ),
  
  #设置垂直线的外观
  vertline_lwd = 1,                           # 垂直线的线宽，可以添加一条额外的垂直线，如果没有就不显示
  vertline_lty = "dashed",                    # 垂直线的线型
  vertline_col = "grey20",                    # 垂直线的颜色
  
  #设置脚注的字体大小、字体样式和颜色
  footnote_gp = gpar(
    # fontfamily = "italic",                  # 脚注文本的字体
    cex = 2,                                  # 脚注字体大小
    fontface = "plain",                       # 脚注字体样式
    col = "black"                             # 脚注文本的颜色
  )
)



# 创建用于森林图的合并表格
# 创建包含标题行的数据框（标题行数值为 NA）
forest_data <- data.frame(
  Vitamin_B2_intake = c(
    "Gastric Cancer",
    "Model 1",
    "Model 2", 
    "Model 3",
    "Colorectal Cancer",
    "Model 1",
    "Model 2",
    "Model 3"
  ),
  N = c(24, "", "", "", 261, "", "", ""),
  Events = c("", 12, 12, 11, "", 109, 109, 104),
  HR = c(NA, 0.117, 0.199, 0.123, NA, 0.680, 0.621, 0.640),  # 标题行为 NA
  CI_lower = c(NA, 0.023, 0.036, 0.012, NA, 0.464, 0.417, 0.421),  # 标题行为 NA
  CI_upper = c(NA, 0.607, 1.103, 1.321, NA, 0.996, 0.924, 0.973),  # 标题行为 NA
  P_value = c("", 0.004, 0.019, 0.006, "", 0.196, 0.129, 0.147),
  Adjustment = c("", "Unadjusted", "+ Gender + Age", "+ Gender + Age +Race +  BMI",
                 "", "Unadjusted", "+ Gender + Age", "+ Gender + Age +Race +  BMI")
)

# 添加用于显示 CI 的文本列
forest_data$`HR (95% CI)` <- ifelse(is.na(forest_data$HR), "",
                                    sprintf("%.3f (%.3f - %.3f)",
                                            forest_data$HR,
                                            forest_data$CI_lower,
                                            forest_data$CI_upper))

# 添加空白列用于森林图（控制宽度）
forest_data$` ` <- paste(rep(" ", 20), collapse = " ")

# # 设置主题
# tm <- forest_theme(base_size = 10,
#                    refline_col = "red",
#                    arrow_type = "closed",
#                    ci_col = "royalblue")

# 绘制森林图（NA 值自动跳过，不显示图形）
p <- forest(forest_data[, c("Vitamin_B2_intake", "N","HR (95% CI)"," ", "P_value", "Adjustment")],
            est = forest_data$HR,
            lower = forest_data$CI_lower,
            upper = forest_data$CI_upper,
            sizes = 1.5,
            ci_column = 4,  # 空白列的位置，用于绘制 CI
            ref_line = 1,
            arrow_lab = c("Low VB2 Better", "High VB2 Better"),
            xlim = c(0, 2),
            ticks_at = c(0, 1,  2),
            theme = tm)
plot(p)


p <- forest(forest_data[, c("Vitamin_B2_intake", "N","HR (95% CI)"," ", "P_value", "Adjustment")],  # 显示的列
            est = res$HR,
            lower = res$CI_lower,
            upper = res$CI_upper,
            sizes = 1.5,
            ci_column = 4,  # CI_text 在第3列
            ref_line = 1,
            xlim = c(0, 2),
            ticks_at = c(0,1,2),
            arrow_lab = c("High VB2 Better", "Low  VB2 Better"),
            footnote = "Cox regression models adjusted for covariates", theme = tm)

plot(p)


ggplot2::ggsave(filename = "forest with subgroups-VB2-High.pdf", plot = p,
                dpi = 300,
                width = 16, height =10, units = "in")

ggplot2::ggsave(filename = "forest with subgroups2-VB2-High.png", plot = p,
                dpi = 300,
                width = 16, height =10, units = "in")

ggplot2::ggsave(filename = "forest with subgroups2-VB2-High.tiff", plot = p,
                dpi = 300,
                width = 16, height =10 , units = "in")
