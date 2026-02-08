# =============================================== #
#                 Compas CI Figure                #
# =============================================== #
rm(list = ls())
library(ggplot2)
library(dplyr)
library(readr)
library(Cairo)
alpha <- 0.95 # 0.9 / 0.95
datanames <- c("African-American","Sex-Age")

# 循环 datanames，为每个绘制一张图
for(dataname in datanames){
  # 读取本地 audit 数据
  audit_data <- read_csv(paste0("data/", dataname, "-Fairness-Audit-", alpha, ".csv"))

  # 添加颜色和线型控制列 (关键修复)
  audit_data <- audit_data %>%
    mutate(
      # 颜色：前4组黑蓝绿红，之后循环
      line_color = c("black", "blue", "green3", "red",
                     rep('black',2),rep('blue',2),rep('green3',2),rep('red',2)),
      # 线型：前4组实线，后续交替点虚线与点线
      line_type = ifelse(row_number() <= 4, "solid",
                         ifelse(row_number() %% 2 == 1, "dotted", "dashed"))
    )

  # 确保分组顺序与数据框行序一致（关键步骤）
  audit_data$group <- factor(audit_data$group, levels = unique(audit_data$group))

  # 画图
  PPV_plot <- ggplot(audit_data, aes(y = group)) +
    geom_segment(aes(x = lb, xend = ub, yend = group,
                 color = line_color, linetype = line_type),
                 lwd = 0.8, alpha = 0.7) +
    geom_point(aes(x = lb, color = line_color), shape = 19, size = 2.5, alpha = 0.9) +
    geom_point(aes(x = ub, color = line_color), shape = 19, size = 2.5, alpha = 0.9) +
    geom_point(aes(x = epsilonG, color = line_color),
               shape = 18, size = 3) +
    scale_color_identity() +
    scale_linetype_identity() +
    geom_vline(xintercept = 0, linetype = "solid",
               color = "gray60", lwd = 0.6) +
    scale_x_continuous(
      limits = c(-0.35, 0.35),
      breaks = seq(-0.3, 0.3, 0.1),
      labels = scales::number_format(accuracy = 0.01),
      expand = c(0, 0)
    ) +
    labs(x = "Postive Predictive Value Disparity", y = NULL,
         title = paste0(dataname, " Subgroups with ", alpha*100, "% Confidence Intervals")) +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA),
      axis.line.x = element_line(color = "black"),
      panel.grid.major.x = element_line(color = "gray90"),
      panel.grid.minor.x = element_blank(),
      axis.text.y = element_text(size = 11, color = "black", margin = margin(r=10)),
      axis.text.x = element_text(size = 10, color = "black"),
      axis.title.x = element_text(size = 12, face = "bold", margin = margin(t=10)),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.margin = unit(c(1,1,1,1), "cm")
    ) +
    scale_y_discrete(limits = rev)

  print(PPV_plot)

  # 保存为 EPS（需要时取消注释）
  # ggsave(
  #   filename = paste0("PPV-", dataname, "-Disparity-", alpha, ".eps"),
  #   plot = PPV_plot,
  #   device = cairo_ps,
  #   width = 8.5, height = 6, units = "in", dpi = 300, bg = "transparent"
  # )
}


