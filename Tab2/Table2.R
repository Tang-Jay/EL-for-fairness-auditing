#==============================================#
#  合并 Model1 运行时间表（EL-EEL + Candes）   #
#==============================================#
# 用法：在 R 中 source 后调用 table2_merge(n)，或通过命令行传参
# 示例：Rscript Table2.R 400

table2_merge <- function(n) {
  data_dir <- "data"
  file_el   <- file.path(data_dir, paste0("Time-Model1-EL-EEL-", n, ".csv"))
  file_candes <- file.path(data_dir, paste0("Time-Model1-Candes-", n, ".csv"))

  if (!file.exists(file_el)) stop("文件不存在: ", file_el)
  if (!file.exists(file_candes)) stop("文件不存在: ", file_candes)

  df_el    <- read.csv(file_el)
  df_candes <- read.csv(file_candes)

  # 合并顺序：先 Candes，再 EL-EEL → 列为 n, m, Candes, EL, EEL
  out <- merge(df_candes, df_el, by = c("n", "m"), sort = FALSE)
  out <- out[order(out$n, out$m), ]

  print(out)
  invisible(out)
}

# 命令行传参：Rscript Table2.R 400
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) {
  n <- as.numeric(args[1])
  if (!is.na(n)) table2_merge(n) else message("请提供有效的 n，例如: Rscript Table2.R 400")
}


# 示例：循环 n = 400, 2000, 4000 输出合并结果
ns <- c(400, 2000, 4000)
for (n in ns) {
  cat("n =", n, "\n")
  table2_merge(n)
  cat("\n")
}
