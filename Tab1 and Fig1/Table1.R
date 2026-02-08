#==============================================#
#  合并 Table1 覆盖率（Bootstrap + EL + EEL）     #
#==============================================#
# 用法：source 后调用 table1_merge() 或指定参数
# 读取 Table1-Bootstrap-CP.R 与 Table1-EL-CP.R 生成的数据，按 Bootstrap、EL、EEL 顺序合并

table1_merge <- function(modelname = "Model1", ns = c(2000, 4000, 8000), ms = c(2, 5, 10), a = 0.95) {
  data_dir <- "data"

  # 1. 读取 Bootstrap 数据（来自 Table1-Bootstrap-CP.R）
  df_bootstrap <- NULL
  for (n in ns) {
    for (m in ms) {
      f <- file.path(data_dir, paste0(modelname, "-Boots-", a, "-", n, "-", m, ".csv"))
      if (file.exists(f)) {
        d <- read.csv(f)
        # 列名通常为 alpha, m, n, CP 或 coverage；取 n, m 及覆盖率列
        cp_col <- if ("CP" %in% names(d)) d$CP[1] else if ("coverage" %in% names(d)) d$coverage[1] else NA
        if (!is.na(cp_col)) {
          row <- data.frame(n = d$n[1], m = d$m[1], Bootstrap = cp_col)
          df_bootstrap <- rbind(df_bootstrap, row)
        }
      }
    }
  }

  # 2. 计算 EL、EEL 覆盖率（与 Table1-EL-CP.R 逻辑一致）
  df_el <- NULL
  for (m in ms) {
    for (n in ns) {
      f_el  <- file.path(data_dir, paste0(modelname, "-EL-",  n, "-", m, ".csv"))
      f_eel <- file.path(data_dir, paste0(modelname, "-EEL-", n, "-", m, ".csv"))
      if (file.exists(f_el) && file.exists(f_eel)) {
        EL  <- read.csv(f_el)
        EEL <- read.csv(f_eel)
        nsim <- nrow(EL)
        q1 <- qchisq(a, m)
        Q1 <- matrix(q1, nsim, 1, byrow = TRUE)
        el_cp  <- colSums(EL < Q1) / nsim
        eel_cp <- colSums(EEL < Q1) / nsim
        df_el <- rbind(df_el, data.frame(n = n, m = m, EL = el_cp, EEL = eel_cp))
      }
    }
  }

  if (is.null(df_bootstrap) && is.null(df_el)) {
    message("没有读取到任何数据")
    return(invisible(NULL))
  }

  # 3. 合并：先 Bootstrap，再 EL-EEL → 列为 n, m, Bootstrap, EL, EEL
  if (!is.null(df_bootstrap) && !is.null(df_el)) {
    out <- merge(df_bootstrap, df_el, by = c("n", "m"), sort = FALSE)
  } else if (!is.null(df_bootstrap)) {
    out <- df_bootstrap
    out$EL <- NA
    out$EEL <- NA
  } else {
    out <- df_el
    out$Bootstrap <- NA
    out <- out[, c("n", "m", "Bootstrap", "EL", "EEL")]
  }

  out <- out[order(out$n, out$m), ]
  rownames(out) <- NULL
  print(out)
  invisible(out)
}

# 示例：Model1，ns = 2000,4000,8000，ms = 2,5,10
ns <- c(2000, 4000, 8000)
ms <- c(2, 5, 10)
for (modelname in c("Model1", "Model2")) {
  cat(modelname, "\n")
  table1_merge(modelname = modelname, ns = ns, ms = ms, a = 0.95)
  cat("\n")
}


