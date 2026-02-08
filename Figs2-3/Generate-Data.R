#==============================================#
#           Compare EL and T-test Power        #  
#==============================================#
# 函数准备
rm(list = ls()) 
source('GlambdaChen.R')
generate_one_hot_matrix <- function(X, m) {
  breaks <- seq(0, 1, length.out = m + 1)
  group_indices <- findInterval(X, breaks, rightmost.closed = TRUE)
  one_hot_matrix <- matrix(0, nrow = length(X), ncol = m)
  one_hot_matrix[cbind(1:length(X), group_indices)] <- 1
  return(one_hot_matrix)
}

# 区间假设检验函数（T-test）
interval_test <- function(sample_data, mu1, mu2, alpha = 0.05){
  test_upper <- t.test(sample_data, mu = mu2, alternative = "greater", conf.level = 1 - alpha)
  test_lower <- t.test(sample_data, mu = mu1, alternative = "less", conf.level = 1 - alpha)
  reject_H0 <- (test_upper$p.value < alpha) || (test_lower$p.value < alpha)
  
  return(list(
    reject = reject_H0,
    p_upper = test_upper$p.value,
    p_lower = test_lower$p.value,
    test_upper = test_upper,
    test_lower = test_lower
  ))
}

#==============================================#
#           Power 计算函数                     #
#==============================================#
calculate_power_EL_Ttest <- function(n = 100,
                                     errorType = "error1",
                                     epsilon0 = 0.05,
                                     alpha = 0.05,
                                     nsim = 1000,
                                     m = 1,
                                     tau_l = -0.5,
                                     tau_u = 0.5,
                                     tau_step = 0.05,
                                     save_csv = FALSE){
  
  # 参数设置（beta0 在函数内部固定为 2）
  beta0 <- 2
  taus <- seq(tau_l, tau_u, tau_step)
  l <- length(taus)
  
  # 临界值
  cut1 <- qchisq(1 - alpha, 1)      # 用于双侧检验
  cut2 <- qchisq(1 - 2*alpha, 1)    # 用于单侧检验
  
  # 初始化power向量（EL方法）
  PowerEL_T04 <- numeric(l)
  PowerEL_T14 <- numeric(l)
  PowerEL_T24 <- numeric(l)
  PowerEL_T34 <- numeric(l)
  
  # 初始化power向量（T-test方法）
  PowerTtest_T04 <- numeric(l)
  PowerTtest_T14 <- numeric(l)
  PowerTtest_T24 <- numeric(l)
  PowerTtest_T34 <- numeric(l)
  
  cat("开始计算 power...\n")
  cat("参数: n =", n, ", nsim =", nsim, ", errorType =", errorType, ", m =", m, "\n\n")
  
  # 对每个tau值计算power
  for(j in 1:l){
    tau <- taus[j]
    beta1 <- beta0 - 2*tau
    epsilon_star <- tau
    
    # 初始化计数器（EL方法）
    f04 <- 0
    f14 <- 0
    f24 <- 0
    f34 <- 0
    
    # 初始化计数器（T-test方法）
    t04_count <- 0
    t14_count <- 0
    t24_count <- 0
    t34_count <- 0
    
    # 进行 nsim 次模拟
    for(i in 1:nsim){
      # 生成样本
      Xi <- runif(n, min = 0, max = 1)
      
      # 根据 errorType 选择误差分布
      if(errorType == "error1"){
        v <- rnorm(n, mean = 0, sd = 1)  # N(0, 1)
      } else if(errorType == "error2"){
        v <- rexp(n, rate = 1) - 1  # Exp(1) - 1
      } else {
        v <- runif(n, min = -1, max = 1 + tau)  # U(-1, 1+tau)
      }
      
      Yi <- beta0 * Xi + v
      X <- matrix(Xi, ncol = 1)
      Y <- matrix(Yi, ncol = 1)
      Y_h <- beta1 * X
      Li <- Y - Y_h
      
      # 设置分组
      S <- generate_one_hot_matrix(X, m)
      L <- matrix(Li, n, m)
      
      # 构造估计函数
      E0 <- matrix(epsilon0, n, m, byrow = TRUE)
      
      #==============================================#
      # T04：双侧检验（检验 H0: μ = μ_0）
      #==============================================#
      # EL方法
      z <- (L - E0) * S
      lam <- lambdaChen(z)
      T04_EL <- 2*sum(log(1 + t(lam) %*% t(z)))
      if(T04_EL > cut1) f04 <- f04 + 1
      
      # T-test方法
      t_test_result <- t.test(Li, mu = epsilon0, alternative = "two.sided")
      if(t_test_result$p.value < alpha) t04_count <- t04_count + 1
      
      #==============================================#
      # T14：检验 H0: μ ≥ μ_0 vs H1: μ < μ_0
      #==============================================#
      # EL方法
      if(mean(L) <= epsilon0){
        z <- (L - E0) * S
        lam <- lambdaChen(z)
        T14_EL <- 2*sum(log(1 + t(lam) %*% t(z)))
      } else {
        T14_EL <- 0
      }
      if(T14_EL > cut2) f14 <- f14 + 1
      
      # T-test方法
      t_test_result <- t.test(Li, mu = epsilon0, alternative = "less")
      if(t_test_result$p.value < alpha) t14_count <- t14_count + 1
      
      #==============================================#
      # T24：检验 H0: μ ≤ μ_0 vs H1: μ > μ_0
      #==============================================#
      # EL方法
      if(mean(L) >= epsilon0){
        z <- (L - E0) * S
        lam <- lambdaChen(z)
        T24_EL <- 2*sum(log(1 + t(lam) %*% t(z)))
      } else {
        T24_EL <- 0
      }
      if(T24_EL > cut2) f24 <- f24 + 1
      
      # T-test方法
      t_test_result <- t.test(Li, mu = epsilon0, alternative = "greater")
      if(t_test_result$p.value < alpha) t24_count <- t24_count + 1
      
      #==============================================#
      # T34：检验 H0: μ₁ ≤ μ ≤ μ₂ vs H1: μ > μ₂ 或 μ < μ₁
      #==============================================#
      # EL方法
      if(mean(L) <= -epsilon0){
        z <- (L + E0) * S  # L - (-E0)
        lam <- lambdaChen(z)
        T34_EL <- 2*sum(log(1 + t(lam) %*% t(z)))
      } else if(mean(L) >= epsilon0){
        z <- (L - E0) * S
        lam <- lambdaChen(z)
        T34_EL <- 2*sum(log(1 + t(lam) %*% t(z)))
      } else {
        T34_EL <- 0
      }
      if(T34_EL >= cut2) f34 <- f34 + 1
      
      # T-test方法
      result <- interval_test(L, -epsilon0, epsilon0, alpha = alpha)
      if(result$reject) t34_count <- t34_count + 1
    }
    
    # 计算power（拒绝H0的比例）
    PowerEL_T04[j] <- f04 / nsim
    PowerEL_T14[j] <- f14 / nsim
    PowerEL_T24[j] <- f24 / nsim
    PowerEL_T34[j] <- f34 / nsim
    
    PowerTtest_T04[j] <- t04_count / nsim
    PowerTtest_T14[j] <- t14_count / nsim
    PowerTtest_T24[j] <- t24_count / nsim
    PowerTtest_T34[j] <- t34_count / nsim
    
    # 显示进度
    if(j %% 5 == 0 || j == l){
      cat(sprintf("tau = %.2f: EL(T04, T14, T24, T34) = (%.4f, %.4f, %.4f, %.4f), ",
                  tau, PowerEL_T04[j], PowerEL_T14[j], PowerEL_T24[j], PowerEL_T34[j]))
      cat(sprintf("Ttest(T04, T14, T24, T34) = (%.4f, %.4f, %.4f, %.4f)\n",
                  PowerTtest_T04[j], PowerTtest_T14[j], PowerTtest_T24[j], PowerTtest_T34[j]))
    }
  }
  
  # 保存CSV数据（循环保存所有检验类型）
  if(save_csv){
    test_types <- c("T04", "T14", "T24", "T34")
    power_el_list <- list(T04 = PowerEL_T04, T14 = PowerEL_T14, T24 = PowerEL_T24, T34 = PowerEL_T34)
    power_ttest_list <- list(T04 = PowerTtest_T04, T14 = PowerTtest_T14, T24 = PowerTtest_T24, T34 = PowerTtest_T34)
    
    for(testType in test_types){
      data_df <- data.frame(
        tau = taus,
        EL = power_el_list[[testType]],
        Ttest = power_ttest_list[[testType]]
      )
      
      filename <- paste0("PowerCompare-", testType, "-", errorType, "-", n, ".csv")
      write.csv(data_df, file = filename, row.names = FALSE)
      cat("已保存:", filename, "\n")
    }
  }


  return(list(
    taus = taus,
    EL = list(T04 = PowerEL_T04, T14 = PowerEL_T14, T24 = PowerEL_T24, T34 = PowerEL_T34),
    Ttest = list(T04 = PowerTtest_T04, T14 = PowerTtest_T14, T24 = PowerTtest_T24, T34 = PowerTtest_T34)
  ))
}

#==============================================#
#           多样本量对比绘图函数               #
# 与 Figs2-3.R 同款图，数据来自内存  #
#==============================================#
plot_multi_n_compare <- function(power_results,
                                 ns,
                                 errorType,
                                 testTypes = c("T04", "T14", "T24", "T34"),
                                 save_eps = FALSE){
  # power_results: list, names = as.character(n)，每个元素为 calculate_power_EL_Ttest 的返回值
  for(testType in testTypes){
    if(save_eps){
      setEPS()
      postscript(paste0('Compare-Multi-', testType, '-', errorType, '.eps'),
                 width = 5.5, height = 5)
    }
    taus <- power_results[[as.character(ns[1])]]$taus
    tau_l <- min(taus)
    tau_u <- max(taus)
    par(oma = c(0, 0, 0, 0), mar = c(2, 3.2, 2, 1))
    plot(taus, taus,
         type = "l",
         xaxs = "i",
         yaxs = "i",
         ylim = c(0, 1),
         xlab = "",
         ylab = "",
         col = "white",
         xaxt = "n",
         yaxt = "n")
    axis(1, seq(tau_l, tau_u, 0.1), mgp = c(4, 0, 0), tcl = 0.5, font = 2, lwd = 2, las = 1)
    axis(2, seq(0, 1, 0.1), mgp = c(4, 0.2, 0), tcl = 0.5, font = 2, lwd = 2, las = 1)
    axis(3, seq(tau_l, tau_u, 0.1), tck = 0.01, labels = F, tcl = 0, lwd = 2)
    axis(4, seq(0, 1, 0.1), labels = F, tcl = 0, lwd = 2)
    hypothesis_map <- list(
      "T04" = expression(paste(H[0], ": ", epsilon[G] == epsilon[0], " vs ", H[1], ": ", epsilon[G] != epsilon[0])),
      "T14" = expression(paste(H[0], ": ", epsilon[G] >= epsilon[0], " vs ", H[1], ": ", epsilon[G] < epsilon[0])),
      "T24" = expression(paste(H[0], ": ", epsilon[G] <= epsilon[0], " vs ", H[1], ": ", epsilon[G] > epsilon[0])),
      "T34" = expression(paste(H[0], ": ", epsilon[1], " < ", epsilon[G], " < ", epsilon[2], " vs ", H[1], ": ", epsilon[G] <= epsilon[1], " or ", epsilon[G] >= epsilon[2]))
    )
    hypothesis_text <- if(testType %in% names(hypothesis_map)) hypothesis_map[[testType]] else paste0("H? vs H", testType)
    title(main = hypothesis_text, cex.main = 1.7)
    title(xlab = "true disparity", line = 1, col.lab = 1, font.lab = 2, cex.lab = 1)
    title(ylab = "power", line = 2, col.lab = 1, font.lab = 2, cex.lab = 1)
    abline(h = 0.05, lwd = 1.5, col = "gray", lty = 2)
    abline(v = 0, lwd = 1.5, col = "gray", lty = 2)
    color_palette_el <- c("coral", "red1", "red3", "red2")
    color_palette_ttest <- c("dodgerblue", "blue1", "blue3", "blue4")
    colors_el <- setNames(color_palette_el[1:length(ns)], as.character(ns))
    colors_ttest <- setNames(color_palette_ttest[1:length(ns)], as.character(ns))
    lty_el <- setNames(rep(1, length(ns)), as.character(ns))
    lty_ttest <- setNames(rep(2, length(ns)), as.character(ns))
    legend_labels <- c()
    legend_cols <- c()
    legend_lty <- c()
    for(n in ns){
      n_str <- as.character(n)
      if(!is.null(power_results[[n_str]]$EL[[testType]])){
        lines(taus, power_results[[n_str]]$EL[[testType]],
              col = colors_el[n_str], lty = lty_el[n_str], lwd = 2)
        legend_labels <- c(legend_labels, paste0("EL(n=", n, ")"))
        legend_cols <- c(legend_cols, colors_el[n_str])
        legend_lty <- c(legend_lty, lty_el[n_str])
      }
    }
    for(n in ns){
      n_str <- as.character(n)
      if(!is.null(power_results[[n_str]]$Ttest[[testType]])){
        lines(taus, power_results[[n_str]]$Ttest[[testType]],
              col = colors_ttest[n_str], lty = lty_ttest[n_str], lwd = 2)
        legend_labels <- c(legend_labels, paste0("T-test(n=", n, ")"))
        legend_cols <- c(legend_cols, colors_ttest[n_str])
        legend_lty <- c(legend_lty, lty_ttest[n_str])
      }
    }
    legend(x = tau_l + 0.03, y = 0.3, inset = 0.02, cex = 0.6,
           legend = legend_labels,
           lty = legend_lty,
           col = legend_cols,
           lwd = 2)
    if(save_eps) dev.off()
  }
}

#==============================================#
#           使用示例                           #
#==============================================#
# 全局参数设置（多样本量对比图）
set.seed(123)
nsim <- 2000   # 我们模拟采用 2000, 数据结果储存在 /data 中
epsilon0 <- 0.05
ns <- c(20, 50, 500)

# 与 Figs2-3 同款：多样本量对比图（数据由本脚本模拟生成，不读 CSV）
for(errorType in c("error1", "error2")){
  power_results <- list()
  for(n in ns){
    power_results[[as.character(n)]] <- calculate_power_EL_Ttest(
      n = n,
      errorType = errorType,
      epsilon0 = epsilon0,
      nsim = nsim,
      save_csv = FALSE
    )
  }
  plot_multi_n_compare(power_results,
                      ns = ns,
                      errorType = errorType,
                      testTypes = c("T04", "T14", "T24", "T34"),
                      save_eps = FALSE)
}



