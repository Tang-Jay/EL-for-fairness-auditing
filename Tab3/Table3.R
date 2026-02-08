#===============================================#
#              EL p-value of FFR                #
#===============================================#
rm(list = ls())
library(emplik)

# 函数准备
# 经验似然比检验函数
el_test <- function(L, epsilon0) {
  # L: 损失向量
  # epsilon0: 零假设阈值
  
  n <- length(L)
  mu_hat <- mean(L)
  
  # 检验 H0: epsilon <= epsilon0
  if (mu_hat <= epsilon0) {
    # 如果估计值小于等于epsilon0，T = 0
    T <- 0
  } else {
    el_result <- el.test(L, mu = epsilon0)
    T <- el_result$`-2LLR`
  }
  
  # 计算p值
  if (T == 0) {
    p_value <- 1
  } else {
    p_value <- 0.5 * (1 - pchisq(T, df = 1))
  }
  
  return(list(p_value = p_value, T = T, mu_hat = mu_hat))
}

# BH程序
bh_procedure <- function(p_values, alpha = 0.05) {
  m <- length(p_values)
  sorted_p <- sort(p_values)
  k <- max(which(sorted_p <= (1:m) * alpha / m), 0)
  
  if (k == 0) {
    return(rep(FALSE, m))
  } else {
    threshold <- sorted_p[k]
    return(p_values <= threshold)
  }
}

# 单次模拟实验
single_experiment <- function(n = 1000, beta0 = 1, tau, epsilon0 = 0.05) {
  # 生成数据
  X <- runif(n)
  Y <- beta0 * X + rnorm(n)
  
  # 模型参数
  beta1 <- beta0 - 2 * tau
  f_X <- beta1 * X
  
  # 计算损失
  L <- Y - f_X
  
  # 分组
  group1 <- which(X < 0.5)
  group2 <- which(X >= 0.5)
  
  # 真实性能差异
  true_epsilon1 <- 2 * tau * mean(X[group1])
  true_epsilon2 <- 2 * tau * mean(X[group2])
  true_epsilons <- c(true_epsilon1, true_epsilon2)
  
  # 零假设：epsilon_G <= epsilon0
  null_true <- c(true_epsilon1 <= epsilon0, true_epsilon2 <= epsilon0)
  
  # 经验似然检验
  p_values <- numeric(2)
  T_stats <- numeric(2)
  mu_hats <- numeric(2)
  
  # 组1检验
  result1 <- el_test(L[group1], epsilon0)
  p_values[1] <- result1$p_value
  T_stats[1] <- result1$T
  mu_hats[1] <- result1$mu_hat
  
  # 组2检验
  result2 <- el_test(L[group2], epsilon0)
  p_values[2] <- result2$p_value
  T_stats[2] <- result2$T
  mu_hats[2] <- result2$mu_hat
  
  # BH程序
  rejections <- bh_procedure(p_values, alpha = 0.05)
  
  return(list(
    p_values = p_values,
    rejections = rejections,
    null_true = null_true,
    true_epsilons = true_epsilons,
    mu_hats = mu_hats,
    T_stats = T_stats
  ))
}

# 主模拟函数
monte_carlo_simulation <- function(nsim = 1000, taus = seq(0, 0.2, 0.02), 
                                   n = 1000, beta0 = 1, epsilon0 = 0.05) {
  results <- list()
  
  for (tau in taus) {
    cat("Simulating tau =", tau, "\n")
    
    ffr <- numeric(nsim)
    power <- numeric(nsim)
    rejection_rates <- matrix(0, nsim, 2)
    
    for (i in 1:nsim) {
      result <- single_experiment(n = n, beta0 = beta0, tau = tau, epsilon0 = epsilon0)
      
      # 计算FFR
      false_rejections <- sum(result$rejections & result$null_true)
      total_rejections <- sum(result$rejections)
      
      if (total_rejections > 0) {
        ffr[i] <- false_rejections / total_rejections
      } else {
        ffr[i] <- 0
      }
      
      # 计算power (只在备择假设为真时)
      alternative_true <- !result$null_true
      if (any(alternative_true)) {
        correct_rejections <- sum(result$rejections & alternative_true)
        power[i] <- correct_rejections / sum(alternative_true)
      } else {
        power[i] <- 0
      }
      
      rejection_rates[i, ] <- result$rejections
    }
    
    # 真实性能差异
    true_epsilon1 <- 2 * tau * 0.25  # E[X|X<0.5] = 0.25
    true_epsilon2 <- 2 * tau * 0.75  # E[X|X>=0.5] = 0.75
    
    results[[as.character(tau)]] <- list(
      tau = tau,
      true_epsilons = c(true_epsilon1, true_epsilon2),
      ffr = mean(ffr),
      ffr_sd = sd(ffr),
      power = mean(power),
      power_sd = sd(power),
      rejection_rate_group1 = mean(rejection_rates[, 1]),
      rejection_rate_group2 = mean(rejection_rates[, 2])
    )
  }
  
  return(results)
}

# 图一：绘制FFR和Power
plot_ffr_power <- function(results_df) {
  par(oma=c(0,0,0,0),mar=c(3,3.5,1.5,1))
  
  # 设置坐标轴范围
  tau_range <- range(results_df$tau)
  y_range <- c(0, 1)
  
  # 创建空图
  plot(results_df$tau, results_df$ffr, type = "n", 
       xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n",
       ylim = y_range, xlim = tau_range,
       xlab = "", ylab = "", main = "ELBH Procedure")
  
  # 添加坐标轴
  axis(1, at = seq(round(tau_range[1], 1), round(tau_range[2], 1), 0.1),
       mgp = c(4, 0  , 0), tcl = 0.5, font = 2, lwd = 2, las = 1)
  axis(2, at = seq(0, 1, 0.2),
       mgp = c(4, 0.2, 0), tcl = 0.5, font = 2, lwd = 2, las = 1)
  axis(3,labels=F,tcl=0,lwd=2)
  axis(4,labels=F,tcl=0,lwd=2)
  # 添加坐标轴标签
  title(xlab = expression(tau), line = 2, font.lab = 2, cex.lab = 1.2)
  title(ylab = "Rate", line = 2.5, font.lab = 2, cex.lab = 1.2)
  
  # 添加参考线
  abline(h = 0.05, lwd = 2, col = "gray", lty = 2)
  
  # 绘制FFR和Power曲线
  lines(results_df$tau, results_df$ffr,   col = "blue",        lwd = 2, lty = 1)
  lines(results_df$tau, results_df$power, col = "aquamarine4", lwd = 2, lty = 2)
  # 添加数据点
  points(results_df$tau, results_df$ffr,   col = "blue",        pch = 16, cex = 1.2)
  points(results_df$tau, results_df$power, col = "aquamarine4", pch = 17, cex = 1.2)
  
  # 添加图例
  legend(x = taus[1]+0.02, y = 0.95, legend = c("Power","FFR"), 
         col = c("aquamarine4", "blue"), lty = c(2, 1), pch = c(17,16), lwd = 1.2,
         bg = "white", cex = 1)
  
}

# 图二：绘制拒绝率
plot_rejection_rates_base <- function(results_df) {
  par(oma=c(0,0,0,0),mar=c(3,3.5,1.5,1))
  
  # 设置坐标轴范围
  tau_range <- range(results_df$tau)
  y_range <- c(0, 1)
  
  # 创建空图
  plot(results_df$tau, results_df$rejection_rate_group1, type = "n", 
       xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n",
       ylim = y_range, xlim = tau_range,
       xlab = "", ylab = "", main = "Rejection Rates for Two Groups")
  
  # 添加坐标轴
  axis(1, at = seq(round(tau_range[1], 1), round(tau_range[2], 1), 0.1),
       mgp = c(4, 0  , 0), tcl = 0.5, font = 2, lwd = 2, las = 1)
  axis(2, at = seq(0, 1, 0.2),
       mgp = c(4, 0.2, 0), tcl = 0.5, font = 2, lwd = 2, las = 1)
  axis(3,labels=F,tcl=0,lwd=2)
  axis(4,labels=F,tcl=0,lwd=2)
  # 添加坐标轴标签
  title(xlab = expression(tau), line = 2, font.lab = 2, cex.lab = 1.2)
  title(ylab = "Rejection Rate", line = 2.5, font.lab = 2, cex.lab = 1.2)
  
  # 绘制拒绝率曲线
  lines(results_df$tau, results_df$rejection_rate_group1, col = "red", lwd = 2, lty = 1)
  lines(results_df$tau, results_df$rejection_rate_group2, col = "blue", lwd = 2, lty = 2)
  
  # 添加数据点
  points(results_df$tau, results_df$rejection_rate_group1, col = "red",  pch = 16, cex = 0.6)
  points(results_df$tau, results_df$rejection_rate_group2, col = "blue", pch = 17, cex = 0.6)
  
  # 标记备择假设为真的情况
  # Group 1: 当 true_epsilon1 > 0.05 时
  alt_true_group1 <- which(results_df$true_epsilon1 > 0.05)
  if (length(alt_true_group1) > 0) {
    points(results_df$tau[alt_true_group1], 
           results_df$rejection_rate_group1[alt_true_group1],
           col = "red", pch = 1, cex = 1.2, lwd = 2)
  }
  
  # Group 2: 当 true_epsilon2 > 0.05 时
  alt_true_group2 <- which(results_df$true_epsilon2 > 0.05)
  if (length(alt_true_group2) > 0) {
    points(results_df$tau[alt_true_group2], 
           results_df$rejection_rate_group2[alt_true_group2],
           col = "blue", pch = 2, cex = 1.2, lwd = 2)
  }
  
  # 添加图例
  legend(x = taus[1]+0.02, y = 0.98, 
         legend = c("Group 1", "Group 2", "Alt Hyp True (Group 1)", "Alt Hyp True (Group 2)"), 
         col = c("red", "blue", "red", "blue"), 
         lty = c(1, 2, NA, NA), 
         pch = c(16, 17, 1, 2),
         lwd = c(2, 2, 2, 2),
         bg = "white", cex = 0.8)
}

#===============================================#
# 运行模拟
set.seed(123)
m = 2
nsim = 2000
n = 1000
taus = seq(-0.2, 0.6, 0.05)
results <- monte_carlo_simulation(nsim, taus, n)

# 结果整理
results_df <- do.call(rbind, lapply(results, function(x) {
  data.frame(
    tau = x$tau,
    true_epsilon1 = x$true_epsilons[1],
    true_epsilon2 = x$true_epsilons[2],
    ffr = x$ffr,
    power = x$power,
    rejection_rate_group1 = x$rejection_rate_group1,
    rejection_rate_group2 = x$rejection_rate_group2
  )
}))

# 可视化结果
plot_ffr_power(results_df)
plot_rejection_rates_base(results_df)

# # 保存图一
# setEPS()
# postscript(paste0('Model3-ELBH-FFR-',m,'.eps'), width=6, height=5)
# plot_ffr_power(results_df)
# dev.off()
# # 保存图二
# setEPS()
# postscript(paste0('Model3-ELBH-Rejection-',m,'.eps'), width=6, height=5)
# plot_rejection_rates_base(results_df)
# dev.off()

# 输出关键结果
cat("ELBH Procedure FFR Control Results:\n")
cat("===================================\n")
print(results_df[, c("tau", "true_epsilon1", "true_epsilon2", "ffr", "power")])

# 检查FFR控制
ffr_controlled <- all(results_df$ffr <= 0.05 ) 
cat("\nFFR Control Check:", ifelse(ffr_controlled, "PASS", "FAIL"), "\n")


