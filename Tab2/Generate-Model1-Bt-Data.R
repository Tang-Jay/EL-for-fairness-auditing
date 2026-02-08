#==============================================#
#           Homoskedastic Linear Model         #  
#==============================================#
# 函数准备
rm(list = ls())
generate_one_hot_matrix <- function(X, m) {
  breaks <- seq(0, 1, length.out = m + 1)
  group_indices <- findInterval(X, breaks, rightmost.closed = TRUE)
  one_hot_matrix <- matrix(0, nrow = length(X), ncol = m)
  # 快速定位需要设为 1 的位置，避免循环
  one_hot_matrix[cbind(1:length(X), group_indices)] <- 1
  return(one_hot_matrix)
}
# 自定义函数统一显示为分钟
display_in_minutes <- function(time_diff) {
  minutes <- as.numeric(time_diff, units = "mins")
  return(round(minutes,4))
}
#==============================================#
set.seed(123)
nsim <- 2000
for(n in c(400,2000,4000)){
  B <- 500
  alpha <- 0.1
  ms <- c(1,2,5,10)
  # 保存运行时间
  TT = matrix(NA,length(ms),3)
  cat('n ',' m ', 'B','1-a','CP','Time','\n')
  colnames(TT) = c('n','m', 'Candes')
  j = 0
  for(m in ms){
    j = j + 1
    counter <- 0
    T1 = 0
    for(i in 1:nsim){
      X <- runif(n, min = 0, max = 1)
      Y <- 2 * X + rnorm(n)  
      D <- data.frame(X = X, Y = Y)
      
      # 训练OLS模型
      model <- lm(Y ~ X - 1, data = D)  
      D$Y_hat <- predict(model)
      D$L <- (D$Y - D$Y_hat)^2 
      D$G <- generate_one_hot_matrix(X, m)
      L <- D$L
      
      # 步骤1：预计算逻辑矩阵
      M <- (D$G == 1)
      
      # 步骤2：计算原始数据集上的统计量
      Pn_G <- colMeans(M)  # 各子总体比例
      epsilon_hat <- colSums(D$L * D$G)/colSums(M) # 各子总体平均损失
      
      t1<-lubridate::now()
      # 步骤3：自助抽样过程
      t_b_vec <- numeric(B)  # 存储t(b)结果
      for (b in 1:B) {
        
        # 3a. 有放回抽样
        idx <- sample(1:n, size = n, replace = TRUE)
        D_star <- D[idx, ]
        L_star <- L[idx]
        M_star <- M[idx, ]
        
        # 3b. 计算每个子总体的项
        if( m == 1){
          num_G_star <- sum(M_star)
        } else {
          num_G_star <- colSums(M_star)
        }
        
        P_b_star <- colSums(D_star$G)/n
        epsilon_b_star <- colSums(D_star$L * D_star$G)/num_G_star
        term_vals <- Pn_G * P_b_star * (epsilon_b_star -  epsilon_hat)
        
        # 3c. 取最大值作为t(b)
        t_b_vec[b] <- max(term_vals, na.rm = TRUE)
      }
      
      # 步骤4：计算(1-alpha)分位数
      t_sorted <- sort(t_b_vec)
      k <- ceiling(B * (1 - alpha))
      t_star <- t_sorted[k]
      
      # 步骤5：计算每个子群体的lb(G)
      lb_G <- epsilon_hat - t_star / (Pn_G^2)
      
      # 步骤6：计算每个子群体的ub(G)
      ub_G <- epsilon_hat + t_star / (Pn_G^2)
      
      # 创建结果数据框
      results <- data.frame(
        Group = paste0("G", seq_len(m)),
        Pn_G = Pn_G,
        epsilon_hat = epsilon_hat,
        lb_G = lb_G,
        ub_G = ub_G
      )
      
      # 在results数据框中添加新列，判断1是否在置信区间内
      contains <- as.integer(all((1 >= results$lb_G) & (1 <= results$ub_G)))
      counter <- counter + contains  # 如果contains为1则加1，否则加0
      
      t2<-lubridate::now()
      T1 = T1 + t2 - t1
      
    }
    
    # 统一转换为分钟
    T1 = display_in_minutes(T1)
    
    # 输出最终结果
    cat(n,m,B,1-alpha,counter/nsim,T1,'\n')
    
    # 保存数据
    TT[j,] = c(n,m,T1)
  }
  
  print(TT)
  # write.csv(TT, file = paste0('Time-Model1-Candes-',n,'.csv'), row.names = FALSE)
}


