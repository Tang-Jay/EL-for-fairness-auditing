#==============================================#
#           Model1 Running time                #
#==============================================#
# 函数准备
rm(list = ls()) 
source('GlambdaChen.R')
generate_one_hot_matrix <- function(X, m){
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
#===============================================
#                   赋值运算            
#===============================================
nsim = 2000
a = 0.95
for(n in c(400,2000,4000)){
  ms = c(1,2,5,10)
  cat('n ',' m ', ' EL ',' EEL','\n')
  
  # 保存运行时间
  TT = matrix(NA,length(ms),4)
  colnames(TT) = c('n','m', 'EL','EEL')
  
  j = 0
  for(m in ms){
    j = j + 1
    set.seed(123)
    beta0 <- 2
    epsilon <- 1
    theta_P <- 0
    f1 = 0 
    f2 = 0
    T1 = 0 
    T2 = 0
    EL = c()
    EEL = c()
    for(i in 1:nsim){ 
      # 估计参数beta
      Xi <- runif(n, min = 0, max = 1)           
      Yi <- rnorm(n, mean = beta0 * Xi, sd = 1)
      X <- matrix(Xi, ncol = 1)  # 转换为 n×1 矩阵
      Y <- matrix(Yi, ncol = 1)  # 转换为 n×1 矩阵
      XtX <- t(X) %*% X          # X转置X
      XtX_inv <- solve(XtX)      # (X转置X)的逆
      XtY <- t(X) %*% Yi         # X转置Y
      beta_h <- XtX_inv %*% XtY 
      beta_h <- as.numeric(beta_h)  # 将矩阵转换为标量
      
      # 计算Loss
      Y_h <- beta_h * X
      Y_h <- beta_h * X
      Li <- (Y - Y_h)^2
      L <- matrix(Li, ncol = 1)
      
      # 设置分组
      S <- generate_one_hot_matrix(X, m)
      
      # 估计方程
      g = ( c(L) - epsilon ) * S
      z = g
      
      # 计算EL值
      t1<-lubridate::now()
      lam = lambdaChen(z)
      el = 2*sum(log(1+t(lam)%*%t(z)))
      if(el < qchisq(a,m)) f1 = f1 + 1
      t2<-lubridate::now()
      T1=T1+t2-t1
      
      # 计算EEL值
      et1<-lubridate::now()
      ez = t(z)
      zbar = rowMeans(ez)
      S0 = (ez-zbar)%*%t(ez-zbar)/(n)
      eel = n*t(zbar) %*% solve(S0) %*% zbar
      if(eel < qchisq(a,m)) f2 = f2 + 1
      et2<-lubridate::now()
      T2=T2+et2-et1
    }
    
    print(T1)
    print(T2)
    
    # 统一转换为分钟
    T1 = display_in_minutes(T1)
    T2 = display_in_minutes(T2)
    
    # 保存数据
    TT[j,] = c(n,m,T1,T2)
    
    # 展示数据
    cat(n, m, f1/nsim, f2/nsim,'\n')
    cat(n, m, T1, T2,'\n\n')
    
  }
  
  print(TT)
  # write.csv(TT, file = paste0('Time-Model1-EL-EEL-', n,'.csv'), row.names = FALSE)
}


