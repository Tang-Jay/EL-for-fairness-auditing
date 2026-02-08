#==============================================#
#           Model Generate Data               #
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
# 主程序准备
generate_Data <- function(n, m, modelname, a=0.95, nsim=2000){
  n <- n
  m <- m
  a <- a
  nsim <- nsim
  
  beta0 <- 2
  theta_P <- 0
  
  # 根据modelname设置epsilon
  if(modelname == "Model1"){
    epsilon <- 1
  } else if(modelname == "Model2"){
    j <- 1:m
    epsilon <- (2 * j - 1) / (2 * m)
  } else {
    stop("modelname must be 'Model1' or 'Model2'")
  }
  
  f1 = 0 
  f2 = 0
  EL = c()
  EEL = c()
  for(i in 1:nsim){ 
    # 估计参数beta
    Xi <- runif(n, min = 0, max = 1)  
    
    # 根据modelname设置Yi的生成方式
    if(modelname == "Model1"){
      Yi <- rnorm(n, mean = beta0 * Xi, sd = 1)
    } else if(modelname == "Model2"){
      Yi <- rnorm(n, mean = beta0 * Xi, sd = sqrt(Xi))
    }
    
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
    
    # 设置分组
    S <- generate_one_hot_matrix(X, m)
    L <- matrix(Li, n, m)
    E <-  matrix(epsilon, n, m, byrow = TRUE)
    
    # 估计方程
    g = (L - E ) * S
    z = g
    
    # 计算EL值
    lam = lambdaChen(z)
    el = 2*sum(log(1+t(lam)%*%t(z)))
    if(el < qchisq(a,m)) f1 = f1 + 1
    
    # 计算EEL值
    ez = t(z)
    zbar = rowMeans(ez)
    S0 = (ez-zbar)%*%t(ez-zbar)/(n)
    eel = n*t(zbar) %*% solve(S0) %*% zbar
    if(eel < qchisq(a,m)) f2 = f2 + 1
    
    # 保存结果
    EL = c(EL,el)
    EEL = c(EEL,eel)
    
    # 检验lam
    # aa <- 1+t(lam)%*%t(z)
    # glam <- rowSums(t(z)/matrix(aa,m,n,byrow = TRUE))
    # cat(i, max(abs(glam)),'\n')
  }
  cat(modelname,n, m, a, f1/nsim, f2/nsim,'\n')
  # 保存结果
  write.csv( EL, file = paste0(modelname, '-EL', '-', n, '-', m,'.csv'), row.names = FALSE)
  write.csv( EEL, file = paste0(modelname, '-EEL', '-', n, '-', m,'.csv'), row.names = FALSE)
}
#===============================================
#                   赋值运算            
#===============================================
modelname <- "Model2"  # 可以改为 "Model1" 或 "Model2"
ns = c(1600)
ms = c(1)
cat('Model',' n ','  m', ' a  ', ' EL ',' EEL','\n')
for(m in ms){
  for(n in ns){
    # set.seed(123)
    generate_Data(n, m, modelname, a=0.95, nsim=2000)
  }
}



