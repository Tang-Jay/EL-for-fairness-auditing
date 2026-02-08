#==============================================#
#           Bootstrap Method CP                #
#==============================================#
# 函数准备
rm(list = ls())
generate_one_hot_matrix <- function(X, m) {
  breaks <- seq(0, 1, length.out = m + 1)
  group_indices <- findInterval(X, breaks, rightmost.closed = TRUE)
  one_hot_matrix <- matrix(0, nrow = length(X), ncol = m)
  one_hot_matrix[cbind(1:length(X), group_indices)] <- 1
  return(one_hot_matrix)
}

# 主程序：单组 (n, m) 的 Bootstrap CP
generate_Boots_CP <- function(n, m, modelname, nsim = 2000, B = 500, alpha = 0.05) {
  a <- alpha
  counter <- 0

  for (i in 1:nsim) {
    X <- runif(n, min = 0, max = 1)

    if (modelname == "Model1") {
      Y <- rnorm(n, mean = 2 * X, sd = 1)
      epsilonG <- rep(1, m)
    } else if (modelname == "Model2") {
      Y <- rnorm(n, mean = 2 * X, sd = sqrt(X))
      j <- 1:m
      epsilonG <- (2 * j - 1) / (2 * m)
    } else {
      stop("modelname must be 'Model1' or 'Model2'")
    }

    D <- data.frame(X = X, Y = Y)
    model <- lm(Y ~ X - 1, data = D)
    D$Y_hat <- predict(model)
    D$L <- (D$Y - D$Y_hat)^2
    D$G <- generate_one_hot_matrix(X, m)
    L <- D$L

    M <- (D$G == 1)
    Pn_G <- colMeans(M)
    epsilon_hat <- colSums(D$L * D$G) / colSums(M)

    t_b_vec <- numeric(B)
    for (b in 1:B) {
      idx <- sample(1:n, size = n, replace = TRUE)
      D_star <- D[idx, ]
      L_star <- L[idx]
      M_star <- M[idx, ]

      if (m == 1) {
        num_G_star <- sum(M_star)
      } else {
        num_G_star <- colSums(M_star)
      }

      P_b_star <- colSums(D_star$G) / n
      epsilon_b_star <- colSums(D_star$L * D_star$G) / num_G_star
      term_vals <- Pn_G * P_b_star * (epsilon_b_star - epsilon_hat)
      t_b_vec[b] <- max(term_vals, na.rm = TRUE)
    }

    t_sorted <- sort(t_b_vec)
    k <- ceiling(B * (1 - alpha))
    t_star <- t_sorted[k]
    lb_G <- epsilon_hat - t_star / (Pn_G^2)
    ub_G <- epsilon_hat + t_star / (Pn_G^2)

    contains <- as.integer(all((epsilonG >= lb_G) & (epsilonG <= ub_G)))
    counter <- counter + contains
  }

  cp <- counter / nsim
  cat(modelname, 1 - alpha, m, n, cp, "\n")

  result_single <- data.frame(alpha = alpha, m = m, n = n, CP = cp)
  write.csv(result_single, file = paste0(modelname, "-Boots-", 1 - a, "-", n, "-", m, ".csv"), row.names = FALSE)
  invisible(result_single)
}

#===============================================
#                   赋值运算
#===============================================
modelname <- "Model1"  # "Model1" 或 "Model2"
nsim <- 2000
B <- 500
alpha <- 0.05
ns <- c(2000, 4000, 8000)
ms <- c(2, 5, 10)

cat("Model", " 1-a ", "m", " n ", " CP", "\n")
for (n in ns) {
  for (m in ms) {
    generate_Boots_CP(n, m, modelname, nsim = nsim, B = B, alpha = alpha)
  }
}
