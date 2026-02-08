# =============================================== #
#    Compas CI (African-American VS Caucasian)            #
# =============================================== #
# A predictive disparity between African-American 
# and Caucasian defendants under the COMPAS risk assessment.
# =============================================== #
rm(list = ls())
source('Functions/GlambdaChen.R')
# 加载必要的包（如未安装需先安装）
# if (!require("dplyr")) install.packages("dplyr")
# if (!require("readr")) install.packages("readr")
# ---------------------获得数据------------------- #
library(dplyr)
library(readr)
# # 从 GitHub 读取 COMPAS 数据集
# compas_url <- "https://raw.githubusercontent.com/propublica/compas-analysis/master/compas-scores-two-years.csv"
# compas_data <- read_csv(compas_url)
# 读取本地 COMPAS 数据集
compas_data <- read_csv("data/compas_data.csv")

# 查看所有列名
# names(compas_data)

# 一、数据处理
# 查看两个 decile_score 列是否相同
identical(compas_data$decile_score...12, compas_data$decile_score...40)
identical(compas_data$priors_count...15, compas_data$priors_count...49)
# 若不同，需根据业务逻辑选择（假设 decile_score...12 是普通再犯评分）
compas_data <- compas_data %>%
  rename(decile_score = decile_score...12) %>%
  select(-decile_score...40)  # 删除重复列
# 同理处理 priors_count
compas_data <- compas_data %>%
  rename(priors_count = priors_count...15) %>%
  select(-priors_count...49)

# 二、数据筛选
# 筛选关键列
clean_data <- compas_data %>%
  select(race,decile_score,two_year_recid
  ) %>%
  filter(
    !is.na(decile_score), 
    !is.na(two_year_recid)
  ) %>%
  rename(Y = two_year_recid)  %>%
  mutate(Y_h = ifelse(decile_score >= 5, 1, 0)) %>% 
  select(-decile_score)

# 生成 L 列（Y == Y_h 时 L=0，否则 L=1）
clean_data <- clean_data %>%
  mutate(L = ifelse(Y == Y_h, 0, 1))
# 查看结果
table(clean_data$L)

# 选出目标族-非裔族
african_american_data <- clean_data %>%
  filter(race == "African-American")

# 选出阳预测非裔族
postive_african_american_data <- african_american_data %>%
  filter(Y_h == 1)
table(postive_african_american_data$L)
# ---------------------估计 theta_P----------------- #

# 选出白人族
caucasian_data <- clean_data %>%
  filter(race == "Caucasian")
# 选出阳预测非裔族
postive_caucasian_data <- caucasian_data %>%
  filter(Y_h == 1)
# 查看数据
theta_P = mean(postive_caucasian_data$Y)
table(postive_caucasian_data$Y)

# ---------------------开始计算------------------- #
# 定断点
a <- -0.2
b <-  0.2
step <- 0.001 # 0.001
alpha <- 0.9 # 0.95

# 计算Loss
Li <- postive_african_american_data$Y
n <- length(Li)
m <- 1
L <- matrix(Li, n, m)
Theta <- matrix(theta_P, n, m)

# 计算cut
cut = exp(-qchisq(alpha,m)/2)

# ---------------------开始模拟------------------- #
lb = 0; 
ub = 0; 
times = 0; 
elrMax_epsilonG = 0; 
elRatio=c(0);

epsilonGs = seq(a+0.01,b,step)
for(epsilonG in epsilonGs){
  # 估计方程
  E <-  matrix(epsilonG, n, m, byrow = TRUE)
  g = L - Theta - E 
  
  # 计算EL值
  z = g
  lam = c(lambdaChen(z))
  npi = 1/(1+lam*z)
  elr = prod(npi)
  # el = 2*sum(log(1+t(lam)%*%t(z)))
  # (-2)*log(elr)==el
  if(elr>max(elRatio)){elrMax_epsilonG=epsilonG}
  elRatio = c(elRatio,elr)
  if(elr>=cut && times==0){lb=epsilonG;times=1}
  if(elr>=cut && times==1){ub=epsilonG}

}

# setEPS()
# postscript(paste0('COMPAS-African-American.eps'), width=5.5, height=4)
# 空白画布
par(mfrow = c(1, 1))
plot(epsilonGs,epsilonGs,
     xlim=c(a, b),ylim=c(0,1),col='white',
     yaxt = 'n', ann = F)
axis(2, las = 1)
title( main=paste0('COMPAS','  African-American ',n),
       xlab='epsilonG',ylab='elr')
legend('topright',legend=c('  ELCP'),
       col=c(col='blue1'),
       lty=c(1),bty='n',lwd=1.5)
# cut图
abline(h=cut,col='red')
text(a+0.01/(b-a),0.2,label = round(cut,3),col='red')

# elr图
points(epsilonGs,elRatio[-1],type='o',lty=1,lwd=0.3,cex=1,pch=20,col='blue1')
abline(v=elrMax_epsilonG,col='blue1')
abline(v=lb,col='blue1')
abline(v=ub,col='blue1')
text(elrMax_epsilonG,0.6,label = elrMax_epsilonG,col='blue1')
text(ub+0.01/(b-a),0.2,label = ub,col='blue1')
text(lb-0.01/(b-a),0.2,label = lb,col='blue1')


# 计算置信区域长度
ELAL = ub-lb
  
cat('n ',' ep_h1',' lb ',' ub ',' ELAL','\n')
cat(n,elrMax_epsilonG,lb,ub,ELAL,'\n')

# P(Y=1|X=aa,Y_h=1) VS P(Y=1|X=ca,Y_h=1)
mu_aa_Y <- mean(postive_african_american_data$Y)
mu_ca_Y <- mean(postive_caucasian_data$Y)
print( mu_aa_Y - mu_ca_Y )



