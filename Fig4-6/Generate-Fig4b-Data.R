# =============================================== #
#                  Compas AL                      #
# =============================================== #
rm(list = ls())
source('Functions/GlambdaChen.R')
source('Functions/epsilonG_CI.R')
# ---------------------获得数据------------------- #
library(dplyr)
library(readr)
# 读取本地 COMPAS 数据集
compas_data <- read_csv("data/compas_data.csv")
table(compas_data$sex)
table(compas_data$age_cat)
table(compas_data$race)

# 一、数据处理
identical(compas_data$decile_score...12, compas_data$decile_score...40)
identical(compas_data$priors_count...15, compas_data$priors_count...49)
compas_data <- compas_data %>%
  rename(decile_score = decile_score...12) %>%
  select(-decile_score...40)  # 删除重复列
compas_data <- compas_data %>%
  rename(priors_count = priors_count...15) %>%
  select(-priors_count...49)  # 删除重复列 

# 二、数据筛选
clean_data <- compas_data %>%
  select(sex,age_cat,race,decile_score,two_year_recid
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

# 选出阳预测
postive_clean_data <- clean_data %>%
  filter(Y_h == 1)
table(postive_clean_data$L)

# ---------------------基准线-------------------- #
# 选出白人族
caucasian_data <- clean_data %>%
  filter(race == "Caucasian")
# 选出阳预测非裔族
postive_caucasian_data <- caucasian_data %>%
  filter(Y_h == 1)
# 查看数据
table(postive_caucasian_data$Y)
table(postive_caucasian_data$L)
# 基准线
theta_P = mean(postive_caucasian_data$Y)
cat('theta_P =', theta_P,'\n')

# ---------------------开始计算------------------- #
# 定区间
a <- -0.4
b <-  0.2
step <- 0.001
alpha <- 0.95

# ALL
# 输出
group_specific_data <- postive_clean_data
cat('\n-------------------------------\n')
cat('ALL,')
cat("分组样本数:", nrow(group_specific_data), "\n")
# 计算置信区间
epsilonG_CI(group_specific_data,theta_P,a,b,step,alpha)


sexs <- unique(postive_clean_data$sex)          
age_cats <- unique(postive_clean_data$age_cat)     

# Sex
for(sex_name in sexs){
  # 选择分组数据
  group_specific_data <- postive_clean_data %>%
    filter(
      sex == sex_name,
    )
  # 输出
  cat('\n-------------------------------\n')
  cat(sex_name,',')
  cat("分组样本数:", nrow(group_specific_data), "\n")
  # 计算置信区间
  epsilonG_CI(group_specific_data,theta_P,a,b,step,alpha)
}

# age_cat
for(age_cat_name in age_cats){
  # 选择分组数据
  group_specific_data <- postive_clean_data %>%
    filter(
      age_cat == age_cat_name
    )
  # 输出
  cat('\n-------------------------------\n')
  cat(age_cat_name,',')
  cat("分组样本数:", nrow(group_specific_data), "\n")
  # 计算置信区间
  epsilonG_CI(group_specific_data,theta_P,a,b,step,alpha)
}

# sex * age_cat
for(sex_name in sexs){
  for(age_cat_name in age_cats){
    # 选择分组数据
    group_specific_data <- postive_clean_data %>%
      filter(
        sex == sex_name,
        age_cat == age_cat_name
      )
    # 输出
    cat('\n-------------------------------\n')
    cat(sex_name,',',age_cat_name,',')
    cat("分组样本数:", nrow(group_specific_data), "\n")
    # 计算置信区间
    epsilonG_CI(group_specific_data,theta_P,a,b,step,alpha)
  }
}






