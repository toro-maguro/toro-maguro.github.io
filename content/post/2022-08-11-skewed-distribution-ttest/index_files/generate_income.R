# incomeデータを国民生活基礎調査に基づいて作成
MakeIncomeData <- function(){
  income_cat <- c(1:25)
  min_income <- c(0, 500000, 1000000, 1500000, 2000000, 2500000, 3000000, 3500000, 
                  4000000, 4500000, 5000000, 5500000, 6000000, 6500000, 7000000, 7500000,
                  8000000, 8500000, 9000000, 9500000, 10000000, 11000000, 12000000, 15000000, 20000000)
  max_income <- c(499999, 999999, 1499999, 1999999, 2499999, 2999999, 3499999, 3999999,
                  4499999, 4999999, 5499999, 5999999, 6499999, 6999999, 7499999, 7999999, 8499999,
                  8999999, 9499999, 9999999, 10999999, 11999999, 14999999, 19999999, 100000000)
  
  df_income <- data.frame(income_cat = income_cat,
                          min_income = min_income, 
                          max_income = max_income)
  return(df_income)
}

# ID, sex, income (confounder), x (binary exposure), y (binary outcome) を含むデータを作成する 
# incomeが少ないほどx and yが発生しやすい-> confounding bias away from the null
# 引数: n_sample=人数、risk_ratio=曝露によるリスク比
MakeData <- function(n_sample=1000, risk_ratio=1.5){
  # income data & income_threshold
  income_cat <- MakeIncomeData()
  income_threshold <- c(0, 0.012, 0.064, 0.127, 0.19, 0.259, 0.326, 0.397, 0.454, 
                        0.51, 0.559, 0.608, 0.646, 0.692, 0.726, 0.759, 0.788,
                        0.814, 0.837, 0.859, 0.878, 0.909, 0.928, 0.966, 0.987, 1)
  
  df <- data.frame(ID = 1:n_sample)
  df <- df %>% 
    mutate(sex = rbinom(nrow(df), 1, 0.5)) %>% 
    mutate(prob_income = runif(n_sample, min=0, max=1)) %>% 
    mutate(income_cat = cut(prob_income,
                            breaks = income_threshold,
                            right = TRUE, 
                            include.lowest = TRUE,
                            labels = c(1:25)
    )) %>% 
    mutate(income_cat = as.integer(income_cat))
  
  df <- left_join(df, income_cat, by = "income_cat") %>% 
    select(ID, sex, income_cat, min_income, max_income) %>% 
    mutate(income = runif(n_sample, min = min_income, max = max_income)) %>% 
    select(ID, sex, income)
  
  df <- df %>% 
    mutate(income_cat = case_when(
      income < 2500000 ~ "<250",
      income >= 2500000 & income < 4000000 ~ "<400",
      income >= 4000000 & income < 7000000 ~ "<700",
      income >= 7000000 ~ "≥700"
    ),
    income_cat = factor(income_cat, levels = c("<250", "<400", "<700", "≥700"))) 
  
  # 曝露の発生確率と曝露(0/1)を定義
  df <- df %>% 
    mutate(x_prob = case_when(
      income < 2500000 ~ 0.25,
      income >= 2500000 & income < 4000000 ~ 0.2,
      income >= 4000000 & income < 7000000 ~ 0.15,
      income >= 7000000 ~ 0.1
    )) %>% 
    mutate(x = rbinom(nrow(df), 1, x_prob)) # 定めた曝露発生確率の二項分布に従って曝露を決定
  
  # 曝露しているとイベント発生確率がrisk_ratio倍になる
  df <- df %>% 
    mutate(y_prob = case_when(
      income < 2500000 ~ 0.25 * (1 + x*(risk_ratio-1)),
      income >= 2500000 & income < 4000000 ~ 0.2 * (1 + x*(risk_ratio-1)),
      income >= 4000000 & income < 7000000 ~ 0.15 * (1 + x*(risk_ratio-1)),
      income >= 7000000 ~ 0.1 * (1 + x*(risk_ratio-1))
    )) %>% 
    mutate(y = rbinom(nrow(df), 1, y_prob))
  
  return(df)
}

# 欠測をMARで発生させる
# MARのベーシックシナリオ: 男性はincomeの欠損が発生しやすい
GenerateMissing_MAR <- function(df, p_missing_men=0.4, p_missing_women=0.2){
  res <- df %>% 
    mutate(p_missing = case_when(
      sex == 0 ~ p_missing_women,
      sex == 1 ~ p_missing_men
    )) %>% 
    mutate(income_missing = rbinom(nrow(df), 1, p_missing)) %>% 
    mutate(income_observed = case_when(
      income_missing == 0 ~ income,
      income_missing == 1 ~ NA_real_
    )) %>% 
    mutate(income_cat_observed = case_when(
      income_observed < 2500000 ~ "<250",
      income_observed >= 2500000 & income_observed < 4000000 ~ "<400",
      income_observed >= 4000000 & income_observed < 7000000 ~ "<700",
      income_observed >= 7000000 ~ "≥700"
    ),
    income_cat_observed = factor(income_cat_observed, levels = c("<250", "<400", "<700", "≥700"))) 
  
  res <- res %>% select(ID, sex, income, income_observed, income_cat, income_cat_observed, x, y)
  
  return(res)
}

# 欠測をMNARで発生させる
# incomeが少ないほど、incomeの欠測が生じる
# p1: income<250万円, p2: 250-400, p3: 400-700, p4: ≥700
GenerateMissing_MNAR <- function(df, p1=0.4, p2=0.3, p3=0.2, p4=0.1){
  res <- df %>% 
    mutate(p_missing = case_when(
      income < 2500000 ~ p1,
      income >= 2500000 & income < 4000000 ~ p2,
      income >= 4000000 & income < 7000000 ~ p3,
      income >= 7000000 ~ p4
    )) %>% 
    mutate(income_missing = rbinom(nrow(df), 1, p_missing)) %>% 
    mutate(income_observed = case_when(
      income_missing == 0 ~ income,
      income_missing == 1 ~ NA_real_
    )) %>% 
    mutate(income_cat_observed = case_when(
      income_observed < 2500000 ~ "<250",
      income_observed >= 2500000 & income_observed < 4000000 ~ "<400",
      income_observed >= 4000000 & income_observed < 7000000 ~ "<700",
      income_observed >= 7000000 ~ "≥700"
    ),
    income_cat_observed = factor(income_cat_observed, levels = c("<250", "<400", "<700", "≥700"))) 
  
  res <- res %>% select(ID, sex, income, income_observed, income_cat, income_cat_observed, x, y)
  
  return(res)
}
