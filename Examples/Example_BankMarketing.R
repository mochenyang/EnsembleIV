# Simulation Example with Bank Marketing Dataset (one iteration)
library(ForestIV)
library(dplyr)
library(randomForest)

Bank = read.csv("bank-full.csv", sep = ";")

# parameters for random forest
ntree = 100
N = nrow(Bank)
N_train = 1500
N_test = 500
N_unlabel = N - N_train - N_test
# parameters for econometric model and simulation
sigma = 2


set.seed(12345)
# make random forest
train = sample(1:nrow(Bank), N_train)
test = sample((1:nrow(Bank))[-train], N_test)
unlabel = sample((1:nrow(Bank))[-c(train, test)], N_unlabel)
Bank.rf=randomForest(y ~ . , data = Bank,
                     mtry = 3, subset = train, ntree = ntree)

# retrieve ground truth and predictions
actual = as.numeric(Bank$y)-1
pred_unlabel = predict(Bank.rf, Bank[unlabel,], predict.all = TRUE)
indiv_pred_unlabel = ifelse(pred_unlabel$individual=="yes",1,0)
aggr_pred_unlabel = as.numeric(pred_unlabel$aggregate)-1

pred_test = predict(Bank.rf, Bank[test,], predict.all = TRUE)
indiv_pred_test = ifelse(pred_test$individual=="yes",1,0)
aggr_pred_test = as.numeric(pred_test$aggregate)-1

# simulate data for econometric model
control1 = runif(N, min = -1, max = 1)
control2 = rnorm(N, sd = 1)
epsilon = rnorm(N, sd = sigma)
control = c("control1", "control2")
Y = 1.0 + 0.5*actual + 2.0*control1 + control2 + epsilon

# prepare various data partitions
data_train = data.frame(Y = Y[train], control1 = control1[train],
                        control2 = control2[train], actual = actual[train])
data_test = data.frame(indiv_pred_test, aggr_pred_test, actual = actual[test])
data_label = data.frame(Y = Y[c(train, test)], control1 = control1[c(train, test)],
                        control2 = control2[c(train, test)], actual = actual[c(train, test)])
data_unlabel = data.frame(Y = Y[unlabel], control1 = control1[unlabel],
                          control2 = control2[unlabel], actual = actual[unlabel],
                          indiv_pred_unlabel, aggr_pred_unlabel)

# biased regression
model_biased = lm(Y~aggr_pred_unlabel+control1+control2, data = data_unlabel)
summary(model_biased)

# unbiased regression
model_unbias = lm(Y~actual+control1+control2, data = data_label)
summary(model_unbias)

# ForestIV estimation
result = ForestIV(data_test, data_unlabel, control, method = "Lasso", iterative = TRUE, ntree, model_unbias, diagnostic = FALSE)
#result = ForestIV(data_test, data_unlabel, control, method = "IIV", ntree = ntree, model_unbias = model_unbias, diagnostic = FALSE)

H_critical = qchisq(0.95, df = length(control)+2)
coef_unbiased = coef(model_unbias)
result %>% rowwise() %>%
  mutate(bias2 = sum((c(beta_1, beta_2, beta_3, beta_4) - coef_unbiased)^2),
         variance = se_1^2 + se_2^2 + se_3^2 + se_4^2,
         mse = bias2+variance) %>% ungroup() %>%
  arrange(mse) %>%
  filter(row_number() == 1 & Hotelling < H_critical)
