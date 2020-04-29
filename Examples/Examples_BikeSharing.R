# Simulation Example with Bike Sharing Dataset (one iteration)
library(ForestIV)
library(dplyr)
library(randomForest)

# import Bike Sharing data
Bike = read.csv("hour.csv", stringsAsFactors = FALSE) %>%
  dplyr::mutate(lnCnt = log(cnt)) %>%
  dplyr::select(-instant, -dteday, -registered, -casual, -cnt)

# parameters for random forest
ntree = 100
N = nrow(Bike)
N_train = 1000
N_test = 200
N_unlabel = N - N_train - N_test
# parameters for econometric model and simulation
sigma = 2


set.seed(123456)
# make random forest
train = sample(1:nrow(Bike), N_train)
test = sample((1:nrow(Bike))[-train], N_test)
unlabel = sample((1:nrow(Bike))[-c(train, test)], N_unlabel)
Bike.rf=randomForest(lnCnt ~ . , data = Bike,
                     mtry = 3, subset = train, ntree = ntree)

# retrieve ground truth and predictions
actual = Bike$lnCnt
pred_unlabel = predict(Bike.rf, Bike[unlabel,], predict.all = TRUE)
indiv_pred_unlabel = pred_unlabel$individual
aggr_pred_unlabel = pred_unlabel$aggregate
pred_test = predict(Bike.rf, Bike[test,], predict.all = TRUE)
indiv_pred_test = pred_test$individual
aggr_pred_test = pred_test$aggregate

# simulate data for econometric model
control1 = runif(N, min = -10, max = 10)
control2 = rnorm(N, sd = 10)
epsilon = rnorm(N, sd = sigma)
control = c("control1", "control2")
Y = 1.0 + 0.5*actual + 2.0*control1 + control2 + epsilon

# prepare various data partitions
data_train = data.frame(Y = Y[train], control1 = control1[train], control2 = control2[train], actual = actual[train])
data_test = data.frame(indiv_pred_test, aggr_pred_test, actual = actual[test])
data_label = data.frame(Y = Y[c(train, test)], control1 = control1[c(train, test)], control2 = control2[c(train, test)], actual = actual[c(train, test)])
data_unlabel = data.frame(Y = Y[unlabel], control1 = control1[unlabel], control2 = control2[unlabel], actual = actual[unlabel], indiv_pred_unlabel, aggr_pred_unlabel)

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
