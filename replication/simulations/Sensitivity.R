#!/usr/bin/env Rscript

# Sensitivity analyses

library(ForestIV)
library(dplyr)
library(randomForest)
library(xgboost)
library(ranger)
library(fastDummies)
library(caret)
source("SimUtility.R")

# IMPORTANT!!! The original script is designed to be run in parallel with different jobIDs on a supercomputing cluster
# If you have access to such a cluster, please set jobID accordingly
# otherwise, you need to run the script in a loop with different jobID values to replicate the full set of results
jobID = 1

# import Bike Sharing data
Bike = read.csv("hour.csv", stringsAsFactors = FALSE) %>%
   dplyr::mutate(lnCnt = log(cnt)) %>%
   dplyr::select(-instant, -dteday, -registered, -casual, -cnt)
dataset = Bike
dataname = "Bike"
actual = Bike$lnCnt
target = "lnCnt"
N = nrow(dataset)

# 4-fold CV to determine how predictive performance change w.r.t. N_label
set.seed(123456)
for (N_label in c(1000, 2000, 3000, 4000, 5000)){
  # cross-validation is only performed on D_train, which is 75% of D_label
  train = sample(1:N, N_label*0.75)
  cf_index = createFolds(y = dataset[train, target], k = 4)
  rmse = c()
  for (fold in 1:4) {
    ml_test = train[cf_index[[fold]]]
    ml_train = train[-cf_index[[fold]]]
    rf_model = ranger(as.formula(paste0(target, "~.")), data = dataset[ml_train,],
                      mtry = 3, num.trees = 100)
    pred_test = predict(rf_model, dataset[ml_test,], predict.all = FALSE)$predictions
    rmse = c(rmse, sqrt(mean((pred_test - actual[ml_test])^2)))
  }
  cat(paste("N_label:", N_label, "RMSE:", mean(rmse)))
  cat("\n")
}

# N_label: 1000 RMSE: 0.70261341689456
# N_label: 2000 RMSE: 0.646193389687828
# N_label: 3000 RMSE: 0.611494496666381
# N_label: 4000 RMSE: 0.571620711371946
# N_label: 5000 RMSE: 0.548108731980563

# 4-fold CV to determine how predictive performance change w.r.t. ntree
set.seed(123456)
# cross-validation is only performed on D_train, which is 75% of D_label
train = sample(1:N, N_label*0.75)
cf_index = createFolds(y = dataset[train, target], k = 4)
for (ntree in c(50, 100, 150, 200, 250)) {
  rmse = c()
  for (fold in 1:4) {
    ml_test = train[cf_index[[fold]]]
    ml_train = train[-cf_index[[fold]]]
    rf_model = ranger(as.formula(paste0(target, "~.")), data = dataset[ml_train,],
                      mtry = 3, num.trees = ntree)
    pred_test = predict(rf_model, dataset[ml_test,], predict.all = FALSE)$predictions
    rmse = c(rmse, sqrt(mean((pred_test - actual[ml_test])^2)))
  }
  cat(paste("Ntree:", ntree, "RMSE:", mean(rmse)))
  cat("\n")
}
# Ntree: 50 RMSE: 0.615893174152542
# Ntree: 100 RMSE: 0.609175510788303
# Ntree: 150 RMSE: 0.606210622731236
# Ntree: 200 RMSE: 0.603136578112415
# Ntree: 250 RMSE: 0.600498870645808


# simulation parameter list - focusing on linear second phase regression and EnsembleIV with PCA selection
link = "linear"
family = switch (link,
                 "linear" = gaussian(link = "identity"),
                 "logit" = binomial(link = "logit")
)
ML_method = "RF"
select_method = "PCA"
# in sensitivity analyses, we vary N_label, ntree, sigma, and cf
# base case: N_label = 3000, ntree = 100, sigma = 2, cf = 4
# vary one parameter at a time
param_list = list(c(N_label = 1000, ntree = 100, sigma = 2, cf = 4),
                  c(N_label = 2000, ntree = 100, sigma = 2, cf = 4),
                  c(N_label = 4000, ntree = 100, sigma = 2, cf = 4),
                  c(N_label = 5000, ntree = 100, sigma = 2, cf = 4),
                  c(N_label = 3000, ntree = 50, sigma = 2, cf = 4),
                  c(N_label = 3000, ntree = 150, sigma = 2, cf = 4),
                  c(N_label = 3000, ntree = 200, sigma = 2, cf = 4),
                  c(N_label = 3000, ntree = 250, sigma = 2, cf = 4),
                  c(N_label = 3000, ntree = 100, sigma = 0.5, cf = 4),
                  c(N_label = 3000, ntree = 100, sigma = 1, cf = 4),
                  c(N_label = 3000, ntree = 100, sigma = 3, cf = 4),
                  c(N_label = 3000, ntree = 100, sigma = 4, cf = 4),
                  c(N_label = 3000, ntree = 100, sigma = 2, cf = 2),
                  c(N_label = 3000, ntree = 100, sigma = 2, cf = 3),
                  c(N_label = 3000, ntree = 100, sigma = 2, cf = 5),
                  c(N_label = 3000, ntree = 100, sigma = 2, cf = 6))

for (param in param_list) {
  N_label = param["N_label"]
  ntree = param["ntree"]
  sigma = param["sigma"]
  cf = param["cf"]
  N_unlabel = N - N_label
  # simulate econometrics data
  reg_data = SimulateData(N, sigma, actual, link)
  vars = reg_data$vars
  control = reg_data$control
  config = paste(dataname, ntree, N_label, cf, sigma, link, ML_method, select_method, sep = "_")
  # results storage
  B_result = data.frame(b_0 = NA, b_1 = NA, b_2 = NA, b_3 = NA,
                        b_0_se = NA, b_1_se = NA, b_2_se = NA, b_3_se = NA,
                        unb_0 = NA, unb_1 = NA, unb_2 = NA, unb_3 = NA,
                        unb_0_se = NA, unb_1_se = NA, unb_2_se = NA, unb_3_se = NA,
                        fold = NA, round_index = NA)
  H_result = data.frame()

  # start
  # label / unlabel partition
  seed = 123455 + as.integer(jobID)
  set.seed(seed)
  label = sample(1:N, N_label)
  unlabel = sample((1:N)[-c(label)], N_unlabel)

  # cross-fitting
  cf_index = createFolds(y = dataset[label, target], k = cf)
  for (fold in 1:cf) {
    test = label[cf_index[[fold]]]
    train = label[-cf_index[[fold]]]
    # build ML model and make predictions on testing and unlabeled data
    pred = BuildModel(dataset, ML_method, target, ntree, train, test, unlabel)
    # get data
    data = GenerateData(train, test, unlabel, pred, vars)
    data_train = data$data_train
    data_test = data$data_test
    data_label = data$data_label
    data_unlabel = data$data_unlabel
    # estimate models
    est = EstimateModels(link, data_unlabel, data_label, data_test, ntree, control, family, select_method)
    # store results
    model_biased = est$model_biased
    B_result[fold,1:4] = summary(model_biased)$coefficients[1:4,1]
    B_result[fold,5:8] = summary(model_biased)$coefficients[1:4,2]
    model_unbias = est$model_unbias
    B_result[fold,9:12] = summary(model_unbias)$coefficients[1:4,1]
    B_result[fold,13:16] = summary(model_unbias)$coefficients[1:4,2]
    B_result[fold, 17] = fold
    B_result[fold, 18] = jobID

    H_temp = est$result
    H_temp$fold = fold
    H_temp$round_index = jobID
    H_result = rbind(H_result, H_temp)
    #print(fold)
  }
  write.csv(B_result, file = paste0("results/", jobID, "_", config, "_B.csv"), row.names = FALSE)
  write.csv(H_result, file = paste0("results/", jobID, "_", config, "_H.csv"), row.names = FALSE)
}

