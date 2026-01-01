#!/usr/bin/env Rscript

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
N_label = 3000
N_unlabel = N - N_label


# import Bank Marketing data
Bank = read.csv("bank-full.csv", sep = ";", stringsAsFactors = TRUE)
dataset = Bank
dataname = "Bank"
actual = as.numeric(Bank$y)-1
target = "y"
N = nrow(dataset)
N_label = 4500
N_unlabel = N - N_label


# other simulation parameter list
ntree = 100
sigma = 2
boot = FALSE
cf = 4

for (link in c("linear", "logit")) {
  family = switch (link,
                   "linear" = gaussian(link = "identity"),
                   "logit" = binomial(link = "logit")
  )
  reg_data = SimulateData(N, sigma, actual, link)
  vars = reg_data$vars
  control = reg_data$control
  for (ML_method in c("RF", "XGB")) {
    for (select_method in c("top3", "PCA", "optimal")) {
      # config params: dataset + ntree + N_label + cf + sigma + link
      config = paste(dataname, ntree, N_label, cf, sigma, link, ML_method, select_method, sep = "_")

      # results storage
      B_result = data.frame(b_0 = NA, b_1 = NA, b_2 = NA, b_3 = NA,
                            b_0_se = NA, b_1_se = NA, b_2_se = NA, b_3_se = NA,
                            unb_0 = NA, unb_1 = NA, unb_2 = NA, unb_3 = NA,
                            unb_0_se = NA, unb_1_se = NA, unb_2_se = NA, unb_3_se = NA,
                            fold = NA, round_index = NA)
      H_result = data.frame()

      # start
      # label / unlabel partition, depending on whether to bootstrap or not
      if (boot) {
        # fixed a single label / unlabel split (assume this is the data you have)
        set.seed(123456)
        label = sample(1:N, N_label)
        unlabel = sample((1:N)[-c(label)], N_unlabel)
        # then sample with replacement different parts
        seed = 123455 + as.integer(jobID)
        set.seed(seed)
        label = sample(label, replace = TRUE)
        unlabel = sample(unlabel, replace = TRUE)
      } else {
        seed = 123455 + as.integer(jobID)
        set.seed(seed)
        label = sample(1:N, N_label)
        unlabel = sample((1:N)[-c(label)], N_unlabel)
      }

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
  }
}

