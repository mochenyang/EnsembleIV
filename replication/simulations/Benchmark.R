# Benchmarking with Regression Calibration and FT-GMM

library(dplyr)
library(ranger)
library(xgboost)
library(fastDummies)
library(caret)
source("SimUtility.R")

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
#p_train = 0.75
N_unlabel = N - N_label


# import Bank Marketing data
Bank = read.csv("bank-full.csv", sep = ";", stringsAsFactors = TRUE)
dataset = Bank
dataname = "Bank"
actual = as.numeric(Bank$y)-1
target = "y"
N = nrow(dataset)
N_label = 4500
#p_train = 0.75
N_unlabel = N - N_label


# other simulation parameter list
ntree = 100
sigma = 2
boot = FALSE
cf = 4
ML_method = "RF"

for (link in c("linear", "logit")) {
  family = switch (link,
                   "linear" = gaussian(link = "identity"),
                   "logit" = binomial(link = "logit")
  )
  reg_data = SimulateData(N, sigma, actual, link)
  vars = reg_data$vars
  control = reg_data$control
  # store results for regression calibration and FongTylor GMM
  results = data.frame()
  # 100 simulation runs
  for (jobID in 1:100) {
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
      rc = RegCal(data_test, data_unlabel, control, link)
      ft = FongTylor(data_train, data_test, data_unlabel, control)
      # store results
      results_temp = data.frame(
        rbind(summary(rc)$coefficients[1:4,1]),
        rbind(summary(rc)$coefficients[1:4,2]),
        rbind(ft$beta[c(4,1,2,3),]),
        rbind(ft$se[c(4,1,2,3)]),
        fold = fold, round_index = jobID
      )
      results = rbind(results, results_temp)
    }
    print(jobID)
  }
  write.csv(results, file = paste(dataname, link, "benchmarking.csv", sep = "_"), row.names = FALSE)
}

