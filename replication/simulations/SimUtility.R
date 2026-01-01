#### Function to build ML model on training data and generate predictions on testing and unlabeled data ####
# because of variations in var names and data types, this function is necessarily ad hoc
# train, test, unlabel: row IDs of training / testing / unlabeled partitions
BuildModel = function(dataset, ML_method, target, ntree, train, test, unlabel) {
  # # train ML model and make predictions
  # if (ML_method == "RF") {
  #   rf_model = randomForest(as.formula(paste0(target, "~.")), data = dataset,
  #                           mtry = 3, subset = train, ntree = ntree)
  #   # retrieve ground truth and predictions
  #   if (target == "lnCnt") {
  #     pred_unlabel = predict(rf_model, dataset[unlabel,], predict.all = TRUE)
  #     indiv_pred_unlabel = pred_unlabel$individual
  #     aggr_pred_unlabel = pred_unlabel$aggregate
  #     pred_test = predict(rf_model, dataset[test,], predict.all = TRUE)
  #     indiv_pred_test = pred_test$individual
  #     aggr_pred_test = pred_test$aggregate
  #   }
  #   if (target == "y") {
  #     pred_unlabel = predict(rf_model, dataset[unlabel,], predict.all = TRUE)
  #     indiv_pred_unlabel = ifelse(pred_unlabel$individual=="yes",1,0)
  #     aggr_pred_unlabel = as.numeric(pred_unlabel$aggregate)-1
  #     pred_test = predict(rf_model, dataset[test,], predict.all = TRUE)
  #     indiv_pred_test = ifelse(pred_test$individual=="yes",1,0)
  #     aggr_pred_test = as.numeric(pred_test$aggregate)-1
  #   }
  # }

  if (ML_method == "RF") {
    # RF with ranger, for pred prob of each individual trees
    if (target == "lnCnt") {
      rf_model = ranger(as.formula(paste0(target, "~.")), data = dataset[train,],
                              mtry = 3, num.trees = ntree)
      indiv_pred_unlabel = predict(rf_model, dataset[unlabel,], predict.all = TRUE)$predictions
      aggr_pred_unlabel = predict(rf_model, dataset[unlabel,], predict.all = FALSE)$predictions
      indiv_pred_test = predict(rf_model, dataset[test,], predict.all = TRUE)$predictions
      aggr_pred_test = predict(rf_model, dataset[test,], predict.all = FALSE)$predictions
    }
    if (target == "y") {
      rf_model = ranger(as.formula(paste0(target, "~.")), data = dataset[train,],
                        mtry = 3, num.trees = ntree, probability = TRUE)
      indiv_pred_unlabel = predict(rf_model, dataset[unlabel,], predict.all = TRUE)$predictions[, 2, ]
      #aggr_pred_unlabel = predict(rf_model, dataset[unlabel,], predict.all = FALSE)$predictions[, 2]
      aggr_pred_unlabel = ifelse(predict(rf_model, dataset[unlabel,], predict.all = FALSE)$predictions[, 2] >= 0.5, 1, 0)
      indiv_pred_test = predict(rf_model, dataset[test,], predict.all = TRUE)$predictions[, 2, ]
      #aggr_pred_test = predict(rf_model, dataset[test,], predict.all = FALSE)$predictions[, 2]
      aggr_pred_test = ifelse(predict(rf_model, dataset[test,], predict.all = FALSE)$predictions[, 2] >= 0.5, 1, 0)
    }
  }

  if (ML_method == "XGB") {
    if (target == "lnCnt") {
      xgboost_model = xgboost(data = as.matrix(dataset[train, names(dataset) != target]),
                              label = dataset[train,target],
                              nthread = 2, nrounds = ntree, objective = "reg:squarederror")
      aggr_pred_test = predict(xgboost_model, as.matrix(dataset[test, names(dataset) != target]))
      aggr_pred_unlabel = predict(xgboost_model, as.matrix(dataset[unlabel, names(dataset) != target]))
      indiv_pred_test = data.frame(matrix(data = NA, nrow = length(test), ncol = ntree))
      indiv_pred_unlabel = data.frame(matrix(data = NA, nrow = length(unlabel), ncol = ntree))
      for (i in 1:ntree) {
        indiv_pred_test[, i] = predict(xgboost_model, as.matrix(dataset[test, names(dataset) != target]), iterationrange = c(1,i+1))
        indiv_pred_unlabel[, i] = predict(xgboost_model, as.matrix(dataset[unlabel, names(dataset) != target]), iterationrange = c(1,i+1))
      }
    }
    if (target == "y") {
      # need to first convert categorical variables into dummies for xgboost to work
      dataset = dummy_cols(dataset,
                           select_columns = c("job", "marital", "education", "default", "housing", "loan", "contact", "month", "poutcome"),
                           remove_selected_columns = TRUE)
      dataset$y = ifelse(dataset$y == "yes", 1, 0)
      # now start training and predictions
      xgboost_model = xgboost(data = as.matrix(dataset[train, names(dataset) != target]),
                              label = dataset[train,target],
                              nthread = 2, nrounds = ntree, objective = "binary:logistic")
      aggr_pred_test = ifelse(predict(xgboost_model, as.matrix(dataset[test, names(dataset) != target])) >= 0.5, 1, 0)
      aggr_pred_unlabel = ifelse(predict(xgboost_model, as.matrix(dataset[unlabel, names(dataset) != target])) >= 0.5, 1, 0)
      indiv_pred_test = data.frame(matrix(data = NA, nrow = length(test), ncol = ntree))
      indiv_pred_unlabel = data.frame(matrix(data = NA, nrow = length(unlabel), ncol = ntree))
      for (i in 1:ntree) {
        indiv_pred_test[, i] = predict(xgboost_model, as.matrix(dataset[test, names(dataset) != target]), iterationrange = c(1,i+1))
        indiv_pred_unlabel[, i] = predict(xgboost_model, as.matrix(dataset[unlabel, names(dataset) != target]), iterationrange = c(1,i+1))
      }
    }
  }
  return(list(indiv_pred_test = indiv_pred_test,
              indiv_pred_unlabel = indiv_pred_unlabel,
              aggr_pred_test = aggr_pred_test,
              aggr_pred_unlabel = aggr_pred_unlabel))
}

#### A function to simulate true data for regressions ####
SimulateData = function(N, sigma, actual, link) {
  # simulate data for econometric model (linear or logit)
  control1 = runif(N, min = -10, max = 10)
  control2 = rnorm(N, sd = 10)
  epsilon = rnorm(N, sd = sigma)
  XB = 1.0 + 0.5*actual + 2.0*control1 + control2
  control = c("control1", "control2")
  Y = switch (link,
              "linear" = XB + epsilon,
              "logit" = rbinom(N, 1, prob = 1/(1+exp(-XB)))
  )
  vars = data.frame(Y = Y, actual = actual, control1 = control1, control2 = control2)
  return(list(vars = vars,
              control = control))
}


#### A function to assemble all data partitions for estimation ####
# train, test, unlabel: row IDs of training / testing / unlabeled partitions
# pred: predictions from ML model
# sigma: sd of error term
# link: model specification (linear, logit, or poisson, etc.)
GenerateData = function(train, test, unlabel, pred, vars) {
  # get predictions from ML model
  indiv_pred_test = pred$indiv_pred_test
  indiv_pred_unlabel = pred$indiv_pred_unlabel
  aggr_pred_test = pred$aggr_pred_test
  aggr_pred_unlabel = pred$aggr_pred_unlabel

  # retrieve regression data
  control_vars = vars[,c("control1", "control2")]
  actual = vars$actual
  Y = vars$Y

  # prepare various data partitions
  data_train = data.frame(Y = Y[train], control_vars[train,], actual = actual[train])
  data_test = data.frame(indiv_pred_test, aggr_pred_test, actual = actual[test], control_vars[test,], Y = Y[test])
  data_label = data.frame(Y = Y[c(train, test)], control_vars[c(train, test),], actual = actual[c(train, test)])
  data_unlabel = data.frame(Y = Y[unlabel], control_vars[unlabel,], actual = actual[unlabel], indiv_pred_unlabel, aggr_pred_unlabel)

  return(list(data_train = data_train,
              data_test = data_test,
              data_label = data_label,
              data_unlabel = data_unlabel))
}

#### A function to estimate biased, unbiased, and ForestIV models ####
EstimateModels = function(link, data_unlabel, data_label, data_test, ntree, control, family, select_method) {
  # estimate biased regression
  model_biased = switch (link,
    "linear" = lm(Y~aggr_pred_unlabel+control1+control2, data = data_unlabel),
    "logit" = glm(Y~aggr_pred_unlabel+control1+control2, data = data_unlabel, family = "binomial")
  )

  # estimate unbiased regression
  model_unbias = switch (link,
    "linear" = lm(Y~actual+control1+control2, data = data_label),
    "logit" = glm(Y~actual+control1+control2, data = data_label, family = "binomial")
  )

  # ForestIV estimation
  result = ForestIV(data_test = data_test, data_unlabel = data_unlabel, control = control, method = "IIV",
                    ntree = ntree, model_unbias = model_unbias, family = family, diagnostic = FALSE, select_method = select_method)

  return(list(model_biased = model_biased,
              model_unbias = model_unbias,
              result = result))
}


#### Regression Calibration as a Benchmark ####
# data_test act as the calibration set
RegCal = function(data_test, data_unlabel, control, link) {
  model_cal = lm(actual~aggr_pred_test+control1+control2, data = data_test)
  pred_cal = predict(model_cal, newdata = data_unlabel %>% select(aggr_pred_test = aggr_pred_unlabel, all_of(control)))
  result = switch (link,
                   "linear" = lm(Y~pred_cal+control1+control2, data = data_unlabel),
                   "logit" = glm(Y~pred_cal+control1+control2, data = data_unlabel, family = "binomial")
  )
  return(result)
}


#### Fong&Tylor GMM as a Benchmark ####
source("../FongTyler Replication/predictionErrorGMM.R")
library(gmm)
FongTylor = function(data_train, data_test, data_unlabel, control) {
  # format data to fit GMM function inputs
  # adapted from "reddit_application.R"
  y = c(data_test$Y, data_unlabel$Y, data_train$Y)
  Xu = matrix(c(data_test$actual, rep(NA, nrow(data_unlabel)), data_train$actual), ncol = 1)
  Xo = as.matrix(rbind(data_test[,control], data_unlabel[,control], data_train[,control]), ncol = 2)
  Zu = matrix(c(data_test$aggr_pred_test, data_unlabel$aggr_pred_unlabel, rep(NA, nrow(data_train))), ncol = 1)
  v <- t <- p <- rep(0, length(y))
  v = c(rep(1, nrow(data_test)), rep(0, nrow(data_unlabel) + nrow(data_train)))
  t = c(rep(0, nrow(data_test) + nrow(data_unlabel)), rep(1, nrow(data_train)))
  p = 1 - v - t

  result = predicted_covariates(y, Xu, Xo, Zu, v, t, p)

}
