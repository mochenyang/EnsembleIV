# Simulation of EnsembleIV extension w.r.t. Peripheral Feature Challenge

# Data Generation Process:
# Ground truth is "actual"
# A pair of prediction variables X1 and X2, each as X plus a prediction error term e1 and e2
# e1 and e2 shares a common term so that prediction variables are not perfect IVs for each other
# A peripheral variable that is the sum of prediction errors in X1 and X2 and some random noise
# Simulate a regression where the peripheral variable is part of the error term
# N: sample size
# e_sd: sd of prediction error, which will determine the strength of correlation between prediction error and peripheral variable
make_data = function(N, e_sd) {
  actual = rnorm(N)
  e = rnorm(N, sd = 0.1)
  e1 = rnorm(N, sd = e_sd)
  e2 = rnorm(N, sd = e_sd)
  X1 = actual + e1 + e
  X2 = actual + e2 + e
  Perf = e1 + e2 + rnorm(N, sd = 0.2)
  # simulate regression data
  W = runif(N)
  eps = Perf + rnorm(N)
  Y = 1 + actual + 0.5*W + eps

  return(data.frame(actual = actual, e1 = e1, e2 = e2, X1 = X1, X2 = X2, Perf = Perf, W = W, Y = Y))
}


# modify the IV transformation step to account for correlation with eps
# need data_train now because need to approximate beta and residual
IIVCreate_Valid = function(data_train, data_test, data_unlabel, regressor, candidates) {
  data_unlabel_new = data_unlabel
  # perform IV creation and transformation
  x_unlabel = data_unlabel[,regressor]
  x_test = data_test[,regressor]
  focal_error = x_test - data_test$actual
  sigma_x = stats::sd(x_test)
  cov_xe = stats::cov(x_test, focal_error)
  # estimate ols on data_train and get residual on data_test
  model_ols = lm(Y ~ actual + W, data = data_train)
  resid_test = data_test$Y - predict(model_ols, data_test)
  beta_est = coef(model_ols)["actual"]
  cov_xeps = stats::cov(x_test, resid_test)
  for (candidate in candidates) {
    z_test = data_test[,candidate]
    z_unlabel = data_unlabel[,candidate]
    sigma_z = stats::sd(z_test)
    cov_ze = stats::cov(z_test, focal_error)
    cov_zeps = stats::cov(z_test, resid_test)
    # updated lambda calculation
    lambda = ((cov_zeps - beta_est*cov_ze) / (cov_xeps - beta_est*cov_xe)) * (sigma_x / sigma_z)
    # transform the corresponding IV in the unlabeled dataset
    data_unlabel_new[,candidate] = sigma_x*z_unlabel - lambda*sigma_z*x_unlabel
  }
  return(data_unlabel_new = data_unlabel_new)
}


# Estimation
# for simplicity, assume X1 as endogenous variable and X2 as candidate IV
# Fix sample size and change strength of correlation
# Suppose N_train = 3000, N_test = 1000, N_unlabel = 10000
# Repeat for varying correlation strength between prediction error and peripheral variable
# for each e_sd, repeat 500 times to obtain coverage statistic of IV estimates
library(AER)
set.seed(123456)
N = 14000
n_iter = 500
grid = 5:20/50  # starting with 0.1 as original IV transformation becomes insufficient
results = data.frame()
for (e_sd in grid) {
  est_beta = c()
  coverage = 0
  for (k in 1:n_iter) {
    data = make_data(N, e_sd)
    data_train = data[1:3000,]
    data_test = data[3001:4000,]
    data_unlabel = data[4001:N,]
    # perform new IV transformation
    data_unlabel_new = IIVCreate_Valid(data_train, data_test, data_unlabel, "X1", "X2")
    # now proceed to IV estimation and store results
    model_iv = ivreg(Y ~ X1 + W | X2 + W, data = data_unlabel_new)
    est_beta = c(est_beta, coef(model_iv)["X1"])
    CI = confint(model_iv)["X1",]
    if (CI[1] <= 1 & CI[2] >= 1) {
      coverage = coverage + 1
    }
  }
  # print results
  cat(paste0("e_sd is ", e_sd, "\n"))
  cat(paste0("IV estimate is ", mean(est_beta), " with sd ", sd(est_beta), "\n"))
  cat(paste0("Coverage rate is ", coverage / n_iter, "\n"))
  # store results for visualization
  results = rbind(results, c(
    e_sd = e_sd,
    beta = mean(est_beta),
    CI_low = quantile(est_beta, 0.025),
    CI_high = quantile(est_beta, 0.975)
  ))
}

# visualization
library(ggplot2)

colnames(results) = c("e_sd", "beta", "CI_low", "CI_high")

ggplot(results, aes(x = e_sd, y = beta)) +
  geom_point(size = 1.5) + geom_line() +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high)) +
  geom_hline(yintercept = 1, color = "red") +
  labs(x = "SD of Common Error", y = "Coefficient Estimate") +
  theme_bw()
