# Simulation of the peripheral feature challenge and EnsembleIV robustness

# Data Generation Process:
# Ground truth is "actual"
# A pair of prediction variables X1 and X2, each as X plus a prediction error termd e1 and e2
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


# Confirm that Perf is correlated with errors, and the correlation can be manipulated by changing e_sd
set.seed(123456)
correlations = c()
N = 5000
for (e_sd in (1:20/50)) {
  data = make_data(N, e_sd)
  correlations = c(correlations, cor(data$e1, data$Perf))
}
plot(x = (1:20/50), y = correlations)


# implement permutation test
# TODO: factor this out into utility
perm_test = function(TS, data_transformed, y_diagnostic, resid_diagnostic, N_perm = 1000) {
  mean_corr_perm = c()
  for (k in 1:N_perm) {
    resid_diagnostic_perm = resid_diagnostic[sample(1:length(resid_diagnostic))]
    corr_perm = abs(cor(data_transformed - y_diagnostic, resid_diagnostic_perm))
    mean_corr_perm = c(mean_corr_perm, mean(corr_perm))
  }
  # get empirical CDF of the test statistic under null
  ecdf_perm = ecdf(mean_corr_perm)
  # then get the p value (prob of observing test statistic or more extreme values given null)
  p_val_perm = 1 - ecdf_perm(TS)

  return(p_val_perm)
}

# Estimation
# for simplicity, assume X1 as endogenous variable and X2 as candidate IV
source("../../R/IIVSelection.R")
library(AER)

# Scenario 1: fix sample size and change strength of correlation
# Suppose N = 5000, and N_train = 3000, N_test = 1000, N_diagnostic = 1000
# Repeat for varying correlation strength between prediction error and peripheral variable
# for each e_sd, repeat 500 times to obtain coverage statistic of IV estimates
set.seed(123456)
N = 5000
n_iter = 500
grid = 1:20/50
ave_corr_before = c()
ave_corr_after = c()
all_p = c()
all_coverage = c()
all_reject = c()
est_results = data.frame()
for (e_sd in grid) {
  corr_before = c()
  corr_after = c()
  est_beta = c()
  coverage = 0
  p = c()
  reject = c()
  for (k in 1:n_iter) {
    data = make_data(N, e_sd)
    data_train = data[1:3000,]
    data_test = data[3001:4000,]
    data_diagnostic = data[4001:5000,]
    # run ols on training + testing
    model_ols = lm(Y ~ actual + W, data = rbind(data_train, data_test))
    # get residual on data_diagnostic
    resid_diagnostic = data_diagnostic$Y - predict(model_ols, data_diagnostic)
    # compute correlation before IV transformation
    corr_before = c(corr_before, cor(data_diagnostic$e2, resid_diagnostic))
    # perform IV transformation using data_test
    data_diagnostic_new = IIVCreate_Valid(data_test, data_diagnostic, "X1", "X2")[[1]]
    # compute correlation after IV transformation
    corr_after = c(corr_after, cor(data_diagnostic_new$X2 - data_diagnostic_new$actual, resid_diagnostic))
    # carry out permutation test
    #pval_resid = cor.test(data_diagnostic_new$X2 - data_diagnostic_new$actual, resid_diagnostic)$p.value
    pval_resid = perm_test(corr_after, data_diagnostic_new$X2, data_diagnostic_new$actual, resid_diagnostic)
    p = c(p, pval_resid)
    reject = c(reject, ifelse(pval_resid <= 0.05, 1, 0))
    # now proceed to IV estimation and store results
    model_iv = ivreg(Y ~ X1 + W | X2 + W, data = data_diagnostic_new)
    est_beta = c(est_beta, coef(model_iv)["X1"])
    CI = confint(model_iv)["X1",]
    if (CI[1] <= 1 & CI[2] >= 1) {
      coverage = coverage + 1
    }
  }
  # before after t-test
  test_result = t.test(abs(corr_before), abs(corr_after), paired = TRUE)
  pval = test_result$p.value
  # print results
  cat(paste0("e_sd is ", e_sd, "\n"))
  cat(paste0("Average correlation before transformation is ", mean(abs(corr_before)), "\n"))
  cat(paste0("Average correlation after transformation is ", mean(abs(corr_after)), "\n"))
  cat(paste0("p value of the change is ", pval, "\n"))
  cat(paste0("IV estimate is ", mean(est_beta), " with sd ", sd(est_beta), "\n"))
  cat(paste0("Coverage rate is ", coverage / n_iter, "\n"))
  # store results for visualization
  ave_corr_before = c(ave_corr_before, mean(abs(corr_before)))
  ave_corr_after = c(ave_corr_after, mean(abs(corr_after)))
  all_coverage = c(all_coverage, coverage / n_iter)
  all_p = c(all_p, mean(p))
  all_reject = c(all_reject, mean(reject))
  est_results = rbind(est_results, c(
    e_sd,
    mean(est_beta),
    quantile(est_beta, 0.025),
    quantile(est_beta, 0.975)
  ))
}

# visualization
library(ggplot2)
library(gridExtra)

plot_corr = data.frame(e_sd = grid,
                       corr = c(ave_corr_before, ave_corr_after),
                       group = rep(c("Before IV Transformation", "After IV Transformation"), each  = 20))
plot_coverage = data.frame(e_sd = grid,
                           coverage = all_coverage)
plot_p = data.frame(e_sd = grid,
                    pval = all_p)
plot_reject = data.frame(e_sd = grid,
                         reject = all_reject)

p1 = ggplot(plot_corr, aes(x = e_sd, y = corr, color = group)) +
  geom_point(size = 1.5) + geom_line() +
  labs(x = "SD of Common Error", y = "Correlation", title = "Correlation between Prediction Error and Peripheral Feature") +
  theme_bw() +
  theme(legend.position = c(0.8, 0.2))

# p2 = ggplot(plot_p, aes(x = e_sd, y = pval)) +
#   geom_point(size = 1.5) + geom_line() +
#   geom_hline(yintercept = 0.05, color = "red") +
#   labs(x = "SD of Common Error", y = "P-Value", title = "P-Value of Correlation Test") +
#   ylim(c(0,1)) +
#   theme_bw()
p2 = ggplot(plot_reject, aes(x = e_sd, y = reject)) +
  geom_point(size = 1.5) + geom_line() +
  geom_hline(yintercept = 0.95, color = "red") +
  labs(x = "SD of Common Error", y = "Null Rejection Rate", title = "Rejection Rate of Correlation Test") +
  ylim(c(0,1)) +
  theme_bw()


p3 = ggplot(plot_coverage, aes(x = e_sd, y = coverage)) +
  geom_point(size = 1.5) + geom_line() +
  geom_hline(yintercept = 0.95, color = "red") +
  labs(x = "SD of Common Error", y = "Coverage", title = "Coverage Rate of IV Estimates") +
  ylim(c(0,1)) +
  theme_bw()

colnames(est_results) = c("e_sd", "beta", "CI_low", "CI_high")
p4 = ggplot(est_results, aes(x = e_sd, y = beta)) +
  geom_point(size = 1.5) + geom_line() +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high)) +
  geom_hline(yintercept = 1, color = "red") +
  labs(x = "SD of Common Error", y = "Coefficient Estimate", title = "IV Estimates and 95% CIs") +
  theme_bw()

grid.arrange(p1, p2, p4, nrow = 1)


# Scenario 2: fix strength of correlation and change sample size
# fix e_sd = 0.16 (because correlation test fails to reject invalid IV)
# also keep the relative proportions of train / test / diagnostic to be 3:1:1
set.seed(123456)
e_sd = 0.24
n_iter = 500
grid = (5:20)*1000
ave_corr_before = c()
ave_corr_after = c()
all_coverage = c()
all_p = c()
all_reject = c()
est_results = data.frame()
for (N in grid) {
  corr_before = c()
  corr_after = c()
  est_beta = c()
  coverage = 0
  p = c()
  for (k in 1:n_iter) {
    data = make_data(N, e_sd)
    data_train = data[1:(N*0.6),]
    data_test = data[(N*0.6):(N*0.8),]
    data_diagnostic = data[(N*0.8):N,]
    # run ols on training + testing
    model_ols = lm(Y ~ actual + W, data = rbind(data_train, data_test))
    # get residual on data_diagnostic
    resid_diagnostic = data_diagnostic$Y - predict(model_ols, data_diagnostic)
    # compute correlation before IV transformation
    corr_before = c(corr_before, cor(data_diagnostic$e2, resid_diagnostic))
    # perform IV transformation using data_test
    data_diagnostic_new = IIVCreate_Valid(data_test, data_diagnostic, "X1", "X2")[[1]]
    # compute correlation after IV transformation
    corr_after = c(corr_after, cor(data_diagnostic_new$X2 - data_diagnostic_new$actual, resid_diagnostic))
    # carry out permutation test
    #pval_resid = cor.test(data_diagnostic_new$X2 - data_diagnostic_new$actual, resid_diagnostic)$p.value
    pval_resid = perm_test(corr_after, data_diagnostic_new$X2, data_diagnostic_new$actual, resid_diagnostic)
    p = c(p, pval_resid)
    reject = c(reject, ifelse(pval_resid <= 0.05, 1, 0))
    # now proceed to IV estimation and store results
    model_iv = ivreg(Y ~ X1 + W | X2 + W, data = data_diagnostic_new)
    est_beta = c(est_beta, coef(model_iv)["X1"])
    CI = confint(model_iv)["X1",]
    if (CI[1] <= 1 & CI[2] >= 1) {
      coverage = coverage + 1
    }
  }
  # before after t-test
  test_result = t.test(abs(corr_before), abs(corr_after), paired = TRUE)
  pval = test_result$p.value
  # print results
  cat(paste0("e_sd is ", e_sd, "\n"))
  cat(paste0("Average correlation before transformation is ", mean(abs(corr_before)), "\n"))
  cat(paste0("Average correlation after transformation is ", mean(abs(corr_after)), "\n"))
  cat(paste0("p value of the change is ", pval, "\n"))
  cat(paste0("IV estimate is ", mean(est_beta), " with sd ", sd(est_beta), "\n"))
  cat(paste0("Coverage rate is ", coverage / n_iter, "\n"))
  # store results for visualization
  ave_corr_before = c(ave_corr_before, mean(abs(corr_before)))
  ave_corr_after = c(ave_corr_after, mean(abs(corr_after)))
  all_coverage = c(all_coverage, coverage / n_iter)
  all_p = c(all_p, mean(p))
  all_reject = c(all_reject, mean(reject))
  est_results = rbind(est_results, c(
    N,
    mean(est_beta),
    quantile(est_beta, 0.025),
    quantile(est_beta, 0.975)
  ))
}

# visualization
plot_corr = data.frame(N = grid,
                       corr = c(ave_corr_before, ave_corr_after),
                       group = rep(c("Before IV Transformation", "After IV Transformation"), each  = length(grid)))
plot_coverage = data.frame(N = grid,
                           coverage = all_coverage)
plot_p = data.frame(N = grid,
                    pval = all_p)
plot_reject = data.frame(e_sd = grid,
                         reject = all_reject)

p1 = ggplot(plot_corr, aes(x = N, y = corr, color = group)) +
  geom_point(size = 1.5) + geom_line() +
  labs(x = "Sample Size", y = "Correlation", title = "Correlation between Prediction Error and Peripheral Feature") +
  theme_bw() +
  theme(legend.position = c(0.8, 0.4))

# p2 = ggplot(plot_p, aes(x = N, y = pval)) +
#   geom_point(size = 1.5) + geom_line() +
#   geom_hline(yintercept = 0.05, color = "red") +
#   labs(x = "Sample Size", y = "P-Value", title = "P-Value of Correlation Test") +
#   ylim(c(0,1)) +
#   theme_bw()
p2 = ggplot(plot_reject, aes(x = e_sd, y = reject)) +
  geom_point(size = 1.5) + geom_line() +
  geom_hline(yintercept = 0.95, color = "red") +
  labs(x = "SD of Common Error", y = "Null Rejection Rate", title = "Rejection Rate of Correlation Test") +
  ylim(c(0,1)) +
  theme_bw()

p3 = ggplot(plot_coverage, aes(x = N, y = coverage)) +
  geom_point(size = 1.5) + geom_line() +
  geom_hline(yintercept = 0.95, color = "red") +
  labs(x = "Sample Size", y = "Coverage", title = "Coverage Rate of IV Estimates") +
  ylim(c(0,1)) +
  theme_bw()

colnames(est_results) = c("N", "beta", "CI_low", "CI_high")
p4 = ggplot(est_results, aes(x = N, y = beta)) +
  geom_point(size = 1.5) + geom_line() +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high)) +
  geom_hline(yintercept = 1, color = "red") +
  labs(x = "SD of Common Error", y = "Coefficient Estimate", title = "IV Estimates and 95% CIs") +
  theme_bw()

grid.arrange(p1, p2, p4, nrow = 1)
