# Post-processing results and make plots
library(dplyr)
library(ggplot2)
library(gridExtra)
library(grid)
setwd("PATH_TO_YOUR_DIR")

#### Helper Function ####
# arrange regression results into latex table
# each input to this function should be a list of coefs, ses, and mse
makeTable = function(...) {
  results = list(...)
  if (length(results) == 0) {
    stop("empty input")
  } else {
    # how many coefs in total
    Ncoef = length(results[[1]]$coefs)
    for (i in 1:Ncoef) {
      all_coef = c()
      all_se = c()
      for (j in 1:length(results)) {
        all_coef = c(all_coef, results[[j]]$coefs[i])
        all_se = c(all_se, paste0("(", results[[j]]$ses[i], ")"))
      }
      cat(paste(all_coef, collapse = " & "))
      cat(" \\\\ \n")
      cat(paste(all_se, collapse = " & "))
      cat(" \\\\ \n")
    }
    # add mse at the end
    all_mse = c()
    for (j in 1:length(results)) {
      all_mse = c(all_mse, results[[j]]$mse[1])
    }
    cat(paste(all_mse, collapse = " & "))
    cat(" \\\\ \n")
  }
}

#### Main Results ####
process_B = function(data) {
  # average over cross-fitting
  data = data %>%
    na.omit() %>%
    group_by(round_index) %>%
    summarise(across(everything(), mean)) %>%
    select(-round_index)

  coefs = apply(data, 2, mean)
  ses = apply(data, 2, sd)
  mse = sum((coefs - c(1, 0.5, 2, 1))^2) + sum(ses^2)
  return(list(coefs = sprintf("%.3f", coefs),
              ses = sprintf("%.3f", ses),
              mse = sprintf("%.3f", mse)))
}

process_H = function(data) {
  data_ave = data %>%
    group_by(round_index, fold) %>%
    summarise(beta_1_ave = mean(beta_1),
              beta_2_ave = mean(beta_2),
              beta_3_ave = mean(beta_3),
              beta_4_ave = mean(beta_4)) %>%
    # further average over cross-fitting
    group_by(round_index) %>%
    summarise(beta_1_ave = mean(beta_1_ave),
              beta_2_ave = mean(beta_2_ave),
              beta_3_ave = mean(beta_3_ave),
              beta_4_ave = mean(beta_4_ave))

  coefs = apply(data_ave, 2, mean)[c("beta_1_ave","beta_2_ave","beta_3_ave","beta_4_ave")]
  ses = apply(data_ave, 2, sd)[c("beta_1_ave","beta_2_ave","beta_3_ave","beta_4_ave")]
  mse = sum((coefs - c(1, 0.5, 2, 1))^2) + sum(ses^2)
  return(list(coefs = sprintf("%.3f", coefs),
              ses = sprintf("%.3f", ses),
              mse = sprintf("%.3f", mse)))
}


biased = list()
unbiased = list()
ensembleIV = list()

for (name in c("Bike_100_3000_4_2", "Bank_100_4500_4_2")) {
  for (link in c("linear", "logit")) {
    for (ML_method in c("RF", "XGB")) {
      config = paste(name, link, ML_method, sep = "_")
      # B_result is the same regardless of selection method
      B_result = read.csv(file = paste0(config, "_PCA_B.csv"))
      biased[[config]] = process_B(B_result[,c(1:4, 18)])
      unbiased[[config]] = process_B(B_result[,c(9:12, 18)])
      for (select_method in c("top3", "PCA", "optimal")) {
        config = paste(name, link, ML_method, select_method, sep = "_")
        H_result = read.csv(file = paste0(config, "_H.csv"))
        ensembleIV[[config]] = process_H(H_result)
      }
    }
  }
}

# Make Tables
# Table 3
makeTable(biased[["Bike_100_3000_4_2_linear_RF"]],
          unbiased[["Bike_100_3000_4_2_linear_RF"]],
          ensembleIV[["Bike_100_3000_4_2_linear_RF_top3"]],
          ensembleIV[["Bike_100_3000_4_2_linear_RF_PCA"]],
          ensembleIV[["Bike_100_3000_4_2_linear_RF_optimal"]],
          biased[["Bike_100_3000_4_2_logit_RF"]],
          unbiased[["Bike_100_3000_4_2_logit_RF"]],
          ensembleIV[["Bike_100_3000_4_2_logit_RF_top3"]],
          ensembleIV[["Bike_100_3000_4_2_logit_RF_PCA"]],
          ensembleIV[["Bike_100_3000_4_2_logit_RF_optimal"]])

# Table 4
makeTable(biased[["Bank_100_4500_4_2_linear_RF"]],
          unbiased[["Bank_100_4500_4_2_linear_RF"]],
          ensembleIV[["Bank_100_4500_4_2_linear_RF_top3"]],
          ensembleIV[["Bank_100_4500_4_2_linear_RF_PCA"]],
          ensembleIV[["Bank_100_4500_4_2_linear_RF_optimal"]],
          biased[["Bank_100_4500_4_2_logit_RF"]],
          unbiased[["Bank_100_4500_4_2_logit_RF"]],
          ensembleIV[["Bank_100_4500_4_2_logit_RF_top3"]],
          ensembleIV[["Bank_100_4500_4_2_logit_RF_PCA"]],
          ensembleIV[["Bank_100_4500_4_2_logit_RF_optimal"]])

# Table 5
makeTable(biased[["Bike_100_3000_4_2_linear_XGB"]],
          unbiased[["Bike_100_3000_4_2_linear_XGB"]],
          ensembleIV[["Bike_100_3000_4_2_linear_XGB_top3"]],
          ensembleIV[["Bike_100_3000_4_2_linear_XGB_PCA"]],
          ensembleIV[["Bike_100_3000_4_2_linear_XGB_optimal"]],
          biased[["Bike_100_3000_4_2_logit_XGB"]],
          unbiased[["Bike_100_3000_4_2_logit_XGB"]],
          ensembleIV[["Bike_100_3000_4_2_logit_XGB_top3"]],
          ensembleIV[["Bike_100_3000_4_2_logit_XGB_PCA"]],
          ensembleIV[["Bike_100_3000_4_2_logit_XGB_optimal"]])

# Table 6
makeTable(biased[["Bank_100_4500_4_2_linear_XGB"]],
          unbiased[["Bank_100_4500_4_2_linear_XGB"]],
          ensembleIV[["Bank_100_4500_4_2_linear_XGB_top3"]],
          ensembleIV[["Bank_100_4500_4_2_linear_XGB_PCA"]],
          ensembleIV[["Bank_100_4500_4_2_linear_XGB_optimal"]],
          biased[["Bank_100_4500_4_2_logit_XGB"]],
          unbiased[["Bank_100_4500_4_2_logit_XGB"]],
          ensembleIV[["Bank_100_4500_4_2_logit_XGB_top3"]],
          ensembleIV[["Bank_100_4500_4_2_logit_XGB_PCA"]],
          ensembleIV[["Bank_100_4500_4_2_logit_XGB_optimal"]])


#### Benchmarking with ForestIV ####
# Function to process ForestIV outputs
process_ForestIV = function(B, H, cf) {
  # process biased and unbiased regressions
  if (cf) {
    # with cross fitting
    B = B %>%
      group_by(round_index) %>%
      summarise(across(everything(), mean)) %>%
      select(-fold)
  } else {
    # without cross fitting, just take fold = 1
    B = B %>%
      filter(fold == 1) %>%
      select(-fold)
  }

  data = left_join(H, B, by = "round_index") %>%
    mutate(bias2 = (beta_2 - unb_1)^2+(beta_1 - unb_0)^2 + (beta_3 - unb_2)^2 + (beta_4 - unb_3)^2,
           variance = se_1^2+se_2^2+se_3^2+se_4^2,
           mse = bias2+variance)

  # select estimates that minimize MSE and pass Hotelling test (with df = 4)
  H_critical = qchisq(0.95, df = 4)
  if (cf) {
    result = data %>%
      group_by(round_index, fold) %>%
      arrange(mse, .by_group = TRUE) %>%
      filter(row_number() == 1 & Hotelling < H_critical) %>%
      ungroup() %>%
      group_by(round_index) %>%
      summarise(across(everything(), mean))
  } else {
    result = data %>%
      filter(fold == 1) %>%
      group_by(round_index) %>%
      arrange(mse, .by_group = TRUE) %>%
      filter(row_number() == 1 & Hotelling < H_critical) %>%
      ungroup()
  }

  est = result %>% select(beta_1:beta_4)
  coefs = apply(est, 2, mean)
  ses = apply(est, 2, sd)
  mse = round(sum((coefs - c(1, 0.5, 2, 1))^2) + sum(ses^2), digits = 3)

  return(list(coefs = sprintf("%.3f", coefs),
              ses = sprintf("%.3f", ses),
              mse = sprintf("%.3f", mse)))
}

forestiv = list()

for (name in c("ForestIV_Bike", "ForestIV_Bank")) {
  for (link in c("linear", "logit")) {
    config = paste(name, link, sep = "_")
    B_result = read.csv(file = paste0(config, "_B.csv"))
    H_result = read.csv(file = paste0(config, "_H.csv"))
    for (cf in c(TRUE, FALSE)) {
      config = paste(name, link, cf, sep = "_")
      forestiv[[config]] = process_ForestIV(B_result, H_result, cf)
    }
  }
}

# Table
makeTable(ensembleIV[["Bike_100_3000_4_2_linear_RF_PCA"]],
          forestiv[["ForestIV_Bike_linear_FALSE"]],
          forestiv[["ForestIV_Bike_linear_TRUE"]],
          ensembleIV[["Bike_100_3000_4_2_logit_RF_PCA"]],
          forestiv[["ForestIV_Bike_logit_FALSE"]],
          forestiv[["ForestIV_Bike_logit_TRUE"]])

makeTable(ensembleIV[["Bank_100_4500_4_2_linear_RF_PCA"]],
          forestiv[["ForestIV_Bank_linear_FALSE"]],
          forestiv[["ForestIV_Bank_linear_TRUE"]],
          ensembleIV[["Bank_100_4500_4_2_logit_RF_PCA"]],
          forestiv[["ForestIV_Bank_logit_FALSE"]],
          forestiv[["ForestIV_Bank_logit_TRUE"]])


#### Benchmarking with Regression Calibration and FT-GMM ####
process_benchmark = function(data, cf) {
  if (cf) {
    # with cross fitting
    data = data %>%
      group_by(round_index) %>%
      summarise(across(everything(), mean)) %>%
      select(-round_index, -fold)
  } else {
    data = data %>%
      filter(fold == 1) %>%
      select(-round_index, -fold)
  }
  coefs = apply(data, 2, mean)
  ses = apply(data, 2, sd)
  mse = sum((coefs - c(1, 0.5, 2, 1))^2) + sum(ses^2)
  return(list(coefs = sprintf("%.3f", coefs),
              ses = sprintf("%.3f", ses),
              mse = sprintf("%.3f", mse)))
}

RC = list()
FT = list()

for (name in c("Bike", "Bank")) {
  for (link in c("linear", "logit")) {
    config = paste(name, link, "benchmarking.csv", sep = "_")
    data = read.csv(config)
    for (cf in c(TRUE, FALSE)) {
      config = paste(name, link, cf, sep = "_")
      RC[[config]] = process_benchmark(data[,c(1:4, 17, 18)], cf)
      FT[[config]] = process_benchmark(data[,c(9:12, 17, 18)], cf)
    }
  }
}

# Table (linear regression)
makeTable(ensembleIV[["Bike_100_3000_4_2_linear_RF_PCA"]],
          RC[["Bike_linear_FALSE"]],
          RC[["Bike_linear_TRUE"]],
          FT[["Bike_linear_FALSE"]],
          FT[["Bike_linear_TRUE"]],
          ensembleIV[["Bike_100_3000_4_2_logit_RF_PCA"]],
          RC[["Bike_logit_FALSE"]],
          RC[["Bike_logit_TRUE"]])

makeTable(ensembleIV[["Bank_100_4500_4_2_linear_RF_PCA"]],
          RC[["Bank_linear_FALSE"]],
          RC[["Bank_linear_TRUE"]],
          FT[["Bank_linear_FALSE"]],
          FT[["Bank_linear_TRUE"]],
          ensembleIV[["Bank_100_4500_4_2_logit_RF_PCA"]],
          RC[["Bank_logit_FALSE"]],
          RC[["Bank_logit_TRUE"]])


#### Plot correlation plots before vs. after transformation ####
basedir = "G:/My Drive/projects/Robust Inference ForestIV/code/IIV simulations/Results for paper/"
H_Bike = read.csv(file = paste0(basedir, "Bike_100_1000_2000_2_linear_RF_PCA_H_diagnostic.csv"))
H_Bank = read.csv(file = paste0(basedir, "Bank_100_1500_3000_2_linear_RF_PCA_H_diagnostic.csv"))

p_Bike_pp = ggplot(data = data.frame(pp = c(H_Bike$pp_abs_before, H_Bike$pp_abs_after),
                                     x = rep(c(" Before", "After"), each = nrow(H_Bike)))) +
  geom_boxplot(aes(x = x, y = pp)) +
  labs(x = "", y = "Relevance Correlation") +
  theme_bw()

p_Bike_pe = ggplot(data = data.frame(pp = c(H_Bike$pe_abs_before, H_Bike$pe_abs_after),
                                     x = rep(c(" Before", "After"), each = nrow(H_Bike)))) +
  geom_boxplot(aes(x = x, y = pp)) +
  labs(x = "", y = "Exclusion Correlation") +
  theme_bw()

p_Bank_pp = ggplot(data = data.frame(pp = c(H_Bank$pp_abs_before, H_Bank$pp_abs_after),
                                     x = rep(c(" Before", "After"), each = nrow(H_Bank)))) +
  geom_boxplot(aes(x = x, y = pp)) +
  labs(x = "", y = "Relevance Correlation") +
  theme_bw()

p_Bank_pe = ggplot(data = data.frame(pp = c(H_Bank$pe_abs_before, H_Bank$pe_abs_after),
                                     x = rep(c(" Before", "After"), each = nrow(H_Bank)))) +
  geom_boxplot(aes(x = x, y = pp)) +
  labs(x = "", y = "Exclusion Correlation") +
  theme_bw()

p_Bike = grid.arrange(p_Bike_pp, p_Bike_pe, nrow = 1, top=textGrob("Bike Sharing Dataset"))
p_Bank = grid.arrange(p_Bank_pp, p_Bank_pe, nrow = 1, top=textGrob("Bank Marketing Dataset"))

# plot dimension: 1000 * 350
grid.arrange(p_Bike, p_Bank, nrow = 1)



#### process facebook study results ####
library(stargazer)

process_B = function(data) {
  data = data %>%
    group_by(round_index) %>%
    summarise(across(everything(), mean)) %>%
    select(-round_index)

  coefs = apply(data, 2, mean)
  ses = apply(data, 2, sd)
  names(coefs) = NULL
  names(ses) = NULL
  return(list(coefs = coefs, ses = ses))
}

process_H = function(data) {
  data_ave = data %>%
    group_by(round_index, fold) %>%
    summarise(across(beta_1:beta_6, mean)) %>%
    group_by(round_index) %>%
    summarise(across(beta_1:beta_6, mean))

  coefs = apply(data_ave, 2, mean)[2:7]
  ses = apply(data_ave, 2, sd)[2:7]
  names(coefs) = NULL
  names(ses) = NULL

  return(list(coefs = coefs, ses = ses))
}

# read results
bow_B = process_B(read.csv("bow_B_LGB.csv"))
bert_B = process_B(read.csv("BERT_B_LGB.csv"))
bow_PCA = process_H(read.csv("bow_PCA_H_LGB.csv"))
bert_PCA = process_H(read.csv("BERT_PCA_H_LGB.csv"))

# use a placeholder regression to apply stargazer
d <- as.data.frame(matrix(rnorm(100 * 6), nc = 6))
names(d) <- c("DV", "sentiment", "wordcount", "type_photo", "type_status", "type_video")
p <- lm(DV ~., d)

stargazer(list(p,p,p,p,p), type = "latex",
          coef = list(bow_B$coefs[13:18], bow_B$coefs[1:6], bow_PCA$coefs, bert_B$coefs[1:6], bert_PCA$coefs),
          se = list(bow_B$coefs[19:24], bow_B$ses[1:6], bow_PCA$ses, bert_B$ses[1:6], bert_PCA$ses),
          omit.stat = "all", digits = 3, star.cutoffs = c(0.05, 0.01, 0.001))



#### Sensitivity analyses ####
# for bike sharing
param_list = list(c(N_label = 1000, ntree = 100, sigma = 2, cf = 4),
                  c(N_label = 2000, ntree = 100, sigma = 2, cf = 4),
                  c(N_label = 3000, ntree = 100, sigma = 2, cf = 4),
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

# for bank marketing
param_list = list(c(N_label = 1500, ntree = 100, sigma = 2, cf = 4),
                  c(N_label = 3000, ntree = 100, sigma = 2, cf = 4),
                  c(N_label = 4500, ntree = 100, sigma = 2, cf = 4),
                  c(N_label = 6000, ntree = 100, sigma = 2, cf = 4),
                  c(N_label = 7500, ntree = 100, sigma = 2, cf = 4),
                  c(N_label = 4500, ntree = 50, sigma = 2, cf = 4),
                  c(N_label = 4500, ntree = 150, sigma = 2, cf = 4),
                  c(N_label = 4500, ntree = 200, sigma = 2, cf = 4),
                  c(N_label = 4500, ntree = 250, sigma = 2, cf = 4),
                  c(N_label = 4500, ntree = 100, sigma = 0.5, cf = 4),
                  c(N_label = 4500, ntree = 100, sigma = 1, cf = 4),
                  c(N_label = 4500, ntree = 100, sigma = 3, cf = 4),
                  c(N_label = 4500, ntree = 100, sigma = 4, cf = 4),
                  c(N_label = 4500, ntree = 100, sigma = 2, cf = 2),
                  c(N_label = 4500, ntree = 100, sigma = 2, cf = 3),
                  c(N_label = 4500, ntree = 100, sigma = 2, cf = 5),
                  c(N_label = 4500, ntree = 100, sigma = 2, cf = 6))

process_result = function(N_label, ntree, sigma, cf) {
  config = paste("Bike", ntree, N_label, cf, sigma, "linear_RF_PCA", sep = "_")
  B = read.csv(file = paste(config, "B.csv", sep = "_"))
  H = read.csv(file = paste(config, "H.csv", sep = "_"))

  # aggregation based on cross-fitting
  B = B %>% group_by(round_index) %>%
    summarise(b_1 = mean(b_1),
              unb_1 = mean(unb_1))
  H = H %>% group_by(round_index, fold) %>%
    summarise(beta_2_ave = mean(beta_2)) %>%
    group_by(round_index) %>%
    summarise(beta_2_ave = mean(beta_2_ave))
  est = left_join(H, B, by = "round_index")
  output = data.frame()
  output = rbind(output, c(N_label, ntree, sigma, cf, "Biased", mean(est$b_1), quantile(est$b_1, 0.025), quantile(est$b_1, 0.975)))
  output = rbind(output, c(N_label, ntree, sigma, cf, "Unbiased", mean(est$unb_1), quantile(est$unb_1, 0.025), quantile(est$unb_1, 0.975)))
  output = rbind(output, c(N_label, ntree, sigma, cf, "EnsembleIV", mean(est$beta_2_ave), quantile(est$beta_2_ave, 0.025), quantile(est$beta_2_ave, 0.975)))
  colnames(output) = c("N_label", "ntree", "sigma", "cf", "estimator", "v_mean", "v_low", "v_high")
  return(output)
}

plot_data = data.frame()
for (param in param_list) {
  N_label = param["N_label"]
  ntree = param["ntree"]
  sigma = param["sigma"]
  cf = param["cf"]
  output = process_result(N_label, ntree, sigma, cf)
  plot_data = rbind(plot_data, output)
}

plot_data = plot_data %>%
  mutate(across(c(1,2,4), as.integer)) %>%
  mutate(across(c(3,6,7,8), as.numeric))

# sensitivity wrt N_label
p_Nlabel = ggplot(plot_data %>% filter(ntree == 100 & sigma == 2 & cf == 4), aes(x = N_label, y = v_mean, color = estimator)) +
  geom_point(size = 2, position = position_dodge(width = 300)) + geom_line(position = position_dodge(width = 300)) +
  geom_errorbar(aes(ymin = v_low, ymax = v_high), width = 300, position = position_dodge(width = 300)) +
  geom_hline(yintercept = 0.5) +
  scale_x_continuous(breaks = c(1000, 2000, 3000, 4000, 5000)) +
  scale_color_manual(values = c("Biased" = "blue", "Unbiased" = "green", "EnsembleIV" = "red")) +
  labs(x = "Size of Labeled Data", y = "Estimates", title = "Varying Size of Labeled Data") +
  theme_bw() + theme(legend.position = "top")

# sensitivity wrt ntree
p_ntree = ggplot(plot_data %>% filter(N_label == 3000 & sigma == 2 & cf == 4), aes(x = ntree, y = v_mean, color = estimator)) +
  geom_point(size = 2, position = position_dodge(width = 15)) + geom_line(position = position_dodge(width = 15)) +
  geom_errorbar(aes(ymin = v_low, ymax = v_high), width = 15, position = position_dodge(width = 15)) +
  geom_hline(yintercept = 0.5) +
  scale_x_continuous(breaks = c(50, 100, 150, 200, 250)) +
  scale_color_manual(values = c("Biased" = "blue", "Unbiased" = "green", "EnsembleIV" = "red")) +
  labs(x = "Number of Weak Learners", y = "Estimates", title = "Varying Number of Weak Learners") +
  theme_bw() + theme(legend.position = "top")

# sensitivity wrt sigma
p_sigma = ggplot(plot_data %>% filter(N_label == 3000 & ntree == 100 & cf == 4), aes(x = sigma, y = v_mean, color = estimator)) +
  geom_point(size = 2, position = position_dodge(width = 0.15)) + geom_line(position = position_dodge(width = 0.15)) +
  geom_errorbar(aes(ymin = v_low, ymax = v_high), width = 0.15, position = position_dodge(width = 0.15)) +
  geom_hline(yintercept = 0.5) +
  scale_x_continuous(breaks = c(0.5, 1.0, 2.0, 3.0, 4.0)) +
  scale_color_manual(values = c("Biased" = "blue", "Unbiased" = "green", "EnsembleIV" = "red")) +
  labs(x = "SD of Error Term", y = "Estimates", title = "Varying Amount of Noise") +
  theme_bw() + theme(legend.position = "top")

# sensitivity wrt cf
p_cf = ggplot(plot_data %>% filter(N_label == 3000 & ntree == 100 & sigma == 2), aes(x = cf, y = v_mean, color = estimator)) +
  geom_point(size = 2, position = position_dodge(width = 0.3)) + geom_line(position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = v_low, ymax = v_high), width = 0.3, position = position_dodge(width = 0.3)) +
  geom_hline(yintercept = 0.5) +
  scale_x_continuous(breaks = c(2,3,4,5,6)) +
  scale_color_manual(values = c("Biased" = "blue", "Unbiased" = "green", "EnsembleIV" = "red")) +
  labs(x = "Number of Folds", y = "Estimates", title = "Varying Number of Folds in Cross-Fitting") +
  theme_bw() + theme(legend.position = "top")

# 1000 * 800
grid.arrange(p_Nlabel, p_ntree, p_sigma, p_cf, nrow = 2)


#### Coverage rate (MS 2nd revision) ####
library(nortest)
getDist = function(data) {
  data_ave = data %>%
    group_by(round_index, fold) %>%
    summarise(beta_1_ave = mean(beta_1),
              beta_2_ave = mean(beta_2),
              beta_3_ave = mean(beta_3),
              beta_4_ave = mean(beta_4)) %>%
    # further average over cross-fitting
    group_by(round_index) %>%
    summarise(beta_1_ave = mean(beta_1_ave),
              beta_2_ave = mean(beta_2_ave),
              beta_3_ave = mean(beta_3_ave),
              beta_4_ave = mean(beta_4_ave))

  beta_MLV = data_ave[,"beta_2_ave"]
  return(beta_MLV$beta_2_ave)
}

ensembleIV = list()

# First, get empirical distributions from non-bootstrapped experiments
setwd("PATH_TO_YOUR_DIR")
for (name in c("Bike_100_3000_4_2", "Bank_100_4500_4_2")) {
  for (link in c("linear", "logit")) {
    for (ML_method in c("RF")) {
      for (select_method in c("top3", "PCA", "optimal")) {
        config = paste(name, link, ML_method, select_method, sep = "_")
        H_result = read.csv(file = paste0(config, "_H.csv"))
        ensembleIV[[config]] = getDist(H_result)
      }
    }
  }
}

# next, get bootstrapped distribution from bootstrapped experiments
setwd("PATH_TO_YOUR_DIR")
for (name in c("Bike_100_3000_4_2", "Bank_100_4500_4_2")) {
  for (link in c("linear", "logit")) {
    for (ML_method in c("RF")) {
      for (select_method in c("top3", "PCA", "optimal")) {
        config = paste(name, link, ML_method, select_method, "bootstrap", sep = "_")
        H_result = read.csv(file = paste0(config, "_H.csv"))
        ensembleIV[[config]] = getDist(H_result)
      }
    }
  }
}

# compare the standard error estimates
for (name in c("Bike_100_3000_4_2", "Bank_100_4500_4_2")) {
  for (link in c("linear", "logit")) {
    for (ML_method in c("RF")) {
      for (select_method in c("top3", "PCA", "optimal")) {
        #config = paste(name, link, ML_method, select_method, sep = "_")
        config = paste(name, link, ML_method, select_method, "bootstrap", sep = "_")
        #cat(paste(sd(ensembleIV[[config]]), sd(ensembleIV[[config_boot]])))

        # normality (K-S test)
        norm_test = lillie.test(ensembleIV[[config]])
        cat(paste(config, norm_test$p.value))

        # coverage rate
        se = sd(ensembleIV[[config]])
        coverage_rate = mean(ifelse((ensembleIV[[config]] - 1.96*se <= 0.5) & (ensembleIV[[config]] + 1.96*se >= 0.5), 1, 0))
        cat(paste(config, coverage_rate))

        cat("\n")
      }
    }
  }
}

#### Variance weighted aggregation of individual ests (MS minor revision) ####

process_H_weighted = function(data) {
  data_ave = data %>%
    group_by(round_index, fold) %>%
    summarise(weight_1 = sum(1/se_1^2),
              weight_2 = sum(1/se_2^2),
              weight_3 = sum(1/se_3^2),
              weight_4 = sum(1/se_4^2),
              beta_1_ave = sum(beta_1/se_1^2)/weight_1,
              beta_2_ave = sum(beta_2/se_2^2)/weight_2,
              beta_3_ave = sum(beta_3/se_3^2)/weight_3,
              beta_4_ave = sum(beta_4/se_4^2)/weight_4) %>%
    # further average over cross-fitting
    group_by(round_index) %>%
    summarise(beta_1_ave = mean(beta_1_ave),
              beta_2_ave = mean(beta_2_ave),
              beta_3_ave = mean(beta_3_ave),
              beta_4_ave = mean(beta_4_ave))

  coefs = apply(data_ave, 2, mean)[c("beta_1_ave","beta_2_ave","beta_3_ave","beta_4_ave")]
  ses = apply(data_ave, 2, sd)[c("beta_1_ave","beta_2_ave","beta_3_ave","beta_4_ave")]
  mse = sum((coefs - c(1, 0.5, 2, 1))^2) + sum(ses^2)
  return(list(coefs = sprintf("%.3f", coefs),
              ses = sprintf("%.3f", ses),
              mse = sprintf("%.3f", mse)))
}

for (name in c("Bike_100_3000_4_2", "Bank_100_4500_4_2")) {
  for (link in c("linear", "logit")) {
    for (ML_method in c("RF")) {
      config = paste(name, link, ML_method, sep = "_")
      # B_result is the same regardless of selection method
      B_result = read.csv(file = paste0(config, "_PCA_B.csv"))
      biased[[config]] = process_B(B_result[,c(1:4, 18)])
      unbiased[[config]] = process_B(B_result[,c(9:12, 18)])
      for (select_method in c("top3", "PCA", "optimal")) {
        config = paste(name, link, ML_method, select_method, sep = "_")
        H_result = read.csv(file = paste0(config, "_H.csv"))
        ensembleIV[[config]] = process_H_weighted(H_result)
      }
    }
  }
}

# Appendix G
makeTable(biased[["Bike_100_3000_4_2_linear_RF"]],
          unbiased[["Bike_100_3000_4_2_linear_RF"]],
          ensembleIV[["Bike_100_3000_4_2_linear_RF_top3"]],
          ensembleIV[["Bike_100_3000_4_2_linear_RF_PCA"]],
          ensembleIV[["Bike_100_3000_4_2_linear_RF_optimal"]],
          biased[["Bike_100_3000_4_2_logit_RF"]],
          unbiased[["Bike_100_3000_4_2_logit_RF"]],
          ensembleIV[["Bike_100_3000_4_2_logit_RF_top3"]],
          ensembleIV[["Bike_100_3000_4_2_logit_RF_PCA"]],
          ensembleIV[["Bike_100_3000_4_2_logit_RF_optimal"]])

makeTable(biased[["Bank_100_4500_4_2_linear_RF"]],
          unbiased[["Bank_100_4500_4_2_linear_RF"]],
          ensembleIV[["Bank_100_4500_4_2_linear_RF_top3"]],
          ensembleIV[["Bank_100_4500_4_2_linear_RF_PCA"]],
          ensembleIV[["Bank_100_4500_4_2_linear_RF_optimal"]],
          biased[["Bank_100_4500_4_2_logit_RF"]],
          unbiased[["Bank_100_4500_4_2_logit_RF"]],
          ensembleIV[["Bank_100_4500_4_2_logit_RF_top3"]],
          ensembleIV[["Bank_100_4500_4_2_logit_RF_PCA"]],
          ensembleIV[["Bank_100_4500_4_2_logit_RF_optimal"]])


#### LIML/Fuller est (MS minor revision) ####
setwd("PATH_TO_YOUR_DIR")

for (name in c("Bike_100_3000_4_2", "Bank_100_4500_4_2")) {
  for (link in c("linear")) {
    for (ML_method in c("RF")) {
      config = paste(name, link, ML_method, sep = "_")
      # B_result is the same regardless of selection method
      B_result = read.csv(file = paste0(config, "_PCA_LIML_B.csv"))
      biased[[config]] = process_B(B_result[,c(1:4, 18)])
      unbiased[[config]] = process_B(B_result[,c(9:12, 18)])
      for (select_method in c("top3", "PCA", "optimal")) {
        for (est_method in c("LIML", "Fuller")){
          config = paste(name, link, ML_method, select_method, est_method, sep = "_")
          H_result = read.csv(file = paste0(config, "_H.csv"))
          ensembleIV[[config]] = process_H(H_result)
        }
      }
    }
  }
}

# Appendix H
makeTable(biased[["Bike_100_3000_4_2_linear_RF"]],
          unbiased[["Bike_100_3000_4_2_linear_RF"]],
          ensembleIV[["Bike_100_3000_4_2_linear_RF_top3_LIML"]],
          ensembleIV[["Bike_100_3000_4_2_linear_RF_PCA_LIML"]],
          ensembleIV[["Bike_100_3000_4_2_linear_RF_optimal_LIML"]],
          ensembleIV[["Bike_100_3000_4_2_linear_RF_top3_Fuller"]],
          ensembleIV[["Bike_100_3000_4_2_linear_RF_PCA_Fuller"]],
          ensembleIV[["Bike_100_3000_4_2_linear_RF_optimal_Fuller"]])

makeTable(biased[["Bank_100_4500_4_2_linear_RF"]],
          unbiased[["Bank_100_4500_4_2_linear_RF"]],
          ensembleIV[["Bank_100_4500_4_2_linear_RF_top3_LIML"]],
          ensembleIV[["Bank_100_4500_4_2_linear_RF_PCA_LIML"]],
          ensembleIV[["Bank_100_4500_4_2_linear_RF_optimal_LIML"]],
          ensembleIV[["Bank_100_4500_4_2_linear_RF_top3_Fuller"]],
          ensembleIV[["Bank_100_4500_4_2_linear_RF_PCA_Fuller"]],
          ensembleIV[["Bank_100_4500_4_2_linear_RF_optimal_Fuller"]])

