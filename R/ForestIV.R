#' ForestIV Main Function
#'
#' This function implements the main ForestIV approach.
#'
#' @param data_test Testing dataframe for random forest, must have a column named "actual" that contains the ground truth, and all trees' predictions.
#' @param data_unlabel Unlabel dataframe for random forest, must have all trees' predictions.
#' @param control A character vector of control variable names. Pass an empty vector if there are no control variables
#' @param method The method for IV selection. Supported values are "Lasso" (lasso-based selection) and "IIV" (imperfect IV method)
#' @param iterative Whether to perform iterative IV selection or not, default to TRUE. Only relevant when method = "Lasso"
#' @param ntree Number of trees in the random forest.
#' @param model_unbias Unbiased estimation.
#' @param family Model specification, same as in the family parameter in glm.
#' @param diagnostic Whether to output diagnostic correlations for instrument validity and strength, default to TRUE.
#' @return ForestIV estimation results
#' @export


ForestIV = function(data_test, data_unlabel, control, method, iterative = TRUE, ntree, model_unbias, family, diagnostic) {
  result = list()
  for (i in 1:ntree) {
    # use i-th tree as the endogenous covariate
    regressor = paste0("X", i)
    if (method == "Lasso") {
      output = LassoSelect(data_test, data_unlabel, iterative, ntree, regressor)
      IVs = output$IVs
      data_unlabel_new = data_unlabel
    }
    if (method == "IIV") {
      output = IIVSelect(data_test, data_unlabel, ntree, regressor)
      IVs = output$IVs
      data_unlabel_new = output$data_unlabel_new
    }

    # Estimate 2SLS
    if (length(IVs) > 0) {
      if (family$family == "gaussian" & family$link == "identity") {
        f = make_formula(regressor, control, IVs, "all")
        model_IV = AER::ivreg(f, data = data_unlabel_new)
        beta_IV = stats::coef(model_IV)
        vcov_IV = stats::vcov(model_IV)
        se_IV = sqrt(diag(vcov_IV))
      }
      else {
        # outdated implementation via OneSampleMR, only supports logit and poisson
        # model_IV = OneSampleMR::tsri(f, data = data_unlabel_new, link = link)
        # varlist = c("(Intercept)", regressor, "rescontrol1", "rescontrol2")
        # beta_IV = stats::coef(model_IV$fit)[varlist]
        # vcov_IV = stats::vcov(model_IV$fit)[varlist,varlist]
        # se_IV = sqrt(diag(vcov_IV))

        # implementation via ivtools, support all glm families
        f1 = make_formula(regressor, control, IVs, "XZ")
        fitX.LZ = stats::glm(f1, data = data_unlabel_new)
        f2 = make_formula(regressor, control, IVs, "YX")
        fitY.LX = stats::glm(f2, data = data_unlabel_new, family = family)
        model_IV = ivtools::ivglm(estmethod = "ts", fitX.LZ = fitX.LZ, fitY.LX = fitY.LX, data = data_unlabel_new, ctrl = TRUE)
        varlist = c("(Intercept)", regressor, control)
        beta_IV = model_IV$est[varlist]
        vcov_IV = model_IV$vcov[varlist,varlist]
        se_IV = sqrt(diag(vcov_IV))
      }
      H_stats = hotelling(beta_IV, vcov_IV, model_unbias)
      corrs = output$correlations
      result[[i]] = c(beta_IV, se_IV, H_stats, corrs)
      #print(i)
    }
  }
  if (diagnostic) {
    result = do.call(rbind.data.frame, result)
    colnames(result) = c(paste0("beta_", c(1:length(beta_IV))), paste0("se_", c(1:length(se_IV))), "Hotelling", "pp_abs_before", "pe_abs_before", "pp_abs_after", "pe_abs_after")
    return(result)
  }
  else {
    result = do.call(rbind.data.frame, result)
    result = result[,1:(ncol(result) - 4)]
    colnames(result) = c(paste0("beta_", c(1:length(beta_IV))), paste0("se_", c(1:length(se_IV))), "Hotelling")
    return(result)
  }
}


