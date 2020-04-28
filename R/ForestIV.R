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
#' @param model_unbias A lm object of unbiased estimation.
#' @param diagnostic Whether to output diagnostic correlations for instrument validity and strength, default to TRUE.
#' @return ForestIV estimation results
#' @export


ForestIV = function(data_test, data_unlabel, control, method, iterative = TRUE, ntree, model_unbias, diagnostic) {
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
      f = make_formula(regressor, control, IVs)
      model_IV = AER::ivreg(f, data = data_unlabel_new)
      H_stats = hotelling(model_IV, model_unbias)
      beta = stats::coef(model_IV)
      se = sqrt(diag(stats::vcov(model_IV)))
      corrs = output$correlations
      result[[i]] = c(beta, se, H_stats, corrs)
    }
  }
  if (diagnostic) {
    result = do.call(rbind.data.frame, result)
    colnames(result) = c(paste0("beta_", c(1:length(beta))), paste0("se_", c(1:length(se))), "Hotelling", "pp_abs_before", "pe_abs_before", "pp_abs_after", "pe_abs_after")
    return(result)
  }
  else {
    result = do.call(rbind.data.frame, result)
    result = result[,1:(ncol(result) - 4)]
    colnames(result) = c(paste0("beta_", c(1:length(beta))), paste0("se_", c(1:length(se))), "Hotelling")
    return(result)
  }
}


