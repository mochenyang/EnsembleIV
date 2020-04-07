#' ForestIV Main Function
#'
#' This function implements the main ForestIV approach.
#'
#' @param data_test Testing dataframe for random forest, must have a column named "actual" that contains the ground truth, and all trees' predictions.
#' @param data_unlabel Unlabel dataframe for random forest, must have all trees' predictions.
#' @param control A character vector of control variable names.
#' @param iterative Whether to perform iterative IV selection or not, default to TRUE.
#' @param ntree Number of trees in the random forest.
#' @param model_unbias A lm object of unbiased estimation.
#' @param diagnostic Whether to output diagnostic correlations for instrument validity and strength, default to TRUE.
#' @return ForestIV estimation results
#' @export


ForestIV = function(data_test, data_unlabel, control, iterative = TRUE, ntree, model_unbias, diagnostic) {
  result = list()
  for (i in 1:ntree) {
    # use i-th tree as the endogenous covariate
    regressor = paste0("X", i)
    output = LassoSelect(data_test, data_unlabel, control, iterative, ntree, regressor)
    IVs = output$IVs
    # Estimate 2SLS
    if (length(IVs) > 0) {
      if (length(control) > 0) {
        f = stats::as.formula(paste0("Y~", regressor, "+", paste0(control, collapse = "+"), "|", paste0(IVs, collapse = "+"), "+", paste0(control, collapse = "+")))
      }
      else {
        f = stats::as.formula(paste0("Y~", regressor, "|", paste0(IVs, collapse = "+")))
      }
        model_IV = AER::ivreg(f, data = data_unlabel)
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


