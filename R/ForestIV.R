#' ForestIV Main Function
#'
#' This function implements the main ForestIV approach.
#'
#' @param data_test Testing dataframe for random forest, must have a column named "actual" that contains the ground truth, and all trees' predictions.
#' @param data_unlabel Unlabel dataframe for random forest, must have all trees' predictions.
#' @param control A character vector of control variable names. Pass an empty vector if there are no control variables
#' @param method The method for IV selection. Supported values are "Lasso" (ForestIV method) and "IIV" (EnsembleIV method)
#' @param iterative Whether to perform iterative IV selection or not, default to TRUE. Only relevant when method = "Lasso"
#' @param ntree Number of trees in the random forest (or number of individual learners in general).
#' @param model_unbias Unbiased estimation model object.
#' @param family Model specification, same as in the family parameter in glm.
#' @param diagnostic Whether to output diagnostic correlations for instrument validity and strength, default to TRUE.
#' @param select_method method of IV selection. Supported values are "optimal", "top3", and "PCA" (see manuscript for details)
#' @return Estimation results
#' @export


ForestIV = function(data_test, data_unlabel, control, method, iterative = TRUE, ntree, model_unbias, family, diagnostic, select_method) {
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
      output = IIVSelect(data_test, data_unlabel, ntree, regressor, select_method)
      IVs = output$IVs
      data_unlabel_new = output$data_unlabel_new
    }

    # Estimate 2SLS
    # TODO: family parameter unnecessarily messy for now. Will clean up after deciding what to do with GLMs
    if (length(IVs) > 0) {
      if (family$family == "gaussian" & family$link == "identity") {
        f = make_formula(regressor = regressor, control = control, IVs = IVs, type = "all")
        model_IV = AER::ivreg(f, data = data_unlabel_new)
        beta_IV = stats::coef(model_IV)
        vcov_IV = stats::vcov(model_IV)
        se_IV = sqrt(diag(vcov_IV))
        # the convergence code is not useful for linear case - 2SLS always converges
        convergence = 0

        # uncomment following lines to use LIML/Fuller fitting using ivmodel package
        # model_IV = ivmodel::ivmodelFormula(f, data = data_unlabel_new)
        # est = model_IV$LIML
        # beta_IV = c(est$point.est, est$point.est.other)[c(4,1,2,3)]
        # se_IV = c(est$std.err, est$std.err.other)[c(4,1,2,3)]
        # vcov_IV = diag(se_IV^2)
        # convergence = 0
      } else {
        # implementation via ivtools package, support all glm families
        f1 = make_formula(regressor = regressor, control = control, IVs = IVs, type = "XZ")
        fitX.LZ = stats::glm(f1, data = data_unlabel_new)
        f2 = make_formula(regressor = regressor, control = control, IVs = IVs, type = "YX")
        fitY.LX = stats::glm(f2, data = data_unlabel_new, family = family)
        model_IV = ivtools::ivglm(estmethod = "ts", fitX.LZ = fitX.LZ, fitY.LX = fitY.LX, data = data_unlabel_new, ctrl = TRUE)
        # for control function estimation, a variable named "R" is attached to the end and need to be removed
        nvar = length(model_IV$est)
        beta_IV = model_IV$est[-nvar]
        vcov_IV = model_IV$vcov[-nvar,-nvar]
        se_IV = sqrt(diag(vcov_IV))
        # convergence label - 1 for converged and 0 otherwise
        convergence = model_IV$converged
      }
      H_stats = hotelling(beta_IV, vcov_IV, model_unbias)
      corrs = output$correlations
      result[[i]] = c(beta_IV, se_IV, H_stats, convergence, corrs)
      #print(i)
    }
  }
  if (diagnostic) {
    result = do.call(rbind.data.frame, result)
    colnames(result) = c(paste0("beta_", c(1:length(beta_IV))), paste0("se_", c(1:length(se_IV))), "Hotelling", "Convergence", "pp_abs_before", "pe_abs_before", "pp_abs_after", "pe_abs_after")
    return(result)
  } else {
    result = do.call(rbind.data.frame, result)
    result = result[,1:(ncol(result) - 4)]
    colnames(result) = c(paste0("beta_", c(1:length(beta_IV))), paste0("se_", c(1:length(se_IV))), "Hotelling", "Convergence")
    return(result)
  }
}


