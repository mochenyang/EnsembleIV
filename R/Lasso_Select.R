# Use Lasso method to select strong IVs
# data_unlabel: unlabeled dataset
# regressor: name of the endogenous tree
# candidates: candidate IVs as a character vector of variable names
IVSelect_Strong = function(data_unlabel, regressor, candidates) {
  # check if there are valid IVs to select from
  if (length(candidates) > 0) {
    f = stats::as.formula(paste0(regressor, "~", paste0(candidates, collapse = "+")))
    lasso.iv = hdm::rlasso(f, data = data_unlabel, post = TRUE)
    # retrieve Lasso selected IVs
    selection = lasso.iv$index
    IVs_strong = candidates[selection]
    #print(length(IVs_strong))
  }
  else {IVs_strong = candidates}
  return(IVs_strong)
}


# Use Lasso to select valid IVs
# data_test: the testing dataset
# ntree: number of trees in the forest
# method: method of IV validity selection, "lasso" is the only one implemented so far
# regressor: name of the endogenous tree
# candidates: candidate IVs as a character vector of variable names
# TODO: allow different names for "actual" column in data_test
IVSelect_Valid = function(data_test, ntree, method, regressor, candidates) {
  # if no IV validity selection is needed, return all others as valid IVs
  if (nrow(data_test) == 0) {
    valid_IVs = candidates
  }
  else if (length(candidates) == 0) {
    return(candidates)
  }
  # otherwise, perform selection
  else {
    focal_pred = data_test[,regressor]
    others_pred = data_test[,candidates]
    actual = data_test$actual
    focal_error = focal_pred - actual

    if (method == "lasso") {
      lasso.error = hdm::rlasso(focal_error ~ as.matrix(others_pred), post = TRUE)
      invalid = lasso.error$index
      valid_IVs = candidates[!invalid]
    }
  }
  return(valid_IVs)
}


# Perform Lasso select for validity and strength for a given endogenous covariate
# data_test: the testing dataset
# data_unlabel: unlabeled dataset
# control: control variables, as a character vector of variable names
# iterative: iterate between IV validity and strength selection? Default to TRUE
# ntree: number of trees in random forest
# regressor: name of the endogenous tree
# method: method of IV validity selection, "lasso" is the only one implemented so far
# TODO: allow different names for trees
LassoSelect = function(data_test, data_unlabel, control, iterative = TRUE, ntree, regressor, method = "lasso"){
  candidates = setdiff(paste0("X", 1:ntree), regressor)
  pp_abs_before = mean(abs(stats::cor(data_unlabel[,regressor], data_unlabel[,candidates])))
  pe_abs_before = mean(abs(stats::cor(data_test[,regressor]-data_test$actual, data_test[,candidates])))
  if (iterative){
    # iterate between validity and strength selection until the selection does not change anymore
    IV_valid = IVSelect_Valid(data_test, ntree, method, regressor, candidates)
    IVs = IVSelect_Strong(data_unlabel, regressor, IV_valid)
    while (length(IVs) != length(candidates)) {
      candidates = IVs
      IV_valid = IVSelect_Valid(data_test, ntree, method, regressor, candidates)
      IVs = IVSelect_Strong(data_unlabel, regressor, IV_valid)
    }
  }
  else {
    # if not iterative selection, perform one pass selection
    IV_valid = IVSelect_Valid(data_test, ntree, method, regressor, candidates)
    IVs = IVSelect_Strong(data_unlabel, regressor, IV_valid)
  }

  if (length(IVs) > 0) {
    pp_abs_after = mean(abs(stats::cor(data_unlabel[,regressor], data_unlabel[,IVs])))
    pe_abs_after = mean(abs(stats::cor(data_test[,regressor]-data_test$actual, data_test[,IVs])))
  }
  else {
    pp_abs_after = NA
    pe_abs_after = NA
  }
  return(list(IVs = IVs, correlations = c(pp_abs_before, pe_abs_before, pp_abs_after, pe_abs_after)))
}
