#### Functions to select strong and valid IVs based on Lasso regression ####


# Use Lasso method to select strong IVs
# data_unlabel: unlabeled dataset
# regressor: name of the endogenous tree
# candidates: candidate IVs as a character vector of variable names
LassoSelect_Strong = function(data_unlabel, regressor, candidates) {
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


# Use Lasso method to select valid IVs
# data_test: the testing dataset
# regressor: name of the endogenous tree
# candidates: candidate IVs as a character vector of variable names
# TODO: allow different names for "actual" column in data_test
LassoSelect_Valid = function(data_test, regressor, candidates) {
  # if no testing data is given, return all others as valid IVs
  if (nrow(data_test) == 0) {
    valid_IVs = candidates
  }
  # skip selection if no candidate is available
  else if (length(candidates) == 0) {
    return(candidates)
  }
  # otherwise, perform selection
  else {
    focal_pred = data_test[,regressor]
    others_pred = data_test[,candidates]
    actual = data_test$actual
    focal_error = focal_pred - actual
    # estimate lasso
    lasso.error = hdm::rlasso(focal_error ~ as.matrix(others_pred), post = TRUE)
    invalid = lasso.error$index
    valid_IVs = candidates[!invalid]
  }
  return(valid_IVs)
}


# Perform Lasso select for validity and strength for a given endogenous covariate
# data_test: the testing dataset
# data_unlabel: unlabeled dataset
# iterative: iterate between IV validity and strength selection? Default to TRUE
# ntree: number of trees in random forest
# regressor: name of the endogenous tree
# TODO: allow different names for trees
LassoSelect = function(data_test, data_unlabel, iterative = TRUE, ntree, regressor){
  candidates = setdiff(paste0("X", 1:ntree), regressor)
  pp_abs_before = get_corrs(lhs = data_unlabel[,regressor], rhs = data_unlabel[,candidates])
  pe_abs_before = get_corrs(lhs = data_test[,regressor]-data_test$actual, rhs = data_test[,candidates])

  if (iterative){
    # iterate between validity and strength selection until the selection does not change anymore
    IV_valid = LassoSelect_Valid(data_test, regressor, candidates)
    IVs = LassoSelect_Strong(data_unlabel, regressor, IV_valid)
    while (length(IVs) != length(candidates)) {
      candidates = IVs
      IV_valid = LassoSelect_Valid(data_test, regressor, candidates)
      IVs = LassoSelect_Strong(data_unlabel, regressor, IV_valid)
    }
  }
  else {
    # if not iterative selection, perform one pass selection
    IV_valid = LassoSelect_Valid(data_test, regressor, candidates)
    IVs = LassoSelect_Strong(data_unlabel, regressor, IV_valid)
  }

  # some post-processing
  if (length(IVs) > 0) {
    pp_abs_after = get_corrs(lhs = data_unlabel[,regressor], rhs = data_unlabel[,IVs])
    pe_abs_after = get_corrs(lhs = data_test[,regressor]-data_test$actual, rhs = data_test[,IVs])
  }
  else {
    pp_abs_after = NA
    pe_abs_after = NA
  }
  return(list(IVs = IVs, correlations = c(pp_abs_before, pe_abs_before, pp_abs_after, pe_abs_after)))
}
