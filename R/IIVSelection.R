#### Functions to create valid IVs from imperfect IVs ####
# Based on imperfect IV work: Nevo and Rosen (2012); Clarke and Matt (2017)


# Create valid IVs from imperfect ones using testing dataset, then transform IVs in unlabeled dataset
# data_test: the testing dataset
# regressor: name of the endogenous tree
# candidates: candidate IVs as a character vector of variable names
# TODO: allow different names for "actual" column in data_test
IIVCreate_Valid = function(data_test, data_unlabel, regressor, candidates) {
  # TODO: add safe checks on data_test and candidates, similar to LassoSelect_Valid

  data_unlabel_new = data_unlabel
  x_unlabel = data_unlabel[,regressor]
  # perform IV creation and transformation
  x = data_test[,regressor]
  focal_error = x - data_test$actual
  sigma_x = stats::sd(x)
  cov_xe = stats::cov(x,focal_error)
  for (candidate in candidates) {
    z = data_test[,candidate]
    sigma_z = stats::sd(z)
    cov_ze = stats::cov(z,focal_error)
    lambda = (cov_ze / cov_xe) * (sigma_x / sigma_z)
    # transform the corresponding IV in the unlabeled dataset
    z_unlabel = data_unlabel[,candidate]
    data_unlabel_new[,candidate] = sigma_x*z_unlabel - lambda*sigma_z*x_unlabel
  }
  return(data_unlabel_new)
}


# Use Lasso method to select strong IVs, only to be used after IIVCreate_Valid
# data_unlabel_new: the transformed unlabeled dataset, output of IIVCreate_Valid function
# regressor: name of the endogenous tree
# candidates: candidate IVs as a character vector of variable names
IIVSelect_Strong = function(data_unlabel_new, regressor, candidates) {
  f = stats::as.formula(paste0(regressor, "~", paste0(candidates, collapse = "+")))
  lasso.iv = hdm::rlasso(f, data = data_unlabel_new, post = TRUE)
  # retrieve Lasso selected IVs
  selection = lasso.iv$index
  IVs_strong = candidates[selection]
  return(IVs_strong)
}


# Perform IIVCreate to get valid IVs, then use Lasso to select strong ones
# data_test: the testing dataset
# data_unlabel: unlabeled dataset
# ntree: number of trees in random forest
# regressor: name of the endogenous tree
# TODO: allow different names for trees
IIVSelect = function(data_test, data_unlabel, ntree, regressor){
  candidates = setdiff(paste0("X", 1:ntree), regressor)
  pp_abs_before = get_corrs(lhs = data_unlabel[,regressor], rhs = data_unlabel[,candidates])
  pe_abs_before = get_corrs(lhs = data_test[,regressor]-data_test$actual, rhs = data_test[,candidates])

  # perform IV creation, then select the strong ones
  data_unlabel_new = IIVCreate_Valid(data_test, data_unlabel, regressor, candidates)
  IVs = IIVSelect_Strong(data_unlabel_new, regressor, candidates)

  # some post-processing
  if (length(IVs) > 0) {
    pp_abs_after = get_corrs(lhs = data_unlabel_new[,regressor], rhs = data_unlabel_new[,IVs])
    pe_abs_after = get_corrs(lhs = data_test[,regressor]-data_test$actual, rhs = data_test[,IVs])
  }
  else {
    pp_abs_after = NA
    pe_abs_after = NA
  }
  return(list(IVs = IVs, data_unlabel_new = data_unlabel_new, correlations = c(pp_abs_before, pe_abs_before, pp_abs_after, pe_abs_after)))
}
