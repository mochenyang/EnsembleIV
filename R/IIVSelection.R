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
  data_test_new = data_test
  # perform IV creation and transformation
  x_unlabel = data_unlabel[,regressor]
  x_test = data_test[,regressor]
  focal_error = x_test - data_test$actual
  sigma_x = stats::sd(x_test)
  cov_xe = stats::cov(x_test, focal_error)
  for (candidate in candidates) {
    z_test = data_test[,candidate]
    z_unlabel = data_unlabel[,candidate]
    sigma_z = stats::sd(z_test)
    cov_ze = stats::cov(z_test, focal_error)
    lambda = (cov_ze / cov_xe) * (sigma_x / sigma_z)
    # transform the corresponding IV in the unlabeled dataset
    # also transform the testing dataset (for diagnostic purpose only)
    data_unlabel_new[,candidate] = sigma_x*z_unlabel - lambda*sigma_z*x_unlabel
    data_test_new[,candidate] = sigma_x*z_test - lambda*sigma_z*x_test
  }
  return(list(data_unlabel_new = data_unlabel_new,
              data_test_new = data_test_new))
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
# select_method: method of IV selection. See manuscript for details.
# TODO: allow different names for trees
IIVSelect = function(data_test, data_unlabel, ntree, regressor, select_method){
  candidates = setdiff(paste0("X", 1:ntree), regressor)
  pp_abs_before = get_corrs(lhs = data_unlabel[,regressor], rhs = data_unlabel[,candidates])
  pe_abs_before = get_corrs(lhs = data_test[,regressor]-data_test$actual, rhs = data_test[,candidates])

  # perform IV creation, then select the strong ones
  data_unlabel_new = IIVCreate_Valid(data_test, data_unlabel, regressor, candidates)[[1]]
  data_test_new = IIVCreate_Valid(data_test, data_unlabel, regressor, candidates)[[2]]
  if (select_method == "optimal") {
    IVs = IIVSelect_Strong(data_unlabel_new, regressor, candidates)
  }
  if (select_method == "strongest one") {
    corrs = stats::cor(data_unlabel_new[,regressor], data_unlabel_new[,candidates])
    IVs = colnames(corrs)[which.max(abs(corrs))]
  }
  if (select_method == "random one") {
    IVs = sample(candidates, size = 1)
  }
  if (select_method == "random 25%") {
    IVs = sample(candidates, size = 0.25*length(candidates))
  }
  if (select_method == "random 50%") {
    IVs = sample(candidates, size = 0.5*length(candidates))
  }
  if (select_method == "random 75%") {
    IVs = sample(candidates, size = 0.75*length(candidates))
  }
  if (select_method == "all other") {
    IVs = candidates
  }
  if (select_method == "top3") {
    corrs = stats::cor(data_unlabel_new[,regressor], data_unlabel_new[,candidates])
    IVs = colnames(corrs)[order(abs(corrs), decreasing = TRUE)[1:3]]
  }
  if (select_method == "PCA") {
    # how many components to take?
    ncomp = 3
    IVs = paste0("PCA_IV", 1:ncomp)
    pca = prcomp(data_unlabel_new[,candidates])$rotation
    data_unlabel_new[IVs] = as.matrix(data_unlabel_new[,candidates]) %*% pca[,1:ncomp]
    data_test_new[IVs] = as.matrix(data_test_new[,candidates]) %*% pca[,1:ncomp]
  }
  if (select_method == "PCA-Lasso") {
    # apply lasso on top of PCA transformed IVs
    pca = prcomp(data_unlabel_new[,candidates])$rotation
    IVs = paste0("PCA_IV", 1:length(candidates))
    data_unlabel_new[IVs] = as.matrix(data_unlabel_new[,candidates]) %*% pca
    data_test_new[IVs] = as.matrix(data_test_new[,candidates]) %*% pca
    IVs = IIVSelect_Strong(data_unlabel_new, regressor, IVs)
  }

  # some post-processing
  if (length(IVs) > 0) {
    pp_abs_after = get_corrs(lhs = data_unlabel_new[,regressor], rhs = data_unlabel_new[,IVs])
    # because ground truth is only observed on data_test, pe_abs_after can only be computed using data_test_new
    # by definition, pe_abs_after will be very close to 0
    pe_abs_after = get_corrs(lhs = data_test_new[,regressor]-data_test_new$actual, rhs = data_test_new[,IVs])
    # in the paper, with synthetic data, pe_abs_after is computed using data_unlabel_new
    # Caution: the following line only make sense in a synthetic context, when "actual" is observed on unlabel data
    #pe_abs_after = get_corrs(lhs = data_unlabel_new[,regressor]-data_unlabel_new$actual, rhs = data_unlabel_new[,IVs])
  }
  else {
    pp_abs_after = NA
    pe_abs_after = NA
  }
  return(list(IVs = IVs, data_unlabel_new = data_unlabel_new, correlations = c(pp_abs_before, pe_abs_before, pp_abs_after, pe_abs_after)))
}
