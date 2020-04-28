#### Some internal utility functions ####

# Compute Hotelling statistic
hotelling = function(model_IV, model_unbias) {
  b_diff = stats::coef(model_IV) - stats::coef(model_unbias)
  var_diff = stats::vcov(model_IV) + stats::vcov(model_unbias)
  H = as.vector(t(b_diff) %*% solve(var_diff) %*% b_diff)
  return(H)
}


# Compute correlations for diagnostics
get_corrs = function(lhs, rhs) {
  c = mean(abs(stats::cor(lhs, rhs)))
  return(c)
}


# make formula for IV regression
make_formula = function(regressor, control, IVs) {
  if (length(control) > 0) {
    f = stats::as.formula(paste0("Y~", regressor, "+", paste0(control, collapse = "+"), "|", paste0(IVs, collapse = "+"), "+", paste0(control, collapse = "+")))
  }
  else {
    f = stats::as.formula(paste0("Y~", regressor, "|", paste0(IVs, collapse = "+")))
  }
  return(f)
}
