#### Some internal utility functions ####

# Compute Hotelling statistic
hotelling = function(beta_IV, vcov_IV, model_unbias) {
  b_diff = beta_IV - stats::coef(model_unbias)
  var_diff = vcov_IV + stats::vcov(model_unbias)
  H = as.vector(t(b_diff) %*% solve(var_diff) %*% b_diff)
  return(H)
}


# Compute correlations for diagnostics
get_corrs = function(lhs, rhs) {
  c = mean(abs(stats::cor(lhs, rhs)))
  return(c)
}


# make formula for IV regression
make_formula = function(regressor, control, IVs, type) {
  if (length(control) > 0) {
    f = switch (type,
      "all" = stats::as.formula(paste0("Y~", regressor, "+", paste0(control, collapse = "+"), "|", paste0(IVs, collapse = "+"), "+", paste0(control, collapse = "+"))),
      "YX" = stats::as.formula(paste0("Y~", regressor, "+", paste(control, collapse = "+"))),
      "XZ" = stats::as.formula(paste0(regressor, "~", paste(IVs, collapse = "+")))
    )
  }
  else {
    f = switch (type,
      "all" = stats::as.formula(paste0("Y~", regressor, "|", paste0(IVs, collapse = "+"))),
      "YX" = stats::as.formula(paste0("Y~", regressor)),
      "XZ" = stats::as.formula(paste0(regressor, "~", paste(IVs, collapse = "+")))
    )
  }
  return(f)
}
