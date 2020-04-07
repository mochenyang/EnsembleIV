# Internal function to compute Hotelling statistic
hotelling = function(model_IV, model_unbias) {
  b_diff = stats::coef(model_IV) - stats::coef(model_unbias)
  var_diff = stats::vcov(model_IV) + stats::vcov(model_unbias)
  H = as.vector(t(b_diff) %*% solve(var_diff) %*% b_diff)
  return(H)
}
