# should works for lmer, lm
# not try with glmer, glm yet
size_default = function(f, d, mf) {
  res = summary(mf(f, d))
  nrow(res$coef) + ifelse(is.null(res$vcov), 0, nrow(res$vcov))
}
