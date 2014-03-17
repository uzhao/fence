# should works for lmer, lm
# not try with glmer, glm yet
size_default = function(res, ...) {
  res = summary(res)
  nrow(res$coef) + ifelse(is.null(res$vcov), 0, nrow(res$vcov))
}
