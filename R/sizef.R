# FIXME: only fixed effect
size_lmerMod = function(f, d) {
  res = lmer(f, d)
  nrow(summary(res)$coef)
}

# FIXME: only fixed effect
#        effect model size is not known
size_glmerMod = function(f, d, family) {
  res = glmer(f, d)
  nrow(summary(res)$coef)
}

size_lm = function(f, d) {
  res = lm(f, d)
  nrow(summary(res)$coef)
}

# FIXME: effect model size is not known
size_glm = function(f, d, family) {
  res = glm(f, d)
  nrow(summary(res)$coef)
}
