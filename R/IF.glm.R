#' Invisible Fence model selection (Linear Mixed Model)
#'
#' Invisible Fence model selection (Linear Mixed Model)
#'
#' @param full formular of full model
#' @param data data
#' @param B number of bootstrap sample, parametric for lmer
#' @param grid grid for c
#' @param fence fence method to be used, (adaptive, nonadaptive, more to add)
#' @param cn cn for nonadaptive
#' @param REML REML
#' @param bandwidth bandwidth for kernel smooth function
#' @return list with whatever
#' @export

IF.glm = function(
  full, data, family = c("b", "p", "g"), B = 100, cpus = 2) {
  
  family = match.arg(family)
  family = switch(family,
                  b = binomial,
                  p = poisson,
                  g = Gamma)
  # model fit function
  mf = function(...) glm(..., family = family)
  # lack of fit function
  lf = function(res) abs(coef(res))[-1]
  # bootstrap sample
  bs = bootstrap.glm(B, full, data, family()$family)
  
  sfInit(parallel = TRUE, cpus = cpus) 
  sfExportAll()
  invisiblefence(mf, full, data, lf, bs)
}

bootstrap.glm = function(B, full, data, family) {
  X = model.matrix(full, data)
  model = glm(full, data, family = family)
  beta = coef(model)
  bootsmp = model.matrix(full, data) %*% beta
  
  ans = replicate(B, data, FALSE)
  if (family == "binomial") {
    for (i in 1:B) {
      ans[[i]][,deparse(full[[2]])] = rbinom(length(bootsmp), 1, exp(bootsmp) / (1 + exp(bootsmp)))
    }
  }
  if (family == "poisson") {
    for (i in 1:B) {
      ans[[i]][,deparse(full[[2]])] = rpois(length(bootsmp), bootsmp)
    }
  }
  if (family == "Gamma") {
    shape = 1 / summary(model)$dispersion
    for (i in 1:B) {
      ans[[i]][,deparse(full[[2]])] = rgamma(length(bootsmp), shape, 1 / bootsmp / shape)
    }
  }
  ans  
}
