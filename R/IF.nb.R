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

IF.nb = function(
  full, data, B = 100, cpus = 2) {
  
  # model fit function
  mf = glm.nb
  # lack of fit function
  lf = function(res) abs(coef(res))[-1]
  # bootstrap sample
  bs = bootstrap.nb(B, full, data)
  
  sfInit(parallel = TRUE, cpus = cpus) 
  sfExportAll()
  invisiblefence(mf, full, data, lf, bs)
}

bootstrap.nb = function(B, full, data) {
  X = model.matrix(full, data)
  model = glm.nb(full, data)
  size = summary(model)$theta
  beta = coef(model)
  bootsmp = model.matrix(full, data) %*% beta
  
  ans = replicate(B, data, FALSE)
  for (i in 1:B) {
    ans[[i]][,deparse(full[[2]])] = rnbinom(length(bootsmp), size, mu = bootsmp)
  }
  ans  
}
