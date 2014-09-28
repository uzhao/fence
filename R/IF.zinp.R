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

IF.zinp = function(
  full, data, B = 100, cpus = 2) {
  
  # model fit function
  mf = function(...) zeroinfl(..., dist = "poisson")
  # lack of fit function
  lf = function(res) abs(coef(res))[-1]
  # bootstrap sample
  bs = bootstrap.zinp(B, full, data)
  
  sfInit(parallel = TRUE, cpus = cpus) 
  sfExportAll()
  sfLibrary(pscl)
  invisiblefence(mf, full, data, lf, bs)
}

bootstrap.zinp = function(B, full, data) {
  X = model.matrix(full, data)
  model = zeroinfl(full, data)
  beta = coef(model)
  alphabeta = beta[(length(beta) / 2 + 1):length(beta)]
  beta = beta[1:(length(beta) / 2)]

  bootsmp = model.matrix(full, data) %*% beta
  zerop = model.matrix(full, data) %*% alphabeta
  zerop = exp(zerop) / (1 + exp(zerop))
  
  ans = replicate(B, data, FALSE)
  
  for (i in 1:B) {
    ans[[i]][,deparse(full[[2]])] = rpois(length(bootsmp), bootsmp)
    ans[[i]][,deparse(full[[2]])][rbinom(length(bootsmp), 1, zerop)] = 0
  }
  ans
}

