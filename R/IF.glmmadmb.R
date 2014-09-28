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

IF.glmmadmb = function(
  full, data, family = c("g", "nb", "p"), zeroInflation = FALSE, B = 100, method = c("marginal", "conditional"), cpus = 2) {
  
  method = match.arg(method)
  family = match.arg(family)

  if (method == "marginal") {
    if (family == "p" & zeroInflation) {
      return(IF.zinp(full, data, B, cpus))
    }
    if (family == "nb") {
      return(IF.nb(full, data, B, cpus))
    }
    return(IF.glm(full, data, family, B, cpus))
  }
  
  family = switch(family,
    g = "gamma",
    nb = "nbinom", 
    p = "poisson")
  # model fit function
  mf = function(...) glmmadmb(..., family = family, zeroInflation = zeroInflation)
  # lack of fit function
  lf = function(res) abs(fixef(res))[-1]
  # bootstrap sample
  bs = bootstrap.glmmadmb(B, full, data, family, zeroInflation)
  
  sfInit(parallel = TRUE, cpus = cpus) 
  sfExportAll()
  sfLibrary(glmmADMB)
  invisiblefence(mf, full, data, lf, bs)
}
