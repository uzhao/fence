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

IF.lmer = function(
  full, data, B = 100, REML = TRUE, method = c("marginal", "conditional")) {
  
  method = match.arg(method)
  if (method == "marginal") {
    return(IF.lm(full, data, B))
  }
  
  # model fit function
  mf = function(...) lmer(..., REML = REML)
  # lack of fit function
  lf = function(res) abs(fixef(res))[-1]
  # bootstrap sample
  bs = bootstrap.lmer(B, full, data, REML)

  IFbase(mf, full, data, lf, bs)
}
