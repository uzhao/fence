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

IF.glmer = function(
  full, data, family = c("b", "p"), B = 100) {
  
  family = match.arg(family)
  family = switch(family,
                  b = binomial,
                  p = poisson)
  # model fit function
  mf = function(...) glmer(..., family = family)
  # lack of fit function
  lf = function(res) abs(fixef(res))[-1]
  # bootstrap sample
  bs = bootstrap.glmer(B, full, data, family)
  
  IFbase(mf, full, data, lf, bs)
}
