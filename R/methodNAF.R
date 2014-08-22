#' Nonadaptive Fence model selection
#'
#' Nonadaptive Fence model selection
#'
#' @param mf function for fitting the model
#' @param f formular of full model
#' @param ms list of formular of candidates models
#' @param d data
#' @param lf function for lack of fit (to minimize)
#' @param pf function for pick measuring (to minimize)
#' @return list with whatever

nonadaptivefence = function(
  mf, f, ms, d, lf, pf, cn) {
  ans = list()

  mf = cmpfun(mf)

  if (missing(ms)) {
    stop("No candidate models specified!")
  }
  ans$models = ms

  eval_models = lapply(ms, function(m) {
    mf(m, d)
  })

  lack_of_fit = sapply(eval_models, lf)
  ans$lack_of_fit = lack_of_fit

  Q = lack_of_fit - min(lack_of_fit)
  pick = sapply(eval_models, pf)

  infence = Q < cn
  pick[!infence] = Inf
  inpick = pick == min(pick)
  Q[!inpick] = Inf
  index = which.min(Q)

  ans$formula = ms[[index]]
  ans$sel_model = mf(ans$formula, d)

  class(ans) = "NAF"
  return(ans)
}

summary.NAF = function(res) {
  print(res$sel_model)
}
