#' Nonparametric Fence
#'
#' Nonparametric Fence
#'
#' @param full formular of full model
#' @param data data
#' @param spline variable need spline terms
#' @param B number of bootstrap sample, parametric for lmer
#' @param grid grid for c
#' @param bandwidth bandwidth for kernel smooth function
#' @return list with whatever
#' @export

NF = function(full, data, spline, B = 100, grid = 101, bandwidth = NA) {
  ps = 1:4
  n = nrow(data)
  if (n < 50) {
    qs = 1:floor(n/4)
  } else {
    qs = ceil(n/5):floor(n/4)
  }

  ms = lapply(ps, function(p) {
    lapply(qs, function(q) {
      generateNFmodels(p, q, full, spline)
    })
  })
    generateNFmodels(ps, qs, full, spline)


  y = data[,resp]

  grouptms = as.list(grouptms)
  grouptms$sep = "+"

  X = model.matrix(as.formula(paste0(resp, "~0+", do.call(paste, grouptms))), data)



}

generateNFmodels = function(p, q, full, spline) {
  resp = as.character(full)[2]
  fulltms = attributes(terms(full))$term.labels
  extratms = fulltms[!(fulltms %in%  spline)]
  extratms = as.list(extratms)
  extratms$sep = "+"

  splinetms = sapply(spline, function(x) {
    xs = paste0(x, "^", 1:p)
    xs = as.list(xs)
    xs$sep = "+"
    do.call(paste, xs)
  })
  splinetms = as.list(splinetms)
  splinetms$sep = "+"
  ans = list()

  ans$x = as.formula(paste0(resp, "~1+", do.call(paste, extratms), "+", do.call(paste, splinetms)))
  ans$z = spline
  ans$q = q
  ans
}

