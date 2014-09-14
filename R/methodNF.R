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
      generateNFmodels(p, q, full, spline, data)
    })
  })

  y = data[,resp]

  grouptms = as.list(grouptms)
  grouptms$sep = "+"

  X = model.matrix(as.formula(paste0(resp, "~0+", do.call(paste, grouptms))), data)
}

generateNFmodels = function(p, q, full, spline, data) {
  resp = as.character(full)[2]
  pred = as.character(full)[3]

  for (sterm in spline) {
    splineterm = addsplineterms(data[,sterm], p, q, sterm)
    pred = paste0(pred, "+", splineterm$pred)
    data = 
  }
}

addsplineterms = function(spvalue, p, q, sterm) {
  knots = seq(min(spvalue), max(spvalue), length.out = q)

}
