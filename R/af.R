#' Adaptive Fence model selection
#'
#' Adaptive Fence model selection
#'
#' @param mf function for fitting the model
#' @param f formular of full model
#' @param ms list of formular of candidates models
#' @param d data
#' @param lf function for lack of fit 
#' @param bootset bootstrap setting
#' B:    number of bootstrap sample
#' mode: parametric("p") or nonparametric("np")
#' @param bs bootstrap sample, default setting is automatically nonparametric bootstrap (B = 100)
#' @param grid grid for c
#' @param pickf function for choosing the model for fixed c
#' @param sizef function for determinant model size
#' @return list with whatever
#' @export

adaptivefence = function(
  # model and lack of fit related
  mf, f, ms, d, lf,
  # bootstrap related
  bootset = list(), bs = NULL,
  # fence related
  grid = 51, pickf = md, sizef = size_default) {

  ans = list(full = f, models = ms, bootset = bootset, pickfunc = pickf, sizefunc = sizef)

  if (missing(ms)) {
    stop("No candidate models specified!")
  }

  # TODO: automatically pick submodel
  #       attention for default max size
  # if (is.null(ms)) {
  #   ms = findsubmodel(f)
  # }

  if (is.null(bs)) {
    bs = bootstrap(bootset, mf, f, d)
  }
  
  # XXX: more elegant implement want
  B = ifelse(is.null(bootset$B), 100, bootset$B)

  model_size = sapply(ms, sizef, d, mf)
  ans$modelsize = model_size

  lack_of_fit_matrix = sapply(ms, function(m)
    sapply(bs, function(b) lf(mf(m, b))))
  ans$lack_of_fit_matrix = lack_of_fit_matrix

  Q_m = sweep(lack_of_fit_matrix, 1, apply(lack_of_fit_matrix, 1, min), '-')
  ans$Qd_matrix = Q_m

  lof_lower = 0
  lof_upper = max(Q_m)

  cs = seq(lof_lower, lof_upper, length.out = grid)

  boot_model = matrix(NA, B, grid)

  ci = 0
  bi = 0
  model_mat = replicate(grid, {
    ci <<- ci + 1
    infence_mat <<- Q_m <= cs[ci]
    bi <<- 0
    replicate(B, {
      bi <<- bi + 1
      candidate_index = which(infence_mat[bi,])
      candidate_index = pickf(candidate_index, model_size, Q_m[bi,])
      mlof(candidate_index, Q_m[bi,])
    })
  })
  force(model_mat)
  colnames(model_mat) = cs
  rm(ci, bi)

  ans$model_mat = model_mat
  
  # if two models have same frequency, this frequency must 
  # be lower than 0.5, so maybe we don't have to worry about 
  # this case too much?
  
  freq_mat = apply(model_mat, 2, function(l) {
    tab = sort(table(l), decreasing = TRUE)
    c(as.numeric(names(tab)[1]), tab[1])
  })
  
  colnames(freq_mat) = cs
  rownames(freq_mat) = c("index", "frequency")
  ans$freq_mat = freq_mat

  cindex = peakw(freq_mat[2,])
  ans$c = cindex

  if (!is.na(cindex)) {
    ans$formula = ms[[freq_mat[1,cindex]]]
    ans$model = mf(ans$formula, d)
  } 
  else {
    ans$formula = NA
    ans$model = NA
  }
  class(ans) = "AF"
  return(ans)
}
