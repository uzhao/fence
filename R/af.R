#' Adaptive Fence model selection
#'
#' Adaptive Fence model selection
#'
#' @param mf function for fitting the model
#' @param f formular of full model
#' @param ms list of formular of candidates models
#' @param d data
#' @param lf function for lack of fit (to minimize)
#' @param bootset bootstrap setting
#' B:    number of bootstrap sample
#' mode: parametric("p") or nonparametric("np")
#' @param bs bootstrap sample, default setting is automatically nonparametric bootstrap (B = 100)
#' @param grid grid for c
#' @param pf function for pick measuring (to minimize)
#' @param bandwidth bandwidth for kernel smooth function
#' @return list with whatever
#' @export

adaptivefence = function(
  # model and lack of fit related
  mf, f, ms, d, lf,
  # bootstrap related
  bootset = list(), bs = NULL,
  # fence related
  grid = 51, pf = size_default, bandwidth) {

  ans = list(full = f, models = ms, bootset = bootset, pickfunc = pf, bandwidth = bandwidth)

  mf = cmpfun(mf)

  if (missing(ms)) {
    stop("No candidate models specified!")
  }

  # TODO: automatically pick submodel
  #       attention for default max size
  # if (is.null(ms)) {
  #   ms = findsubmodel(f)
  # }

  if (!is.null(bs)) {
    B = length(bs)
  }
  else {
    B = ifelse((missing(bootset) | is.null(bootset$B)), 100, bootset$B)
    bs = bootstrap(bootset, mf, f, d)
  }
  ans$B = B
  
  eval_models = lapply(ms, function(m) {
    lapply(bs, function(b) {
      mf(m, b)
    })
  })

  mi = 0
  bi = 0
  lack_of_fit_matrix = replicate(length(ms), {
    mi <<- mi + 1
    bi <<- 0
    replicate(B, {
      bi <<- bi + 1
      lf(eval_models[[mi]][[bi]])
    })
  })
  ans$lack_of_fit_matrix = lack_of_fit_matrix

  mi = 0
  bi = 0
  pick_matrix = replicate(length(ms), {
    mi <<- mi + 1
    bi <<- 0
    replicate(B, {
      bi <<- bi + 1
      pf(eval_models[[mi]][[bi]], bs[[bi]])
    })
  })
  ans$pick_matrix = pick_matrix

  rm(mi, bi)
  
  Q_m = sweep(lack_of_fit_matrix, 1, apply(lack_of_fit_matrix, 1, min), '-')
  ans$Qd_matrix = Q_m
  lof_lower = 0
  lof_upper = max(Q_m)
  cs = seq(lof_lower, lof_upper, length.out = grid)

  model_mat = matrix(NA, nrow = B, ncol = grid)
  for (i in 1:length(cs)) {
    infence_matrix = Q_m <= cs[i]
    for (bi in 1:B) {
      b_infence = infence_matrix[bi,]
      b_lack = lack_of_fit_matrix[bi,]
      b_pick = pick_matrix[bi,]
      b_pick[!b_infence] = Inf
      b_pick = which(b_pick == min(b_pick))
      model_mat[bi, i] = b_pick[which.min(b_lack[b_pick])]
    }
  }
  ans$model_mat = model_mat

  # if two models have same frequency, this frequency must 
  # be lower than 0.5, so maybe we don't have to worry about 
  # this case too much?
  
  freq_mat = apply(model_mat, 2, function(l) {
    tab = sort(table(l), decreasing = TRUE)
    c(as.numeric(names(tab)[1]), tab[1])
  })
  freq_mat[2,] = freq_mat[2,] / B
  freq_mat = rbind(freq_mat, ksmooth(cs, freq_mat[2,], kernel = "normal", bandwidth = bandwidth, x.point = cs)$y)
  
  colnames(freq_mat) = cs
  rownames(freq_mat) = c("index", "frequency", "smooth_frequency")
  ans$freq_mat = freq_mat

  cindex = peakw(cs, freq_mat[3,], 2)
  ans$c = cs[cindex]

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
