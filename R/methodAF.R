#' Adaptive Fence model selection
#'
#' Adaptive Fence model selection
#'
#' @param mf function for fitting the model
#' @param f formular of full model
#' @param ms list of formular of candidates models
#' @param d data
#' @param lf function for lack of fit (to minimize)
#' @param pf function for pick measuring (to minimize)
#' @param bs bootstrap sample
#' @param grid grid for c
#' @param bandwidth bandwidth for kernel smooth function
#' @return list with whatever
#' @export

adaptivefence = function(
  # model and lack of fit related
  mf, f, ms, d, lf, pf,
  # bootstrap sample
  bs,
  # fence related
  grid = 101, bandwidth) {

  ans = list(full = f, models = ms, pickfunc = pf)
  mf = cmpfun(mf)

  if (missing(ms)) {
    stop("No candidate models specified!")
  }

  if (missing(bs)) {
    stop("No bootstrap sample specified!")
  }

  eval_models = lapply(ms, function(m) {
    lapply(bs, function(b) {
      try(mf(m, b), silent = TRUE)
    })
  })

  em = sapply(eval_models, function(eval_model) sapply(eval_model, class))
  eb = rowSums(em == "try-error") == 0
  if (sum(eb) != length(bs)) {
    warning(paste0("Some bootstrap sample are not avaiable, new bootstrap size is ", sum(eb)))
  }
  B = sum(eb)
  ans$B = sum(eb)
  for (i in 1:length(ms)) {
    eval_models[[i]] = eval_models[[i]][eb]
  }

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
      pf(eval_models[[mi]][[bi]])
    })
  })
  ans$pick_matrix = pick_matrix

  rm(mi, bi)

  Q_m = sweep(lack_of_fit_matrix, 1, apply(lack_of_fit_matrix, 1, min), '-')
  ans$Qd_matrix = Q_m
  lof_lower = 0
  lof_upper = max(Q_m)
  cs = seq(lof_lower, lof_upper, length.out = grid)

  if (is.na(bandwidth)) {
    bandwidth = (cs[2] - cs[1]) * 3
  }
  ans$bandwidth = bandwidth * 1

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
    ans$sel_model = mf(ans$formula, d)
  }
  else {
    ans$formula = NA
    ans$sel_model = NA
  }
  class(ans) = "AF"
  return(ans)
}

#' Plot Adaptive Fence model selection
#'
#' @export
plot.AF = function(res) {
  tmp = data.frame(c = as.numeric(colnames(res$freq_mat)),
                   p = res$freq_mat[2, ],
                   sp = res$freq_mat[3, ],
                   m = as.factor(res$freq_mat[1, ]))
  p = ggplot(tmp) +
    geom_point(aes(x = c, y = p, colour = m)) +
    geom_line(aes(x = c, y = sp), linetype="dashed") +
    geom_line(aes(x = c, y = sp, colour = m)) +
    ylim(0, 1)
  p
}

#' Summary Adaptive Fence model selection
#'
#' @export
summary.AF = function(res) {
  print(res$sel_model)
}
