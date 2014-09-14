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

NF = function(full, data, spline, ps = 1:4, qs = NA, B = 1000, grid = 101, bandwidth = NA) {
  n = nrow(data)
  if (is.na(qs)) {
    if (n < 50) {
      qs = 2:floor(n/4)
    } else {
      qs = ceiling(n/5):floor(n/4)
    }    
  }

  y = matrix(data[,as.character(full)[2]], ncol = 1)
  baseX = model.matrix(full, data)
  fulladdtionX = genaddX(data[,spline], ps, max(qs))  
  bs = genNFbs(B, y, baseX, fulladdtionX)
  
  addtionX = genaddX(data[,spline], ps, qs)
  pind = 0
  lacks = replicate(length(ps), {
    pind <<- pind + 1
    p = ps[pind]
    if (p == 1) {
      curX = baseX
    }
    else {
      curX = cbind(baseX, addtionX$pmatrix[, 2:p])
    }
    qind = 0
    replicate(length(qs), {
      qind <<- qind + 1
      X = cbind(curX, addtionX$qmatrix[[pind]][[qind]])
      colSums((X %*% ginv(t(X) %*% X) %*% t(X) %*% bs - bs)^2)
    }, FALSE)
  }, FALSE)
  lacks = matrix(unlist(lacks), B)
  Qdiff = sweep(lacks, 1, apply(lacks, 1, min), '-')
  cs = seq(0, max(Qdiff), length.out = grid)
  ms = sapply(cs, function(cv) as.integer(names(sort(table(apply(Qdiff <= cv, 1, function(x) min(which(x)))), decreasing = TRUE)[1])))
  fs = sapply(cs, function(cv) max(colSums(Qdiff <= cv))) / B
  
}

genaddX = function(svalue, ps, qs) {
  mins = min(svalue)
  maxs = max(svalue)
  psm = outer(svalue, ps, '^')
  qsm = lapply(ps, function(p) {
    res = lapply(qs, function(q) {
      knots = seq(mins, maxs, length.out = q + 1)[-(q+1)]
      outer(svalue, knots, function(x, y) ifelse(x < y, 0, x - y)) ^ p
    })
    names(res) = paste0("q", qs)
    res
  })
  names(qsm) = paste0("p", ps)
  list(pmatrix = psm, qmatrix = qsm)
}

genNFbs = function(B, y, baseX, fulladdtionX) {
  X = cbind(baseX, fulladdtionX$pmatrix, fulladdtionX$qmatrix[[length(fulladdtionX$qmatrix)]][[1]])
  betaest = ginv(t(X) %*% X) %*% t(X) %*% y
  sigma = sqrt(sum((X %*% betaest - y)^2) / (nrow(X) - ncol(X)))
  ybase = X %*% betaest
  matrix(as.vector(ybase) + rnorm(B * nrow(X), 0, sigma), nrow = nrow(X))
}
