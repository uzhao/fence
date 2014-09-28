#' Invisible Fence model selection
#'
#' Invisible Fence model selection
#'
#' @param mf function for fitting the model
#' @param f formular of full model
#' @param d data
#' @param lf a function provides lack of fit for all terms, e.g. function(x) abs(coef(x))[-1]
#' @param bs bootstrap sample
#' @return list with whatever
#' @export

invisiblefence = function(
  # model and lack of fit related
  mf, f, d, lf, 
  # bootstrap sample
  bs) {
  
  ans = list(full = f)
  mf = cmpfun(mf)
  
  if (missing(bs)) {
    stop("No bootstrap sample specified!")
  }
  
  boot_evaluations = sfClusterApplyLB(bs, function(b) {
    try(mf(f, b), silent = TRUE)
  })
  sfStop()
  
  bm = sapply(boot_evaluations, function(x) class(x)[1]) 
  bb = sum(bm != "try-error")
  if (bb != length(bs)) {
    warning(paste0("Some bootstrap sample are not avaiable, new bootstrap size is ", sum(bb)))
  }
  B = bb
  ans$B = bb
  
  boot_evaluations = boot_evaluations[bm != "try-error"]
  
  s_mat = sapply(boot_evaluations, lf)
  orders_mat = apply(s_mat, 2, function(x) order(x, decreasing = TRUE))
    
  freq = sapply(1:nrow(orders_mat), function(size) {
    max(table(apply(orders_mat, 2, function(x) do.call(paste, as.list(sort(x[1:size]))))))
  })
  freq = freq / B
  
  size = peakglobal(freq[1:(nrow(orders_mat) - 1)])
  orlof = lf(mf(f, d))
  terms = names(sort(orlof, decreasing = TRUE))
  model = terms[1:size]
  ans$modeltest = names(sort(table(apply(orders_mat, 2, function(x) do.call(paste, as.list(sort(x[1:size]))))))[size])

  ans$freq = freq
  ans$size = size
  ans$terms = terms
  ans$model = model
  class(ans) = "IF"
  return(ans)
}


