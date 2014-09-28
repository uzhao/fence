#' Fence model selection (Small Area Estmation)
#'
#' Fence model selection (Small Area Estmation)
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

fence.sae = function(
  full, data, B = 100, grid = 101, fence = c("adaptive", "nonadaptive"),
  cn = NA, method = c("F-H", "NER"), D = NA, REML = FALSE, bandwidth = NA, cpus = 2) {

  fence = match.arg(fence)
  if (fence == "adaptive" & !is.na(cn) |
      fence == "nonadaptive" & is.na(cn)) {
    stop("Adaptive agreement doesn't match!")
  }

  method = match.arg(method)
  if (method == "NER") {
    return(fence.lmer(full, data, B, grid, fence, cn, REML, bandwidth, cpus))
  }
  # if (method == "F-H") {
  #   return(fence.fh(full, data, B, grid, fence, cn, D, bandwidth, cpus))
  # }
  # find all candidate submodels
  ms = findsubmodel.fh(full)
  # model fit function
  mf = function(m, b) eblupFH(formula = m, vardir = D, data = b, method = "FH")
  # lack of fit function
  lf = function(res) -res$fit$goodness[1]
  # pick up function
  pf = function(res) nrow(res$fit$estcoef) - 1

  if (fence == "nonadaptive") {
    return(nonadaptivefence(mf = mf, f = full, ms = ms, d = data, lf = lf, pf = pf,
      cn = cn))
  }

  if (fence == "adaptive")    {
    sfInit(parallel = TRUE, cpus = cpus) 
    sfExportAll()
    sfLibrary(sae)
    return(   adaptivefence(mf = mf, f = full, ms = ms, d = data, lf = lf, pf = pf,
      bs = bootstrap.fh(B, full, data, D), grid = grid, bandwidth = bandwidth))
  }
}

findsubmodel.fh = function(full) {
  resp = as.character(full)[2]
  tms = attributes(terms(full))$term.labels
  res = paste(resp, "~", sep = "")
  for (tm in tms) {
    res = as.vector(sapply(res, function(x) paste(x, c("", paste("+", tm, sep = "")), sep = "")))
  }
  res = gsub("~ +", "~", res)
  res = res[res != "y~"]
  lapply(res, as.formula)
}

bootstrap.fh = function(B, full, data, D) {
  X = model.matrix(full, data)
  model = eblupFH(formula = full, vardir = D, data = data, method = "FH")
  beta = model$fit$estcoef[,1]
  tau = model$fit$refvar

  ans = replicate(B, data, FALSE)
  bootsmp = as.vector(X %*% beta) + replicate(B, rnorm(nrow(X), 0, sqrt(tau + data[,deparse(substitute(D))])))
  for (i in 1:B) {
    ans[[i]][,deparse(full[[2]])] = bootsmp[,i]
  }
  ans
}