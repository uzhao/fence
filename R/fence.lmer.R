#' Fence model selection (Linear Mixed Model)
#'
#' Fence model selection (Linear Mixed Model)
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

fence.lmer = function(
  full, data, B = 100, grid = 101, fence = c("adaptive", "nonadaptive"),
  cn = NA, REML = TRUE, bandwidth = NA, cpus = 2) {

  fence = match.arg(fence)
  if (fence == "adaptive" & !is.na(cn) |
      fence == "nonadaptive" & is.na(cn)) {
    stop("Adaptive agreement doesn't match!")
  }

  # find all candidate submodels
  ms = findsubmodel.lmer(full)
  # model fit function
  mf = function(...) lmer(..., REML = REML)
  # lack of fit function
  lf = function(res) -logLik(res)
  # pick up function
  pf = size.lmer

  if (fence == "nonadaptive") {
    return(nonadaptivefence(mf = mf, f = full, ms = ms, d = data, lf = lf, pf = pf,
      cn = cn))
  }
  if (fence == "adaptive")    {
    sfInit(parallel = TRUE, cpus = cpus) 
    sfExportAll()
    sfLibrary(lme4)
    return(   adaptivefence(mf = mf, f = full, ms = ms, d = data, lf = lf, pf = pf,
      bs = bootstrap.lmer(B, full, data, REML), grid = grid, bandwidth = bandwidth))
  }
}

findsubmodel.lmer = function(full) {
  resp = as.character(full)[2]
  tms = attributes(terms(full))$term.labels
  fr = grepl("\\|", tms)
  fixs = tms[!fr]
  rans = tms[fr]
  rans = paste("(", rans, ")", sep = "")
  res = paste(resp, "~0+", rans, sep = "")
  for (fix in fixs) {
    res = as.vector(sapply(res, function(x) paste(x, c("", paste("+", fix, sep = "")), sep = "")))
  }
  lapply(res, as.formula)
}

bootstrap.lmer = function (B, f, data, REML) {
  full = lmer(formula = f, data = data, REML = REML)
  ans = replicate(B, data, FALSE)
  X = full@pp$X
  beta = fixef(full)
  Z = t(as.matrix(full@pp$Zt))
  # random intercept only
  # can't handle random slope or interactive term for now
  # TODO: sigma = full@pp$Lambdat * tau
  #       sigma = sigma %*% t(sigma)
  tau = attr(VarCorr(full), "sc")
  sigmas = diag(as.matrix(full@pp$Lambdat) * tau)
  n = nrow(X)

  # generate
  fe = X %*% beta
  a = matrix(rnorm(B * length(sigmas), 0, sigmas), nrow = length(sigmas))
  re = Z %*% a
  bootsmp = as.matrix(as.vector(fe) + re + rnorm(n * B, 0, tau))
  for (i in 1:B) {
    ans[[i]][,deparse(f[[2]])] = bootsmp[,i]
  }
  ans
}

size.lmer = function(res) {
  length(fixef(res))
}
