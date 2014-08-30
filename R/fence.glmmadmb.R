#' Fence model selection (Generalized Linear Mixed Model by glmmadmb)
#'
#' Fence model selection (Generalized Linear Mixed Model by glmmadmb)
#'
#' @param full formular of full model
#' @param data data
#' @param family gamma("g") or negative binomial("n")
#' @param B number of bootstrap sample, parametric for lmer
#' @param grid grid for c
#' @param fence fence method to be used, (adaptive, nonadaptive, more to add)
#' @param cn cn for nonadaptive
#' @param bandwidth bandwidth for kernel smooth function
#' @return list with whatever
#' @note glmmadmb has more families but slower, more link need to be added and may fail in adaptive
#' @export

fence.glmmadmb = function(
  full, data, family = c("g", "n", "p"), zeroInflation = FALSE, B = 100, grid = 101, fence = c("adaptive", "nonadaptive"),
  cn = NA, bandwidth = NA) {

  family = match.arg(family)
  family = switch(family,
    g = "gamma",
    n = "nbinom", 
    p = "poisson")

  fence = match.arg(fence)
  if (fence == "adaptive" & !is.na(cn) |
      fence == "nonadaptive" & is.na(cn)) {
    stop("Adaptive agreement doesn't match!")
  }

  # find all candidate submodels
  ms = findsubmodel.glmmadmb(full)
  # model fit function
  mf = function(...) glmmadmb(..., family = family, zeroInflation = zeroInflation)
  # lack of fit function
  lf = function(res) -logLik(res)
  # pick up function
  pf = size.glmmadmb

  if (fence == "nonadaptive") {
    return(nonadaptivefence(mf = mf, f = full, ms = ms, d = data, lf = lf, pf = pf,
      cn = cn))
  }
  if (fence == "adaptive")    {
    return(   adaptivefence(mf = mf, f = full, ms = ms, d = data, lf = lf, pf = pf,
      bs = bootstrap.glmmadmb(B, full, data, family, zeroInflation), grid = grid, bandwidth = bandwidth))
  }
}

# same to lmer
findsubmodel.glmmadmb = function(full) {
  resp = as.character(full)[2]
  tms = attributes(terms(full))$term.labels
  fr = grepl("\\|", tms)
  fixs = tms[!fr]
  rans = tms[fr]
  rans = paste("(", rans, ")", sep = "")
  res = paste(resp, "~", rans, sep = "")
  for (fix in fixs) {
    res = as.vector(sapply(res, function(x) paste(x, c("", paste("+", fix, sep = "")), sep = "")))
  }
  lapply(res, as.formula)
}

bootstrap.glmmadmb = function (B, f, data, family, zeroInflation) {
  full = glmmadmb(formula = f, data = data, family = family, zeroInflation = zeroInflation)
  support = lmer(formula = f, data = data)
  ans = replicate(B, data, FALSE)

  X = support@pp$X
  beta = fixef(full)
  Z = t(as.matrix(support@pp$Zt))
  # random intercept only
  # could be more complex than lmer
  # FIXME: only work for one random effect
  sigmas = matrix(unlist(VarCorr(full)), nrow = ncol(Z))
  n = nrow(X)

  # generate
  fe = X %*% beta
  a = matrix(rnorm(B * length(sigmas), 0, sigmas), nrow = length(sigmas))
  re = Z %*% a

  if (family == "gamma") {
    bootsmp = exp(as.matrix(as.vector(fe) + re))
    bootsmp = matrix(rgamma(length(bootsmp), full$alpha, scale = bootsmp / full$alpha), nrow = nrow(bootsmp))
  }

  if (family == "nbinom") {
    bootsmp = exp(as.matrix(as.vector(fe) + re))
    bootsmp = matrix(rnbinom(length(bootsmp), full$alpha, mu = bootsmp), nrow = nrow(bootsmp))
  }

  if (family == "poisson") {
    bootsmp = exp(as.matrix(as.vector(fe) + re))
    bootsmp = matrix(rpois(length(bootsmp), lambda = bootsmp), nrow = nrow(bootsmp))
  }

  if (zeroInflation) {
    bootsmp[rbinom(length(bootsmp), 1, full$pz)] = 0
  }

  for (i in 1:B) {
    ans[[i]][,deparse(f[[2]])] = bootsmp[,i]
  }
  ans
}

# remove NA check to avoid warning
size.glmmadmb = function(res) {
  length(fixef(res))
}
