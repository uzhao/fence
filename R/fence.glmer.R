#' Fence model selection (Generalized Linear Mixed Model)
#'
#' Fence model selection (Generalized Linear Mixed Model)
#'
#' @param full formular of full model
#' @param data data
#' @param family binary("b") or poisson("p")
#' @param B number of bootstrap sample, parametric for lmer
#' @param grid grid for c
#' @param fence fence method to be used, (adaptive, nonadaptive, more to add)
#' @param cn cn for nonadaptive
#' @param bandwidth bandwidth for kernel smooth function
#' @return list with whatever
#' @note more link to be added
#' @author Jianyang Zhao
#' @export

fence.glmer = function(
  full, data, family = c("b", "p"), B = 100, grid = 101, fence = c("adaptive", "nonadaptive"),
  cn = NA, bandwidth = NA, cpus = 2) {

  family = match.arg(family)
  family = switch(family,
    b = binomial,
    p = poisson)

  fence = match.arg(fence)
  if (fence == "adaptive" & !is.na(cn) |
      fence == "nonadaptive" & is.na(cn)) {
    stop("Adaptive agreement doesn't match!")
  }

  # find all candidate submodels
  ms = findsubmodel.glmer(full)
  # model fit function
  mf = function(...) glmer(..., family = family)
  # lack of fit function
  lf = function(res) -logLik(res)
  # pick up function
  pf = size.glmer

  if (fence == "nonadaptive") {
    return(nonadaptivefence(mf = mf, f = full, ms = ms, d = data, lf = lf, pf = pf,
      cn = cn))
  }
  if (fence == "adaptive")    {
    sfInit(parallel = TRUE, cpus = cpus) 
    sfExportAll()
    sfLibrary(lme4)
    return(   adaptivefence(mf = mf, f = full, ms = ms, d = data, lf = lf, pf = pf,
      bs = bootstrap.glmer(B, full, data, family), grid = grid, bandwidth = bandwidth))
  }
}

# same to lmer
findsubmodel.glmer = function(full) {
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

bootstrap.glmer = function (B, f, data, family) {
  full = glmer(formula = f, data = data, family = family)
  ans = replicate(B, data, FALSE)

  X = full@pp$X
  beta = fixef(full)
  Z = t(as.matrix(full@pp$Zt))
  # random intercept only
  # could be more complex than lmer
  # FIXME: need more careful check
  #        perhaps only work for one random effect
  sigmas = diag(as.matrix(full@pp$Lambdat))
  n = nrow(X)

  # generate
  fe = X %*% beta
  a = matrix(rnorm(B * length(sigmas), 0, sigmas), nrow = length(sigmas))
  re = Z %*% a
  bootsmp = family()$linkinv(as.matrix(as.vector(fe) + re))

  if (family()$family == "binomial") {
    bootsmp = matrix(rbinom(length(bootsmp), 1, bootsmp), nrow = nrow(bootsmp))
  }

  if (family()$family == "poisson") {
    bootsmp = matrix(rpois(length(bootsmp), bootsmp), nrow = nrow(bootsmp))
  }

  for (i in 1:B) {
    ans[[i]][,deparse(f[[2]])] = bootsmp[,i]
  }
  ans
}

# same to lmer
size.glmer = function(res) {
  length(fixef(res))
}
