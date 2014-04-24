bootstrap = function(set, mf, f, d) {
  mode = ifelse(is.null(set$mode), "np", match.arg(set$mode, c("np", "p")))
  B = ifelse(is.null(set$B), 100, set$B)

  # FIXME: bootstrap residual, not subject
  if (mode == "np") {
    if (class(d) == "data.frame") {
      return(replicate(B, d[sample(1:nrow(d), nrow(d), TRUE),], FALSE))
    }
  }
  if (mode == "p") {
    # TODO: only data.frame type support, add more type in future
    full = mf(f, d)
    ans = replicate(B, d, FALSE)
    bootsmp = switch(class(full)[1],
      lm = pboot_lm(B, full),
      glm = pboot_glm(B, full),
      lmerMod = pboot_lmerMod(B, full),
      glmerMod = pboot_lmerMod(B, full),
      NULL)

    if (is.null(bootsmp)) {
      warning("Switch to nonparametric bootstrap due to model is unsupported!")
      return(bootstrap(list(B = B), mf, f, d))
    } else {
      for (i in 1:B) {
        ans[[i]][,deparse(f[[2]])] = bootsmp[,i]
      }
      return(ans)
    }
  }
}

pboot_lmerMod = function(B, full) {
  # obtain parameters
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
  as.matrix(as.vector(fe) + re + rnorm(n * B, 0, tau))
}
