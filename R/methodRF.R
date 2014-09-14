#' Restricted Fence
#'
#' Restricted Fence
#'
#' @param full formular of full model
#' @param data data
#' @param B number of bootstrap sample, parametric for lmer
#' @param grid grid for c
#' @param bandwidth bandwidth for kernel smooth function
#' @return list with whatever
#' @export

RF = function(full, data, groups, B = 100, grid = 101, bandwidth = NA, plot = FALSE, method = c("marginal", "conditional"), id = "id", cpus = 2) {
  if (is.na(bandwidth)) {
    stop("Must assign a bandwidth.")
  }


  method = match.arg(method)

  # 1st step
  cands = character(0)
  
  resp = as.character(full)[2]
  y = as.matrix(data[,resp])

  for (group in groups) {
    bs = bootstrap.RF(B, full, group, data)
    ms = findsubmodel.RF(full, group, data) 
    mf = function(model, bootsmp) {
      if (is.data.frame(bootsmp)) {
        return(model$model)
      }
      list(Q = t(bootsmp) %*% model$meat %*% bootsmp, 
           size = length(model$model))
    }
    lf = function(x) {
      x$Q
    }
    pf = function(x) {
      x$size
    }

    sfInit(parallel = TRUE, cpus = cpus) 
    sfExportAll()
    groupaf = adaptivefence(mf = mf, f = group, ms = ms, d = data, lf = lf, pf = pf, bs = bs, grid = grid, bandwidth = bandwidth)
    if (plot) {
      print(plot(groupaf))
      dev.off()
    }    
    groupc = groupaf$c

    Qs = sapply(ms, function(m) t(y) %*% m$meat %*% y)
    Qs = Qs - min(Qs)
    infence = Qs < groupc
    sizes = sapply(ms, function(x) length(x$model))
    sizes[!infence] = Inf
    Qs[sizes != min(sizes)] = Inf
    cands = c(cands, ms[[which.min(Qs)]]$model)
  }
  cands = unique(cands)

  # 2nd step
  if (method == "marginal") {
    resp = as.character(full)[2]
    full2nd = as.formula(paste0(resp, "~", gsub(" ", "+", do.call(paste, as.list(cands)))))
    bs = bootstrap.RF2(B, full2nd, data)
    ms = findsubmodel.RF2(resp, cands)
    sfInit(parallel = TRUE, cpus = cpus) 
    sfExportAll()
    res = adaptivefence(mf = lm, f = full2nd, ms = ms, d = data, lf = function(x) -logLik(x), 
      pf = function(x) length(attributes(terms(x))$term.labels), bs = bs, bandwidth = NA)
  }

  if (method == "conditional") {
    resp = as.character(full)[2]
    full2nd = as.formula(paste0(resp, "~", gsub(" ", "+", do.call(paste, as.list(cands))), "+(1|", id, ")"))
    res = fence.lmer(full2nd, data, B = B, grid = grid)
  }

  class(res) = "RF"
  res
}

bootstrap.RF = function(B, full, group, data) {
  resp = as.character(full)[2]
  fulltms = attributes(terms(full))$term.labels
  grouptms = attributes(terms(group))$term.labels
  groupxtms = fulltms[!(fulltms %in%  grouptms)]

  y = data[,resp]

  grouptms = as.list(grouptms)
  grouptms$sep = "+"
  groupxtms = as.list(groupxtms)
  groupxtms$sep = "+"

  X1 = model.matrix(as.formula(paste0(resp, "~0+", do.call(paste, grouptms))), data)
  X2 = model.matrix(as.formula(paste0(resp, "~0+", do.call(paste, groupxtms))), data)

  PX2O = diag(nrow(X2)) - X2 %*% ginv(t(X2) %*% X2) %*% t(X2)
  PX2OMX1 = PX2O - PX2O %*% X1 %*% ginv(t(X1) %*% PX2O %*% X1) %*% t(X1) %*% PX2O

  df = nrow(data) - length(fulltms)
  beta1 = ginv(t(X1) %*% PX2O %*% X1) %*% t(X1) %*% PX2O %*% matrix(y, ncol = 1)

  sigma = matrix(y, nrow = 1) %*% PX2OMX1 %*% matrix(y, ncol = 1) / df
  base = X1 %*% beta1
  lapply(1:B, function(i) matrix(base + rnorm(nrow(data), 0, sqrt(sigma)), ncol = 1))
}

findsubmodel.RF = function(full, group, data) {
  resp = as.character(full)[2]
  fulltms = attributes(terms(full))$term.labels
  grouptms = attributes(terms(group))$term.labels
  groupxtms = fulltms[!(fulltms %in%  grouptms)]

  y = data[,resp]

  grouptms = as.list(grouptms)
  grouptms$sep = "+"
  groupxtms = as.list(groupxtms)
  groupxtms$sep = "+"

  X1 = model.matrix(as.formula(paste0(resp, "~0+", do.call(paste, grouptms))), data)
  X2 = model.matrix(as.formula(paste0(resp, "~0+", do.call(paste, groupxtms))), data)

  PX2O = diag(nrow(X2)) - X2 %*% ginv(t(X2) %*% X2) %*% t(X2)

  res = ""
  for (i in 1:length(grouptms)) {
    res = append(lapply(res, function(x) c(x, grouptms[i])), lapply(res, function(x) x))
  }
  res = lapply(res, function(x) x[-1])[-length(res)]
  for (i in 1:length(res)) {
    X1 = as.matrix(data[,unlist(res[[i]][-length(res[[i]])])])
    res[[i]] = list(model = res[[i]], meat = PX2O - PX2O %*% X1 %*% ginv(t(X1) %*% PX2O %*% X1) %*% t(X1) %*% PX2O)
  }
  res
}

bootstrap.RF2 = function(B, full2nd, data) {
  resp = as.character(full2nd)[2]
  m = lm(full2nd, data)
  base = fitted(m)
  sigma = summary(m)$sigma
  newresp = lapply(1:B, function(i) base + rnorm(length(base), 0, sigma))
  ans = list()
  for (i in 1:B) {
    data[,resp] = newresp[[i]]
    ans[[i]] = data
  }
  ans
}

findsubmodel.RF2 = function(resp, tms) {
  res = paste(resp, "~1", sep = "")
  for (fix in tms) {
    res = as.vector(sapply(res, function(x) paste(x, c("", paste("+", fix, sep = "")), sep = "")))
  }
  lapply(res, as.formula)
}

