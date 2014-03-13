## model related
# mf:   model function
# f:    formular for full model
# ms:   the set of candidates model
# d:    data
# lf:   lack of fit function to MINIMIZE

## bootstrap related
# bootset: bootstrap setting
#          B:    number of bootstrap sample
#          mode: parametric("p") or nonparametric("np")
# bs:      bootstrap sample, default setting is nonparametric bootstrap

adaptivefence = function(
  # model and lack of fit related
  mf, f, ms = NULL, d, lf,
  # bootstrap related
  bootset = list(), bs = NULL,
  # fence related
  grid = 51, pickf = md, sizef = NULL) {

  ans = list(full = f, ms = ms, bootset = bootset, pickf = pickf, sizef = sizef)

  if (is.null(ms)) {
    ms = findallsubmodel(f)
  }

  if (is.null(bs)) {
    bs = bootstrap(bootset, mf, f, d)
  }

  B = ifelse(is.null(bootset$B), 100, bootset$B)

  # Is this one necessary?
  # lofbaseline = min(sapply(bs, function(b) lf(mf(f, b))))

  model_size = sapply(ms, sizef, d)
  ans$modelsize = model_size

  lack_of_fit_matrix = sapply(ms, function(m)
    sapply(bs, function(b) lf(mf(m, b))))

  lof_lower = min(lack_of_fit_matrix)
  lof_upper = max(lack_of_fit_matrix)

  total_lack_of_fit = colSums(lack_of_fit_matrix)
  force(total_lack_of_fit)

  cs = seq(lof_lower, lof_upper, length.out = grid)

  freq_mat = sapply(cs, function(c) {
    infence_freq= colSums(lack_of_fit_matrix <= c)
    candidate_index = which(infence_freq== max(infence_freq))
    champion = pickf(candidate_index, model_size, total_lack_of_fit)
    c(champion, infence_freq[champion] / B)
  })
  colnames(freq_mat) = cs
  rownames(freq_mat) = c("index", "frequency")
  ans$freq_mat = freq_mat

  cindex = peakw(freq_mat[2,])

  if (!is.na(cindex)) {
    ans$formula = ms[[freq_mat[1,cindex]]]
    ans$model = mf(ans$formula, d)
  }
  return(ans)
}
