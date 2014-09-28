#' Invisible Fence gene selection
#'
#' Invisible Fence gene selection
#'
#' @param gene     a matrix, each column for a subject
#' @param response a factor vector, each element for a subject
#' @param genename a character vector
#' @param geneset  a character list
#' @param minsize  a integer, only use geneset which size larger or equal to this
#' @param B        a integer, number of bootstrap samples
#' @return list with whatever
#' @export

# Two class unpaired type in GSA
IF.GSA = function(
  gene, response,
  genenames, genesets, minsize = 0,
  B = 100, limitboot = NA) {
    ans = IFbase.GSA(gene, response, genenames, genesets, minsize, B, limitboot)
    ans$RIF = NULL
    ans
}

#' Relative Invisible Fence gene selection
#'
#' Relative Invisible Fence gene selection
#'
#' @param gene     a matrix, each column for a subject
#' @param response a factor vector, each element for a subject
#' @param genename a character vector
#' @param geneset  a character list
#' @param minsize  a integer, only use geneset which size larger or equal to this
#' @param B        a integer, number of bootstrap samples
#' @return list with whatever
#' @export

# Two class unpaired type in GSA
RIF.GSA = function(
  ){
  "Hello"
  # gene, response,
  # genenames, genesets, minsize = 0,
  # B = 100, limitboot = NA) {
  #   ans = IFbase.GSA(gene, response, genenames, genesets, minsize, B, limitboot)
  #   ans$IF = NULL
  #   ans
}

IFbase.GSA = function(gene, response, genenames, genesets, minsize, B, limitboot) {
  response = factor(response)
  groupa = response == levels(response)[1]
  groupb = !groupa
  na = sum(groupa)
  nb = sum(groupb)

  # remove too small genesets definition
  genesets = genesets[sapply(genesets, length) >= minsize]
  # remove duplicated genesets definition
  dup = duplicated(sapply(genesets, function(geneset) {
    do.call(paste, as.list(sort(geneset)))
  }))
  if (sum(dup) > 1) {
    cat("These genesets are removed due to duplicate!\n")
    print(names(genesets)[dup])
  }
  genesets = genesets[!dup]

  genesetsrow = lapply(genesets, function(geneset) match(geneset, genenames))
  anymissgene = sapply(genesetsrow, function(genesetrow) any(is.na(genesetrow)))
  allmissgene = sapply(genesetsrow, function(genesetrow) all(is.na(genesetrow)))
  if (any(anymissgene != allmissgene)) {
    cat("Some genesets missing part of gene(s) data, please check the data.\n")
  }
  if (any(anymissgene)) {
    cat("Some genesets missing all genes data, will be removed.\n")
  }

  genesets = genesets[!allmissgene]
  genesetsrow = genesetsrow[!allmissgene]
  cat("Totally ", length(genesets)," gene sets remain.\n")

  rawscores = get_geneset_scores(gene, response, genesetsrow)
  genesetsneedboot = rank(-rawscores$geneset_scores) <= limitboot
  rowsneedboot = na.omit(unique(unlist(genesetsrow[genesetsneedboot])))

  genesets = genesets[genesetsneedboot]
  genesetsrow = genesetsrow[genesetsneedboot]
  nogs = length(genesets)

  Bsample = gene
  Bsample_genesetorder = replicate(B, {
    Bsample[,groupa] = gene[,sample(which(groupa), na, TRUE)]
    Bsample[,groupb] = gene[,sample(which(groupb), nb, TRUE)]
    get_geneset_scores(Bsample, response, genesetsrow)$order
  })

  if (is.na(limitboot)) {
    limitboot = nogs
  }

  value = matrix(NA, 2, limitboot)
  for (size in 1:limitboot) {
    cand = apply(matrix(Bsample_genesetorder[1:size,], nrow = size), 2, function(x) do.call(paste, as.list(sort(x))))
    raw = sort(table(cand), decreasing = TRUE)[1]

    value[1, size] = as.numeric(raw)
    # freq[2, size] = relative(as.numeric(raw))
  }

  IF = peakglobal(value[1,1:(limitboot-1)])
  # RIF = freq[[which.min(sapply(freq, function(x) x$relative))]]

  ans = list(genesets = genesets, value = value, IF = IF)
  # ans = list(genesets = genesets, freq = freq, IF = IF, RIF = RIF)
  class(ans) = "IF.GSA"
  ans
}

get_gene_scores = function(gene, response) {
  ni = table(response)
  groupa = gene[,response == levels(response)[1]]
  groupb = gene[,response == levels(response)[2]]

  meana = rowMeans(groupa)
  meanb = rowMeans(groupb)

  svara = rowSums(sweep(groupa, 1, meana, '-') ^ 2)
  svarb = rowSums(sweep(groupb, 1, meanb, '-') ^ 2)

  sd = sqrt((svara + svarb) * sum(1 / ni) / sum(ni - 2))

  (meana - meanb) / sd
}

maxmean = function (gene_scores, genesetsrow) {
  i = 0
  replicate(length(genesetsrow), {
    i <<- i + 1
    cur_score = gene_scores[genesetsrow[[i]]]
    cur_score[is.na(cur_score)] = 0
    pos = mean((cur_score + abs(cur_score)) / 2)
    neg = mean((cur_score - abs(cur_score)) / 2)
    max(pos, -neg)
    })
}

get_geneset_scores = function(gene, response, genesetsrow) {
  gene_scores = get_gene_scores(gene, response)
  geneset_scores = maxmean(gene_scores, genesetsrow)
  names(geneset_scores) = names(genesetsrow)
  list(gene_scores = gene_scores,
       geneset_scores = geneset_scores,
       order = order(-geneset_scores))
}

# hugepbirthday = function(n, classes, coincident) {
#   k = coincident
#   c = classes    # could be huge
#   if (k < 2) 
#     return(1)
#   if (k == 2)  {
#     ans = as.bigq(1)
#     for (i in 1:(k - 1)) {
#       ans = ans * (c - i) / c 
#     }
#     return(ans)
#   }
#   if (k > n) 
#     return(0)
# }

# relative = function(B, times, nogs, size) {
#   ans = hugepbirthday(B, chooseZ(nogs, size), times)
#   names(ans) = NULL
#   ans
# }
