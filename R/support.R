findsubmodel = function(full, maxsize = 5) {
  
}

peakw = function(freq, req = 2) {
  of = order(freq, decreasing = TRUE)
  for (i in 1:length(of)) {
    index = of[i]
    if ((index - req < 1) | (index + req > length(freq))) {
      next
    }
    start = index - req
    end = index + req
    # FIXME: can not identify more than one peak with same freq
    if (all(freq[start:end] <= freq[index]) && any(freq[index:length(freq)] < freq[index])) {
      return(index)
    }
  }
  warning("No peak identified!")
  return(NA)
}

plot.AF = function(res) {
  tmp = data.frame(c = as.numeric(colnames(res$freq_mat)), 
                   p = res$freq_mat[2, ],
                   m = as.factor(res$freq_mat[1, ]))
  p = ggplot(tmp) + geom_point(aes(x = c, y = p, colour = m)) + geom_line(aes(x = c, y = p), linetype="dashed") + geom_line(aes(x = c, y = p, colour = m)) + ylim(0, 1)
  p
}

summary.AF = function(res) {
  print(res$formula)
}
