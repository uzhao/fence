findsubmodel = function(full, maxsize = 5) {
  
}

pickbycn = function(res, cn) {
  model_vector = rep(NA, res$B)
  infence_matrix = res$Qd_matrix <= cn
  for (bi in 1:res$B) {
    b_infence = infence_matrix[bi,]
    b_lack = res$lack_of_fit_matrix[bi,]
    b_pick = res$pick_matrix[bi,]
    b_pick[!b_infence] = Inf
    b_pick = which(b_pick == min(b_pick))
    model_vector[bi] = b_pick[which.min(b_lack[b_pick])]
  }
  res$models[[as.integer(names(sort(table(model_vector), decreasing = TRUE))[1])]]
}

peakw = function(cs, freq, req = 2) {
  for (index in 1:length(freq)) {
    start = index - req
    end = index + req
    if (start < 1 | end > length(freq)) {
      next
    }
    if (all(freq[start:end] <= freq[index])) {
      return(index)
    }
  }
  warning("No peak identified!")
  return(NA)
}

plot.AF = function(res) {
  tmp = data.frame(c = as.numeric(colnames(res$freq_mat)), 
                   p = res$freq_mat[2, ],
                   sp = res$freq_mat[3, ],
                   m = as.factor(res$freq_mat[1, ]))
  p = ggplot(tmp) + 
    geom_point(aes(x = c, y = p, colour = m)) + 
    geom_line(aes(x = c, y = sp), linetype="dashed") + 
    geom_line(aes(x = c, y = sp, colour = m)) + 
    ylim(0, 1)
  p
}

summary.AF = function(res) {
  print(res$formula)
}
