findallsubmodel = function(full) {
  
}

peakw = function(freq, req = 2) {
  of = order(freq, decreasing = TRUE)
  for (i in 1:length(of)) {
    index = of[i]
    start = max(1, index - req)
    end = min(length(freq), index + req)
    # FIXME: can not identify more than one peak with same freq
    if (all(freq[start:end] <= freq[index])) {
      return(index)
    }
  }
  warning("No peak identified!")
  return(NA)
}
