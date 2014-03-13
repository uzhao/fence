# minimum sum of bootstrap lack of fit criterion
msblof = function(candidate_index, model_size, total_lack_of_fit) {
  candidate_index[which.min(total_lack_of_fit[candidate_index])]
}

# minimum-dimension criterion
md = function(candidate_index, model_size, total_lack_of_fit) {
  # pick the model with minimum-dimension
  candidate_size = model_size[candidate_index]
  candidate_index = candidate_index[candidate_size == min(candidate_size)]
  # within the model with minimum-dimension pick the one with minimum sum of bootstrap lack of fit
  msblof(candidate_index, model_size, total_lack_of_fit)
}
