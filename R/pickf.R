# minimum lack of fit criterion
# should not use directly, just as backup plan for tie
mlof = function(candidate_index, lack_of_fit) {
  candidate_index[which.min(lack_of_fit[candidate_index])]
}

# minimum-dimension criterion
md = function(candidate_index, model_size, lack_of_fit) {
  # pick the model with minimum-dimension
  candidate_size = model_size[candidate_index]
  candidate_index[candidate_size == min(candidate_size)]
}
