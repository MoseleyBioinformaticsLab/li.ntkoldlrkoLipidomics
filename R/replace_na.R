replace_na = function(data_matrix){
  zero_matrix = data_matrix
  zero_matrix[is.na(data_matrix)] = 0
  zero_matrix
}

