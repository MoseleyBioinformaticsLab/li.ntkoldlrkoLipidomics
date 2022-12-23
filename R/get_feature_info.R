get_feature_info = function(concentration_data){
  feature_df = concentration_data[, c(1:3)]
  feature_df$feature_id = paste0("f_", stringr::str_pad(seq(1, nrow(feature_df)), width = 4, pad = "0"))
  feature_df = janitor::clean_names(feature_df)
  feature_df
}

df_to_matrix = function(concentration_data, feature_info, sample_info){
  concentration_matrix = as.matrix(concentration_data[, seq(4, ncol(concentration_data))])
  colnames(concentration_matrix)  = sample_info$client_identifier
  rownames(concentration_matrix) = feature_info$feature_id
  
  concentration_matrix
}
