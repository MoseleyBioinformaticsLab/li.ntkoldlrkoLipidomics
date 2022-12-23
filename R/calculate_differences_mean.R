#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param ko_normalized
#' @return
#' @author rmflight
#' @export
calculate_differences_mean <- function(ko_normalized) {

  log_norm = log(ko_normalized$normalized_data)
  log_norm[log_norm == log(ko_normalized$imputed_value)] = NA
  
  row_mean = rowMeans(log_norm, na.rm = TRUE)
  row_mean_matrix = matrix(row_mean, nrow = nrow(log_norm), ncol = ncol(log_norm), byrow = FALSE)
  diff_mean = log_norm - row_mean_matrix
  diff_mean
}
