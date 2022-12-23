#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param ko_info
#' @param ko_matrix_single
#' @return
#' @author rmflight
#' @export
remove_outliers <- function(ko_info, ko_matrix_single) {

  ko_cor = ici_kendalltau(t(ko_matrix_single))$cor
  ko_cor = ko_cor[ko_info$sample_id, ko_info$sample_id]
  ko_med = median_correlations(ko_cor, ko_info$group_short)
  
  ko_out = determine_outliers(ko_med, outlier_fraction = NULL)
  
  dplyr::left_join(ko_info, ko_out[, c("sample_id", "med_cor", "outlier")], by = "sample_id")

}
