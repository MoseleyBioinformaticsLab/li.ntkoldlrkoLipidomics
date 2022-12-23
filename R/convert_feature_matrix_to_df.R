#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param ko_normalized
#' @param sample_info
#' @return
#' @author rmflight
#' @export
convert_feature_matrix_to_df <- function(ko_normalized, sample_info) {

  t_norm = t(ko_normalized$normalized_data) %>%
    as.data.frame()
  t_norm$sample_id = rownames(t_norm)
  t_norm = dplyr::left_join(t_norm, sample_info, by = "sample_id")
  t_norm = t_norm %>%
    dplyr::arrange(group_short)
  t_norm
}
