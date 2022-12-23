#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param ko_normalized
#' @param ko_info_nooutlier
#' @return
#' @author rmflight
#' @export
calculate_summaries = function(ko_normalized, ko_info_nooutlier) {

  split_group = split(ko_info_nooutlier$sample_id, ko_info_nooutlier$group_number)
  ko_matrix = ko_normalized$normalized_data
  ko_matrix[ko_matrix == ko_normalized$imputed_value] = NA
  
  group_summary = purrr::imap_dfr(split_group, function(in_group, group_number){
    group_matrix = ko_matrix[, in_group]
    group_means = rowMeans(group_matrix, na.rm = TRUE)
    group_sd = apply(group_matrix, 1, sd, na.rm = TRUE)
    group_rsd = group_sd / group_means
    group_n = apply(group_matrix, 1, function(.x){sum(!is.na(.x))})
    
    out_data = data.frame(feature_id = rownames(ko_matrix),
                          mean = group_means,
                          sd = group_sd,
                          rsd = group_rsd,
                          n = group_n,
                          group = group_number)
    out_data
  })
  group_summary

}
