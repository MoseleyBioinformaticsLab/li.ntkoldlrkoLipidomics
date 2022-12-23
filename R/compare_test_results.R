#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param differential_results
#' @param genotype_results
#' @return
#' @author rmflight
#' @export
compare_test_results <- function(differential_results, genotype_results) {

  all_diff = purrr::map_df(differential_results, ~ .x)
  all_diff = all_diff %>%
    dplyr::mutate(feature_num_denom = paste0(feature, ":", genotype_numerator, ":", genotype_denominator))
  
  all_genotype = purrr::map_df(genotype_results, ~ .x)
  all_genotype = all_genotype %>%
    dplyr::mutate(feature_num_denom = paste0(feature_id, ":", genotype_numerator, ":", genotype_denominator))
  
  combine_results = dplyr::left_join(all_genotype, all_diff, by = "feature_num_denom", suffix = c(".genotype", ".diff"))
  
  combine_results %>%
    ggplot(aes(x = logFC.genotype, y = logFC.diff)) +
    geom_point() +
    geom_abline(slope = 1, color = "red") +
    coord_equal()
}
