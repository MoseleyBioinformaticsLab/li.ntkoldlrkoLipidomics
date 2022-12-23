#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param ko_matrix_single
#' @return
#' @author rmflight
#' @export
normalize_data = function(ko_matrix_single) {

  na_matrix = ko_matrix_single
  na_matrix[ko_matrix_single == 0] = NA
  na_df = as.data.frame(na_matrix)
  na_long = na_df %>%
    tidyr::pivot_longer(cols = everything(),
                        names_to = "sample",
                        values_to = "intensity")
  norm_factors = na_long %>%
    dplyr::mutate(intensity = log(intensity)) %>%
    mu_calc_normalization() %>%
    exp()
  norm_matrix = mu_apply_normalization(na_matrix, norm_factors)
  
  norm_long = as.data.frame(norm_matrix) %>%
    tidyr::pivot_longer(cols = everything(),
                        names_to = "sample_id",
                        values_to = "intensity") %>%
    dplyr::filter(!is.na(intensity)) %>%
    dplyr::mutate(intensity = log(intensity))
  
  stats_norm = boxplot.stats(norm_long$intensity)
  min_val = min(norm_long$intensity)
  thresh_value = min_val + (min_val * 0.01)
  exp_value = exp(thresh_value)
  norm_matrix[is.na(norm_matrix)] = exp_value
  
  list(normalized_data = norm_matrix, imputed_value = exp_value)
}
