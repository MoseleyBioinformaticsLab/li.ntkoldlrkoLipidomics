#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param ko_normalized
#' @param ko_info_nooutlier
#' @param pairwise_tests
#' @return
#' @author rmflight
#' @export
dasev_differential_test <- function(ko_normalized, ko_info_nooutlier,
                                    pairwise_tests) {

  ko_data = ko_normalized$normalized_data
  ko_data[ko_data == ko_normalized$imputed_value] = 0
  detection_limit = log(ko_normalized$imputed_value)
  
  dasev_results = purrr::imap(pairwise_tests, function(in_test, test_id){
    #message(test_id)
    tmp_info = ko_info_nooutlier %>%
      dplyr::filter(group_number %in% in_test)
    matrix_full = ko_data[, tmp_info$sample_id]
    matrix_data = t(keep_non_zero_percentage(t(matrix_full), tmp_info$group_number,
                                            keep_num = 0.33))
    tmp_info = tmp_info %>%
      dplyr::mutate(main_factor = dplyr::case_when(
        group_number == in_test[1] ~ 1,
        group_number == in_test[2] ~ 0
      ))
    matrix_status = as.matrix(tmp_info$main_factor)
    matrix_status = cbind(rep(1, nrow(tmp_info)), matrix_status)
    tmp_dasev = try(dasev_test(matrix_data, matrix_status, detection_limit))
    if (inherits(tmp_dasev, "try-error")) {
      return(NULL)
    } else {
      tmp_dasev$test = test_id
      return(tmp_dasev)
    }
    
  })
  dasev_results
}

dasev_test = function(matrix_data, matrix_status, detection_limit){
  tmp_res = DASEV(matrix_data, cov.matrix = matrix_status, test_cov = 2,
                  DL_method = "Fixed Value", DL_value = detection_limit)
  as.data.frame(tmp_res)
}
