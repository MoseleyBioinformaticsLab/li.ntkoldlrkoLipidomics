#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param ko_normalized
#' @return
#' @author rmflight
#' @export
basic_differential_test <- function(ko_normalized, ko_info, feature_info, tests) {

  norm_matrix = log2(ko_normalized$normalized_data)
  ko_info = as.data.frame(ko_info)
  
  test_res = purrr::imap(tests, function(in_groups, test_id){
    numerator_samples = ko_info %>%
      dplyr::filter(group_number %in% in_groups[1]) %>%
        dplyr::pull(sample_id)
    denominator_samples = ko_info %>%
      dplyr::filter(group_number %in% in_groups[2]) %>%
      dplyr::pull(sample_id)
    sample_info = data.frame(sample_id = c(numerator_samples,
                                             denominator_samples),
                               disease = c(rep(1, length(numerator_samples)),
                                           rep(0, length(denominator_samples))))
    tmp_matrix = norm_matrix[, sample_info$sample_id]
    
    out_data = run_diff_test(tmp_matrix, sample_info)
    out_data$test_id = gsub("x", "", test_id)
    out_data$genotype_numerator = ko_info[ko_info$group_number == in_groups[1], "group_short"][1]
    out_data$genotype_denominator = ko_info[ko_info$group_number == in_groups[2], "group_short"][1]
    out_data = dplyr::left_join(out_data, dplyr::select(feature_info, -group_units), by = c("feature" = "feature_id"))
    out_data
  })
  
  test_res
  

}

run_diff_test = function(data_matrix, sample_info){
  
  disease_status = factor(sample_info$disease)
  
  design = model.matrix(~disease_status)
  
  log_fit = lmFit(data_matrix, design)
  log_fit = eBayes(log_fit)
  log_limma = topTable(log_fit, coef = "disease_status1", number = Inf, p.value = 1)
  log_limma$feature = rownames(log_limma)
  rownames(log_limma) = NULL
  log_limma
}

write_basic_stats_results = function(results_list, feature_intensities, file_location){
  run_comparisons = purrr::map_chr(results_list, function(.x){
    paste0(.x$genotype_numerator[1], "__", .x$genotype_denominator[1])
  })
  
  names(results_list) = run_comparisons
  results_list$feature_intensities = feature_intensities
  write.xlsx(results_list, file = file_location, asTable = TRUE)
  file_location
}

basic_ttest = function(ko_normalized, ko_info, feature_info, tests){
  norm_matrix = log2(ko_normalized$normalized_data)
  ko_info = as.data.frame(ko_info)
  
  test_res = purrr::imap(tests, function(in_groups, test_id){
    numerator_samples = ko_info %>%
      dplyr::filter(group_short %in% in_groups[1]) %>%
      dplyr::pull(sample_id)
    denominator_samples = ko_info %>%
      dplyr::filter(group_short %in% in_groups[2]) %>%
      dplyr::pull(sample_id)
    t_res = purrr::map_dfr(seq(1, nrow(norm_matrix)), function(in_row){
      #message(in_row)
      x = norm_matrix[in_row, numerator_samples]
      y = norm_matrix[in_row, denominator_samples]
      if (sum(c(x, y) == min(norm_matrix)) == (length(x) + length(y))) {
        return(NULL)
      }
      tmp_res = broom::tidy(t.test(x, y))
      tmp_res$feature_id = rownames(norm_matrix)[in_row]
      names(tmp_res)[1:3] = c("logFC", "mean_x", "mean_y")
      tmp_res
    })
    t_res$test_id = test_id
    t_res$genotype_numerator = in_groups[1]
    t_res$genotype_denominator = in_groups[2]
    t_res
  })
  
  test_res
}
