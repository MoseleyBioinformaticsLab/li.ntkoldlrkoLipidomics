#' Calculating our own log-fold-changes because we don't trust what {limma} is
#' spitting out because it's not necessarily making sure a certain number of things
#' are non-zero.
#'
#' .. content for \details{} ..
#'
#' @title
#' @param ko_normalized
#' @param ko_info_nooutlier
#' @return
#' @author rmflight
#' @export
calculate_lfc <- function(ko_normalized, ko_info_nooutlier, pairwise_tests, feature_info) {

  # tar_load(ko_normalized)
  # tar_load(ko_info_nooutlier)
  # tar_load(pairwise_tests)
  # # 
  ko_values = ko_normalized$normalized_data
  
  out_lfc = purrr::map(pairwise_tests, function(in_test){
    sample_numerator = ko_info_nooutlier %>%
      dplyr::filter(group_number %in% in_test[1]) %>%
      dplyr::pull(sample_id)
    sample_denominator = ko_info_nooutlier %>%
      dplyr::filter(group_number %in% in_test[2]) %>%
      dplyr::pull(sample_id)
    sample_df = data.frame(sample_id = c(sample_numerator, sample_denominator)) %>%
      dplyr::left_join(., ko_info_nooutlier[, c("sample_id", "group_number", "group_short")],
                       by = "sample_id")
    tmp_vals = ko_values[, sample_df$sample_id]
    tmp_vals_3 = t(visualizationQualityControl::keep_non_zero_percentage(t(tmp_vals),
                    sample_classes = sample_df$group_number,
                    keep_num = 0.33,
                    zero_value = ko_normalized$imputed_value))
    tmp_lfc = full_lfc(tmp_vals_3[, sample_numerator], tmp_vals_3[, sample_denominator])
    tmp_lfc$group_numerator = in_test[1]
    tmp_lfc$group_denominator = in_test[2]
    tmp_lfc$which = "all"
    
    other_lfc = alt_lfc(tmp_vals_3[, sample_numerator], tmp_vals_3[, sample_denominator], ko_normalized$imputed_value)
    other_lfc$p.adjust = stats::p.adjust(other_lfc$p.value, method = "BH")
    other_lfc$q.value = qvalue::qvalue(other_lfc$p.value)$qvalues
    
    other_lfc$group_numerator = in_test[1]
    other_lfc$group_denominator = in_test[2]
    other_lfc$which = "trimmed"
    
    tmp_lfc = dplyr::left_join(tmp_lfc, dplyr::select(feature_info, -group_units), by = "feature_id")
    other_lfc = dplyr::left_join(other_lfc, dplyr::select(feature_info, -group_units), by = "feature_id")
    
    list(all = tmp_lfc, trimmed = other_lfc)
  })

}

full_lfc = function(numerator, denominator){
  mean_numerator = rowMeans(log2(numerator))
  mean_denominator = rowMeans(log2(denominator))
  sample_df = data.frame(sample_id = c(colnames(numerator),
                                       colnames(denominator)),
                         group = rep(c(1, 2),
                                     c(ncol(numerator),
                                       ncol(denominator))))
  t_res = ttest_matrix(log(cbind(numerator, denominator)), sample_df$group)
  names(t_res)[1:3] = c("diff", "numerator", "denominator")
  lfc_res = data.frame(feature_id = rownames(numerator),
             logFC = mean_numerator - mean_denominator)
  dplyr::left_join(lfc_res, t_res, by = c("feature_id" = "feature"))
  
}

alt_lfc = function(numerator, denominator, imputed_value){
  all_lfc = purrr::map_dfr(seq(1, nrow(numerator)), function(irow){
    #message(irow)
    alt_single_lfc(numerator[irow, , drop = FALSE], denominator[irow, , drop = FALSE], imputed_value)
  })
}

alt_single_lfc = function(num_vector, den_vector, imputed_value){
  null_test = structure(list(feature_id = "NA", logFC = NA, 
                             diff = NA, 
                             mean_numerator = NA,
                             mean_denominator = NA,
                             numerator = NA, 
                             denominator = NA, 
                             statistic = c(t = NA), 
                             p.value = NA, 
                             parameter = c(df = NA), 
                             conf.low = NA, 
                             conf.high = NA, 
                             method = "Welch Two Sample t-test", 
                             alternative = "two.sided"), 
                        class = "data.frame", 
                        row.names = c(1L))
  if (sum(num_vector == imputed_value) == length(num_vector)) {
    num_use = as.vector(num_vector)
  } else {
    num_use = num_vector[!(num_vector == imputed_value), drop = TRUE]
  }
  if (sum(den_vector == imputed_value) == length(den_vector)) {
    den_use = as.vector(den_vector)
  } else {
    den_use = den_vector[!(den_vector == imputed_value), drop = TRUE]
  }
  
  t_res = try(broom::tidy(stats::t.test(log(num_use), log(den_use))))
  if (!inherits(t_res, "try-error")) {
    names(t_res)[1:3] = c("diff", "numerator", "denominator")
    t_res$feature = rownames(num_vector)[1]
    lfc = data.frame(feature_id = rownames(num_vector)[1],
                     logFC = mean(log2(num_use)) - mean(log2(den_use)),
                     mean_numerator = mean(log2(num_use)),
                     mean_denominator = mean(log2(den_use)))
    out_res = dplyr::left_join(lfc, t_res, by = c("feature_id" = "feature"))
  } else {
    #message(rownames(num_vector)[1])
    out_res = null_test
    out_res$feature_id = rownames(num_vector)[1]
    out_res$logFC = mean(log2(num_use)) - mean(log2(den_use))
    out_res$mean_numerator = mean(log2(num_use))
    out_res$mean_denominator = mean(log2(den_use))
  }
  out_res = tibble::as_tibble(out_res)
  out_res$numerator = list(log2(num_use))
  out_res$denominator = list(log2(den_use))
  
  out_res
}
