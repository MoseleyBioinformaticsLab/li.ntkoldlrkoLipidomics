#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param ko_noormalized
#' @param ko_info
#' @param feature_info
#' @return
#' @author rmflight
#' @export
sdams_differential_test <- function(ko_normalized, ko_info, feature_info, pairwise_tests) {
  #tar_load(ko_normalized)
  #tar_load(ko_info)
  #tar_load(feature_info)
  #tar_load(pairwise_tests)
  norm_matrix = ko_normalized$normalized_data[, ko_info$sample_id]
  impute_value = ko_normalized$imputed_value
  norm_matrix[norm_matrix == impute_value] = 0
  ko_info = as.data.frame(ko_info)
  
  ko_seo = createSEFromMatrix(norm_matrix, ko_info)
  
  test_res = purrr::imap(pairwise_tests, function(in_groups, test_id){
    message(test_id)
    match_groups = ko_seo$group_number %in% in_groups
    tmp_seo = ko_seo[, match_groups]
    
    filter_assay = t(visualizationQualityControl::keep_non_zero_percentage(t(assay(tmp_seo)), sample_classes = tmp_seo$group_number, keep_num = 3))
    
    tmp_seo = tmp_seo[rownames(filter_assay), ]
    group_match = vector("character", ncol(tmp_seo))
    group_match = rep("numerator", length(group_match))
    group_match[tmp_seo$group_number %in% in_groups[2]] = "other"
    colData(tmp_seo) = NULL
    tmp_seo$test = as.factor(group_match)
    all_zero = apply(assay(tmp_seo), 1, function(.x){sum(.x == 0) == length(match_groups)})
    tmp_seo = tmp_seo[!all_zero, ]
    
    sda_res = SDA(tmp_seo, cleanData = FALSE, correctPValues = FALSE)
    sda_res2 = run_qvalues(sda_res)
    sda_df = as.data.frame(sda_res2)
    sda_df$test = test_id
    sda_df
  })
  test_res
}

run_qvalues = function(sda_list){
  match_p_q = c("pv_gamma" = "qv_gamma",
                "pv_beta" = "qv_beta",
                "pv_2part" = "qv_2part")
  for (ip in names(match_p_q)) {
    tmp_q = try(qvalue::qvalue(sda_list[[ip]]))
    if (inherits(tmp_q, "try-error")) {
      tmp_q = try(qvalue::qvalue(sda_list[[ip]], pi0 = 1))
      if (inherits(tmp_q, "qvalue")) {
        sda_list[[match_p_q[ip]]] = tmp_q$qvalues
      } 
    } else {
      sda_list[[match_p_q[ip]]] = tmp_q$qvalues
    }
  }
  sda_list
}
