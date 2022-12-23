#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param ko_normalized
#' @return
#' @author rmflight
#' @export
anova_add_differential_test <- function(ko_normalized, ko_info, feature_info) {

  norm_matrix = log2(ko_normalized$normalized_data[, ko_info$client_identifier])
  ko_info = as.data.frame(ko_info)
  
  #ko_info$group_short = as.factor(ko_info$group_short)
  ko_info$nt_status = factor(ko_info$nt_status, levels = c("wt", "ko"))
  ko_info$ldlr_status = factor(ko_info$ldlr_status, levels = c("wt", "ko"))
  
  design = model.matrix(~ nt_status + ldlr_status, ko_info)
  
  log_fit = lmFit(norm_matrix, design)
  log_fit = eBayes(log_fit)
  log_limma = topTable(log_fit, coef = "nt_statusko", number = Inf, p.value = 1)
  
  log_limma$feature = rownames(log_limma)
  rownames(log_limma) = NULL
  out_data = dplyr::left_join(log_limma, dplyr::select(feature_info, -group_units), by = c("feature" = "feature_id"))
  design2 = cbind(ko_info$group_short, as.data.frame(design))
  names(design2)[1] = "group_short"
  list(stats = out_data, design = design2)
}


anova_group_differential_test <- function(ko_normalized, ko_info, feature_info) {
  
  norm_matrix = log2(ko_normalized$normalized_data[,  ko_info$client_identifier])
  ko_info = as.data.frame(ko_info)
  ko_info$group_short = as.factor(ko_info$group_short)
  design = model.matrix(~ 0 + group_short, ko_info)
  colnames(design) = gsub("group_short", "", colnames(design))
  
  contrasts = makeContrasts(
    NtKOvsWT = (ntKO_ldKO + ntKO_ldWT) / 2 - (ntWT_ldKO + ntWT_ldWT) / 2, levels = colnames(design))
  
  log_fit = lmFit(norm_matrix, design)
  cont_fit = contrasts.fit(log_fit, contrasts)
  cont_fit = eBayes(cont_fit)
  cont_limma = topTable(cont_fit, number = Inf, p.value = 1)
  
  cont_limma$feature = rownames(cont_limma)
  rownames(cont_limma) = NULL
  out_data = dplyr::left_join(cont_limma, dplyr::select(feature_info, -group_units), by = c("feature" = "feature_id"))
  design2 = cbind(ko_info$group_short, as.data.frame(design))
  names(design2)[1] = "group_short"
  list(stats = out_data, design = design2, contrasts = contrasts)
}
