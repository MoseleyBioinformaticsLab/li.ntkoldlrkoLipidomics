#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param differential_results
#' @param sig_cutoff
#' @return
#' @author rmflight
#' @export
filter_sig <- function(differential_results, sig_cutoff = 0.01) {

  differential_results %>%
    dplyr::filter(adj.P.Val <= sig_cutoff)

}
