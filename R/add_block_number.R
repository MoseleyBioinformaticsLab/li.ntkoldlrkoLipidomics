#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param ko_info
#' @return
#' @author rmflight
#' @export
add_block_number <- function(ko_info) {

  group_run = rle(ko_info$group_number)
  block_groups = rep(seq(1, length(group_run$values)), group_run$lengths)
  ko_info$block = block_groups
  ko_info

}
