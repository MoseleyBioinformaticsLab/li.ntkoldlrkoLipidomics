#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param manifest_file
#' @return
#' @author rmflight
#' @export
get_manifest_data <- function(manifest_file = here::here("data/Metabolon_Manifest_UKY Jing Li-081821.xlsx")) {

  manifest_data = as.data.frame(readxl::read_xlsx(manifest_file,
                                                  sheet = "Sample Manifest",
                                                  range = readxl::anchored("A14", dim = c(50, 18)),
                                                  col_names = TRUE)) %>%
    janitor::clean_names() %>%
      dplyr::mutate(unique_tube_label_id = substr(unique_tube_label_id, 3, 8),
                    box_number = substr(box_number, 3, 7))
  manifest_data
}
