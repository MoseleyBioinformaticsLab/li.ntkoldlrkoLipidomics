get_concentrations = function(){
  readxl::read_xlsx(here::here("data/UNKY-03-20VWCLP+ PLASMA CLP 6-TAB_working.xlsx"), 
                    sheet = "Species Concentrations", skip = 29)
}


get_sample_info = function(){
  tmp_frame = as.data.frame(readxl::read_xlsx(here::here("data/UNKY-03-20VWCLP+ PLASMA CLP 6-TAB_working.xlsx"),
                    sheet = "Species Concentrations",
                    range = readxl::anchored("C6", dim = c(23, 19)),
                    col_names = FALSE))
  frame_names = tmp_frame[[1]]
  ncol_frame = ncol(tmp_frame)
  nrow_frame = nrow(tmp_frame)
  frame_classed = purrr::map_dfc(seq(1, nrow_frame), function(in_row){
    row_data = unlist(tmp_frame[in_row, seq(2, ncol_frame), drop = TRUE])
    names(row_data) = NULL
    tmp_num = as.numeric(row_data)
    if (sum(is.na(tmp_num)) == 0) {
      return(tmp_num)
    } else {
      return(row_data)
    }
  })
  names(frame_classed) = janitor::make_clean_names(frame_names)
  frame_classed$sample_id = frame_classed$client_identifier
  frame_classed = frame_classed %>%
    dplyr::mutate(nt_status = dplyr::case_when(
      grepl("Nt-/-", group_name) ~ "ko",
      TRUE ~ "wt"
    ),
    ldlr_status = dplyr::case_when(
      grepl("Ldlr-/-", group_name) ~ "ko",
      TRUE ~ "wt"
    ),
    group_short = dplyr::case_when(
      group_name %in% "Nt+/+; Ldlr+/+" ~ "ntWT_ldWT",
      group_name %in% "Nt-/-; Ldlr+/+" ~ "ntKO_ldWT",
      group_name %in% "Nt+/+; Ldlr-/-" ~ "ntWT_ldKO",
      group_name %in% "Nt-/-; Ldlr-/-" ~ "ntKO_ldKO"
    ))
  frame_classed
}
