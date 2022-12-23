print_table = function(comparison_table){
  comparison_table %>%
    dplyr::select(class, p.value, adj.P.Val, n_pos, n_neg, direction) %>%
    dplyr::arrange(adj.P.Val) %>%
    dplyr::mutate(p.value = format(p.value, digits = 2, scientific = TRUE, justify = "right"),
                  adj.P.Val = format(adj.P.Val, digits = 2, scientific = TRUE, justify = "right")) %>%
    knitr::kable(digits = 2) %>%
    print()
}
