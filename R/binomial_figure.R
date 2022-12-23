create_binomial_plot = function(binomial_results){
  binomial_results = binomial_results %>%
    dplyr::arrange(padjust) %>%
    dplyr::mutate(label = format(padjust, digits = 2))
  long_binomial = binomial_results %>%
    dplyr::select(annotation, num_positive, num_negative, padjust) %>%
    dplyr::mutate(positive = num_positive,
                  negative = num_negative,
                  num_positive = NULL,
                  num_negative = NULL) %>%
    tidyr::pivot_longer(cols = !c(annotation, padjust), values_to = "count",
                        names_to = "direction") %>%
    dplyr::mutate(count = dplyr::case_when(
      direction == "positive" ~ 1 * count,
      direction == "negative" ~ -1 * count
    ))
  long_binomial$annotation = factor(long_binomial$annotation, levels = binomial_results$annotation, ordered = TRUE)
  long_binomial$direction = factor(long_binomial$direction,
                                   levels = c("positive", "negative"),
                                   ordered = TRUE)
  
  binomial_sig = sum(binomial_results$padjust <= 0.05) + 0.5
  
  label_df = binomial_results %>%
    dplyr::select(annotation, label, num_positive) %>%
    mutate(count = num_positive + 10,
           direction = "positive")
  
  out_plot = long_binomial %>%
    ggplot(aes(x = annotation, y = count, fill = direction)) +
    scale_fill_discrete() +
    geom_bar(stat = "identity") +
    geom_hline(color = "black", yintercept = 0) +
    geom_text(data = label_df, aes(x = annotation, y = count, label = label), size = 3) +
    geom_vline(color = "blue", xintercept = binomial_sig) +
    theme(axis.text.x = element_text(angle = 90),
          legend.position = "none")
  out_plot
}

breakdown_lipid = function(lipid_str){
  lipid_str2 = gsub("CER|CE|DAG|DCER|HCER|LCER|LPC|LPE|MAG|PC|PE|PI|SM|TAG|FA|O\\-|P\\-|\\(|\\)", "", lipid_str)
  lipid_split = strsplit(lipid_str2, ":|\\/|\\-")
  carbon_df = purrr::map_dfr(lipid_split, function(.x){
    if (length(.x) > 3) {
      out_frame = data.frame(C_total = .x[1], Dbl_total = .x[2],
                             C_chain = .x[3], Dbl_chain = .x[4])
    } else {
      out_frame = data.frame(C_total = .x[1], Dbl_total = .x[2],
                             C_chain = "0", Dbl_chain = "0")
    }
    out_frame
  })
  carbon_df = purrr::map_dfc(carbon_df, ~ as.numeric(.x))
  carbon_df$name = lipid_str
  carbon_long = carbon_df %>%
    tidyr::pivot_longer(cols = !name, values_to = "value",
                        names_to = "what")
  carbon_long
}
