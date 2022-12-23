lyso_multiple_test = function(stats_list, classes = c("LPC", "LPE"), id = class, fc = "logFC"){
  NULL
}

annotation_direction_test = function(stats_df, annotations, min_features = 6){
  # stats_df = tar_read(pairwise_trim)$x4_3
  # annotations = tar_read(feature_annotations)
  stats_pos = stats_df |>
    dplyr::filter(logFC > 0) |>
    dplyr::pull(name)
  stats_neg = stats_df |>
    dplyr::filter(logFC < 0) |>
    dplyr::pull(name)
  
  binom_results = binomial_feature_enrichment(
    new("binomial_features",
        positivefc = stats_pos,
        negativefc = stats_neg,
        annotation = annotations),
    min_features = min_features
  )
  
  binom_stats = tibble::as_tibble(binom_results@statistics@statistic_data)
  binom_stats$annotation = names(binom_results@statistics@statistic_data$statistic)
  
  binom_stats = binom_stats |>
    dplyr::mutate(overall_fc = dplyr::case_when(
      num_positive > num_negative ~ "positive",
      num_negative > num_positive ~ "negative"
    ))
  
  
  binom_stats
}

chain_total_len_db = function(lipid_annotation){
  # lipid_annotation = tar_read(feature_annotations)
  annotation_list = lipid_annotation@annotation_features
  len_db_entries = grep("^tot_len_db|^chain_len_db", names(annotation_list), value = TRUE)
  split_entries = strsplit(len_db_entries, ":", fixed = TRUE)
  split_df = purrr::map_df(split_entries, function(.x){
    tibble(id_1 = .x[1],
           n_c = as.numeric(.x[2]),
           n_dbl = as.numeric(.x[3])) %>%
      dplyr::mutate(type = dplyr::case_when(
        grepl("^tot", id_1) ~ "total",
        grepl("^chain", id_1) ~ "chain"
      ))
  })
  split_df$full_len_db = len_db_entries
  
  class_entries = grep("^class:", names(annotation_list), value = TRUE)
  split_class = strsplit(class_entries, ":", fixed = TRUE)
  classes = tibble(class_id = class_entries, class = purrr::map_chr(split_class, ~ .x[2]))
  
  classes_df = purrr::imap_dfr(annotation_list[class_entries], function(.x, .y){
    tibble(name = .x, class_id = .y)
  })
  classes_df = dplyr::left_join(classes_df, classes, by = "class_id")
  
  annot_df = purrr::imap_dfr(annotation_list[len_db_entries], function(.x, .y){
    tibble(full_len_db = .y,
           name = .x)
  })
  annot_df = dplyr::left_join(annot_df, split_df, by = "full_len_db")
  annot_df = dplyr::left_join(annot_df, classes_df, by = "name")
  annot_df
}

up_down_len = function(lipid_len_sat2){
  
  summarized_lists = list("C_total" = lipid_len_sat2 %>%
                            dplyr::filter(type %in% "total") %>%
                            dplyr::group_by(n_c) %>%
                            dplyr::summarize(n_pos = sum(logFC > 0), n_neg = -1 * sum(logFC < 0)) %>%
                            tidyr::pivot_longer(cols = c(n_pos, n_neg), names_to = "direction", values_to = "count") %>%
                            dplyr::mutate(x = n_c,
                                          direction = factor(direction, levels = c("n_pos", "n_neg"), ordered = TRUE)),
                          "Dbl_total" = lipid_len_sat2 %>%
                            dplyr::filter(type %in% "total") %>%
                            dplyr::group_by(n_dbl) %>%
                            dplyr::summarize(n_pos = sum(logFC > 0), n_neg = -1 * sum(logFC < 0)) %>%
                            tidyr::pivot_longer(cols = c(n_pos, n_neg), names_to = "direction", values_to = "count") %>%
                            dplyr::mutate(x = n_dbl,
                                          direction = factor(direction, levels = c("n_pos", "n_neg"), ordered = TRUE)),
                          "C_chain" = lipid_len_sat2 %>%
                            dplyr::filter(type %in% "chain") %>%
                            dplyr::group_by(n_c) %>%
                            dplyr::summarize(n_pos = sum(logFC > 0), n_neg = -1 * sum(logFC < 0)) %>%
                            tidyr::pivot_longer(cols = c(n_pos, n_neg), names_to = "direction", values_to = "count") %>%
                            dplyr::mutate(x = n_c,
                                          direction = factor(direction, levels = c("n_pos", "n_neg"), ordered = TRUE)),
                          "Dbl_chain" = lipid_len_sat2 %>%
                            dplyr::filter(type %in% "chain") %>%
                            dplyr::group_by(n_dbl) %>%
                            dplyr::summarize(n_pos = sum(logFC > 0), n_neg = -1 * sum(logFC < 0)) %>%
                            tidyr::pivot_longer(cols = c(n_pos, n_neg), names_to = "direction", values_to = "count") %>%
                            dplyr::mutate(x = n_dbl,
                                          direction = factor(direction, levels = c("n_pos", "n_neg"), ordered = TRUE)))
  
  plot_lists = purrr::imap(summarized_lists, function(in_summary, plot_id){
    if (grepl("Dbl", plot_id)) {
      x_lab = "# of Double Bonds"
    } else {
      x_lab = "# of Carbons"
    }
    if (grepl("total", plot_id)) {
      sub_title = "Total"  
    } else {
      sub_title = "Chain"
    }
    if (plot_id %in% "C_total") {
      y_lab = "# of Positive & Negative Features"
    } else {
      y_lab = NULL
    }
    
    if (plot_id %in% "C_total") {
      use_breaks = seq(min(in_summary$x), max(in_summary$x), 4)
    } else if (plot_id %in% "Dbl_total") {
      use_breaks = seq(min(in_summary$x), max(in_summary$x), 4)
    } else if (plot_id %in% "C_chain") {
      use_breaks = seq(min(in_summary$x), max(in_summary$x), 2)
    } else if (plot_id %in% "Dbl_chain") {
      use_breaks = seq(min(in_summary$x), max(in_summary$x), 2)
    }
    
    max_val = max(abs(in_summary$count))
    ylim = c(-1 * max_val, max_val)
    tmp_plot = in_summary %>%
      ggplot(aes(x = x, y = count, fill = direction)) +
      scale_fill_discrete() +
      geom_bar(stat = "identity") +
      geom_hline(color = "black", yintercept = 0) +
      scale_x_continuous(breaks = use_breaks) +
      coord_cartesian(ylim = ylim) +
      labs(subtitle = sub_title, x = x_lab, y = y_lab) +
      theme(legend.position = "none")
  })
  
  plot_lists
}


binomial_test = function(fc_values){
  n_up = sum(fc_values > 0)
  n_dn = sum(fc_values < 0)
  b_res = binom.test(c(n_up, n_dn))
  b_tidy = broom::tidy(b_res)
  b_tidy = b_tidy %>%
    dplyr::mutate(direction = dplyr::case_when(
      n_up >= n_dn ~ "pos",
      n_up < n_dn ~ "neg"
    ))
  b_tidy$n_pos = n_up
  b_tidy$n_neg = n_dn
  b_tidy
}
