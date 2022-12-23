calculate_icikt_species = function(lipotype_pnt, lipotype_all_species)
{
  # tar_load(lipotype_pnt)
  # tar_load(lipotype_all_species)
  all_data = dplyr::bind_rows(lipotype_pnt |>
                                dplyr::select(USUBJID, PNTS) |>
                                dplyr::mutate(feature2 = "PNTS",
                                              amount = PNTS) |>
                              dplyr::select(USUBJID, feature2, amount),
                              lipotype_all_species |>
                                dplyr::select(USUBJID, feature2, amount))
  wide_data = all_data |>
    tidyr::pivot_wider(id_cols = USUBJID, names_from = feature2, values_from = amount) |>
    dplyr::select(-USUBJID) |>
    as.matrix()
  
  ici_cor = ici_kendalltau(t(wide_data), perspective = "global", scale_max = TRUE, diag_good = FALSE, include_only = "PNTS")
  ici_res = tibble(correlation = ici_cor$cor[1, -1], p_value = ici_cor$pval[1, -1]) |>
    dplyr::mutate(padjust = p.adjust(p_value, method = "BH"),
                  feature2 = colnames(ici_cor$cor)[-1])
  ici_res
}

check_fa = function(lipotype_fa)
{
  fa_each = lipotype_fa |>
    dplyr::group_by(USUBJID, feature) |>
    dplyr::summarise(n_fa = dplyr::n(), n_unique = length(unique(amount)))
  split_feature = strsplit(fa_each$feature, " ", fixed = TRUE) %>% 
    purrr::map_chr(., ~ .x[1])
  fa_each$lipid_class = split_feature
  fa_each |>
    dplyr::filter(n_fa > 1, n_unique > 1) %>%
    dplyr::pull(lipid_class) %>% 
    unique()
}

extract_tagfa = function(lipotype_fa)
{
  tag_features = lipotype_fa |>
    dplyr::filter(grepl("^TAG", feature))
  tag_features = tag_features |>
    dplyr::mutate(feature2 = paste0(feature, "_", FA)) |>
    dplyr::select(USUBJID, feature, amount, feature2)
  
  tag_features
}

annotate_tag = function(lipotype_all_species)
{
  # tar_load(lipotype_all_species)
  just_species = lipotype_all_species |>
    dplyr::select(feature2) |>
    dplyr::distinct() |>
    dplyr::filter(grepl("^TAG", feature2)) |>
    dplyr::mutate(tmp1 = gsub(";[0-9]+", "", feature2),
                  lipid_id = gsub("_", "/", tmp1))
  species_annote = lipidr::annotate_lipids(just_species$lipid_id)
  species_annote = species_annote |>
    dplyr::mutate(total_cl = l_1,
                  total_cs = s_1,
                  chain1 = chain2,
                  l_1 = l_2,
                  s_1 = s_2,
                  chain2 = "",
                  l_2 = NA,
                  s_2 = NA)
  
  lt_annote = convert_lipidr_lipotype(species_annote)
  out_annote = dplyr::bind_cols(just_species[, "feature2"], lt_annote)
  tag_f1_f2 = lipotype_all_species |>
    dplyr::filter(grepl("^TAG", feature)) |>
    dplyr::select(feature, feature2) |>
    dplyr::distinct()
  out_annote = dplyr::left_join(out_annote, tag_f1_f2, by = "feature2")
  out_annote
}

convert_lipidr_lipotype = function(tag_lipidr)
{
  # tar_load(tag_lipidr)
  new_names = list(class = "Class",
                   totallength = "total_cl",
                   totaldb = "total_cs",
                   totaloh = as.numeric(NA),
                   FA1kind = "FA",	
                   FA1length = "l_1",
                   FA1db = "s_1",
                   FA1oh = as.numeric(NA),
                   FA2kind = as.character(NA),
                   FA2length = as.numeric(NA),
                   FA2db = as.numeric(NA),
                   FA2oh = as.numeric(NA),
                   shorthand_notation = as.character(NA))
  
  n_entry = nrow(tag_lipidr)
  
  renamed_data = purrr::map_dfc(new_names, function(lipid_match){
    #message(lipid_match)
    if (lipid_match %in% names(tag_lipidr)) {
      out_var = tag_lipidr[[lipid_match]]
    } else {
      out_var = rep(lipid_match, n_entry)
    }
    out_var
  })
  renamed_data

}

combine_lt_annotations = function(lipotype_annotations, tag_lipotype_annotation)
{
  # tar_load(lipotype_annotations)
  # tar_load(tag_lipotype_annotation)
  tmp_lt = lipotype_annotations |>
    dplyr::filter(!grepl("^TAG", feature)) |>
    dplyr::mutate(feature2 = feature)
  tag_lt = lipotype_annotations |>
    dplyr::filter(grepl("^TAG", feature)) |>
    dplyr::select(feature, swissname, swissid, swissrank)
  tag2_lt = dplyr::left_join(tag_lt, tag_lipotype_annotation, by = "feature") |>
    dplyr::select(dplyr::all_of(names(tmp_lt)))
  
  complete_annotations = dplyr::bind_rows(tmp_lt,
                                          tag2_lt)
  
  complete_annotations
}

create_lt_feature_annotations = function(lt_all_annotation, remove_regex = c("NA$|oh:0$"))
{
  lt_all_annotation = lt_all_annotation |>
    dplyr::mutate(class_totallength = paste0(class, "_", totallength),
                  class_totaldb = paste0(class, "_", totaldb),
                  class_fa1length = paste0(class, "_", FA1length),
                  class_fa2length = paste0(class, "_", FA2length))
  not_annotation = c("feature", "shorthand_notation", "swissname", "swissid", "swissrank", "feature2")
  annotation_cols = names(lt_all_annotation)[!(names(lt_all_annotation) %in% not_annotation)]
  
  
  feature_annotations = purrr::map(annotation_cols, function(use_annotation){
    out_annot = split(lt_all_annotation$feature2, lt_all_annotation[[use_annotation]])
    names(out_annot) = paste0(use_annotation, ":", names(out_annot))
    out_annot
  })
  feature_annotations = unlist(feature_annotations, recursive = FALSE)
  
  remove_annotations = grepl(remove_regex, names(feature_annotations))
  
  lt_annotation = categoryCompare2::annotation(feature_annotations[!remove_annotations], annotation_type = "lipotype", feature_type = "lipids")
  lt_annotation
}

run_binomial_enrichment = function(lt_sig_cor, median_cor, lt_feature_annotations)
{
  # tar_load(lt_sig_cor)
  # tar_load(lt_feature_annotations)
  
  up_feat = lt_sig_cor |>
    dplyr::filter(significant, correlation > median_cor) |>
    dplyr::pull(feature2)
  down_feat = lt_sig_cor |>
    dplyr::filter(significant, correlation < median_cor) |>
    dplyr::pull(feature2)
  
  binom_res = binomial_feature_enrichment(
    new("binomial_features",
        positivefc = up_feat,
        negativefc = down_feat,
        annotation = lt_feature_annotations),
    min_features = 6
  )
  binom_stats = as_tibble(binom_res@statistics@statistic_data)
  binom_stats$annotation = binom_res@statistics@annotation_id
  binom_stats = binom_stats |>
    dplyr::mutate(overall_cor = dplyr::case_when(
      num_positive > num_negative ~ "positive",
      num_positive < num_negative ~ "negative"
    ))
  binom_stats
}

format_p_values = function(in_df, p_columns = c("p", "padjust"), digits = 2)
{
  # in_df = pathway_up_down
  for (icol in p_columns) {
    in_df[[icol]] = purrr::map_chr(in_df[[icol]], format, digits = digits, scientific = -3)
  }
  in_df
}

calculate_median_sig_cor_heatmaps = function(lt_sig_cor, lt_all_annotation)
{
  # tar_load(lt_sig_cor)
  # tar_load(lt_all_annotation)
  
  lipid_class_order = tibble(class = c("CE", 
                                       "LPC", "PC", "PC O-", 
                                       "LPE", "PE", "PE O-",
                                       "MAG", "DAG", "TAG"))
  
  lt_sig_cor = lt_sig_cor |>
    dplyr::filter(significant)
  
  is_fa = lt_all_annotation |>
    dplyr::filter((FA1kind %in% "FA") | (FA2kind %in% "FA") | (class %in% lipid_class_order$class))
  total_db = is_fa |>
    dplyr::filter(!(class %in% "CE") & !(class %in% "TAG")) |>
    dplyr::select(feature2, class, totallength, totaldb) |>
    dplyr::mutate(len_db = paste0(totallength, ":", totaldb))
  sig_total = dplyr::inner_join(lt_sig_cor, total_db, by = "feature2")
  
  fa1_db = is_fa |>
    dplyr::filter(!(class %in% "CE")) |>
    dplyr::select(feature2, class, FA1length, FA1db) |>
    dplyr::mutate(len_db = paste0(FA1length, ":", FA1db)) |>
    dplyr::select(feature2, class, len_db)
  fa2_db = is_fa |>
    dplyr::select(feature2, class, FA2length, FA2db) |>
    dplyr::mutate(len_db = paste0(FA2length, ":", FA2db)) |>
    dplyr::select(feature2, class, len_db)
  fa_db = dplyr::bind_rows(fa1_db,
                           fa2_db) |>
    dplyr::filter(!grepl("NA", len_db))
  
  sig_fa = dplyr::inner_join(lt_sig_cor, fa_db, by = "feature2")
  
  tag_db = is_fa |>
    dplyr::filter(class %in% "TAG") |>
    dplyr::select(feature2, class, totallength, totaldb) |>
    dplyr::mutate(len_db = paste0(totallength, ":", totaldb))
  
  sig_tag = dplyr::inner_join(lt_sig_cor, tag_db, by = "feature2")
  
  fa_matrix = create_median_matrix(sig_fa, lipid_class_order)
  total_matrix = create_median_matrix(sig_total, lipid_class_order)
  tag_matrix = create_median_matrix(sig_tag, lipid_class_order)
  
  
  
  list(total = total_matrix,
       tag = tag_matrix,
       fa = fa_matrix)
}

create_median_matrix = function(sig_df, 
                                lipid_classes)
{
  # sig_df = sig_fa
  # lipid_classes = lipid_class_order
  
  median_df = sig_df |>
    dplyr::mutate(id = paste0(class, " ", len_db)) |>
    dplyr::group_by(id) |>
    dplyr::summarise(med_cor = median(correlation),
                     class = class[1],
                     len_db = len_db[1])
  
  in_classes = lipid_classes |>
    dplyr::filter(class %in% median_df$class) |>
    dplyr::pull(class)
  
  len_dbl = sort(unique(median_df$len_db))
  
  lipid_matrix = matrix(NA, nrow = length(in_classes),
                     ncol = length(len_dbl))
  rownames(lipid_matrix) = in_classes
  colnames(lipid_matrix) = len_dbl
  
  lipid_tibble = tidyr::expand_grid(class = rownames(lipid_matrix),
                                 len_db = colnames(lipid_matrix)) |>
    dplyr::mutate(id = paste0(class, " ", len_db))
  median_df_lipid = dplyr::inner_join(median_df[, c("id", "med_cor")], lipid_tibble, by = "id")
  
  for (irow in seq_len(nrow(median_df_lipid))) {
    lipid_matrix[median_df_lipid$class[irow], median_df_lipid$len_db[irow]] = median_df_lipid$med_cor[irow]
  }
  lipid_matrix
}

calc_lt_updown = function(lt_sig_cor, lt_all_annotation)
{
  tar_load(lt_sig_cor)
  tar_load(lt_all_annotation)
  
  lt_sig_cor = lt_sig_cor |>
    dplyr::filter(significant)
  
  is_fa = lt_all_annotation |>
    dplyr::filter((FA1kind %in% "FA") | (FA2kind %in% "FA"))
  total_db = is_fa |>
    dplyr::select(feature2, class, totallength, totaldb) |>
    dplyr::mutate(length = totallength,
                  db = totaldb,
                  source = "total",
                  totallength = NULL,
                  totaldb = NULL)
  fa1_db = is_fa |>
    dplyr::filter(!(class %in% "CE")) |>
    dplyr::select(feature2, class, FA1length, FA1db) |>
    dplyr::mutate(length = FA1length,
                  db = FA1db,
                  FA1length = NULL,
                  FA1db = NULL,
                  source = "fa")
  fa2_db = is_fa |>
    dplyr::select(feature2, class, FA2length, FA2db) |>
    dplyr::mutate(length = FA2length,
                  db = FA2db,
                  FA2length = NULL,
                  FA2db = NULL,
                  source = "fa")
  all_len_db = dplyr::bind_rows(
    total_db,
    fa1_db,
    fa2_db
  )
  split_len_db = split(all_len_db, all_len_db$class)
  
  up_down_cor_by_class = purrr::map(split_len_db, function(class_len_db){
    # class_len_db = split_len_db[[1]]
    class_len_db = dplyr::inner_join(class_len_db, lt_sig_cor[, c("feature2", "correlation")], by = "feature2")
    
    length_count = class_len_db |>
      dplyr::select(-db) |>
      dplyr::filter(!is.na(length)) |>
      dplyr::group_by(length, source) |>
      dplyr::summarise(positive = sum(correlation > 0),
                       negative = -1 * sum(correlation < 0)) |>
      tidyr::pivot_longer(cols = c(positive, negative), names_to = "direction",
                          values_to = "count")
    db_count = class_len_db |>
      dplyr::select(-length) |>
      dplyr::filter(!is.na(db)) |>
      dplyr::group_by(db, source) |>
      dplyr::summarise(positive = sum(correlation > 0),
                       negative = -1 * sum(correlation < 0)) |>
      tidyr::pivot_longer(cols = c(positive, negative), names_to = "direction",
                          values_to = "count")
    list(length = length_count,
         db = db_count)
    
  })
  up_down_cor_by_class
}

plot_lt_updown = function(lt_sig_updown)
{
  
  calc_breaks = function(length_values){
    range_lengths = range(length_values)
    if ((range_lengths[1] != 0) && ((range_lengths[1] %% 2) != 0)) {
      break_start = range_lengths[1] - 1
    } else {
      break_start = range_lengths[1]
    }
    if (diff(range_lengths) < 16) {
      out_breaks = seq(break_start, range_lengths[2], 2)
    } else {
      out_breaks = seq(break_start, range_lengths[2], 4)
    }
    out_breaks
  }
  # tar_load(lt_sig_updown)
  all_plots = purrr::imap(lt_sig_updown, function(in_data, in_id){
    message(in_id)
    len_total = in_data$length |>
      dplyr::filter(source %in% "total") |>
      dplyr::mutate(direction = factor(direction, levels = c("positive", "negative"), ordered = TRUE))
    len_max_total = max(abs(len_total$count))
    len_total_plot = len_total |>
      ggplot(aes(x = length, y = count, fill = direction)) +
      geom_col() +
      scale_fill_discrete() +
      scale_x_continuous(breaks = calc_breaks(len_total$length)) +
      geom_hline(yintercept = 0, color = "black") +
      theme(legend.position = "none") +
      coord_cartesian(ylim = c(-1 * len_max_total, len_max_total)) +
      labs(subtitle = "Total", x = "# of Carbons",
           y = "# of Positive & Negative Features")
    len_fa = in_data$length |>
      dplyr::filter(source %in% "fa") |>
      dplyr::mutate(direction = factor(direction, levels = c("positive", "negative"), ordered = TRUE))
    
    len_max_fa = max(abs(len_fa$count))
    len_fa_plot = len_fa |>
      ggplot(aes(x = length, y = count, fill = direction)) +
      geom_col() +
      scale_fill_discrete() +
      geom_hline(yintercept = 0, color = "black") +
      theme(legend.position = "none") +
      scale_x_continuous(breaks = calc_breaks(len_fa$length)) +
      coord_cartesian(ylim = c(-1 * len_max_fa, len_max_fa)) +
      labs(subtitle = "Chain", x = "# of Carbons",
           y = NULL)
    
    db_total = in_data$db |>
      dplyr::filter(source %in% "total") |>
      dplyr::mutate(direction = factor(direction, levels = c("positive", "negative"), ordered = TRUE))
    db_max_total = max(abs(db_total$count))
    db_total_plot = db_total |>
      ggplot(aes(x = db, y = count, fill = direction)) +
      geom_col() +
      scale_fill_discrete() +
      geom_hline(yintercept = 0, color = "black") +
      theme(legend.position = "none") +
      scale_x_continuous(breaks = calc_breaks(db_total$db)) +
      coord_cartesian(ylim = c(-1 * db_max_total, db_max_total)) +
      labs(subtitle = "Total", x = "# of Double Bonds",
           y = NULL)
    db_fa = in_data$db |>
      dplyr::filter(source %in% "fa") |>
      dplyr::mutate(direction = factor(direction, levels = c("positive", "negative"), ordered = TRUE))
    db_max_fa = max(abs(db_fa$count))
    db_fa_plot = db_fa |>
      ggplot(aes(x = db, y = count, fill = direction)) +
      geom_col() +
      scale_fill_discrete() +
      geom_hline(yintercept = 0, color = "black") +
      scale_x_continuous(breaks = calc_breaks(db_fa$db)) +
      theme(legend.position = "none") +
      coord_cartesian(ylim = c(-1 * db_max_fa, db_max_fa)) +
      labs(subtitle = "Chain", x = "# of Double Bonds", y = NULL)
    
  list(c_total = len_total_plot,
       db_total = db_total_plot,
       c_fa = len_fa_plot,
       db_fa = db_fa_plot)
  })
  all_plots
}

double_check_cor = function(lt_sig_cor, lipotype_pnt, lipotype_all_species)
{
  big_cor = lt_sig_cor |>
    dplyr::filter(significant) |>
    dplyr::arrange(dplyr::desc(abs(correlation))) |>
    dplyr::slice_head(n = 10)
  
  all_data = dplyr::bind_rows(lipotype_pnt |>
                                dplyr::select(USUBJID, PNTS) |>
                                dplyr::mutate(feature2 = "PNTS",
                                              amount = PNTS) |>
                                dplyr::select(USUBJID, feature2, amount),
                              lipotype_all_species |>
                                dplyr::select(USUBJID, feature2, amount))
  wide_data = all_data |>
    tidyr::pivot_wider(id_cols = USUBJID, names_from = feature2, values_from = amount) |>
    dplyr::select(-USUBJID)
  
  use_features = c("PNTS", big_cor$feature2)
  wide_data = wide_data |>
    dplyr::select(all_of(use_features))
  
  
}
