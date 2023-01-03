## Load your packages, e.g. library(targets).
source("./packages.R")

## Load your R files
lapply(list.files("./R", full.names = TRUE), source)

## tar_plan supports drake-style targets and also tar_target()
tar_plan(
  
  # loading metabolon data --------
  concentration_data = get_concentrations(),
  sample_info = get_sample_info(),
  feature_info = get_feature_info(concentration_data),
  feature_annotation_df = metabolonDataProcessing::mdp_annotate_lipids(feature_info, lipid_maps = metabolonDataProcessing::mdp_get_lipid_matrices()) |>
    dplyr::mutate(annotation2 = dplyr::case_when(
      grepl("PE\\(P", name) & (type %in% "class") ~ "PE P-",
      grepl("PE\\(O", name) & (type %in% "class") ~ "PE O-",
      TRUE ~ annotation
    )) |>
    dplyr::mutate(full_annotation = paste0(type, ":", annotation2),
                  annotation = annotation2,
                  annotation2 = NULL),
  feature_annotations = cc2_feature_annotations(feature_annotation_df),
  concentration_matrix = df_to_matrix(concentration_data, feature_info, sample_info),
  
  # qc-qa metabolon ------------
  ko_info = sample_info %>%
    dplyr::filter(age %in% "25 weeks"),
  ko_info_block = add_block_number(ko_info),
  ko_matrix = concentration_matrix[, ko_info$sample_id],
  ko_matrix_0 = replace_na(ko_matrix),
  ko_matrix_three = t(visualizationQualityControl::keep_non_zero_percentage(t(ko_matrix_0), sample_classes = ko_info$group_name, keep_num = 3)),
  
  tar_render(quality_control, "doc/quality_control.Rmd"),
  
  ko_outinfo = remove_outliers(ko_info, ko_matrix_three),
  ko_info_nooutlier = ko_outinfo %>%
    dplyr::filter(!outlier),
  ko_matrix_nooutlier = ko_matrix_three[, ko_info_nooutlier$sample_id],
  ko_matrix_3_nooutlier = t(visualizationQualityControl::keep_non_zero_percentage(t(ko_matrix_nooutlier), sample_classes = ko_info_nooutlier$group_name, keep_num = 0.33)),
  
  ko_normalized = normalize_data(ko_matrix_3_nooutlier),
  
  # differential analysis metabolon -------
  
  ko_diffs_mean = calculate_differences_mean(ko_normalized),
  pairwise_tests = list(x4_3 = c(4, 3)),
  differential_results = basic_differential_test(ko_normalized, ko_info_nooutlier, feature_info, pairwise_tests),
  
  pairwise_lfc = calculate_lfc(ko_normalized, ko_info_nooutlier, pairwise_tests, feature_info),
  
  ko_norm_df = convert_feature_matrix_to_df(ko_normalized, sample_info),
  
  #tar_render(basic_differential_results, "doc/basic_differential_results.Rmd"),
  basic_results_file = write_basic_stats_results(differential_results, ko_norm_df, here::here("doc/pairwise_comparisons_results.xlsx")),
  
  #tar_render(reduce_ldlr_variance, "doc/reduce_ldlr_variance.Rmd"),
  
  ko_anova_add = anova_add_differential_test(ko_normalized, ko_info_nooutlier, feature_info),
  #tar_render(differential_report, "doc/differential_results.Rmd"),
  
  pairwise_trim = purrr::map(pairwise_lfc, ~ .x$trimmed),
  feature_binom_diff = purrr::map(pairwise_trim["x4_3"], annotation_direction_test, feature_annotations),
  
  genotype_tests = list(ntkoldko_ntwtldko = c("ntKO_ldKO", "ntWT_ldKO")),
  genotype_results = basic_ttest(ko_normalized, ko_info_nooutlier, feature_info, genotype_tests),
  genotype_vs_pairwise = compare_test_results(differential_results, genotype_results),
  
  ko_sdams = sdams_differential_test(ko_normalized, ko_info_nooutlier, feature_info, pairwise_tests),
  
  tar_render(check_fold_changes, "doc/check_fold_changes.Rmd"),
  
  ko_dasev = dasev_differential_test(ko_normalized, ko_info_nooutlier, pairwise_tests),
  
  tar_render(pvalue_weirdness, "doc/pvalue_weirdness.Rmd"),
  
  ko_summary = calculate_summaries(ko_normalized, ko_info_nooutlier),
  
  tar_render(example_binomial, "doc/example_binomial.Rmd"),
  
  tar_render(split_classes, "doc/split_classes.Rmd"),
  
  tar_render(nt_in_ldlrko, "doc/nt_in_ldlrko.Rmd"),
  
  # lipotype data --------
  lipotype_pnt = readr::read_csv("data/lipotype_fattyacids/unfiltered/subject_data.csv") |>
    dplyr::filter(!is.na(PNTS)),
  lipotype_species = readr::read_csv("data/lipotype_fattyacids/unfiltered/species_nf.zip") |>
    dplyr::filter(USUBJID %in% lipotype_pnt$USUBJID),
  lt_species_filtered = filter_lt_species(lipotype_species, 0.7),
  lipotype_fa = readr::read_csv("data/lipotype_fattyacids/unfiltered/fa_data_nf.zip"),
  
  lipotype_tag = extract_tagfa(lipotype_fa),
  lipotype_all_species = dplyr::bind_rows(lipotype_species |>
                                            dplyr::filter(!grepl("^TAG", feature)) |>
                                            dplyr::mutate(feature2 = feature),
                                          lipotype_tag),
  
  lipotype_annotations = readr::read_csv("data/lipotype_fattyacids/unfiltered/lipid_data_nf.csv"),
  
  #tar_quarto(lipotype_evaluation, "doc/evaluate_lipotype.qmd"),
  
  lt_pnt_species_cor = calculate_icikt_species(lipotype_pnt, lipotype_all_species),
  # lt_median_cor = lt_pnt_species_cor |>
  #   dplyr::filter(padjust >= 0.5) |>
  #   dplyr::pull(correlation) |>
  #   median(),
  lt_median_cor = 0,
  
  padjust_limit = 0.05,
  lt_sig_cor = lt_pnt_species_cor |>
    dplyr::mutate(significant = padjust <= padjust_limit),
  
  tag_lipotype_annotation = annotate_tag(lipotype_all_species),
  lt_all_annotation = combine_lt_annotations(lipotype_annotations, tag_lipotype_annotation),
  
  lt_feature_annotations = create_lt_feature_annotations(lt_all_annotation),
  lt_correlation_enrichment = run_binomial_enrichment(lt_sig_cor, lt_median_cor, lt_feature_annotations),
  lt_correlation_sig = lt_correlation_enrichment |>
    dplyr::filter(padjust <= padjust_limit),
  
  lt_sig_heatmaps = calculate_median_sig_cor_heatmaps(lt_sig_cor, lt_all_annotation),
  lt_sig_updown = calc_lt_updown(lt_sig_cor, lt_all_annotation),
  lt_updown_plots = plot_lt_updown(lt_sig_updown),
  
  tar_target(eval_lipotype_file,
             "doc/evaluate_lipotype.qmd",
             format = "file"),
  
  evaluate_lipotype_hash = digest::digest(eval_lipotype_file, algo = "sha256", file = TRUE),
  
  tar_target(tag_xlsx,
             "data/TG_input_data.xlsx",
             format = "file"),
  
  tar_render(tag_results, "doc/tag_lipid_analysis.rmd"),
  
  tar_quarto(lipotype_results, "doc/lipotype_results.qmd"),
  
  tar_target(renv_lockfile,
             "renv.lock",
             format = "file"),
  
  tar_target(dependencies,
             {
               renv_lockfile
               devtools::session_info(pkgs = "loaded", to_file = "session-info.txt")
             }),
  
  tar_quarto(software_code_data_availability, "doc/software_code_data_availability.qmd"),
  
  tar_target(beep_end,
             {
               #differential_report
               example_binomial
               split_classes
               quality_control
               lipotype_results
               nt_in_ldlrko
               tag_results
               software_code_data_availability
               beepr::beep(4)
             })
  
# target = function_to_make(arg), ## drake style

# tar_target(target2, function_to_make2(arg)) ## targets style

)
