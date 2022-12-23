filter_lt_species = function(lipotype_species, min_presence = 0.7)
{
  n_sub = length(unique(lipotype_species$USUBJID))
  lt_presence = lipotype_species |>
    dplyr::group_by(feature) %>%
    dplyr::summarise(n_present = dplyr::n(),
                     perc_present = n_present / n_sub)
  lt_keep = lt_presence |>
    dplyr::filter(perc_present >= min_presence)
  lt_out = lipotype_species |>
    dplyr::filter(feature %in% lt_keep$feature)
  lt_out
}
