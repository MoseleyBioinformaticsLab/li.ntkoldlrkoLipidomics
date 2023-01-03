lipotype_pnt = readr::read_csv("data/lipotype_fattyacids/unfiltered/subject_data.csv")
set.seed(20230103)

# create sequential IDs, and then randomly shuffle them around
# so that original data order doesn't match them.
n_id = nrow(lipotype_pnt)
max_width = nchar(n_id)
out_format = paste0("%0", max_width, "d")
new_id = sprintf(out_format, seq(1000, 1000 + n_id))
new_id = sample(new_id, n_id)

lipotype_pnt$NEWID = new_id
lipotype_id_map = lipotype_pnt |>
  dplyr::select(USUBJID, NEWID)

lipotype_pnt_anon = lipotype_pnt |>
  dplyr::select(NEWID, PNTS) |>
  dplyr::transmute(USUBJID = NEWID,
                   PNTS = PNTS)

lipotype_species = readr::read_csv("data/lipotype_fattyacids/unfiltered/species_nf.zip")

lipotype_species_anon = dplyr::left_join(lipotype_species, lipotype_id_map, by = "USUBJID")
lipotype_species_anon = lipotype_species_anon |>
  dplyr::transmute(USUBJID = NEWID,
                   feature = feature,
                   amount = amount)

lipotype_fa = readr::read_csv("data/lipotype_fattyacids/unfiltered/fa_data_nf.zip")
lipotype_fa_anon = dplyr::left_join(lipotype_fa, lipotype_id_map, by = "USUBJID")
lipotype_fa_anon = lipotype_fa_anon |>
  dplyr::transmute(USUBJID = NEWID,
                   feature = feature,
                   class = class,
                   FA = FA,
                   FAKind = FAKind,
                   amount = amount)

write.table(lipotype_pnt_anon, file = "data/lipotype_fattyacids/anonymized/subject_data.csv",
            row.names = FALSE, col.names = TRUE, sep = ",",
            )
write.table(lipotype_species_anon, file = "data/lipotype_fattyacids/anonymized/species_nf.csv",
            row.names = FALSE, col.names = TRUE, sep = ",")
write.table(lipotype_fa_anon, file = "data/lipotype_fattyacids/anonymized/fa_data_nf.csv",
            row.names = FALSE, col.names = TRUE, sep = ",")

# update the column names while we are at it
lipid_annotation = readr::read_csv("data/lipotype_fattyacids/unfiltered/lipid_data_nf.csv")
lipid_annotation2 = lipid_annotation |>
  dplyr::mutate(swisslipids_name = swissname,
                swisslipids_id = swissid,
                swisslipids_rank = swissrank,
                swissname = NULL,
                swissid = NULL,
                swissrank = NULL)
write.table(lipid_annotation2, file = "data/lipotype_fattyacids/anonymized/lipid_data_nf.csv",
            row.names = FALSE, col.names = TRUE, sep = ",")
