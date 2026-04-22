#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

derived_dir <- file.path("data", "derived")
dir.create(derived_dir, recursive = TRUE, showWarnings = FALSE)

collapse_unique <- function(x) {
  x <- unique(as.character(x[!is.na(x) & x != ""]))
  if (length(x) == 0) {
    return(NA_character_)
  }
  paste(x, collapse = "|")
}

full_extension <- read_csv(
  file.path(derived_dir, "extension_au_non_au_matched_paired_labelled.csv"),
  show_col_types = FALSE
)

full_block_folds <- read_csv(
  file.path(derived_dir, "extension_au_non_au_site_block_folds.csv"),
  show_col_types = FALSE
)

full_loso_folds <- read_csv(
  file.path(derived_dir, "extension_au_non_au_site_loso_folds.csv"),
  show_col_types = FALSE
)

south_island_site_lookup <- tribble(
  ~sampling_site_label, ~sampling_location_specific, ~is_south_island_breeder, ~classification_basis, ~classification_note,
  "Cass River", "Cass River, Canterbury", TRUE, "site_label_and_specific_location", "South Island site retained in South-Island-only branch.",
  "Kaikoura", "Kaikoura, Canterbury", TRUE, "site_label_and_specific_location", "South Island site retained in South-Island-only branch.",
  "Kaitorete", "Kaitorete Spit, Canterbury", TRUE, "site_label_and_specific_location", "South Island site retained in South-Island-only branch.",
  "Motueka Spit", "Motueka Spit, Tasman Bay", TRUE, "site_label_and_specific_location", "South Island site retained in South-Island-only branch.",
  "Napier", "Napier, Hawkes Bay", FALSE, "site_label_and_specific_location", "North Island breeding site excluded from the South-Island-only branch.",
  "Pisa Range", "Pisa Range, Otago", TRUE, "site_label_and_specific_location", "South Island site retained in South-Island-only branch.",
  "Porters Beach", "Porters Beach, Abel Tasman National Park", TRUE, "site_label_and_specific_location", "South Island site retained in South-Island-only branch.",
  "Rangipo Desert", "Rangipo Desert, Tongariro National Park", FALSE, "site_label_and_specific_location", "North Island breeding site excluded from the South-Island-only branch.",
  "Tasman River", "Tasman River, Mt. Cook National Park", TRUE, "site_label_and_specific_location", "South Island site retained in South-Island-only branch.",
  "Te Anau", "Te Anau, Southland", TRUE, "site_label_and_specific_location", "South Island site retained in South-Island-only branch.",
  "Tiwai Point", "Tiwai Point, Canterbury", TRUE, "site_label_and_specific_location", "South Island site retained in South-Island-only branch."
)

write_csv(
  south_island_site_lookup,
  file.path(derived_dir, "extension_au_non_au_south_island_site_eligibility.csv")
)

south_island_extension <- full_extension %>%
  left_join(
    south_island_site_lookup,
    by = c("sampling_site_label", "sampling_location_specific")
  ) %>%
  mutate(
    south_island_rule = "Retain rows whose sampling_site_label and sampling_location_specific map to a South Island breeding site in the explicit site eligibility table; do not use sampling_location_general because it is non-discriminative in this dataset.",
    comparison_subset_id = "au_non_au_south_island_matched_paired_complete_cnh",
    extension_target = "AU_vs_non_AU_SouthIslandOnly"
  ) %>%
  filter(is_south_island_breeder %in% TRUE) %>%
  arrange(ring)

write_csv(
  south_island_extension,
  file.path(derived_dir, "extension_au_non_au_south_island_matched_paired_labelled.csv")
)

south_island_breast_features <- south_island_extension %>%
  transmute(
    comparison_subset_id,
    extension_target,
    south_island_rule,
    ring,
    sampling_site_label,
    sampling_location_specific,
    sampling_location_general,
    longitude_capture,
    latitude_capture,
    status,
    status_bin,
    status_rung2_3class,
    status_au_non_au,
    sex,
    tag_id_adjusted,
    migratory_status_comment,
    is_south_island_breeder,
    classification_basis,
    normalised_d13c_Breast,
    normalised_d15n_Breast,
    normalised_d2h_Breast,
    source_phenotype_row_ids,
    source_measurement_ids_Breast,
    source_measurement_ids_model_Breast
  )

write_csv(
  south_island_breast_features,
  file.path(derived_dir, "extension_au_non_au_south_island_breast_features.csv")
)

south_island_paired_contrast_features <- south_island_extension %>%
  transmute(
    comparison_subset_id,
    extension_target,
    south_island_rule,
    ring,
    sampling_site_label,
    sampling_location_specific,
    sampling_location_general,
    longitude_capture,
    latitude_capture,
    status,
    status_bin,
    status_rung2_3class,
    status_au_non_au,
    sex,
    tag_id_adjusted,
    migratory_status_comment,
    is_south_island_breeder,
    classification_basis,
    normalised_d13c_Breast,
    normalised_d13c_Primary,
    normalised_d15n_Breast,
    normalised_d15n_Primary,
    normalised_d2h_Breast,
    normalised_d2h_Primary,
    delta_d13c_breast_minus_primary,
    delta_d15n_breast_minus_primary,
    delta_d2h_breast_minus_primary,
    abs_delta_d13c,
    abs_delta_d15n,
    abs_delta_d2h,
    source_phenotype_row_ids,
    source_measurement_ids_Breast,
    source_measurement_ids_Primary,
    source_measurement_ids_model_Breast,
    source_measurement_ids_model_Primary
  )

write_csv(
  south_island_paired_contrast_features,
  file.path(derived_dir, "extension_au_non_au_south_island_paired_contrast_features.csv")
)

south_island_candidate_definitions <- tibble(
  extension_target = "AU_vs_non_AU_SouthIslandOnly",
  candidate_id = c(
    "au_non_au_south_island_breast_cnh",
    "au_non_au_south_island_paired_contrast_cnh"
  ),
  candidate_label = c(
    "South-Island-only AU vs non_AU: Breast-only C/N/H",
    "South-Island-only AU vs non_AU: Paired-contrast C/N/H"
  ),
  predictor_columns = c(
    "normalised_d13c_Breast|normalised_d15n_Breast|normalised_d2h_Breast",
    "abs_delta_d13c|abs_delta_d15n|abs_delta_d2h"
  ),
  source_table = c(
    "extension_au_non_au_south_island_breast_features.csv",
    "extension_au_non_au_south_island_paired_contrast_features.csv"
  ),
  matched_subset_required = TRUE,
  benchmark_table = "extension_au_non_au_south_island_matched_paired_labelled.csv",
  note = c(
    "Binary AU-vs-non_AU breast-only comparison on the frozen South-Island-only paired complete-C/N/H subset.",
    "Binary AU-vs-non_AU paired-contrast comparison on the same frozen South-Island-only paired complete-C/N/H subset. Paired-distance transforms, if used later, should be recomputed within each training split."
  )
)

write_csv(
  south_island_candidate_definitions,
  file.path(derived_dir, "extension_au_non_au_south_island_candidate_definitions.csv")
)

south_island_class_counts <- bind_rows(
  south_island_extension %>%
    count(status_au_non_au, name = "n") %>%
    transmute(
      summary_type = "south_island_binary_target",
      class = status_au_non_au,
      n
    ),
  south_island_extension %>%
    count(status_rung2_3class, name = "n") %>%
    transmute(
      summary_type = "south_island_3class_benchmark",
      class = status_rung2_3class,
      n
    ),
  full_extension %>%
    count(status_au_non_au, name = "n") %>%
    transmute(
      summary_type = "full_extension_binary_reference",
      class = status_au_non_au,
      n
    )
)

write_csv(
  south_island_class_counts,
  file.path(derived_dir, "extension_au_non_au_south_island_class_counts.csv")
)

south_island_site_counts <- south_island_extension %>%
  group_by(
    sampling_site_label,
    sampling_location_specific
  ) %>%
  summarise(
    n_total = n(),
    n_AU = sum(status_au_non_au == "AU", na.rm = TRUE),
    n_non_AU = sum(status_au_non_au == "non_AU", na.rm = TRUE),
    n_resident = sum(status_rung2_3class == "resident", na.rm = TRUE),
    n_NZ_migrant = sum(status_rung2_3class == "NZ migrant", na.rm = TRUE),
    n_AU_migrant = sum(status_rung2_3class == "AU migrant", na.rm = TRUE),
    source_rings = collapse_unique(ring),
    .groups = "drop"
  ) %>%
  arrange(desc(n_total), sampling_site_label)

write_csv(
  south_island_site_counts,
  file.path(derived_dir, "extension_au_non_au_south_island_site_class_counts.csv")
)

south_island_site_block_folds <- full_block_folds %>%
  semi_join(
    south_island_extension %>%
      select(ring),
    by = "ring"
  ) %>%
  mutate(
    comparison_subset_id = "au_non_au_south_island_matched_paired_complete_cnh",
    extension_target = "AU_vs_non_AU_SouthIslandOnly"
  )

write_csv(
  south_island_site_block_folds,
  file.path(derived_dir, "extension_au_non_au_south_island_site_block_folds.csv")
)

south_island_site_loso_folds <- full_loso_folds %>%
  semi_join(
    south_island_extension %>%
      select(ring),
    by = "ring"
  ) %>%
  mutate(
    comparison_subset_id = "au_non_au_south_island_matched_paired_complete_cnh",
    extension_target = "AU_vs_non_AU_SouthIslandOnly"
  )

write_csv(
  south_island_site_loso_folds,
  file.path(derived_dir, "extension_au_non_au_south_island_site_loso_folds.csv")
)

fold_summary <- bind_rows(
  south_island_site_block_folds %>%
    group_by(validation_layer, split_id) %>%
    summarise(
      n_rows = n(),
      n_sites = n_distinct(sampling_site_label),
      n_AU = sum(status_au_non_au == "AU", na.rm = TRUE),
      n_non_AU = sum(status_au_non_au == "non_AU", na.rm = TRUE),
      sites = collapse_unique(sampling_site_label),
      .groups = "drop"
    ),
  south_island_site_loso_folds %>%
    group_by(validation_layer, split_id) %>%
    summarise(
      n_rows = n(),
      n_sites = n_distinct(sampling_site_label),
      n_AU = sum(status_au_non_au == "AU", na.rm = TRUE),
      n_non_AU = sum(status_au_non_au == "non_AU", na.rm = TRUE),
      sites = collapse_unique(sampling_site_label),
      .groups = "drop"
    )
)

write_csv(
  fold_summary,
  file.path(derived_dir, "extension_au_non_au_south_island_fold_summary.csv")
)

subset_ring_key <- sort(unique(south_island_extension$ring))

subset_check <- tibble(
  table_name = c(
    "extension_au_non_au_south_island_matched_paired_labelled.csv",
    "extension_au_non_au_south_island_breast_features.csv",
    "extension_au_non_au_south_island_paired_contrast_features.csv",
    "extension_au_non_au_south_island_site_block_folds.csv",
    "extension_au_non_au_south_island_site_loso_folds.csv"
  ),
  n_rows = c(
    nrow(south_island_extension),
    nrow(south_island_breast_features),
    nrow(south_island_paired_contrast_features),
    nrow(south_island_site_block_folds),
    nrow(south_island_site_loso_folds)
  ),
  n_unique_rings = c(
    n_distinct(south_island_extension$ring),
    n_distinct(south_island_breast_features$ring),
    n_distinct(south_island_paired_contrast_features$ring),
    n_distinct(south_island_site_block_folds$ring),
    n_distinct(south_island_site_loso_folds$ring)
  ),
  rings_match_frozen_subset = c(
    identical(sort(unique(south_island_extension$ring)), subset_ring_key),
    identical(sort(unique(south_island_breast_features$ring)), subset_ring_key),
    identical(sort(unique(south_island_paired_contrast_features$ring)), subset_ring_key),
    identical(sort(unique(south_island_site_block_folds$ring)), subset_ring_key),
    identical(sort(unique(south_island_site_loso_folds$ring)), subset_ring_key)
  ),
  note = c(
    "Frozen South-Island-only binary extension benchmark table.",
    "Breast-only South-Island-only AU-vs-non_AU modelling-ready feature table.",
    "Paired-contrast South-Island-only AU-vs-non_AU modelling-ready feature table.",
    "Reuses accepted grouped site-block definitions after South Island filtering.",
    "Reuses accepted LOSO site definitions after South Island filtering."
  )
)

write_csv(
  subset_check,
  file.path(derived_dir, "extension_au_non_au_south_island_subset_check.csv")
)

summarise_site_imbalance <- function(site_counts, dataset_scope, removed_sites) {
  total_au <- sum(site_counts$n_AU)
  total_non_au <- sum(site_counts$n_non_AU)
  max_index <- which.max(site_counts$n_AU)

  tibble(
    dataset_scope = dataset_scope,
    n_rows = sum(site_counts$n_total),
    n_sites = nrow(site_counts),
    n_AU = total_au,
    n_non_AU = total_non_au,
    prop_AU = total_au / sum(site_counts$n_total),
    n_mixed_sites = sum(site_counts$n_AU > 0 & site_counts$n_non_AU > 0),
    n_AU_only_sites = sum(site_counts$n_AU > 0 & site_counts$n_non_AU == 0),
    n_non_AU_only_sites = sum(site_counts$n_non_AU > 0 & site_counts$n_AU == 0),
    largest_AU_site = site_counts$sampling_site_label[max_index],
    largest_AU_site_n_AU = site_counts$n_AU[max_index],
    largest_AU_site_share_of_AU = site_counts$n_AU[max_index] / total_au,
    removed_sites = removed_sites
  )
}

full_site_imbalance <- full_extension %>%
  group_by(sampling_site_label) %>%
  summarise(
    n_total = n(),
    n_AU = sum(status_au_non_au == "AU", na.rm = TRUE),
    n_non_AU = sum(status_au_non_au == "non_AU", na.rm = TRUE),
    .groups = "drop"
  )

south_island_site_imbalance <- south_island_extension %>%
  group_by(sampling_site_label) %>%
  summarise(
    n_total = n(),
    n_AU = sum(status_au_non_au == "AU", na.rm = TRUE),
    n_non_AU = sum(status_au_non_au == "non_AU", na.rm = TRUE),
    .groups = "drop"
  )

imbalance_comparison <- bind_rows(
  summarise_site_imbalance(full_site_imbalance, "full_extension", NA_character_),
  summarise_site_imbalance(south_island_site_imbalance, "south_island_only_extension", "Napier|Rangipo Desert")
)

write_csv(
  imbalance_comparison,
  file.path(derived_dir, "extension_au_non_au_south_island_imbalance_comparison.csv")
)
