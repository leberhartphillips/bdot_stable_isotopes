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

base_live <- read_csv(
  file.path(derived_dir, "live_screening_ready_paired_labelled.csv"),
  show_col_types = FALSE
)

site_block_folds <- read_csv(
  file.path(derived_dir, "model_stage2_site_block_folds.csv"),
  show_col_types = FALSE
)

site_loso_folds <- read_csv(
  file.path(derived_dir, "model_stage2_site_loso_folds.csv"),
  show_col_types = FALSE
)

extension_base <- base_live %>%
  mutate(
    comparison_subset_id = "au_non_au_matched_paired_complete_cnh",
    extension_target = "AU_vs_non_AU",
    status_rung2_3class = case_when(
      status == "resident" ~ "resident",
      status == "AU migrant" ~ "AU migrant",
      status %in% c("SI migrant", "NI migrant") ~ "NZ migrant",
      TRUE ~ NA_character_
    ),
    status_au_non_au = case_when(
      status == "AU migrant" ~ "AU",
      status %in% c("resident", "SI migrant", "NI migrant") ~ "non_AU",
      TRUE ~ NA_character_
    ),
    delta_d13c_breast_minus_primary = normalised_d13c_Breast - normalised_d13c_Primary,
    delta_d15n_breast_minus_primary = normalised_d15n_Breast - normalised_d15n_Primary,
    delta_d2h_breast_minus_primary = normalised_d2h_Breast - normalised_d2h_Primary,
    abs_delta_d13c = abs(delta_d13c_breast_minus_primary),
    abs_delta_d15n = abs(delta_d15n_breast_minus_primary),
    abs_delta_d2h = abs(delta_d2h_breast_minus_primary)
  ) %>%
  filter(
    has_complete_cnh_pair %in% TRUE,
    !is.na(status_rung2_3class),
    !is.na(status_au_non_au)
  ) %>%
  arrange(ring)

write_csv(
  extension_base,
  file.path(derived_dir, "extension_au_non_au_matched_paired_labelled.csv")
)

extension_breast_features <- extension_base %>%
  transmute(
    comparison_subset_id,
    extension_target,
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
    normalised_d13c_Breast,
    normalised_d15n_Breast,
    normalised_d2h_Breast,
    source_phenotype_row_ids,
    source_measurement_ids_Breast,
    source_measurement_ids_model_Breast
  )

write_csv(
  extension_breast_features,
  file.path(derived_dir, "extension_au_non_au_breast_features.csv")
)

extension_paired_contrast_features <- extension_base %>%
  transmute(
    comparison_subset_id,
    extension_target,
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
  extension_paired_contrast_features,
  file.path(derived_dir, "extension_au_non_au_paired_contrast_features.csv")
)

candidate_definitions <- tibble(
  extension_target = "AU_vs_non_AU",
  candidate_id = c(
    "au_non_au_breast_cnh",
    "au_non_au_paired_contrast_cnh"
  ),
  candidate_label = c(
    "AU vs non_AU: Breast-only C/N/H",
    "AU vs non_AU: Paired-contrast C/N/H"
  ),
  predictor_columns = c(
    "normalised_d13c_Breast|normalised_d15n_Breast|normalised_d2h_Breast",
    "abs_delta_d13c|abs_delta_d15n|abs_delta_d2h"
  ),
  source_table = c(
    "extension_au_non_au_breast_features.csv",
    "extension_au_non_au_paired_contrast_features.csv"
  ),
  matched_subset_required = TRUE,
  benchmark_table = "extension_au_non_au_matched_paired_labelled.csv",
  note = c(
    "Binary AU-vs-non_AU breast-only comparison on the frozen paired complete-C/N/H subset.",
    "Binary AU-vs-non_AU paired-contrast comparison on the same frozen paired complete-C/N/H subset. Paired-distance transforms, if used later, should be recomputed within each training split."
  )
)

write_csv(
  candidate_definitions,
  file.path(derived_dir, "extension_au_non_au_candidate_definitions.csv")
)

overall_counts <- bind_rows(
  extension_base %>%
    count(status_au_non_au, name = "n") %>%
    transmute(
      summary_type = "overall_binary_target",
      class = status_au_non_au,
      n
    ),
  extension_base %>%
    count(status_rung2_3class, name = "n") %>%
    transmute(
      summary_type = "matched_3class_benchmark",
      class = status_rung2_3class,
      n
    )
)

write_csv(
  overall_counts,
  file.path(derived_dir, "extension_au_non_au_class_counts.csv")
)

site_counts <- extension_base %>%
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
  site_counts,
  file.path(derived_dir, "extension_au_non_au_site_class_counts.csv")
)

extension_site_block_folds <- site_block_folds %>%
  inner_join(
    extension_base %>%
      select(
        ring,
        comparison_subset_id,
        extension_target,
        sampling_site_label,
        sampling_location_specific,
        status,
        status_rung2_3class,
        status_au_non_au
      ),
    by = c("ring", "sampling_site_label", "sampling_location_specific")
  )

write_csv(
  extension_site_block_folds,
  file.path(derived_dir, "extension_au_non_au_site_block_folds.csv")
)

extension_site_loso_folds <- site_loso_folds %>%
  inner_join(
    extension_base %>%
      select(
        ring,
        comparison_subset_id,
        extension_target,
        sampling_site_label,
        sampling_location_specific,
        status,
        status_rung2_3class,
        status_au_non_au
      ),
    by = c("ring", "sampling_site_label", "sampling_location_specific")
  )

write_csv(
  extension_site_loso_folds,
  file.path(derived_dir, "extension_au_non_au_site_loso_folds.csv")
)

subset_ring_key <- sort(unique(extension_base$ring))

subset_check <- tibble(
  table_name = c(
    "extension_au_non_au_matched_paired_labelled.csv",
    "extension_au_non_au_breast_features.csv",
    "extension_au_non_au_paired_contrast_features.csv",
    "extension_au_non_au_site_block_folds.csv",
    "extension_au_non_au_site_loso_folds.csv"
  ),
  n_rows = c(
    nrow(extension_base),
    nrow(extension_breast_features),
    nrow(extension_paired_contrast_features),
    nrow(extension_site_block_folds),
    nrow(extension_site_loso_folds)
  ),
  n_unique_rings = c(
    dplyr::n_distinct(extension_base$ring),
    dplyr::n_distinct(extension_breast_features$ring),
    dplyr::n_distinct(extension_paired_contrast_features$ring),
    dplyr::n_distinct(extension_site_block_folds$ring),
    dplyr::n_distinct(extension_site_loso_folds$ring)
  ),
  rings_match_frozen_subset = c(
    identical(sort(unique(extension_base$ring)), subset_ring_key),
    identical(sort(unique(extension_breast_features$ring)), subset_ring_key),
    identical(sort(unique(extension_paired_contrast_features$ring)), subset_ring_key),
    identical(sort(unique(extension_site_block_folds$ring)), subset_ring_key),
    identical(sort(unique(extension_site_loso_folds$ring)), subset_ring_key)
  ),
  note = c(
    "Frozen binary extension benchmark table.",
    "Breast-only AU-vs-non_AU modelling-ready feature table.",
    "Paired-contrast AU-vs-non_AU modelling-ready feature table.",
    "Reuses accepted grouped site-block definitions with the binary target attached.",
    "Reuses accepted LOSO site definitions with the binary target attached."
  )
)

write_csv(
  subset_check,
  file.path(derived_dir, "extension_au_non_au_subset_check.csv")
)
