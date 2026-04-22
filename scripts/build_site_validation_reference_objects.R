#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
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

safe_mean <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (all(is.na(x))) {
    return(NA_real_)
  }
  mean(x, na.rm = TRUE)
}

safe_var <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[!is.na(x)]
  if (length(x) <= 1) {
    return(NA_real_)
  }
  stats::var(x)
}

safe_sd <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[!is.na(x)]
  if (length(x) <= 1) {
    return(NA_real_)
  }
  stats::sd(x)
}

pairwise_n <- function(x, y) {
  sum(stats::complete.cases(x, y))
}

safe_cov <- function(x, y) {
  keep <- stats::complete.cases(x, y)
  if (sum(keep) <= 1) {
    return(NA_real_)
  }
  stats::cov(x[keep], y[keep])
}

safe_cor <- function(x, y) {
  keep <- stats::complete.cases(x, y)
  if (sum(keep) <= 1) {
    return(NA_real_)
  }
  stats::cor(x[keep], y[keep])
}

derive_stage2_features_local <- function(data) {
  data %>%
    mutate(
      delta_d13c_breast_minus_primary = normalised_d13c_Breast - normalised_d13c_Primary,
      delta_d15n_breast_minus_primary = normalised_d15n_Breast - normalised_d15n_Primary,
      delta_d2h_breast_minus_primary = normalised_d2h_Breast - normalised_d2h_Primary,
      abs_delta_d13c = abs(delta_d13c_breast_minus_primary),
      abs_delta_d15n = abs(delta_d15n_breast_minus_primary),
      abs_delta_d2h = abs(delta_d2h_breast_minus_primary),
      status_rung1 = factor(status_bin, levels = c("resident", "migrant")),
      status_rung2 = factor(
        case_when(
          status == "resident" ~ "resident",
          status == "AU migrant" ~ "AU migrant",
          status %in% c("SI migrant", "NI migrant") ~ "NZ migrant",
          TRUE ~ NA_character_
        ),
        levels = c("resident", "NZ migrant", "AU migrant")
      ),
      coord_group = paste(longitude_capture, latitude_capture, sep = ",")
    )
}

make_balanced_group_folds <- function(data, outcome_col, group_col, v) {
  outcome_sym <- rlang::sym(outcome_col)
  group_sym <- rlang::sym(group_col)

  group_counts <- data %>%
    count(!!group_sym, !!outcome_sym, name = "n") %>%
    pivot_wider(
      names_from = !!outcome_sym,
      values_from = n,
      values_fill = 0
    ) %>%
    rename(group_id = !!rlang::sym(group_col))

  class_cols <- setdiff(colnames(group_counts), "group_id")
  group_counts$total_n <- rowSums(group_counts[, class_cols, drop = FALSE])
  group_counts <- group_counts %>%
    arrange(desc(total_n))

  target_class <- colSums(group_counts[, class_cols, drop = FALSE]) / v
  target_total <- sum(group_counts$total_n) / v
  target_mat <- matrix(rep(target_class, each = v), nrow = v)

  fold_class <- matrix(0, nrow = v, ncol = length(class_cols))
  colnames(fold_class) <- class_cols
  fold_total <- rep(0, v)
  assigned_fold <- integer(nrow(group_counts))

  for (i in seq_len(nrow(group_counts))) {
    row_counts <- as.numeric(group_counts[i, class_cols, drop = FALSE])
    row_total <- group_counts$total_n[i]

    scores <- vapply(
      seq_len(v),
      function(fold_id) {
        new_fold_class <- fold_class
        new_fold_class[fold_id, ] <- new_fold_class[fold_id, ] + row_counts

        new_fold_total <- fold_total
        new_fold_total[fold_id] <- new_fold_total[fold_id] + row_total

        class_score <- sum(((new_fold_class - target_mat) / pmax(target_mat, 1))^2)
        total_score <- sum(((new_fold_total - target_total) / pmax(target_total, 1))^2)
        empty_penalty <- sum(new_fold_total == 0) * 100

        class_score + total_score + empty_penalty
      },
      numeric(1)
    )

    best_folds <- which(scores == min(scores))
    if (length(best_folds) > 1) {
      best_fold <- best_folds[which.min(fold_total[best_folds])]
    } else {
      best_fold <- best_folds
    }

    assigned_fold[i] <- best_fold
    fold_class[best_fold, ] <- fold_class[best_fold, ] + row_counts
    fold_total[best_fold] <- fold_total[best_fold] + row_total
  }

  tibble(
    group_id = group_counts$group_id,
    split_id = paste0("Fold", assigned_fold)
  )
}

build_covariance_long <- function(data, group_cols, feature_cols, reference_object, grouping_variable) {
  feature_pairs <- combn(feature_cols, 2, simplify = FALSE)

  bind_rows(
    lapply(
      feature_pairs,
      function(pair) {
        feature_x <- pair[[1]]
        feature_y <- pair[[2]]

        data %>%
          group_by(across(all_of(group_cols))) %>%
          summarise(
            reference_object = reference_object,
            grouping_variable = grouping_variable,
            feature_x = feature_x,
            feature_y = feature_y,
            n_complete_pairs = pairwise_n(.data[[feature_x]], .data[[feature_y]]),
            covariance = safe_cov(.data[[feature_x]], .data[[feature_y]]),
            correlation = safe_cor(.data[[feature_x]], .data[[feature_y]]),
            .groups = "drop"
          )
      }
    )
  )
}

feature_summary_cols <- function(prefix, feature_names, stats = c("mean", "var", "sd")) {
  cols <- character(0)
  for (feature in feature_names) {
    for (stat in stats) {
      cols <- c(cols, paste0(prefix, stat, "_", feature))
    }
  }
  cols
}

live_phenotypes <- read_csv(
  file.path(derived_dir, "live_phenotypes_canonical.csv"),
  show_col_types = FALSE
)

live_paired <- read_csv(
  file.path(derived_dir, "live_cn_paired_by_ring.csv"),
  show_col_types = FALSE
)

live_labelled <- read_csv(
  file.path(derived_dir, "live_screening_ready_paired_labelled.csv"),
  show_col_types = FALSE
) %>%
  derive_stage2_features_local()

live_tissue <- read_csv(
  file.path(derived_dir, "live_cn_tissue_summary.csv"),
  show_col_types = FALSE
)

museum_winter <- read_csv(
  file.path(derived_dir, "museum_screening_ready_breast_homogenate_winter.csv"),
  show_col_types = FALSE
)

site_lookup <- live_phenotypes %>%
  group_by(
    sampling_location_general,
    sampling_location_specific,
    sampling_site_label
  ) %>%
  summarise(
    n_live_phenotype_rings = n(),
    n_status_known = sum(status_known %in% TRUE, na.rm = TRUE),
    n_status_unknown = sum(!(status_known %in% TRUE), na.rm = TRUE),
    source_phenotype_row_ids = collapse_unique(source_phenotype_row_ids),
    sample_identifier_examples = collapse_unique(sample_identifier_examples),
    .groups = "drop"
  ) %>%
  arrange(sampling_site_label)

site_label_mapping_check <- site_lookup %>%
  count(sampling_site_label, name = "n_specific_locations_for_site_label")

specific_mapping_check <- site_lookup %>%
  count(sampling_location_specific, name = "n_site_labels_for_specific_location")

site_lookup <- site_lookup %>%
  left_join(site_label_mapping_check, by = "sampling_site_label") %>%
  left_join(specific_mapping_check, by = "sampling_location_specific") %>%
  mutate(
    grouping_ambiguity = n_specific_locations_for_site_label > 1 | n_site_labels_for_specific_location > 1,
    harmonization_note = case_when(
      sampling_location_specific == "Kaitorete Spit, Canterbury" ~
        "sampling_site_label harmonized from 'Kaitorete Spit' to 'Kaitorete'",
      TRUE ~ NA_character_
    )
  )

write_csv(
  site_lookup,
  file.path(derived_dir, "model_stage2_site_lookup.csv")
)

paired_site_counts <- live_paired %>%
  group_by(sampling_site_label) %>%
  summarise(
    n_paired_rings_total = n(),
    n_paired_rings_status_known = sum(status_known %in% TRUE, na.rm = TRUE),
    n_paired_rings_status_unknown = sum(!(status_known %in% TRUE), na.rm = TRUE),
    n_paired_complete_cnh = sum(has_complete_cnh_pair %in% TRUE, na.rm = TRUE),
    n_missing_site_fields = sum(
      is.na(sampling_site_label) |
        is.na(sampling_location_specific) |
        is.na(longitude_capture) |
        is.na(latitude_capture)
    ),
    source_phenotype_row_ids_paired = collapse_unique(source_phenotype_row_ids),
    .groups = "drop"
  )

rung1_site_summary <- live_labelled %>%
  group_by(sampling_site_label) %>%
  summarise(
    n_labelled_rung1 = n(),
    rung1_n_resident = sum(status_rung1 == "resident", na.rm = TRUE),
    rung1_n_migrant = sum(status_rung1 == "migrant", na.rm = TRUE),
    rung1_n_classes_present = sum(c(rung1_n_resident, rung1_n_migrant) > 0),
    rung1_max_class_prop = max(c(rung1_n_resident, rung1_n_migrant)) / n_labelled_rung1,
    .groups = "drop"
  )

rung2_site_summary <- live_labelled %>%
  group_by(sampling_site_label) %>%
  summarise(
    n_labelled_rung2 = n(),
    rung2_n_resident = sum(status_rung2 == "resident", na.rm = TRUE),
    rung2_n_nz_migrant = sum(status_rung2 == "NZ migrant", na.rm = TRUE),
    rung2_n_au_migrant = sum(status_rung2 == "AU migrant", na.rm = TRUE),
    rung2_n_classes_present = sum(c(rung2_n_resident, rung2_n_nz_migrant, rung2_n_au_migrant) > 0),
    rung2_max_class_prop = max(c(rung2_n_resident, rung2_n_nz_migrant, rung2_n_au_migrant)) / n_labelled_rung2,
    .groups = "drop"
  )

primary_reference_input <- live_tissue %>%
  filter(feather_type == "Primary", !is.na(sampling_site_label))

primary_site_counts <- primary_reference_input %>%
  group_by(sampling_site_label) %>%
  summarise(
    n_primary_reference_rings = n(),
    n_primary_labelled_rings = sum(status_known %in% TRUE, na.rm = TRUE),
    n_primary_unknown_status_rings = sum(!(status_known %in% TRUE), na.rm = TRUE),
    .groups = "drop"
  )

site_unit_summary <- site_lookup %>%
  select(
    sampling_site_label,
    sampling_location_general,
    sampling_location_specific,
    n_live_phenotype_rings,
    n_status_known,
    n_status_unknown,
    grouping_ambiguity,
    harmonization_note
  ) %>%
  left_join(paired_site_counts, by = "sampling_site_label") %>%
  left_join(rung1_site_summary, by = "sampling_site_label") %>%
  left_join(rung2_site_summary, by = "sampling_site_label") %>%
  left_join(primary_site_counts, by = "sampling_site_label") %>%
  mutate(
    n_paired_rings_total = coalesce(n_paired_rings_total, 0L),
    n_paired_rings_status_known = coalesce(n_paired_rings_status_known, 0L),
    n_paired_rings_status_unknown = coalesce(n_paired_rings_status_unknown, 0L),
    n_paired_complete_cnh = coalesce(n_paired_complete_cnh, 0L),
    n_missing_site_fields = coalesce(n_missing_site_fields, 0L),
    n_labelled_rung1 = coalesce(n_labelled_rung1, 0L),
    rung1_n_resident = coalesce(rung1_n_resident, 0L),
    rung1_n_migrant = coalesce(rung1_n_migrant, 0L),
    rung1_n_classes_present = coalesce(rung1_n_classes_present, 0L),
    n_labelled_rung2 = coalesce(n_labelled_rung2, 0L),
    rung2_n_resident = coalesce(rung2_n_resident, 0L),
    rung2_n_nz_migrant = coalesce(rung2_n_nz_migrant, 0L),
    rung2_n_au_migrant = coalesce(rung2_n_au_migrant, 0L),
    rung2_n_classes_present = coalesce(rung2_n_classes_present, 0L),
    n_primary_reference_rings = coalesce(n_primary_reference_rings, 0L),
    n_primary_labelled_rings = coalesce(n_primary_labelled_rings, 0L),
    n_primary_unknown_status_rings = coalesce(n_primary_unknown_status_rings, 0L),
    loso_flag_small_assessment = n_labelled_rung1 > 0 & n_labelled_rung1 < 3,
    loso_flag_single_class_rung1 = n_labelled_rung1 > 0 & rung1_n_classes_present < 2,
    loso_flag_missing_rung2_class = n_labelled_rung2 > 0 & rung2_n_classes_present < 3,
    validation_unit_note = case_when(
      n_labelled_rung1 == 0 ~ "Not represented in the current labelled paired live modelling table.",
      loso_flag_small_assessment & loso_flag_missing_rung2_class ~
        "Very small labelled site and incomplete Rung 2 class coverage; LOSO assessment would be unstable.",
      loso_flag_small_assessment ~
        "Very small labelled site; LOSO assessment would be unstable.",
      loso_flag_single_class_rung1 | loso_flag_missing_rung2_class ~
        "Labelled site is class-imbalanced; LOSO assessment would not represent all outcome classes.",
      TRUE ~ "Usable as a site validation unit, but still contributes to grouped blocked folds because site sizes remain uneven."
    )
  ) %>%
  arrange(desc(n_labelled_rung1), sampling_site_label)

write_csv(
  site_unit_summary,
  file.path(derived_dir, "model_stage2_site_unit_summary.csv")
)

site_validation_completeness <- bind_rows(
  tibble(
    dataset = "live_phenotypes_canonical",
    n_rows = nrow(live_phenotypes),
    n_distinct_sites = dplyr::n_distinct(live_phenotypes$sampling_site_label),
    n_missing_sampling_site_label = sum(is.na(live_phenotypes$sampling_site_label)),
    n_missing_sampling_location_specific = sum(is.na(live_phenotypes$sampling_location_specific)),
    n_missing_longitude = sum(is.na(live_phenotypes$longitude_capture)),
    n_missing_latitude = sum(is.na(live_phenotypes$latitude_capture)),
    n_missing_any_validation_field = sum(
      is.na(live_phenotypes$sampling_site_label) |
        is.na(live_phenotypes$sampling_location_specific) |
        is.na(live_phenotypes$longitude_capture) |
        is.na(live_phenotypes$latitude_capture)
    ),
    missing_ring_ids = collapse_unique(live_phenotypes$ring[
      is.na(live_phenotypes$sampling_site_label) |
        is.na(live_phenotypes$sampling_location_specific) |
        is.na(live_phenotypes$longitude_capture) |
        is.na(live_phenotypes$latitude_capture)
    ])
  ),
  tibble(
    dataset = "live_cn_paired_by_ring",
    n_rows = nrow(live_paired),
    n_distinct_sites = dplyr::n_distinct(live_paired$sampling_site_label),
    n_missing_sampling_site_label = sum(is.na(live_paired$sampling_site_label)),
    n_missing_sampling_location_specific = sum(is.na(live_paired$sampling_location_specific)),
    n_missing_longitude = sum(is.na(live_paired$longitude_capture)),
    n_missing_latitude = sum(is.na(live_paired$latitude_capture)),
    n_missing_any_validation_field = sum(
      is.na(live_paired$sampling_site_label) |
        is.na(live_paired$sampling_location_specific) |
        is.na(live_paired$longitude_capture) |
        is.na(live_paired$latitude_capture)
    ),
    missing_ring_ids = collapse_unique(live_paired$ring[
      is.na(live_paired$sampling_site_label) |
        is.na(live_paired$sampling_location_specific) |
        is.na(live_paired$longitude_capture) |
        is.na(live_paired$latitude_capture)
    ])
  ),
  tibble(
    dataset = "live_screening_ready_paired_labelled",
    n_rows = nrow(live_labelled),
    n_distinct_sites = dplyr::n_distinct(live_labelled$sampling_site_label),
    n_missing_sampling_site_label = sum(is.na(live_labelled$sampling_site_label)),
    n_missing_sampling_location_specific = sum(is.na(live_labelled$sampling_location_specific)),
    n_missing_longitude = sum(is.na(live_labelled$longitude_capture)),
    n_missing_latitude = sum(is.na(live_labelled$latitude_capture)),
    n_missing_any_validation_field = sum(
      is.na(live_labelled$sampling_site_label) |
        is.na(live_labelled$sampling_location_specific) |
        is.na(live_labelled$longitude_capture) |
        is.na(live_labelled$latitude_capture)
    ),
    missing_ring_ids = collapse_unique(live_labelled$ring[
      is.na(live_labelled$sampling_site_label) |
        is.na(live_labelled$sampling_location_specific) |
        is.na(live_labelled$longitude_capture) |
        is.na(live_labelled$latitude_capture)
    ])
  ),
  tibble(
    dataset = "live_primary_reference_input",
    n_rows = nrow(primary_reference_input),
    n_distinct_sites = dplyr::n_distinct(primary_reference_input$sampling_site_label),
    n_missing_sampling_site_label = sum(is.na(primary_reference_input$sampling_site_label)),
    n_missing_sampling_location_specific = sum(is.na(primary_reference_input$sampling_location_specific)),
    n_missing_longitude = sum(is.na(primary_reference_input$longitude_capture)),
    n_missing_latitude = sum(is.na(primary_reference_input$latitude_capture)),
    n_missing_any_validation_field = sum(
      is.na(primary_reference_input$sampling_site_label) |
        is.na(primary_reference_input$sampling_location_specific) |
        is.na(primary_reference_input$longitude_capture) |
        is.na(primary_reference_input$latitude_capture)
    ),
    missing_ring_ids = collapse_unique(primary_reference_input$ring[
      is.na(primary_reference_input$sampling_site_label) |
        is.na(primary_reference_input$sampling_location_specific) |
        is.na(primary_reference_input$longitude_capture) |
        is.na(primary_reference_input$latitude_capture)
    ])
  )
) %>%
  mutate(
    validation_ready = n_missing_any_validation_field == 0
  )

write_csv(
  site_validation_completeness,
  file.path(derived_dir, "model_stage2_site_validation_completeness.csv")
)

site_fold_map <- make_balanced_group_folds(
  data = live_labelled,
  outcome_col = "status_rung2",
  group_col = "sampling_site_label",
  v = 4L
)

site_block_folds <- live_labelled %>%
  left_join(site_fold_map, by = c("sampling_site_label" = "group_id")) %>%
  transmute(
    validation_layer = "blocked_site_cv",
    split_id,
    ring,
    sampling_site_label,
    sampling_location_specific,
    coord_group,
    status_rung1 = as.character(status_rung1),
    status_rung2 = as.character(status_rung2)
  )

write_csv(
  site_block_folds,
  file.path(derived_dir, "model_stage2_site_block_folds.csv")
)

site_loso_folds <- live_labelled %>%
  transmute(
    validation_layer = "leave_one_site_out_cv",
    split_id = paste0("Site__", make.names(sampling_site_label)),
    ring,
    sampling_site_label,
    sampling_location_specific,
    coord_group,
    status_rung1 = as.character(status_rung1),
    status_rung2 = as.character(status_rung2)
  )

write_csv(
  site_loso_folds,
  file.path(derived_dir, "model_stage2_site_loso_folds.csv")
)

fold_summary <- bind_rows(site_block_folds, site_loso_folds) %>%
  group_by(validation_layer, split_id) %>%
  summarise(
    n_assessment = n(),
    n_sites_assessment = dplyr::n_distinct(sampling_site_label),
    assessment_sites = collapse_unique(sampling_site_label),
    assessment_specific_locations = collapse_unique(sampling_location_specific),
    rung1_n_resident = sum(status_rung1 == "resident", na.rm = TRUE),
    rung1_n_migrant = sum(status_rung1 == "migrant", na.rm = TRUE),
    rung1_all_classes_present = rung1_n_resident > 0 & rung1_n_migrant > 0,
    rung2_n_resident = sum(status_rung2 == "resident", na.rm = TRUE),
    rung2_n_nz_migrant = sum(status_rung2 == "NZ migrant", na.rm = TRUE),
    rung2_n_au_migrant = sum(status_rung2 == "AU migrant", na.rm = TRUE),
    rung2_all_classes_present = rung2_n_resident > 0 &
      rung2_n_nz_migrant > 0 &
      rung2_n_au_migrant > 0,
    .groups = "drop"
  ) %>%
  mutate(
    fold_design_note = case_when(
      validation_layer == "blocked_site_cv" & rung2_all_classes_present ~
        "Recommended grouped site-block design. Every assessment fold retains all three Rung 2 classes.",
      validation_layer == "blocked_site_cv" ~
        "Grouped site-block design, but at least one assessment fold is missing a Rung 2 class.",
      validation_layer == "leave_one_site_out_cv" & n_assessment < 3 ~
        "LOSO split is very small and should be treated as a stress test rather than a stable performance estimate.",
      validation_layer == "leave_one_site_out_cv" & !rung2_all_classes_present ~
        "LOSO split is site-pure or site-imbalanced and does not represent all Rung 2 classes.",
      TRUE ~ NA_character_
    )
  ) %>%
  arrange(validation_layer, split_id)

write_csv(
  fold_summary,
  file.path(derived_dir, "model_stage2_site_fold_summary.csv")
)

primary_reference_summary <- primary_reference_input %>%
  group_by(
    sampling_site_label,
    sampling_location_specific,
    sampling_location_general
  ) %>%
  summarise(
    reference_object = "live_primary_breeding_site",
    grouping_variable = "sampling_site_label",
    contributor_source_table = "live_cn_tissue_summary.csv",
    n_rings_total = n(),
    n_labelled_rings = sum(status_known %in% TRUE, na.rm = TRUE),
    n_unknown_status_rings = sum(!(status_known %in% TRUE), na.rm = TRUE),
    n_with_d13c = sum(!is.na(normalised_d13c)),
    n_with_d15n = sum(!is.na(normalised_d15n)),
    n_with_d2h = sum(!is.na(normalised_d2h)),
    mean_longitude_capture = safe_mean(longitude_capture),
    mean_latitude_capture = safe_mean(latitude_capture),
    mean_normalised_d13c = safe_mean(normalised_d13c),
    mean_normalised_d15n = safe_mean(normalised_d15n),
    mean_normalised_d2h = safe_mean(normalised_d2h),
    var_normalised_d13c = safe_var(normalised_d13c),
    var_normalised_d15n = safe_var(normalised_d15n),
    var_normalised_d2h = safe_var(normalised_d2h),
    sd_normalised_d13c = safe_sd(normalised_d13c),
    sd_normalised_d15n = safe_sd(normalised_d15n),
    sd_normalised_d2h = safe_sd(normalised_d2h),
    source_rings = collapse_unique(ring),
    source_phenotype_row_ids = collapse_unique(source_phenotype_row_ids),
    source_measurement_ids = collapse_unique(source_measurement_ids),
    contributor_status_values = collapse_unique(status),
    .groups = "drop"
  ) %>%
  arrange(desc(n_rings_total), sampling_site_label)

write_csv(
  primary_reference_summary,
  file.path(derived_dir, "reference_live_primary_site_summary.csv")
)

primary_reference_covariance <- build_covariance_long(
  data = primary_reference_input,
  group_cols = c("sampling_site_label", "sampling_location_specific", "sampling_location_general"),
  feature_cols = c("normalised_d13c", "normalised_d15n", "normalised_d2h"),
  reference_object = "live_primary_breeding_site",
  grouping_variable = "sampling_site_label"
) %>%
  arrange(sampling_site_label, feature_x, feature_y)

write_csv(
  primary_reference_covariance,
  file.path(derived_dir, "reference_live_primary_site_covariance_long.csv")
)

winter_reference_summary <- museum_winter %>%
  group_by(geo_region) %>%
  summarise(
    reference_object = "museum_winter_breast_region",
    grouping_variable = "geo_region",
    contributor_source_table = "museum_screening_ready_breast_homogenate_winter.csv",
    n_specimens_total = n(),
    n_with_d13c = sum(!is.na(value_mean_normalised_d13c)),
    n_with_d15n = sum(!is.na(value_mean_normalised_d15n)),
    n_with_d18o = sum(!is.na(value_mean_normalised_d18o)),
    n_with_d2h = sum(!is.na(value_mean_normalised_d2h)),
    mean_normalised_d13c = safe_mean(value_mean_normalised_d13c),
    mean_normalised_d15n = safe_mean(value_mean_normalised_d15n),
    mean_normalised_d18o = safe_mean(value_mean_normalised_d18o),
    mean_normalised_d2h = safe_mean(value_mean_normalised_d2h),
    var_normalised_d13c = safe_var(value_mean_normalised_d13c),
    var_normalised_d15n = safe_var(value_mean_normalised_d15n),
    var_normalised_d18o = safe_var(value_mean_normalised_d18o),
    var_normalised_d2h = safe_var(value_mean_normalised_d2h),
    sd_normalised_d13c = safe_sd(value_mean_normalised_d13c),
    sd_normalised_d15n = safe_sd(value_mean_normalised_d15n),
    sd_normalised_d18o = safe_sd(value_mean_normalised_d18o),
    sd_normalised_d2h = safe_sd(value_mean_normalised_d2h),
    source_specimen_ids = collapse_unique(specimen_id),
    source_metadata_row_ids = collapse_unique(source_metadata_row_ids),
    source_branches_all = collapse_unique(source_branches_all),
    source_measurement_ids_normalised_d13c = collapse_unique(source_measurement_ids_normalised_d13c),
    source_measurement_ids_normalised_d15n = collapse_unique(source_measurement_ids_normalised_d15n),
    source_measurement_ids_normalised_d18o = collapse_unique(source_measurement_ids_normalised_d18o),
    source_measurement_ids_normalised_d2h = collapse_unique(source_measurement_ids_normalised_d2h),
    .groups = "drop"
  ) %>%
  arrange(geo_region)

write_csv(
  winter_reference_summary,
  file.path(derived_dir, "reference_museum_winter_breast_region_summary.csv")
)

winter_reference_covariance <- build_covariance_long(
  data = museum_winter %>%
    transmute(
      geo_region,
      normalised_d13c = value_mean_normalised_d13c,
      normalised_d15n = value_mean_normalised_d15n,
      normalised_d18o = value_mean_normalised_d18o,
      normalised_d2h = value_mean_normalised_d2h
    ),
  group_cols = c("geo_region"),
  feature_cols = c("normalised_d13c", "normalised_d15n", "normalised_d18o", "normalised_d2h"),
  reference_object = "museum_winter_breast_region",
  grouping_variable = "geo_region"
) %>%
  arrange(geo_region, feature_x, feature_y)

write_csv(
  winter_reference_covariance,
  file.path(derived_dir, "reference_museum_winter_breast_region_covariance_long.csv")
)

reference_manifest <- tibble(
  table_name = c(
    "reference_live_primary_site_summary.csv",
    "reference_live_primary_site_covariance_long.csv",
    "reference_museum_winter_breast_region_summary.csv",
    "reference_museum_winter_breast_region_covariance_long.csv"
  ),
  contributor_source_table = c(
    "live_cn_tissue_summary.csv",
    "live_cn_tissue_summary.csv",
    "museum_screening_ready_breast_homogenate_winter.csv",
    "museum_screening_ready_breast_homogenate_winter.csv"
  ),
  contributor_filter = c(
    "Primary tissue rows with non-missing sampling_site_label; includes status-known and unknown-status live birds, but only rows with primary tissue can contribute.",
    "Same contributors as reference_live_primary_site_summary.csv.",
    "Museum breast homogenate winter subset grouped only by geo_region; New Zealand specimens are not relabelled as residents or NZ migrants.",
    "Same contributors as reference_museum_winter_breast_region_summary.csv."
  ),
  intended_representation = c(
    "Breeding-site primary-feather reference centroids and within-site spread for site-aware modelling extensions.",
    "Pairwise primary-feather covariance structure within live breeding sites.",
    "Broad winter breast-feather reference structure for New Zealand versus Australia museum winter specimens.",
    "Pairwise covariance structure for the broad museum winter breast reference groups."
  ),
  leakage_safe_use = c(
    "Use only as a full-data descriptive reference outside model fitting. If site summaries are used as training features inside resampling, rebuild them on the training fold only so assessment-site information does not leak.",
    "Same leakage rule as the primary site summary: recompute within the training partition whenever covariance-derived features are used inside resampling.",
    "Use as descriptive regional reference summaries. If converted into reference-distance features inside cross-validation, recompute from the training reference subset only.",
    "Same leakage rule as the winter region summary: do not compute test-aware covariance features from the full table during resampling."
  )
)

write_csv(
  reference_manifest,
  file.path(derived_dir, "reference_object_manifest.csv")
)
