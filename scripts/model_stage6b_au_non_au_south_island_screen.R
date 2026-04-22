#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(tidyr)
})

source(file.path("R", "model_stage_repeatability_helpers.R"))
source(file.path("R", "model_stage2_helpers.R"))
source(file.path("R", "model_stage6_au_non_au_helpers.R"))

derived_dir <- file.path("data", "derived")
ensure_dir(derived_dir)

positive_class <- "AU"
indeterminate_threshold <- 0.8

format_num <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits))
}

evaluate_direct_binary_split <- function(data, fold_rows, candidate, seed) {
  test_ids <- fold_rows$ring
  train_data <- data %>% filter(!(ring %in% test_ids))
  test_data <- data %>% filter(ring %in% test_ids)

  fit <- fit_predict_ridge(
    train_data = train_data,
    test_data = test_data,
    candidate = list(
      predictors = candidate$predictors,
      distance_vars = candidate$distance_vars,
      distance_name = candidate$distance_name
    ),
    outcome_col = "status_au_non_au",
    family = "binomial",
    seed = seed
  )

  if (!identical(fit$status, "ok")) {
    return(list(status = fit$status, message = fit$message %||% NA_character_))
  }

  pred <- augment_binary_predictions(
    predictions = fit$predictions %>%
      mutate(truth = factor(truth, levels = c("non_AU", "AU"))),
    positive_class = positive_class,
    indeterminate_threshold = indeterminate_threshold
  )

  cal <- binary_calibration_summary(pred, positive_class = positive_class)

  metrics <- binary_au_prediction_metrics(
    prediction_tbl = pred,
    positive_class = positive_class,
    indeterminate_threshold = indeterminate_threshold
  ) %>%
    mutate(
      selected_lambda = fit$metrics$selected_lambda[[1]],
      ece = cal %>% filter(metric == "ece") %>% pull(value),
      calibration_intercept = cal %>% filter(metric == "calibration_intercept") %>% pull(value),
      calibration_slope = cal %>% filter(metric == "calibration_slope") %>% pull(value)
    ) %>%
    relocate(selected_lambda, ece, calibration_intercept, calibration_slope, .after = brier)

  list(
    status = "ok",
    predictions = pred,
    metrics = metrics,
    predictors_used = paste(fit$predictors_used, collapse = "|"),
    n_analysis = nrow(train_data),
    n_assessment = nrow(test_data)
  )
}

summarise_binary_split_from_predictions <- function(prediction_tbl) {
  pred <- augment_binary_predictions(
    predictions = prediction_tbl,
    positive_class = positive_class,
    indeterminate_threshold = indeterminate_threshold
  )

  cal <- binary_calibration_summary(pred, positive_class = positive_class)

  binary_au_prediction_metrics(
    prediction_tbl = pred,
    positive_class = positive_class,
    indeterminate_threshold = indeterminate_threshold
  ) %>%
    mutate(
      selected_lambda = NA_real_,
      ece = cal %>% filter(metric == "ece") %>% pull(value),
      calibration_intercept = cal %>% filter(metric == "calibration_intercept") %>% pull(value),
      calibration_slope = cal %>% filter(metric == "calibration_slope") %>% pull(value)
    ) %>%
    relocate(selected_lambda, ece, calibration_intercept, calibration_slope, .after = brier)
}

matched_tbl <- read_csv(
  file.path(derived_dir, "extension_au_non_au_south_island_matched_paired_labelled.csv"),
  show_col_types = FALSE
)

breast_tbl <- read_csv(
  file.path(derived_dir, "extension_au_non_au_south_island_breast_features.csv"),
  show_col_types = FALSE
)

paired_tbl <- read_csv(
  file.path(derived_dir, "extension_au_non_au_south_island_paired_contrast_features.csv"),
  show_col_types = FALSE
)

candidate_input_tbl <- read_csv(
  file.path(derived_dir, "extension_au_non_au_south_island_candidate_definitions.csv"),
  show_col_types = FALSE
)

site_block_folds <- read_csv(
  file.path(derived_dir, "extension_au_non_au_south_island_site_block_folds.csv"),
  show_col_types = FALSE
)

site_loso_folds <- read_csv(
  file.path(derived_dir, "extension_au_non_au_south_island_site_loso_folds.csv"),
  show_col_types = FALSE
)

subset_check <- read_csv(
  file.path(derived_dir, "extension_au_non_au_south_island_subset_check.csv"),
  show_col_types = FALSE
)

fold_summary <- read_csv(
  file.path(derived_dir, "extension_au_non_au_south_island_fold_summary.csv"),
  show_col_types = FALSE
)

imbalance_comparison <- read_csv(
  file.path(derived_dir, "extension_au_non_au_south_island_imbalance_comparison.csv"),
  show_col_types = FALSE
)

stage5_predictions <- read_csv(
  file.path(derived_dir, "model_stage5_siteaware_assessment_predictions.csv"),
  show_col_types = FALSE
)

stage6_full_summary <- read_csv(
  file.path(derived_dir, "model_stage6_au_non_au_model_summary.csv"),
  show_col_types = FALSE
)

if (!all(subset_check$rings_match_frozen_subset)) {
  stop("The South-Island-only AU vs non_AU extension tables do not all match the frozen 50-ring subset.")
}

expected_rings <- sort(unique(matched_tbl$ring))

if (!identical(sort(unique(breast_tbl$ring)), expected_rings)) {
  stop("South-Island-only breast feature table does not match the frozen ring set.")
}

if (!identical(sort(unique(paired_tbl$ring)), expected_rings)) {
  stop("South-Island-only paired-contrast feature table does not match the frozen ring set.")
}

fold_lookup <- bind_rows(site_block_folds, site_loso_folds) %>%
  distinct(ring, coord_group)

prepare_binary_table <- function(tbl) {
  tbl %>%
    left_join(fold_lookup, by = "ring") %>%
    mutate(
      status_au_non_au = factor(status_au_non_au, levels = c("non_AU", "AU"))
    )
}

breast_data <- prepare_binary_table(breast_tbl)
paired_data <- prepare_binary_table(paired_tbl)

all_folds <- bind_rows(site_block_folds, site_loso_folds) %>%
  arrange(validation_layer, split_id, ring)

write_csv(
  tibble(
    parameter = c("positive_class", "indeterminate_threshold"),
    value = c(positive_class, as.character(indeterminate_threshold))
  ),
  file.path(derived_dir, "model_stage6b_au_non_au_south_island_config.csv")
)

candidate_defs <- tibble(
  candidate_id = c(
    "au_non_au_si_a_breast_cnh",
    "au_non_au_si_b_breast_cnh_plus_contrast_cnh",
    "benchmark_si_r2s2_collapsed_au_non_au",
    "benchmark_si_r2c_collapsed_au_non_au"
  ),
  candidate_label = c(
    "South-Island-only AU vs non_AU A: Breast-only C/N/H",
    "South-Island-only AU vs non_AU B: Breast C/N/H + paired contrasts",
    "Benchmark SI: R2 S2 collapsed to AU vs non_AU",
    "Benchmark SI: R2 C collapsed to AU vs non_AU"
  ),
  model_type = c(
    "ridge_binomial",
    "ridge_binomial",
    "collapsed_multiclass_benchmark",
    "collapsed_multiclass_benchmark"
  ),
  candidate_group = c(
    "direct_binary_baseline",
    "direct_binary_extension",
    "collapsed_3class_benchmark",
    "collapsed_3class_benchmark"
  ),
  caution_level = c(
    "south_island_breast_cnh_binary_baseline",
    "south_island_breast_cnh_plus_paired_h_high_caution",
    "collapsed_from_siteaware_staged_rung2_on_si_subset",
    "collapsed_from_siteaware_direct_rung2_on_si_subset"
  ),
  source_table = c(
    "extension_au_non_au_south_island_breast_features.csv",
    "extension_au_non_au_south_island_paired_contrast_features.csv",
    "model_stage5_siteaware_assessment_predictions.csv",
    "model_stage5_siteaware_assessment_predictions.csv"
  ),
  predictors = c(
    "normalised_d13c_Breast|normalised_d15n_Breast|normalised_d2h_Breast",
    "normalised_d13c_Breast|normalised_d15n_Breast|normalised_d2h_Breast|abs_delta_d13c|abs_delta_d15n|abs_delta_d2h",
    "collapsed_3class_probabilities_to_binary",
    "collapsed_3class_probabilities_to_binary"
  )
)

write_csv(
  candidate_defs,
  file.path(derived_dir, "model_stage6b_au_non_au_south_island_candidate_definitions.csv")
)

direct_candidates <- list(
  list(
    candidate_id = "au_non_au_si_a_breast_cnh",
    candidate_label = "South-Island-only AU vs non_AU A: Breast-only C/N/H",
    model_type = "ridge_binomial",
    candidate_group = "direct_binary_baseline",
    caution_level = "south_island_breast_cnh_binary_baseline",
    predictors = c(
      "normalised_d13c_Breast",
      "normalised_d15n_Breast",
      "normalised_d2h_Breast"
    ),
    distance_vars = character(0),
    distance_name = NA_character_,
    source_name = "breast_data"
  ),
  list(
    candidate_id = "au_non_au_si_b_breast_cnh_plus_contrast_cnh",
    candidate_label = "South-Island-only AU vs non_AU B: Breast C/N/H + paired contrasts",
    model_type = "ridge_binomial",
    candidate_group = "direct_binary_extension",
    caution_level = "south_island_breast_cnh_plus_paired_h_high_caution",
    predictors = c(
      "normalised_d13c_Breast",
      "normalised_d15n_Breast",
      "normalised_d2h_Breast",
      "abs_delta_d13c",
      "abs_delta_d15n",
      "abs_delta_d2h"
    ),
    distance_vars = character(0),
    distance_name = NA_character_,
    source_name = "paired_data"
  )
)

split_keys <- all_folds %>%
  distinct(validation_layer, split_id) %>%
  arrange(validation_layer, split_id)

direct_split_rows <- list()
direct_prediction_rows <- list()
row_idx <- 1L

for (cand in direct_candidates) {
  data_obj <- if (identical(cand$source_name, "breast_data")) breast_data else paired_data

  for (split_idx in seq_len(nrow(split_keys))) {
    split_meta <- split_keys[split_idx, ]
    fold_rows <- all_folds %>%
      filter(
        validation_layer == split_meta$validation_layer[[1]],
        split_id == split_meta$split_id[[1]]
      )

    split_eval <- evaluate_direct_binary_split(
      data = data_obj,
      fold_rows = fold_rows,
      candidate = cand,
      seed = 970000L + row_idx
    )

    split_row <- fold_rows %>%
      summarise(
        validation_layer = first(validation_layer),
        split_id = first(split_id),
        n_assessment = n(),
        n_sites_assessment = n_distinct(sampling_site_label),
        assessment_sites = paste(sort(unique(sampling_site_label)), collapse = "|"),
        .groups = "drop"
      ) %>%
      mutate(
        candidate_id = cand$candidate_id,
        candidate_label = cand$candidate_label,
        model_type = cand$model_type,
        candidate_group = cand$candidate_group,
        caution_level = cand$caution_level,
        fit_status = split_eval$status
      )

    if (!identical(split_eval$status, "ok")) {
      split_row <- split_row %>%
        mutate(
          n_analysis = nrow(data_obj) - n_assessment,
          log_loss = NA_real_,
          balanced_accuracy = NA_real_,
          accuracy = NA_real_,
          brier = NA_real_,
          selected_lambda = NA_real_,
          ece = NA_real_,
          calibration_intercept = NA_real_,
          calibration_slope = NA_real_,
          mean_top_prob = NA_real_,
          weak_support_rate = NA_real_,
          indeterminate_rate = NA_real_,
          au_sensitivity = NA_real_,
          au_specificity = NA_real_,
          au_ppv = NA_real_,
          au_npv = NA_real_,
          predictors_used = paste(cand$predictors, collapse = "|"),
          fit_message = split_eval$message %||% NA_character_
        )
    } else {
      split_row <- split_row %>%
        mutate(
          n_analysis = split_eval$n_analysis,
          log_loss = split_eval$metrics$log_loss[[1]],
          balanced_accuracy = split_eval$metrics$balanced_accuracy[[1]],
          accuracy = split_eval$metrics$accuracy[[1]],
          brier = split_eval$metrics$brier[[1]],
          selected_lambda = split_eval$metrics$selected_lambda[[1]],
          ece = split_eval$metrics$ece[[1]],
          calibration_intercept = split_eval$metrics$calibration_intercept[[1]],
          calibration_slope = split_eval$metrics$calibration_slope[[1]],
          mean_top_prob = split_eval$metrics$mean_top_prob[[1]],
          weak_support_rate = split_eval$metrics$weak_support_rate[[1]],
          indeterminate_rate = split_eval$metrics$indeterminate_rate[[1]],
          au_sensitivity = split_eval$metrics$au_sensitivity[[1]],
          au_specificity = split_eval$metrics$au_specificity[[1]],
          au_ppv = split_eval$metrics$au_ppv[[1]],
          au_npv = split_eval$metrics$au_npv[[1]],
          predictors_used = split_eval$predictors_used,
          fit_message = NA_character_
        )

      direct_prediction_rows[[length(direct_prediction_rows) + 1L]] <- split_eval$predictions %>%
        mutate(
          candidate_id = cand$candidate_id,
          candidate_label = cand$candidate_label,
          model_type = cand$model_type,
          candidate_group = cand$candidate_group,
          caution_level = cand$caution_level,
          validation_layer = split_meta$validation_layer[[1]],
          split_id = split_meta$split_id[[1]]
        ) %>%
        relocate(candidate_id, candidate_label, model_type, candidate_group, caution_level, validation_layer, split_id)
    }

    direct_split_rows[[length(direct_split_rows) + 1L]] <- split_row
    row_idx <- row_idx + 1L
  }
}

direct_split_metrics <- bind_rows(direct_split_rows)
direct_predictions <- bind_rows(direct_prediction_rows)

benchmark_map <- tibble(
  stage5_candidate_id = c(
    "r2s_b_soft_hierarchical_r1c_then_breast_cnh",
    "r2_c_breast_cnh_plus_paired_contrast"
  ),
  candidate_id = c(
    "benchmark_si_r2s2_collapsed_au_non_au",
    "benchmark_si_r2c_collapsed_au_non_au"
  ),
  candidate_label = c(
    "Benchmark SI: R2 S2 collapsed to AU vs non_AU",
    "Benchmark SI: R2 C collapsed to AU vs non_AU"
  ),
  model_type = "collapsed_multiclass_benchmark",
  candidate_group = "collapsed_3class_benchmark",
  caution_level = c(
    "collapsed_from_siteaware_staged_rung2_on_si_subset",
    "collapsed_from_siteaware_direct_rung2_on_si_subset"
  )
)

benchmark_source <- stage5_predictions %>%
  filter(
    ring %in% expected_rings,
    candidate_id %in% benchmark_map$stage5_candidate_id,
    validation_layer %in% c("blocked_site_cv", "leave_one_site_out_cv")
  ) %>%
  inner_join(benchmark_map, by = c("candidate_id" = "stage5_candidate_id")) %>%
  transmute(
    candidate_id = candidate_id.y,
    candidate_label = candidate_label.y,
    model_type = model_type.y,
    candidate_group = candidate_group.y,
    caution_level = caution_level,
    validation_layer,
    split_id,
    ring,
    coord_group,
    longitude_capture,
    latitude_capture,
    truth,
    resident,
    `NZ migrant`,
    `AU migrant`
  )

benchmark_predictions <- benchmark_source %>%
  collapse_multiclass_predictions_to_au_binary(
    indeterminate_threshold = indeterminate_threshold
  ) %>%
  bind_cols(
    benchmark_source %>%
      select(candidate_id, candidate_label, model_type, candidate_group, caution_level, validation_layer, split_id)
  ) %>%
  relocate(candidate_id, candidate_label, model_type, candidate_group, caution_level, validation_layer, split_id)

benchmark_split_metrics <- benchmark_predictions %>%
  group_by(candidate_id, candidate_label, model_type, candidate_group, caution_level, validation_layer, split_id) %>%
  group_modify(
    ~ {
      fold_rows <- all_folds %>%
        filter(
          validation_layer == .y$validation_layer[[1]],
          split_id == .y$split_id[[1]]
        )

      metrics <- summarise_binary_split_from_predictions(.x)

      tibble(
        n_assessment = nrow(.x),
        n_sites_assessment = n_distinct(fold_rows$sampling_site_label),
        assessment_sites = paste(sort(unique(fold_rows$sampling_site_label)), collapse = "|"),
        n_analysis = length(expected_rings) - nrow(.x),
        fit_status = "ok",
        log_loss = metrics$log_loss[[1]],
        balanced_accuracy = metrics$balanced_accuracy[[1]],
        accuracy = metrics$accuracy[[1]],
        brier = metrics$brier[[1]],
        selected_lambda = NA_real_,
        ece = metrics$ece[[1]],
        calibration_intercept = metrics$calibration_intercept[[1]],
        calibration_slope = metrics$calibration_slope[[1]],
        mean_top_prob = metrics$mean_top_prob[[1]],
        weak_support_rate = metrics$weak_support_rate[[1]],
        indeterminate_rate = metrics$indeterminate_rate[[1]],
        au_sensitivity = metrics$au_sensitivity[[1]],
        au_specificity = metrics$au_specificity[[1]],
        au_ppv = metrics$au_ppv[[1]],
        au_npv = metrics$au_npv[[1]],
        predictors_used = "collapsed_3class_probabilities_to_binary",
        fit_message = NA_character_
      )
    }
  ) %>%
  ungroup()

all_predictions <- bind_rows(direct_predictions, benchmark_predictions)
all_split_metrics <- bind_rows(direct_split_metrics, benchmark_split_metrics)

model_summary <- summarise_binary_au_models(all_split_metrics) %>%
  left_join(
    all_predictions %>%
      group_by(candidate_id, candidate_label, model_type, candidate_group, caution_level, validation_layer) %>%
      group_modify(
        ~ binary_au_pooled_summary(
          prediction_tbl = .x,
          positive_class = positive_class,
          indeterminate_threshold = indeterminate_threshold
        )
      ) %>%
      ungroup(),
    by = c("candidate_id", "candidate_label", "model_type", "candidate_group", "caution_level", "validation_layer")
  ) %>%
  group_by(validation_layer) %>%
  mutate(
    pooled_log_loss_rank = min_rank(pooled_log_loss),
    pooled_bal_accuracy_rank = min_rank(desc(pooled_bal_accuracy))
  ) %>%
  ungroup() %>%
  arrange(validation_layer, pooled_log_loss_rank, pooled_bal_accuracy_rank)

calibration_summary <- all_predictions %>%
  group_by(candidate_id, candidate_label, model_type, candidate_group, caution_level, validation_layer) %>%
  group_modify(
    ~ binary_calibration_summary(
      predictions = .x %>% select(truth, .pred_positive),
      positive_class = positive_class
    )
  ) %>%
  ungroup()

comparison_summary <- model_summary %>%
  select(
    validation_layer,
    candidate_id,
    candidate_label,
    candidate_group,
    caution_level,
    pooled_log_loss,
    pooled_bal_accuracy,
    pooled_au_sensitivity,
    pooled_au_specificity,
    pooled_au_ppv,
    pooled_au_npv,
    pooled_brier,
    pooled_ece,
    pooled_indeterminate_rate,
    mean_split_log_loss,
    sd_split_log_loss,
    mean_split_bal_accuracy,
    sd_split_bal_accuracy,
    mean_split_indeterminate_rate
  ) %>%
  arrange(validation_layer, pooled_log_loss, desc(pooled_bal_accuracy))

full_map <- tibble(
  south_candidate_id = c(
    "au_non_au_si_a_breast_cnh",
    "au_non_au_si_b_breast_cnh_plus_contrast_cnh",
    "benchmark_si_r2s2_collapsed_au_non_au",
    "benchmark_si_r2c_collapsed_au_non_au"
  ),
  full_candidate_id = c(
    "au_non_au_a_breast_cnh",
    "au_non_au_b_breast_cnh_plus_contrast_cnh",
    "benchmark_r2s2_collapsed_au_non_au",
    "benchmark_r2c_collapsed_au_non_au"
  )
)

vs_full <- comparison_summary %>%
  left_join(full_map, by = c("candidate_id" = "south_candidate_id")) %>%
  left_join(
    stage6_full_summary %>%
      select(
        validation_layer,
        candidate_id,
        full_pooled_log_loss = pooled_log_loss,
        full_pooled_bal_accuracy = pooled_bal_accuracy,
        full_pooled_au_sensitivity = pooled_au_sensitivity,
        full_pooled_au_specificity = pooled_au_specificity,
        full_pooled_brier = pooled_brier,
        full_pooled_ece = pooled_ece,
        full_pooled_indeterminate_rate = pooled_indeterminate_rate
      ),
    by = c("validation_layer", "full_candidate_id" = "candidate_id")
  ) %>%
  mutate(
    delta_log_loss_vs_full = pooled_log_loss - full_pooled_log_loss,
    delta_bal_accuracy_vs_full = pooled_bal_accuracy - full_pooled_bal_accuracy,
    delta_au_sensitivity_vs_full = pooled_au_sensitivity - full_pooled_au_sensitivity,
    delta_au_specificity_vs_full = pooled_au_specificity - full_pooled_au_specificity,
    delta_brier_vs_full = pooled_brier - full_pooled_brier,
    delta_ece_vs_full = pooled_ece - full_pooled_ece,
    delta_indeterminate_rate_vs_full = pooled_indeterminate_rate - full_pooled_indeterminate_rate
  ) %>%
  select(
    validation_layer,
    candidate_id,
    candidate_label,
    full_candidate_id,
    pooled_log_loss,
    full_pooled_log_loss,
    delta_log_loss_vs_full,
    pooled_bal_accuracy,
    full_pooled_bal_accuracy,
    delta_bal_accuracy_vs_full,
    pooled_au_sensitivity,
    full_pooled_au_sensitivity,
    delta_au_sensitivity_vs_full,
    pooled_au_specificity,
    full_pooled_au_specificity,
    delta_au_specificity_vs_full,
    pooled_brier,
    full_pooled_brier,
    delta_brier_vs_full,
    pooled_ece,
    full_pooled_ece,
    delta_ece_vs_full,
    pooled_indeterminate_rate,
    full_pooled_indeterminate_rate,
    delta_indeterminate_rate_vs_full
  ) %>%
  arrange(validation_layer, candidate_id)

write_csv(
  all_folds,
  file.path(derived_dir, "model_stage6b_au_non_au_south_island_outer_folds.csv")
)

write_csv(
  all_split_metrics,
  file.path(derived_dir, "model_stage6b_au_non_au_south_island_split_metrics.csv")
)

write_csv(
  all_predictions,
  file.path(derived_dir, "model_stage6b_au_non_au_south_island_assessment_predictions.csv")
)

write_csv(
  model_summary,
  file.path(derived_dir, "model_stage6b_au_non_au_south_island_model_summary.csv")
)

write_csv(
  calibration_summary,
  file.path(derived_dir, "model_stage6b_au_non_au_south_island_calibration_summary.csv")
)

write_csv(
  comparison_summary,
  file.path(derived_dir, "model_stage6b_au_non_au_south_island_vs_benchmarks.csv")
)

write_csv(
  vs_full,
  file.path(derived_dir, "model_stage6b_au_non_au_south_island_vs_full_extension.csv")
)

blocked_top <- comparison_summary %>%
  filter(validation_layer == "blocked_site_cv")

loso_top <- comparison_summary %>%
  filter(validation_layer == "leave_one_site_out_cv")

si_imbalance_row <- imbalance_comparison %>%
  filter(dataset_scope == "south_island_only_extension")

summary_lines <- c(
  "# South-Island-only AU vs non_AU exploratory extension",
  "",
  paste0(
    "Main grouped site-block layer: best pooled log loss was ",
    blocked_top$candidate_label[[1]],
    " (",
    format_num(blocked_top$pooled_log_loss[[1]]),
    "), with pooled balanced accuracy ",
    format_num(blocked_top$pooled_bal_accuracy[[1]]),
    " and indeterminate-rate proxy ",
    format_num(blocked_top$pooled_indeterminate_rate[[1]])
  ),
  paste0(
    "LOSO stress test: best pooled log loss was ",
    loso_top$candidate_label[[1]],
    " (",
    format_num(loso_top$pooled_log_loss[[1]]),
    "), with pooled balanced accuracy ",
    format_num(loso_top$pooled_bal_accuracy[[1]]),
    " and indeterminate-rate proxy ",
    format_num(loso_top$pooled_indeterminate_rate[[1]])
  ),
  "",
  paste0(
    "South-Island-only subset size: ",
    nrow(matched_tbl),
    " rings; AU = ",
    sum(matched_tbl$status_au_non_au == "AU", na.rm = TRUE),
    ", non_AU = ",
    sum(matched_tbl$status_au_non_au == "non_AU", na.rm = TRUE),
    "."
  ),
  paste0(
    "Imbalance still remains concentrated at Tasman River: ",
    si_imbalance_row$largest_AU_site_n_AU[[1]],
    " of ",
    si_imbalance_row$n_AU[[1]],
    " AU birds."
  ),
  "",
  "Indeterminate burden here is a deterministic confidence proxy based on binary top-class probability < 0.80.",
  "It is reported on the same footing for the South-Island-only direct models and the collapsed 3-class benchmarks, but it is not a Stage 4-style uncertainty-aware operational rule."
)

writeLines(
  summary_lines,
  file.path(derived_dir, "model_stage6b_au_non_au_south_island_summary.md")
)
