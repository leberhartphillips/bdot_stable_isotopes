#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(tidyr)
})

source(file.path("R", "model_stage_repeatability_helpers.R"))
source(file.path("R", "model_stage2_helpers.R"))
source(file.path("R", "model_stage3_helpers.R"))
source(file.path("R", "model_stage4_helpers.R"))
source(file.path("R", "model_stage2b_rung2_staged_helpers.R"))
source(file.path("R", "model_stage4b_rung2_staged_helpers.R"))
source(file.path("R", "model_stage5_siteaware_helpers.R"))

derived_dir <- file.path("data", "derived")
ensure_dir(derived_dir)

format_num <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits))
}

binary_pooled_summary <- function(prediction_tbl) {
  cal <- binary_calibration_summary(
    predictions = prediction_tbl %>%
      transmute(
        truth = factor(truth, levels = c("resident", "migrant")),
        .pred_positive = .pred_positive
      ),
    positive_class = "migrant"
  )

  tibble(
    pooled_log_loss = binary_log_loss(prediction_tbl$truth, prediction_tbl$.pred_positive, positive_class = "migrant"),
    pooled_bal_accuracy = macro_bal_accuracy(
      factor(prediction_tbl$truth, levels = c("resident", "migrant")),
      factor(prediction_tbl$.pred_class, levels = c("resident", "migrant"))
    ),
    pooled_accuracy = mean(prediction_tbl$.pred_class == prediction_tbl$truth),
    pooled_brier = binary_brier(prediction_tbl$truth, prediction_tbl$.pred_positive, positive_class = "migrant"),
    pooled_ece = cal %>% filter(metric == "ece") %>% pull(value),
    pooled_mean_top_prob = mean(prediction_tbl$top_prob, na.rm = TRUE),
    pooled_weak_support_rate = mean(prediction_tbl$weak_support, na.rm = TRUE)
  )
}

multiclass_pooled_summary <- function(prediction_tbl, class_levels) {
  prob_mat <- as.matrix(prediction_tbl[, class_levels, drop = FALSE])
  cal <- multiclass_calibration_summary(
    predictions = prediction_tbl %>%
      transmute(
        truth = factor(truth, levels = class_levels),
        resident = resident,
        `NZ migrant` = `NZ migrant`,
        `AU migrant` = `AU migrant`
      ),
    class_levels = class_levels
  )

  tibble(
    pooled_log_loss = multiclass_log_loss(prediction_tbl$truth, prob_mat),
    pooled_bal_accuracy = macro_bal_accuracy(
      factor(prediction_tbl$truth, levels = class_levels),
      factor(prediction_tbl$.pred_class, levels = class_levels)
    ),
    pooled_accuracy = mean(prediction_tbl$.pred_class == prediction_tbl$truth),
    pooled_brier = multiclass_brier(prediction_tbl$truth, prob_mat),
    pooled_ece = cal %>% filter(metric == "confidence_ece", class == "overall") %>% pull(value),
    pooled_mean_top_prob = mean(prediction_tbl$top_prob, na.rm = TRUE),
    pooled_weak_support_rate = mean(prediction_tbl$weak_support, na.rm = TRUE)
  )
}

evaluate_binary_candidate_split <- function(train_data, test_data, candidate, seed) {
  fit <- fit_predict_ridge(
    train_data = train_data,
    test_data = test_data,
    candidate = list(
      predictors = candidate$predictors,
      distance_vars = candidate$distance_vars,
      distance_name = candidate$distance_name
    ),
    outcome_col = "status_rung1",
    family = "binomial",
    seed = seed
  )

  if (!identical(fit$status, "ok")) {
    return(list(status = fit$status, message = fit$message %||% NA_character_))
  }

  pred <- fit$predictions %>%
    mutate(
      top_prob = pmax(.pred_negative, .pred_positive),
      weak_support = top_prob < 0.5
    )

  cal <- binary_calibration_summary(
    predictions = pred %>%
      transmute(
        truth = factor(truth, levels = c("resident", "migrant")),
        .pred_positive = .pred_positive
      ),
    positive_class = "migrant"
  )

  metrics <- fit$metrics %>%
    transmute(
      log_loss,
      balanced_accuracy,
      accuracy,
      brier,
      selected_lambda,
      ece = cal %>% filter(metric == "ece") %>% pull(value),
      calibration_intercept = cal %>% filter(metric == "calibration_intercept") %>% pull(value),
      calibration_slope = cal %>% filter(metric == "calibration_slope") %>% pull(value),
      mean_top_prob = mean(pred$top_prob, na.rm = TRUE),
      weak_support_rate = mean(pred$weak_support, na.rm = TRUE)
    )

  list(
    status = "ok",
    predictions = pred,
    metrics = metrics,
    predictors_used = paste(fit$predictors_used, collapse = "|")
  )
}

evaluate_multiclass_candidate_split <- function(train_data, test_data, candidate, seed, class_levels) {
  fit <- fit_predict_ridge(
    train_data = train_data,
    test_data = test_data,
    candidate = list(
      predictors = candidate$predictors,
      distance_vars = candidate$distance_vars,
      distance_name = candidate$distance_name
    ),
    outcome_col = "status_rung2",
    family = "multinomial",
    seed = seed
  )

  if (!identical(fit$status, "ok")) {
    return(list(status = fit$status, message = fit$message %||% NA_character_))
  }

  pred <- fit$predictions %>%
    mutate(
      top_prob = pmax(resident, `NZ migrant`, `AU migrant`),
      weak_support = top_prob < 0.5
    )

  cal <- multiclass_calibration_summary(
    predictions = pred %>% select(truth, all_of(class_levels)),
    class_levels = class_levels
  )

  metrics <- fit$metrics %>%
    transmute(
      log_loss,
      balanced_accuracy,
      accuracy,
      brier,
      selected_lambda,
      ece = cal %>% filter(metric == "confidence_ece", class == "overall") %>% pull(value),
      mean_top_prob = mean(pred$top_prob, na.rm = TRUE),
      weak_support_rate = mean(pred$weak_support, na.rm = TRUE)
    )

  list(
    status = "ok",
    predictions = pred,
    metrics = metrics,
    predictors_used = paste(fit$predictors_used, collapse = "|")
  )
}

evaluate_staged_candidate_split <- function(train_data, test_data, r1_candidate, migrant_candidate, seed_r1, seed_sub, class_levels) {
  fit <- fit_stage4b_r2s2(
    train_data = train_data,
    r1_candidate = r1_candidate,
    migrant_candidate = migrant_candidate,
    seed_r1 = seed_r1,
    seed_submodel = seed_sub
  )

  if (!identical(fit$status, "ok")) {
    return(list(status = fit$status, message = fit$detail$message %||% NA_character_))
  }

  pred <- predict_stage4b_r2s2(fit, test_data) %>%
    bind_cols(
      test_data %>% select(ring, coord_group, longitude_capture, latitude_capture, truth = status_rung2)
    ) %>%
    relocate(ring, coord_group, longitude_capture, latitude_capture, truth) %>%
    rename(.pred_class = base_class) %>%
    mutate(
      top_prob = pmax(resident, `NZ migrant`, `AU migrant`),
      weak_support = top_prob < 0.5
    )

  metrics <- stage2b_prediction_metrics(pred, class_levels) %>%
    rename(
      log_loss = log_loss,
      balanced_accuracy = balanced_accuracy,
      accuracy = accuracy,
      brier = brier,
      mean_top_prob = mean_top_prob,
      weak_support_rate = weak_support_rate
    )

  cal <- multiclass_calibration_summary(
    predictions = pred %>% select(truth, all_of(class_levels)),
    class_levels = class_levels
  )

  metrics <- metrics %>%
    mutate(
      ece = cal %>% filter(metric == "confidence_ece", class == "overall") %>% pull(value),
      selected_lambda = NA_real_
    )

  list(
    status = "ok",
    predictions = pred,
    metrics = metrics,
    predictors_used = paste(c(r1_candidate$predictors, migrant_candidate$predictors), collapse = "|")
  )
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) {
    y
  } else {
    x
  }
}

class_levels_r2 <- c("resident", "NZ migrant", "AU migrant")

stage_data <- read_csv(
  file.path(derived_dir, "live_screening_ready_paired_labelled.csv"),
  show_col_types = FALSE
) %>%
  derive_stage2_features() %>%
  mutate(
    status_nzau = factor(
      ifelse(status_rung2 == "AU migrant", "AU migrant", "NZ migrant"),
      levels = c("NZ migrant", "AU migrant")
    )
  )

live_primary_reference_all <- read_csv(
  file.path(derived_dir, "live_cn_tissue_summary.csv"),
  show_col_types = FALSE
) %>%
  filter(feather_type == "Primary", !is.na(sampling_site_label)) %>%
  transmute(
    ring,
    status_known,
    sampling_site_label,
    sampling_location_specific,
    normalised_d13c,
    normalised_d15n,
    normalised_d2h
  )

museum_reference_all <- read_csv(
  file.path(derived_dir, "museum_screening_ready_breast_homogenate_winter.csv"),
  show_col_types = FALSE
) %>%
  transmute(
    specimen_id,
    geo_region,
    normalised_d13c_Breast = value_mean_normalised_d13c,
    normalised_d15n_Breast = value_mean_normalised_d15n,
    normalised_d2h_Breast = value_mean_normalised_d2h
  )

site_block_folds <- read_csv(
  file.path(derived_dir, "model_stage2_site_block_folds.csv"),
  show_col_types = FALSE
)

site_loso_folds <- read_csv(
  file.path(derived_dir, "model_stage2_site_loso_folds.csv"),
  show_col_types = FALSE
)

site_fold_summary <- read_csv(
  file.path(derived_dir, "model_stage2_site_fold_summary.csv"),
  show_col_types = FALSE
)

site_unit_summary <- read_csv(
  file.path(derived_dir, "model_stage2_site_unit_summary.csv"),
  show_col_types = FALSE
)

site_lookup <- read_csv(
  file.path(derived_dir, "model_stage2_site_lookup.csv"),
  show_col_types = FALSE
)

site_validation_completeness <- read_csv(
  file.path(derived_dir, "model_stage2_site_validation_completeness.csv"),
  show_col_types = FALSE
)

reference_manifest <- read_csv(
  file.path(derived_dir, "reference_object_manifest.csv"),
  show_col_types = FALSE
)

previous_stage2_rung1 <- read_csv(
  file.path(derived_dir, "model_stage2_rung1_model_comparison.csv"),
  show_col_types = FALSE
)

previous_stage2_rung2 <- read_csv(
  file.path(derived_dir, "model_stage2_rung2_model_comparison.csv"),
  show_col_types = FALSE
)

previous_stage2b_rung2 <- read_csv(
  file.path(derived_dir, "model_stage2b_rung2_staged_model_comparison.csv"),
  show_col_types = FALSE
)

fold_definitions <- bind_rows(site_block_folds, site_loso_folds) %>%
  arrange(validation_layer, split_id, ring)

write_csv(
  fold_definitions,
  file.path(derived_dir, "model_stage5_siteaware_fold_definitions.csv")
)

write_csv(
  site_fold_summary,
  file.path(derived_dir, "model_stage5_siteaware_fold_summary_input.csv")
)

write_csv(
  site_unit_summary,
  file.path(derived_dir, "model_stage5_siteaware_site_unit_summary_input.csv")
)

write_csv(
  site_lookup,
  file.path(derived_dir, "model_stage5_siteaware_site_lookup_input.csv")
)

write_csv(
  site_validation_completeness,
  file.path(derived_dir, "model_stage5_siteaware_validation_completeness_input.csv")
)

write_csv(
  reference_manifest,
  file.path(derived_dir, "model_stage5_siteaware_reference_manifest_input.csv")
)

r1_candidate_c <- list(
  candidate_id = "r1_c_paired_contrast_cn",
  candidate_label = "Rung 1 C: Paired-contrast C/N",
  rung = "Rung 1",
  model_type = "ridge_binomial",
  candidate_group = "accepted_baseline",
  caution_level = "paired_cn_benchmark",
  predictors = c("abs_delta_d13c", "abs_delta_d15n", "paired_distance_cn"),
  distance_vars = c("abs_delta_d13c", "abs_delta_d15n"),
  distance_name = "paired_distance_cn",
  uses_primary_site_reference = FALSE,
  uses_winter_reference = FALSE,
  uses_r1_upstream = FALSE,
  uses_breast_h = FALSE,
  uses_primary_h = FALSE
)

r1_candidate_d <- list(
  candidate_id = "r1_d_structured_paired_cn",
  candidate_label = "Rung 1 D: Structured paired C/N",
  rung = "Rung 1",
  model_type = "ridge_binomial",
  candidate_group = "accepted_comparator",
  caution_level = "paired_cn_benchmark",
  predictors = c("normalised_d13c_Breast", "normalised_d15n_Breast", "abs_delta_d13c", "abs_delta_d15n"),
  distance_vars = character(0),
  distance_name = NA_character_,
  uses_primary_site_reference = FALSE,
  uses_winter_reference = FALSE,
  uses_r1_upstream = FALSE,
  uses_breast_h = FALSE,
  uses_primary_h = FALSE
)

r1_candidate_x <- list(
  candidate_id = "r1_x_paired_contrast_cn_plus_primary_site_reference",
  candidate_label = "Rung 1 X: Paired-contrast C/N + primary-site reference distance",
  rung = "Rung 1",
  model_type = "ridge_binomial",
  candidate_group = "exploratory_extension",
  caution_level = "training_fold_primary_site_reference_exploratory",
  predictors = c("abs_delta_d13c", "abs_delta_d15n", "paired_distance_cn", "primary_ref_min_site_dist_cn"),
  distance_vars = c("abs_delta_d13c", "abs_delta_d15n"),
  distance_name = "paired_distance_cn",
  uses_primary_site_reference = TRUE,
  uses_winter_reference = FALSE,
  uses_r1_upstream = FALSE,
  uses_breast_h = FALSE,
  uses_primary_h = FALSE
)

r2_candidate_c <- list(
  candidate_id = "r2_c_breast_cnh_plus_paired_contrast",
  candidate_label = "Rung 2 C: Breast C/N/H + paired contrasts",
  rung = "Rung 2",
  model_type = "ridge_multinomial",
  candidate_group = "co_leading_direct",
  caution_level = "primary_h_proxy_variance_high_caution",
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
  uses_primary_site_reference = FALSE,
  uses_winter_reference = FALSE,
  uses_r1_upstream = FALSE,
  uses_breast_h = TRUE,
  uses_primary_h = TRUE
)

r2_candidate_x <- list(
  candidate_id = "r2_x_direct_cnh_plus_paired_contrast_plus_winter_reference",
  candidate_label = "Rung 2 X: Direct C/N/H + paired contrasts + winter reference distances",
  rung = "Rung 2",
  model_type = "ridge_multinomial",
  candidate_group = "exploratory_extension",
  caution_level = "direct_cnh_plus_training_fold_winter_reference",
  predictors = c(
    "normalised_d13c_Breast",
    "normalised_d15n_Breast",
    "normalised_d2h_Breast",
    "abs_delta_d13c",
    "abs_delta_d15n",
    "abs_delta_d2h",
    "winter_ref_dist_nz_cnh",
    "winter_ref_dist_au_cnh"
  ),
  distance_vars = character(0),
  distance_name = NA_character_,
  uses_primary_site_reference = FALSE,
  uses_winter_reference = TRUE,
  uses_r1_upstream = FALSE,
  uses_breast_h = TRUE,
  uses_primary_h = TRUE
)

r1_upstream_candidate <- r1_candidate_c

r2s2_migrant_base <- list(
  predictors = c("normalised_d13c_Breast", "normalised_d15n_Breast", "normalised_d2h_Breast"),
  distance_vars = character(0),
  distance_name = NA_character_
)

r2s2_migrant_winter_ref <- list(
  predictors = c(
    "normalised_d13c_Breast",
    "normalised_d15n_Breast",
    "normalised_d2h_Breast",
    "winter_ref_dist_nz_cnh",
    "winter_ref_dist_au_cnh"
  ),
  distance_vars = character(0),
  distance_name = NA_character_
)

staged_candidate_s2 <- list(
  candidate_id = "r2s_b_soft_hierarchical_r1c_then_breast_cnh",
  candidate_label = "Rung 2 S2: Soft hierarchical R1 C then breast C/N/H NZ-vs-AU",
  rung = "Rung 2",
  model_type = "staged_soft_hierarchical",
  candidate_group = "co_leading_staged",
  caution_level = "staged_r1c_then_breast_cnh_no_primary_h",
  migrant_candidate = r2s2_migrant_base,
  uses_primary_site_reference = FALSE,
  uses_winter_reference = FALSE,
  uses_r1_upstream = TRUE,
  uses_breast_h = TRUE,
  uses_primary_h = FALSE
)

staged_candidate_x <- list(
  candidate_id = "r2x_b_staged_r1c_then_breast_cnh_plus_winter_reference",
  candidate_label = "Rung 2 X2: Soft hierarchical R1 C then breast C/N/H + winter reference distances",
  rung = "Rung 2",
  model_type = "staged_soft_hierarchical",
  candidate_group = "exploratory_extension",
  caution_level = "staged_r1c_then_breast_cnh_plus_training_fold_winter_reference",
  migrant_candidate = r2s2_migrant_winter_ref,
  uses_primary_site_reference = FALSE,
  uses_winter_reference = TRUE,
  uses_r1_upstream = TRUE,
  uses_breast_h = TRUE,
  uses_primary_h = FALSE
)

candidate_definitions <- bind_rows(
  lapply(
    list(r1_candidate_c, r1_candidate_d, r1_candidate_x, r2_candidate_c, r2_candidate_x, staged_candidate_s2, staged_candidate_x),
    function(candidate) {
      tibble(
        rung = candidate$rung,
        candidate_id = candidate$candidate_id,
        candidate_label = candidate$candidate_label,
        model_type = candidate$model_type,
        candidate_group = candidate$candidate_group,
        caution_level = candidate$caution_level,
        uses_primary_site_reference = candidate$uses_primary_site_reference,
        uses_winter_reference = candidate$uses_winter_reference,
        uses_r1_upstream = candidate$uses_r1_upstream,
        uses_breast_h = candidate$uses_breast_h,
        uses_primary_h = candidate$uses_primary_h,
        predictors = if (!is.null(candidate$predictors)) paste(candidate$predictors, collapse = "|") else "",
        distance_vars = if (!is.null(candidate$distance_vars)) paste(candidate$distance_vars, collapse = "|") else "",
        migrant_submodel_predictors = if (!is.null(candidate$migrant_candidate)) paste(candidate$migrant_candidate$predictors, collapse = "|") else ""
      )
    }
  )
)

write_csv(
  candidate_definitions,
  file.path(derived_dir, "model_stage5_siteaware_candidate_definitions.csv")
)

split_keys <- fold_definitions %>%
  distinct(validation_layer, split_id) %>%
  arrange(validation_layer, split_id)

split_metric_rows <- list()
prediction_rows <- list()
row_idx <- 1L

for (split_idx in seq_len(nrow(split_keys))) {
  split_meta <- split_keys[split_idx, ]

  split_fold_rows <- fold_definitions %>%
    filter(
      validation_layer == split_meta$validation_layer[[1]],
      split_id == split_meta$split_id[[1]]
    )

  test_ids <- split_fold_rows$ring
  assessment_sites <- unique(split_fold_rows$sampling_site_label)

  train_data <- stage_data %>% filter(!(ring %in% test_ids))
  test_data <- stage_data %>% filter(ring %in% test_ids)

  primary_reference_train <- live_primary_reference_all %>%
    filter(!(sampling_site_label %in% assessment_sites))

  train_aug <- train_data %>%
    add_primary_site_reference_features(primary_reference_train) %>%
    add_winter_region_reference_features(museum_reference_all)

  test_aug <- test_data %>%
    add_primary_site_reference_features(primary_reference_train) %>%
    add_winter_region_reference_features(museum_reference_all)

  split_context <- site_fold_summary %>%
    filter(
      validation_layer == split_meta$validation_layer[[1]],
      split_id == split_meta$split_id[[1]]
    ) %>%
    transmute(
      validation_layer,
      split_id,
      n_assessment,
      n_sites_assessment,
      assessment_sites,
      fold_design_note,
      rung1_all_classes_present,
      rung2_all_classes_present
    )

  if (nrow(split_context) != 1) {
    stop(paste("Could not match site fold summary for", split_meta$validation_layer[[1]], split_meta$split_id[[1]]))
  }

  primary_reference_note <- tibble(
    primary_reference_rows = nrow(primary_reference_train),
    primary_reference_sites = n_distinct(primary_reference_train$sampling_site_label),
    museum_reference_rows_core3 = sum(stats::complete.cases(
      museum_reference_all[, c("normalised_d13c_Breast", "normalised_d15n_Breast", "normalised_d2h_Breast")]
    ))
  )

  binary_candidates <- list(r1_candidate_c, r1_candidate_d, r1_candidate_x)

  for (candidate in binary_candidates) {
    res <- evaluate_binary_candidate_split(
      train_data = train_aug,
      test_data = test_aug,
      candidate = candidate,
      seed = 910000L + split_idx * 100L + match(candidate$candidate_id, c(r1_candidate_c$candidate_id, r1_candidate_d$candidate_id, r1_candidate_x$candidate_id))
    )

    split_metric_rows[[row_idx]] <- bind_cols(
      tibble(
        rung = candidate$rung,
        candidate_id = candidate$candidate_id,
        candidate_label = candidate$candidate_label,
        model_type = candidate$model_type,
        candidate_group = candidate$candidate_group,
        caution_level = candidate$caution_level,
        uses_primary_site_reference = candidate$uses_primary_site_reference,
        uses_winter_reference = candidate$uses_winter_reference,
        uses_r1_upstream = candidate$uses_r1_upstream,
        uses_breast_h = candidate$uses_breast_h,
        uses_primary_h = candidate$uses_primary_h,
        fit_status = res$status
      ),
      split_context,
      tibble(n_analysis = nrow(train_data)),
      primary_reference_note,
      if (identical(res$status, "ok")) {
        bind_cols(
          res$metrics,
          tibble(predictors_used = res$predictors_used)
        )
      } else {
        tibble(
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
          predictors_used = NA_character_
        )
      }
    )

    if (identical(res$status, "ok")) {
      prediction_rows[[row_idx]] <- bind_cols(
        tibble(
          rung = candidate$rung,
          candidate_id = candidate$candidate_id,
          candidate_label = candidate$candidate_label,
          model_type = candidate$model_type,
          candidate_group = candidate$candidate_group,
          validation_layer = split_context$validation_layer,
          split_id = split_context$split_id
        ),
        res$predictions
      )
    }

    row_idx <- row_idx + 1L
  }

  direct_candidates <- list(r2_candidate_c, r2_candidate_x)

  for (candidate in direct_candidates) {
    res <- evaluate_multiclass_candidate_split(
      train_data = train_aug,
      test_data = test_aug,
      candidate = candidate,
      seed = 920000L + split_idx * 100L + match(candidate$candidate_id, c(r2_candidate_c$candidate_id, r2_candidate_x$candidate_id)),
      class_levels = class_levels_r2
    )

    split_metric_rows[[row_idx]] <- bind_cols(
      tibble(
        rung = candidate$rung,
        candidate_id = candidate$candidate_id,
        candidate_label = candidate$candidate_label,
        model_type = candidate$model_type,
        candidate_group = candidate$candidate_group,
        caution_level = candidate$caution_level,
        uses_primary_site_reference = candidate$uses_primary_site_reference,
        uses_winter_reference = candidate$uses_winter_reference,
        uses_r1_upstream = candidate$uses_r1_upstream,
        uses_breast_h = candidate$uses_breast_h,
        uses_primary_h = candidate$uses_primary_h,
        fit_status = res$status
      ),
      split_context,
      tibble(n_analysis = nrow(train_data)),
      primary_reference_note,
      if (identical(res$status, "ok")) {
        bind_cols(
          res$metrics,
          tibble(predictors_used = res$predictors_used)
        )
      } else {
        tibble(
          log_loss = NA_real_,
          balanced_accuracy = NA_real_,
          accuracy = NA_real_,
          brier = NA_real_,
          selected_lambda = NA_real_,
          ece = NA_real_,
          mean_top_prob = NA_real_,
          weak_support_rate = NA_real_,
          predictors_used = NA_character_
        )
      }
    )

    if (identical(res$status, "ok")) {
      prediction_rows[[row_idx]] <- bind_cols(
        tibble(
          rung = candidate$rung,
          candidate_id = candidate$candidate_id,
          candidate_label = candidate$candidate_label,
          model_type = candidate$model_type,
          candidate_group = candidate$candidate_group,
          validation_layer = split_context$validation_layer,
          split_id = split_context$split_id
        ),
        res$predictions
      )
    }

    row_idx <- row_idx + 1L
  }

  staged_candidates <- list(staged_candidate_s2, staged_candidate_x)

  for (candidate in staged_candidates) {
    res <- evaluate_staged_candidate_split(
      train_data = train_aug,
      test_data = test_aug,
      r1_candidate = r1_upstream_candidate,
      migrant_candidate = candidate$migrant_candidate,
      seed_r1 = 930000L + split_idx * 100L + ifelse(candidate$candidate_id == staged_candidate_s2$candidate_id, 1L, 2L),
      seed_sub = 940000L + split_idx * 100L + ifelse(candidate$candidate_id == staged_candidate_s2$candidate_id, 1L, 2L),
      class_levels = class_levels_r2
    )

    split_metric_rows[[row_idx]] <- bind_cols(
      tibble(
        rung = candidate$rung,
        candidate_id = candidate$candidate_id,
        candidate_label = candidate$candidate_label,
        model_type = candidate$model_type,
        candidate_group = candidate$candidate_group,
        caution_level = candidate$caution_level,
        uses_primary_site_reference = candidate$uses_primary_site_reference,
        uses_winter_reference = candidate$uses_winter_reference,
        uses_r1_upstream = candidate$uses_r1_upstream,
        uses_breast_h = candidate$uses_breast_h,
        uses_primary_h = candidate$uses_primary_h,
        fit_status = res$status
      ),
      split_context,
      tibble(n_analysis = nrow(train_data)),
      primary_reference_note,
      if (identical(res$status, "ok")) {
        bind_cols(
          res$metrics,
          tibble(predictors_used = res$predictors_used)
        )
      } else {
        tibble(
          log_loss = NA_real_,
          balanced_accuracy = NA_real_,
          accuracy = NA_real_,
          brier = NA_real_,
          selected_lambda = NA_real_,
          ece = NA_real_,
          mean_top_prob = NA_real_,
          weak_support_rate = NA_real_,
          predictors_used = NA_character_
        )
      }
    )

    if (identical(res$status, "ok")) {
      prediction_rows[[row_idx]] <- bind_cols(
        tibble(
          rung = candidate$rung,
          candidate_id = candidate$candidate_id,
          candidate_label = candidate$candidate_label,
          model_type = candidate$model_type,
          candidate_group = candidate$candidate_group,
          validation_layer = split_context$validation_layer,
          split_id = split_context$split_id
        ),
        res$predictions
      )
    }

    row_idx <- row_idx + 1L
  }
}

split_metrics <- bind_rows(split_metric_rows) %>%
  arrange(rung, validation_layer, candidate_id, split_id)

assessment_predictions <- bind_rows(prediction_rows) %>%
  arrange(rung, validation_layer, candidate_id, split_id, ring)

write_csv(
  split_metrics,
  file.path(derived_dir, "model_stage5_siteaware_split_metrics.csv")
)

write_csv(
  assessment_predictions,
  file.path(derived_dir, "model_stage5_siteaware_assessment_predictions.csv")
)

pooled_summary <- bind_rows(
  assessment_predictions %>%
    filter(rung == "Rung 1") %>%
    group_by(rung, candidate_id, candidate_label, model_type, candidate_group, validation_layer) %>%
    group_modify(~ {
      binary_pooled_summary(.x)
    }) %>%
    ungroup(),
  assessment_predictions %>%
    filter(rung == "Rung 2") %>%
    group_by(rung, candidate_id, candidate_label, model_type, candidate_group, validation_layer) %>%
    group_modify(~ {
      multiclass_pooled_summary(.x, class_levels = class_levels_r2)
    }) %>%
    ungroup()
)

model_summary <- split_metrics %>%
  group_by(
    rung,
    candidate_id,
    candidate_label,
    model_type,
    candidate_group,
    caution_level,
    uses_primary_site_reference,
    uses_winter_reference,
    uses_r1_upstream,
    uses_breast_h,
    uses_primary_h,
    validation_layer
  ) %>%
  summarise(
    n_splits = n(),
    n_successful_splits = sum(fit_status == "ok", na.rm = TRUE),
    mean_split_log_loss = mean(log_loss, na.rm = TRUE),
    sd_split_log_loss = stats::sd(log_loss, na.rm = TRUE),
    mean_split_bal_accuracy = mean(balanced_accuracy, na.rm = TRUE),
    sd_split_bal_accuracy = stats::sd(balanced_accuracy, na.rm = TRUE),
    mean_split_accuracy = mean(accuracy, na.rm = TRUE),
    mean_split_brier = mean(brier, na.rm = TRUE),
    mean_split_ece = mean(ece, na.rm = TRUE),
    mean_split_top_prob = mean(mean_top_prob, na.rm = TRUE),
    mean_split_weak_support_rate = mean(weak_support_rate, na.rm = TRUE),
    mean_primary_reference_sites = mean(primary_reference_sites, na.rm = TRUE),
    mean_primary_reference_rows = mean(primary_reference_rows, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    pooled_summary,
    by = c("rung", "candidate_id", "candidate_label", "model_type", "candidate_group", "validation_layer")
  ) %>%
  group_by(rung, validation_layer) %>%
  mutate(
    pooled_log_loss_rank = min_rank(pooled_log_loss),
    pooled_bal_accuracy_rank = min_rank(desc(pooled_bal_accuracy))
  ) %>%
  ungroup() %>%
  arrange(rung, validation_layer, pooled_log_loss_rank, pooled_bal_accuracy_rank)

write_csv(
  model_summary,
  file.path(derived_dir, "model_stage5_siteaware_model_summary.csv")
)

confusion_summary <- bind_rows(
  assessment_predictions %>%
    filter(rung == "Rung 1") %>%
    transmute(rung, validation_layer, candidate_id, candidate_label, truth = as.character(truth), predicted = as.character(.pred_class)),
  assessment_predictions %>%
    filter(rung == "Rung 2") %>%
    transmute(rung, validation_layer, candidate_id, candidate_label, truth = as.character(truth), predicted = as.character(.pred_class))
) %>%
  count(rung, validation_layer, candidate_id, candidate_label, truth, predicted, name = "n")

write_csv(
  confusion_summary,
  file.path(derived_dir, "model_stage5_siteaware_confusion_summary.csv")
)

calibration_summary <- bind_rows(
  assessment_predictions %>%
    filter(rung == "Rung 1") %>%
    group_by(rung, candidate_id, candidate_label, validation_layer) %>%
    group_modify(~ {
      binary_calibration_summary(
        predictions = .x %>%
          transmute(
            truth = factor(truth, levels = c("resident", "migrant")),
            .pred_positive = .pred_positive
          ),
        positive_class = "migrant"
      )
    }) %>%
    ungroup(),
  assessment_predictions %>%
    filter(rung == "Rung 2") %>%
    group_by(rung, candidate_id, candidate_label, validation_layer) %>%
    group_modify(~ {
      multiclass_calibration_summary(
        predictions = .x %>%
          transmute(
            truth = factor(truth, levels = class_levels_r2),
            resident = resident,
            `NZ migrant` = `NZ migrant`,
            `AU migrant` = `AU migrant`
          ),
        class_levels = class_levels_r2
      )
    }) %>%
    ungroup()
)

write_csv(
  calibration_summary,
  file.path(derived_dir, "model_stage5_siteaware_calibration_summary.csv")
)

previous_baseline_summary <- bind_rows(
  previous_stage2_rung1 %>%
    filter(candidate_id %in% c("r1_c_paired_contrast_cn", "r1_d_structured_paired_cn"),
           validation_layer %in% c("blocked_coordinate_cv", "repeated_stratified_cv")) %>%
    transmute(
      rung = "Rung 1",
      candidate_id,
      previous_validation_layer = validation_layer,
      previous_mean_log_loss = mean_log_loss,
      previous_mean_bal_accuracy = mean_bal_accuracy,
      previous_mean_brier = mean_brier
    ),
  previous_stage2_rung2 %>%
    filter(candidate_id == "r2_c_breast_cnh_plus_paired_contrast",
           validation_layer %in% c("blocked_coordinate_cv", "repeated_stratified_cv")) %>%
    transmute(
      rung = "Rung 2",
      candidate_id,
      previous_validation_layer = validation_layer,
      previous_mean_log_loss = mean_log_loss,
      previous_mean_bal_accuracy = mean_bal_accuracy,
      previous_mean_brier = mean_brier
    ),
  previous_stage2b_rung2 %>%
    filter(candidate_id == "r2s_b_soft_hierarchical_r1c_then_breast_cnh",
           validation_layer %in% c("blocked_coordinate_cv", "repeated_stratified_cv")) %>%
    transmute(
      rung = "Rung 2",
      candidate_id,
      previous_validation_layer = validation_layer,
      previous_mean_log_loss = mean_log_loss,
      previous_mean_bal_accuracy = mean_bal_accuracy,
      previous_mean_brier = mean_brier
    )
)

vs_previous <- model_summary %>%
  filter(candidate_id %in% c(
    "r1_c_paired_contrast_cn",
    "r1_d_structured_paired_cn",
    "r2_c_breast_cnh_plus_paired_contrast",
    "r2s_b_soft_hierarchical_r1c_then_breast_cnh"
  )) %>%
  transmute(
    rung,
    candidate_id,
    candidate_label,
    siteaware_validation_layer = validation_layer,
    pooled_log_loss,
    pooled_bal_accuracy,
    pooled_brier,
    pooled_ece,
    pooled_mean_top_prob,
    pooled_weak_support_rate
  ) %>%
  crossing(
    tibble(previous_validation_layer = c("blocked_coordinate_cv", "repeated_stratified_cv"))
  ) %>%
  left_join(
    previous_baseline_summary,
    by = c("rung", "candidate_id", "previous_validation_layer")
  ) %>%
  mutate(
    delta_log_loss_vs_previous = pooled_log_loss - previous_mean_log_loss,
    delta_bal_accuracy_vs_previous = pooled_bal_accuracy - previous_mean_bal_accuracy,
    delta_brier_vs_previous = pooled_brier - previous_mean_brier
  )

write_csv(
  vs_previous,
  file.path(derived_dir, "model_stage5_siteaware_vs_previous_validation.csv")
)

blocked_r1 <- model_summary %>%
  filter(rung == "Rung 1", validation_layer == "blocked_site_cv")

blocked_r2 <- model_summary %>%
  filter(rung == "Rung 2", validation_layer == "blocked_site_cv")

loso_r2 <- model_summary %>%
  filter(rung == "Rung 2", validation_layer == "leave_one_site_out_cv")

summary_lines <- c(
  "# Stage 5 Site-Aware Validation and Reference Extension Summary",
  "",
  "This extension branch keeps all previously accepted results intact and adds a separate site-aware validation pass.",
  "",
  "## Validation layers",
  "",
  "- Primary new evaluation layer: grouped 4-fold site-block validation (`sampling_site_label`).",
  "- Secondary stress-test layer: leave-one-site-out validation.",
  "- Formal site-level partial-pooling models were not fitted here because the held-out-site validation target and the small number of labelled birds per site make that extension unstable and difficult to interpret.",
  "",
  "## Grouped site-block highlights",
  "",
  paste0(
    "- Best grouped-site Rung 1 pooled log loss: ",
    blocked_r1$candidate_id[[which.min(blocked_r1$pooled_log_loss)]],
    " = ",
    format_num(min(blocked_r1$pooled_log_loss, na.rm = TRUE)),
    "."
  ),
  paste0(
    "- Best grouped-site Rung 2 pooled log loss: ",
    blocked_r2$candidate_id[[which.min(blocked_r2$pooled_log_loss)]],
    " = ",
    format_num(min(blocked_r2$pooled_log_loss, na.rm = TRUE)),
    "."
  ),
  paste0(
    "- Best grouped-site Rung 2 pooled balanced accuracy: ",
    blocked_r2$candidate_id[[which.max(blocked_r2$pooled_bal_accuracy)]],
    " = ",
    format_num(max(blocked_r2$pooled_bal_accuracy, na.rm = TRUE)),
    "."
  ),
  "",
  "## LOSO stress-test highlights",
  "",
  "- LOSO should be read cautiously because most held-out sites do not contain all Rung 2 classes.",
  paste0(
    "- Best LOSO Rung 2 pooled log loss: ",
    loso_r2$candidate_id[[which.min(loso_r2$pooled_log_loss)]],
    " = ",
    format_num(min(loso_r2$pooled_log_loss, na.rm = TRUE)),
    "."
  ),
  "",
  "## Extension scope",
  "",
  "- `R1 X` adds a training-fold-only primary-site reference distance built from live primary feathers. This is where the unknown-status live birds contribute: they help characterize breeding-site primary structure without being given phenotype labels.",
  "- `R2 X` adds museum winter breast reference distances to the direct Rung 2 model.",
  "- `R2 X2` adds the same museum winter reference distances to the staged `R2 S2` migrant submodel.",
  "",
  "See the CSV outputs in `data/derived/` for pooled summaries, split metrics, confusion tables, and the comparison against previous validation layers."
)

writeLines(
  summary_lines,
  file.path(derived_dir, "model_stage5_siteaware_summary.md")
)
