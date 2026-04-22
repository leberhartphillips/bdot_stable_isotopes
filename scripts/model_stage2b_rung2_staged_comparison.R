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
source(file.path("R", "model_stage2b_rung2_staged_helpers.R"))

derived_dir <- file.path("data", "derived")
ensure_dir(derived_dir)

format_num <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits))
}

class_levels <- c("resident", "NZ migrant", "AU migrant")

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

r1_candidate <- list(
  candidate_id = "r1_c_paired_contrast_cn",
  candidate_label = "Rung 1 C: Paired-contrast C/N",
  predictors = c("abs_delta_d13c", "abs_delta_d15n", "paired_distance_cn"),
  distance_vars = c("abs_delta_d13c", "abs_delta_d15n"),
  distance_name = "paired_distance_cn"
)

staged_candidates <- tribble(
  ~candidate_id, ~candidate_label, ~strategy_family, ~caution_level, ~uses_r1_upstream, ~uses_breast_h, ~uses_primary_h, ~uses_hard_or_soft_staging,
  "r2s_a_breast_cnh_plus_r1c_prob", "Rung 2 S1: Breast C/N/H + out-of-fold R1 migrant probability", "probability_augmented_direct", "staged_r1c_prob_plus_breast_cnh_no_primary_h", TRUE, TRUE, FALSE, "stacked_direct",
  "r2s_b_soft_hierarchical_r1c_then_breast_cnh", "Rung 2 S2: Soft hierarchical R1 C then breast C/N/H NZ-vs-AU", "soft_hierarchical_two_step", "staged_r1c_then_breast_cnh_no_primary_h", TRUE, TRUE, FALSE, "soft_hierarchical"
)

write_csv(
  staged_candidates,
  file.path(derived_dir, "model_stage2b_rung2_staged_candidate_definitions.csv")
)

repeated_folds <- read_csv(
  file.path(derived_dir, "model_stage2_rung2_repeated_cv_folds.csv"),
  show_col_types = FALSE
)

blocked_fold_file <- read_csv(
  file.path(derived_dir, "model_stage2_coordinate_block_folds.csv"),
  show_col_types = FALSE
)

blocked_folds <- blocked_fold_file %>%
  transmute(
    rung = "rung2",
    validation_layer = "blocked_coordinate_cv",
    split_id,
    repeat_id = NA_character_,
    fold_id = split_id,
    ring,
    coord_group,
    outcome = status_rung2
  )

all_fold_definitions <- bind_rows(repeated_folds, blocked_folds)

write_csv(
  all_fold_definitions,
  file.path(derived_dir, "model_stage2b_rung2_staged_outer_folds.csv")
)

split_keys <- all_fold_definitions %>%
  distinct(validation_layer, split_id, repeat_id, fold_id) %>%
  arrange(validation_layer, split_id)

augmented_candidate <- list(
  predictors = c(
    "normalised_d13c_Breast",
    "normalised_d15n_Breast",
    "normalised_d2h_Breast",
    "r1_prob_migrant"
  ),
  distance_vars = character(0),
  distance_name = NA_character_
)

migrant_submodel_candidate <- list(
  predictors = c(
    "normalised_d13c_Breast",
    "normalised_d15n_Breast",
    "normalised_d2h_Breast"
  ),
  distance_vars = character(0),
  distance_name = NA_character_
)

resample_rows <- vector("list", nrow(split_keys) * nrow(staged_candidates))
prediction_rows <- vector("list", nrow(split_keys) * nrow(staged_candidates))
row_idx <- 1L

for (split_idx in seq_len(nrow(split_keys))) {
  split_meta <- split_keys[split_idx, ]

  test_ids <- all_fold_definitions %>%
    filter(
      validation_layer == split_meta$validation_layer[[1]],
      split_id == split_meta$split_id[[1]]
    ) %>%
    pull(ring)

  train_data <- stage_data %>% filter(!(ring %in% test_ids))
  test_data <- stage_data %>% filter(ring %in% test_ids)

  r1_train_oof <- crossfit_stage2b_r1_probabilities(
    train_data = train_data,
    r1_candidate = r1_candidate,
    seed = 200000L + split_idx * 100L
  )

  r1_full_fit <- fit_stage3_ridge_binomial(
    train_data = train_data,
    candidate = r1_candidate,
    outcome_col = "status_rung1",
    seed = 300000L + split_idx
  )

  if (!identical(r1_full_fit$status, "ok")) {
    stop(paste("Outer R1 fit failed for split", split_meta$split_id[[1]], "with status", r1_full_fit$status))
  }

  r1_test_pred <- predict_stage3_ridge_binomial(r1_full_fit, test_data)

  train_aug <- train_data %>%
    left_join(r1_train_oof, by = "ring")

  if (any(is.na(train_aug$r1_prob_migrant))) {
    stop(paste("Missing out-of-fold R1 probabilities in staged Rung 2 training data for split", split_meta$split_id[[1]]))
  }

  test_aug <- test_data %>%
    mutate(r1_prob_migrant = r1_test_pred$.pred_migrant)

  fit_s1 <- fit_predict_ridge(
    train_data = train_aug,
    test_data = test_aug,
    candidate = augmented_candidate,
    outcome_col = "status_rung2",
    family = "multinomial",
    seed = 400000L + split_idx
  )

  if (!identical(fit_s1$status, "ok")) {
    stop(paste("Staged candidate S1 failed for split", split_meta$split_id[[1]], "with status", fit_s1$status))
  }

  pred_s1 <- fit_s1$predictions %>%
    mutate(
      top_prob = pmax(`resident`, `NZ migrant`, `AU migrant`),
      weak_support = top_prob < 0.5,
      r1_prob_migrant = r1_test_pred$.pred_migrant,
      p_au_given_migrant = NA_real_
    )

  metrics_s1 <- stage2b_prediction_metrics(pred_s1, class_levels)

  migrant_train <- train_data %>% filter(status_rung2 != "resident")
  fit_s2_sub <- fit_stage2b_ridge_binomial_generic(
    train_data = migrant_train,
    candidate = migrant_submodel_candidate,
    outcome_col = "status_nzau",
    seed = 500000L + split_idx
  )

  if (!identical(fit_s2_sub$status, "ok")) {
    stop(paste("Staged candidate S2 migrant submodel failed for split", split_meta$split_id[[1]], "with status", fit_s2_sub$status))
  }

  pred_s2_sub <- predict_stage2b_ridge_binomial_generic(fit_s2_sub, test_data)
  prob_s2 <- build_stage2b_hierarchical_probabilities(
    p_migrant = r1_test_pred$.pred_migrant,
    p_au_given_migrant = pred_s2_sub$prob_positive
  )

  pred_s2 <- stage2b_prediction_table(test_data, prob_s2) %>%
    mutate(
      r1_prob_migrant = r1_test_pred$.pred_migrant,
      p_au_given_migrant = pred_s2_sub$prob_positive
    )

  metrics_s2 <- stage2b_prediction_metrics(pred_s2, class_levels)

  split_info <- tibble(
    validation_layer = split_meta$validation_layer[[1]],
    split_id = split_meta$split_id[[1]],
    repeat_id = split_meta$repeat_id[[1]],
    fold_id = split_meta$fold_id[[1]],
    n_analysis = nrow(train_data),
    n_assessment = nrow(test_data)
  )

  resample_rows[[row_idx]] <- bind_cols(
    staged_candidates %>% slice(1),
    split_info,
    metrics_s1
  )

  prediction_rows[[row_idx]] <- bind_cols(
    staged_candidates %>% slice(1),
    split_info %>% select(validation_layer, split_id, repeat_id, fold_id),
    pred_s1
  )

  row_idx <- row_idx + 1L

  resample_rows[[row_idx]] <- bind_cols(
    staged_candidates %>% slice(2),
    split_info,
    metrics_s2
  )

  prediction_rows[[row_idx]] <- bind_cols(
    staged_candidates %>% slice(2),
    split_info %>% select(validation_layer, split_id, repeat_id, fold_id),
    pred_s2
  )

  row_idx <- row_idx + 1L
}

resample_tbl <- bind_rows(resample_rows) %>%
  arrange(validation_layer, candidate_id, split_id)

prediction_tbl <- bind_rows(prediction_rows) %>%
  arrange(validation_layer, candidate_id, split_id, ring)

write_csv(
  resample_tbl,
  file.path(derived_dir, "model_stage2b_rung2_staged_resample_metrics.csv")
)

write_csv(
  prediction_tbl,
  file.path(derived_dir, "model_stage2b_rung2_staged_assessment_predictions.csv")
)

staged_model_comparison <- resample_tbl %>%
  group_by(
    candidate_id,
    candidate_label,
    strategy_family,
    caution_level,
    uses_r1_upstream,
    uses_breast_h,
    uses_primary_h,
    uses_hard_or_soft_staging,
    validation_layer
  ) %>%
  summarise(
    n_successful_splits = n(),
    mean_log_loss = mean(log_loss, na.rm = TRUE),
    sd_log_loss = stats::sd(log_loss, na.rm = TRUE),
    mean_bal_accuracy = mean(balanced_accuracy, na.rm = TRUE),
    sd_bal_accuracy = stats::sd(balanced_accuracy, na.rm = TRUE),
    mean_accuracy = mean(accuracy, na.rm = TRUE),
    mean_brier = mean(brier, na.rm = TRUE),
    mean_top_prob = mean(mean_top_prob, na.rm = TRUE),
    mean_weak_support_rate = mean(weak_support_rate, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(validation_layer) %>%
  mutate(
    log_loss_rank = rank(mean_log_loss, ties.method = "min"),
    bal_accuracy_rank = rank(-mean_bal_accuracy, ties.method = "min")
  ) %>%
  ungroup() %>%
  arrange(validation_layer, log_loss_rank, bal_accuracy_rank)

write_csv(
  staged_model_comparison,
  file.path(derived_dir, "model_stage2b_rung2_staged_model_comparison.csv")
)

staged_calibration <- stage2b_calibration_summary(
  prediction_tbl = prediction_tbl,
  class_levels = class_levels
)

write_csv(
  staged_calibration,
  file.path(derived_dir, "model_stage2b_rung2_staged_calibration_summary.csv")
)

staged_support_summary <- prediction_tbl %>%
  group_by(candidate_id, candidate_label, validation_layer, strategy_family, caution_level) %>%
  summarise(
    n_assessment = n(),
    mean_top_prob = mean(top_prob, na.rm = TRUE),
    weak_support_rate = mean(weak_support, na.rm = TRUE),
    mean_r1_prob_migrant = mean(r1_prob_migrant, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(validation_layer, candidate_id)

write_csv(
  staged_support_summary,
  file.path(derived_dir, "model_stage2b_rung2_staged_support_summary.csv")
)

direct_comparison <- read_csv(
  file.path(derived_dir, "model_stage2_rung2_model_comparison.csv"),
  show_col_types = FALSE
) %>%
  filter(candidate_id %in% c("r2_a_breast_cn", "r2_c_breast_cnh_plus_paired_contrast")) %>%
  mutate(
    strategy_family = "accepted_direct",
    uses_r1_upstream = FALSE,
    uses_hard_or_soft_staging = "none"
  )

direct_support <- read_csv(
  file.path(derived_dir, "model_stage2_rung2_assessment_predictions.csv"),
  show_col_types = FALSE
) %>%
  filter(candidate_id %in% c("r2_a_breast_cn", "r2_c_breast_cnh_plus_paired_contrast")) %>%
  mutate(
    top_prob = pmax(`resident`, `NZ migrant`, `AU migrant`),
    weak_support = top_prob < 0.5
  ) %>%
  group_by(candidate_id, candidate_label, validation_layer) %>%
  summarise(
    mean_top_prob = mean(top_prob, na.rm = TRUE),
    mean_weak_support_rate = mean(weak_support, na.rm = TRUE),
    .groups = "drop"
  )

comparison_against_direct <- bind_rows(
  direct_comparison %>%
    select(
      candidate_id,
      candidate_label,
      validation_layer,
      strategy_family,
      caution_level,
      uses_r1_upstream,
      uses_breast_h,
      uses_primary_h,
      uses_hard_or_soft_staging,
      mean_log_loss,
      mean_bal_accuracy,
      mean_accuracy,
      mean_brier
    ) %>%
    left_join(
      direct_support,
      by = c("candidate_id", "candidate_label", "validation_layer")
    ),
  staged_model_comparison %>%
    select(
      candidate_id,
      candidate_label,
      validation_layer,
      strategy_family,
      caution_level,
      uses_r1_upstream,
      uses_breast_h,
      uses_primary_h,
      uses_hard_or_soft_staging,
      mean_log_loss,
      mean_bal_accuracy,
      mean_accuracy,
      mean_brier,
      mean_top_prob,
      mean_weak_support_rate
    )
) %>%
  group_by(validation_layer) %>%
  mutate(
    combined_log_loss_rank = rank(mean_log_loss, ties.method = "min"),
    combined_bal_accuracy_rank = rank(-mean_bal_accuracy, ties.method = "min")
  ) %>%
  ungroup() %>%
  arrange(validation_layer, combined_log_loss_rank, combined_bal_accuracy_rank)

write_csv(
  comparison_against_direct,
  file.path(derived_dir, "model_stage2b_rung2_staged_vs_direct_comparison.csv")
)

repeated_best <- comparison_against_direct %>%
  filter(validation_layer == "repeated_stratified_cv") %>%
  slice_min(order_by = combined_log_loss_rank, n = 1)

blocked_best <- comparison_against_direct %>%
  filter(validation_layer == "blocked_coordinate_cv") %>%
  slice_min(order_by = combined_log_loss_rank, n = 1)

summary_lines <- c(
  "# Staged Rung 2 Summary",
  "",
  "## Leakage control",
  "",
  "- The probability-augmented direct model used cross-fitted Rung 1 migrant probabilities inside each outer training split and a separate outer-trained Rung 1 fit for each assessment split.",
  "- The soft hierarchical model fit the accepted Rung 1 model and the NZ-vs-AU migrant submodel on the outer analysis set only; no assessment rows contributed to either fit.",
  "",
  "## Staged candidates",
  "",
  paste0(
    "- ",
    staged_candidates$candidate_label,
    ": strategy family `",
    staged_candidates$strategy_family,
    "`, caution `",
    staged_candidates$caution_level,
    "`."
  ),
  "",
  "## Headline comparison versus accepted direct Rung 2 models",
  "",
  paste0(
    "- Best repeated-CV log-loss rank: ",
    repeated_best$candidate_label[[1]],
    " (log loss ",
    format_num(repeated_best$mean_log_loss[[1]]),
    ", balanced accuracy ",
    format_num(repeated_best$mean_bal_accuracy[[1]]),
    ")."
  ),
  paste0(
    "- Best blocked-CV log-loss rank: ",
    blocked_best$candidate_label[[1]],
    " (log loss ",
    format_num(blocked_best$mean_log_loss[[1]]),
    ", balanced accuracy ",
    format_num(blocked_best$mean_bal_accuracy[[1]]),
    ")."
  ),
  "",
  "## Interpretation notes",
  "",
  "- These staged approaches were evaluated as additional Stage 2-style candidates only. They do not replace the accepted Stage 4 uncertainty-aware Rung 2 results.",
  "- Deterministic weak-support rates are reported here from top-class probabilities < 0.5, but Stage 4-style uncertainty-aware indeterminate behavior was not re-estimated for these staged variants.",
  "- The soft hierarchical model avoids primary H entirely; any gain from that variant should therefore be interpreted as improved biological alignment rather than as proof that the accepted direct R2 C model should be displaced."
)

writeLines(
  summary_lines,
  file.path(derived_dir, "model_stage2b_rung2_staged_summary.md")
)
