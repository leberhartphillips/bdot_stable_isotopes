#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(rsample)
  library(tidyr)
  library(tibble)
  library(stringr)
})

source(file.path("R", "model_stage_repeatability_helpers.R"))
source(file.path("R", "model_stage2_helpers.R"))

derived_dir <- file.path("data", "derived")
ensure_dir(derived_dir)

format_num <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits))
}

build_repeated_cv_fold_file <- function(rs_obj, outcome_col, rung_name, layer_name) {
  rs_tbl <- as_tibble(rs_obj) %>%
    mutate(
      repeat_id = if ("id" %in% colnames(.)) id else NA_character_,
      fold_id = if ("id2" %in% colnames(.)) id2 else id,
      split_id = if ("id2" %in% colnames(.)) {
        paste(id, id2, sep = "__")
      } else {
        as.character(id)
      }
    )

  bind_rows(
    lapply(
      seq_len(nrow(rs_tbl)),
      function(i) {
        assessment_dat <- assessment(rs_tbl$splits[[i]])
        tibble(
          rung = rung_name,
          validation_layer = layer_name,
          split_id = rs_tbl$split_id[[i]],
          repeat_id = rs_tbl$repeat_id[[i]],
          fold_id = rs_tbl$fold_id[[i]],
          ring = assessment_dat$ring,
          coord_group = assessment_dat$coord_group,
          outcome = as.character(assessment_dat[[outcome_col]])
        )
      }
    )
  )
}

evaluate_candidates <- function(data, candidates, fold_definitions, outcome_col, family, rung_name) {
  split_ids <- unique(fold_definitions$split_id)

  metric_rows <- vector("list", length(split_ids) * nrow(candidates))
  prediction_rows <- vector("list", length(split_ids) * nrow(candidates))
  row_index <- 1L

  for (split_idx in seq_along(split_ids)) {
    split_id <- split_ids[[split_idx]]
    test_ids <- fold_definitions %>%
      filter(split_id == !!split_id) %>%
      pull(ring)

    train_data <- data %>% filter(!(ring %in% test_ids))
    test_data <- data %>% filter(ring %in% test_ids)

    for (candidate_idx in seq_len(nrow(candidates))) {
      candidate <- candidates[candidate_idx, ]
      fit <- fit_predict_ridge(
        train_data = train_data,
        test_data = test_data,
        candidate = list(
          predictors = candidate$predictors[[1]],
          distance_vars = candidate$distance_vars[[1]],
          distance_name = candidate$distance_name[[1]]
        ),
        outcome_col = outcome_col,
        family = family,
        seed = 1000L + split_idx * 100L + candidate_idx
      )

      metric_rows[[row_index]] <- bind_cols(
        candidate %>%
          select(
            rung,
            candidate_id,
            candidate_label,
            model_family,
            caution_level,
            preferred_default,
            uses_breast_h,
            uses_primary_h,
            uses_paired_distance
          ),
        fold_definitions %>%
          filter(split_id == !!split_id) %>%
          summarise(
            validation_layer = first(validation_layer),
            split_id = first(split_id),
            repeat_id = first(repeat_id),
            fold_id = first(fold_id),
            n_assessment = n()
          ),
        tibble(
          n_analysis = nrow(train_data),
          fit_status = fit$status
        )
      )

      if (identical(fit$status, "ok")) {
        metric_rows[[row_index]] <- bind_cols(metric_rows[[row_index]], fit$metrics) %>%
          mutate(
            predictors_used = paste(fit$predictors_used, collapse = "|")
          )

        prediction_rows[[row_index]] <- bind_cols(
          candidate %>%
            select(
              rung,
              candidate_id,
              candidate_label,
              model_family,
              caution_level,
              preferred_default,
              uses_breast_h,
              uses_primary_h,
              uses_paired_distance
            ),
          fold_definitions %>%
            filter(split_id == !!split_id) %>%
            summarise(
              validation_layer = first(validation_layer),
              split_id = first(split_id),
              repeat_id = first(repeat_id),
              fold_id = first(fold_id)
            ),
          fit$predictions
        )
      } else {
        metric_rows[[row_index]] <- metric_rows[[row_index]] %>%
          mutate(
            log_loss = NA_real_,
            balanced_accuracy = NA_real_,
            accuracy = NA_real_,
            brier = NA_real_,
            selected_lambda = NA_real_,
            predictors_used = NA_character_
          )
      }

      row_index <- row_index + 1L
    }
  }

  metric_tbl <- bind_rows(metric_rows)
  prediction_tbl <- bind_rows(prediction_rows)

  list(metrics = metric_tbl, predictions = prediction_tbl)
}

summarise_model_comparison <- function(metric_tbl) {
  metric_tbl %>%
    group_by(
      rung,
      validation_layer,
      candidate_id,
      candidate_label,
      model_family,
      caution_level,
      preferred_default,
      uses_breast_h,
      uses_primary_h,
      uses_paired_distance
    ) %>%
    summarise(
      n_successful_splits = sum(fit_status == "ok", na.rm = TRUE),
      mean_log_loss = mean(log_loss, na.rm = TRUE),
      sd_log_loss = stats::sd(log_loss, na.rm = TRUE),
      mean_bal_accuracy = mean(balanced_accuracy, na.rm = TRUE),
      sd_bal_accuracy = stats::sd(balanced_accuracy, na.rm = TRUE),
      mean_accuracy = mean(accuracy, na.rm = TRUE),
      mean_brier = mean(brier, na.rm = TRUE),
      median_lambda = stats::median(selected_lambda, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    group_by(rung, validation_layer) %>%
    mutate(
      log_loss_rank = rank(mean_log_loss, ties.method = "min"),
      bal_accuracy_rank = rank(-mean_bal_accuracy, ties.method = "min")
    ) %>%
    ungroup() %>%
    arrange(validation_layer, log_loss_rank, bal_accuracy_rank)
}

summarise_ranking_stability <- function(metric_tbl) {
  split_ranks <- metric_tbl %>%
    filter(fit_status == "ok") %>%
    group_by(rung, validation_layer, split_id) %>%
    mutate(
      split_log_loss_rank = rank(log_loss, ties.method = "average"),
      split_bal_accuracy_rank = rank(-balanced_accuracy, ties.method = "average")
    ) %>%
    ungroup()

  layer_summary <- split_ranks %>%
    group_by(rung, validation_layer, candidate_id, candidate_label) %>%
    summarise(
      mean_split_log_loss_rank = mean(split_log_loss_rank, na.rm = TRUE),
      sd_split_log_loss_rank = stats::sd(split_log_loss_rank, na.rm = TRUE),
      mean_split_bal_accuracy_rank = mean(split_bal_accuracy_rank, na.rm = TRUE),
      sd_split_bal_accuracy_rank = stats::sd(split_bal_accuracy_rank, na.rm = TRUE),
      .groups = "drop"
    )

  primary_tbl <- layer_summary %>%
    filter(validation_layer == "repeated_stratified_cv") %>%
    rename_with(~ paste0("primary_", .x), c(mean_split_log_loss_rank, sd_split_log_loss_rank, mean_split_bal_accuracy_rank, sd_split_bal_accuracy_rank)) %>%
    select(-validation_layer)

  blocked_tbl <- layer_summary %>%
    filter(validation_layer == "blocked_coordinate_cv") %>%
    rename_with(~ paste0("blocked_", .x), c(mean_split_log_loss_rank, sd_split_log_loss_rank, mean_split_bal_accuracy_rank, sd_split_bal_accuracy_rank)) %>%
    select(-validation_layer)

  full_join(primary_tbl, blocked_tbl, by = c("rung", "candidate_id", "candidate_label")) %>%
    mutate(
      log_loss_rank_shift = blocked_mean_split_log_loss_rank - primary_mean_split_log_loss_rank,
      bal_accuracy_rank_shift = blocked_mean_split_bal_accuracy_rank - primary_mean_split_bal_accuracy_rank
    ) %>%
    arrange(rung, primary_mean_split_log_loss_rank)
}

summarise_calibration <- function(prediction_tbl, rung_name, class_levels) {
  grouped <- split(prediction_tbl, interaction(prediction_tbl$candidate_id, prediction_tbl$validation_layer, drop = TRUE))

  bind_rows(
    lapply(
      grouped,
      function(df) {
        base_info <- df %>%
          summarise(
            rung = first(rung),
            candidate_id = first(candidate_id),
            candidate_label = first(candidate_label),
            validation_layer = first(validation_layer)
          )

        if (rung_name == "rung1") {
          cal <- binary_calibration_summary(df, positive_class = "migrant")
        } else {
          cal <- multiclass_calibration_summary(df, class_levels = class_levels)
        }

        bind_cols(base_info[rep(1, nrow(cal)), ], cal)
      }
    )
  )
}

stage2_data <- read_csv(
  file.path(derived_dir, "live_screening_ready_paired_labelled.csv"),
  show_col_types = FALSE
) %>%
  derive_stage2_features()

adequacy <- read_csv(
  file.path(derived_dir, "model_stage1_tissue_summary_adequacy.csv"),
  show_col_types = FALSE
)

breast_h_status <- adequacy %>%
  filter(feather_type == "Breast", isotope == "normalised_d2h") %>%
  pull(adequacy_for_stage2)

primary_h_status <- adequacy %>%
  filter(feather_type == "Primary", isotope == "normalised_d2h") %>%
  pull(adequacy_for_stage2)

if (!identical(breast_h_status, "adequate_with_technical_caution")) {
  stop("Stage 1 breast H adequacy no longer matches the expected Stage 2 constraint.")
}

if (!identical(primary_h_status, "provisionally_adequate_with_high_caution")) {
  stop("Stage 1 primary H adequacy no longer matches the expected Stage 2 constraint.")
}

if (any(!stage2_data$has_complete_cnh_pair)) {
  stop("Stage 2 expected a complete paired C/N/H labelled set, but missing C/N/H pairs were found.")
}

rung1_candidates <- tibble(
  rung = "rung1",
  candidate_id = c(
    "r1_a_breast_cn",
    "r1_b_primary_cn",
    "r1_c_paired_contrast_cn",
    "r1_d_structured_paired_cn",
    "r1_e_breast_cnh",
    "r1_f_paired_contrast_cnh",
    "r1_g_structured_paired_cnh",
    "r1_h_primary_cnh"
  ),
  candidate_label = c(
    "Rung 1 A: Breast-only C/N",
    "Rung 1 B: Primary-only C/N",
    "Rung 1 C: Paired-contrast C/N",
    "Rung 1 D: Structured paired C/N",
    "Rung 1 E: Breast-only C/N/H",
    "Rung 1 F: Paired-contrast C/N/H",
    "Rung 1 G: Structured paired C/N/H",
    "Rung 1 H: Primary-only C/N/H"
  ),
  model_family = "ridge_glmnet_binomial",
  predictors = list(
    c("normalised_d13c_Breast", "normalised_d15n_Breast"),
    c("normalised_d13c_Primary", "normalised_d15n_Primary"),
    c("abs_delta_d13c", "abs_delta_d15n", "paired_distance_cn"),
    c("normalised_d13c_Breast", "normalised_d15n_Breast", "abs_delta_d13c", "abs_delta_d15n"),
    c("normalised_d13c_Breast", "normalised_d15n_Breast", "normalised_d2h_Breast"),
    c("abs_delta_d13c", "abs_delta_d15n", "abs_delta_d2h", "paired_distance_cnh"),
    c("normalised_d13c_Breast", "normalised_d15n_Breast", "normalised_d2h_Breast", "abs_delta_d13c", "abs_delta_d15n", "abs_delta_d2h"),
    c("normalised_d13c_Primary", "normalised_d15n_Primary", "normalised_d2h_Primary")
  ),
  distance_vars = list(
    character(0),
    character(0),
    c("abs_delta_d13c", "abs_delta_d15n"),
    character(0),
    character(0),
    c("abs_delta_d13c", "abs_delta_d15n", "abs_delta_d2h"),
    character(0),
    character(0)
  ),
  distance_name = list(
    NA_character_,
    NA_character_,
    "paired_distance_cn",
    NA_character_,
    NA_character_,
    "paired_distance_cnh",
    NA_character_,
    NA_character_
  ),
  uses_breast_h = c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE),
  uses_primary_h = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE),
  uses_paired_distance = c(FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE),
  preferred_default = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE),
  caution_level = c(
    "baseline_cn_benchmark",
    "primary_cn_limited_biological_support",
    "paired_cn_benchmark",
    "paired_cn_benchmark",
    "breast_h_provisional_pending_triplicate_check",
    "primary_h_provisional_pending_triplicate_check",
    "primary_h_provisional_pending_triplicate_check",
    "primary_h_provisional_pending_triplicate_check"
  ),
  constraint_note = c(
    "C/N benchmark.",
    "Primary C/N retained as benchmark but with limited direct biological support.",
    "C/N paired contrast benchmark with within-fold standardized distance.",
    "C/N structured paired benchmark.",
    "Includes breast H only; interpret results as provisional pending confirmation that the museum H triplicate workbook is fully incorporated.",
    "Uses paired H contrasts and therefore imports primary H caution; interpret results as provisional pending the museum H triplicate check.",
    "Uses breast H plus paired H contrasts; interpret gains as provisional pending the museum H triplicate check.",
    "Primary H comparator only; not a preferred operational default and still provisional pending the museum H triplicate check."
  )
)

rung2_candidates <- tibble(
  rung = "rung2",
  candidate_id = c(
    "r2_a_breast_cn",
    "r2_b_breast_cnh",
    "r2_c_breast_cnh_plus_paired_contrast",
    "r2_d_breast_cnh_plus_paired_distance"
  ),
  candidate_label = c(
    "Rung 2 A: Breast-only C/N",
    "Rung 2 B: Breast-only C/N/H",
    "Rung 2 C: Breast C/N/H + paired contrasts",
    "Rung 2 D: Breast C/N/H + paired distance"
  ),
  model_family = "ridge_glmnet_multinomial",
  predictors = list(
    c("normalised_d13c_Breast", "normalised_d15n_Breast"),
    c("normalised_d13c_Breast", "normalised_d15n_Breast", "normalised_d2h_Breast"),
    c("normalised_d13c_Breast", "normalised_d15n_Breast", "normalised_d2h_Breast", "abs_delta_d13c", "abs_delta_d15n", "abs_delta_d2h"),
    c("normalised_d13c_Breast", "normalised_d15n_Breast", "normalised_d2h_Breast", "paired_distance_cnh")
  ),
  distance_vars = list(
    character(0),
    character(0),
    character(0),
    c("abs_delta_d13c", "abs_delta_d15n", "abs_delta_d2h")
  ),
  distance_name = list(
    NA_character_,
    NA_character_,
    NA_character_,
    "paired_distance_cnh"
  ),
  uses_breast_h = c(FALSE, TRUE, TRUE, TRUE),
  uses_primary_h = c(FALSE, FALSE, TRUE, TRUE),
  uses_paired_distance = c(FALSE, FALSE, FALSE, TRUE),
  preferred_default = c(TRUE, TRUE, TRUE, FALSE),
  caution_level = c(
    "baseline_cn_benchmark",
    "breast_h_provisional_pending_triplicate_check",
    "primary_h_provisional_pending_triplicate_check",
    "primary_h_provisional_pending_triplicate_check"
  ),
  constraint_note = c(
    "C/N benchmark.",
    "Includes breast H only; interpret results as provisional pending confirmation that the museum H triplicate workbook is fully incorporated.",
    "Paired H contrasts import primary H caution and require clear validated gains; treat any gains as provisional pending the museum H triplicate check.",
    "Paired distance uses C/N/H and therefore imports primary H caution; treat any gains as provisional pending the museum H triplicate check."
  )
)

feature_definitions <- bind_rows(rung1_candidates, rung2_candidates) %>%
  mutate(
    predictors = vapply(predictors, paste, collapse = "|", character(1)),
    distance_vars = vapply(distance_vars, paste, collapse = "|", character(1)),
    distance_name = vapply(distance_name, function(x) ifelse(all(is.na(x)), "", x), character(1))
  )

write_csv(
  feature_definitions,
  file.path(derived_dir, "model_stage2_feature_set_definitions.csv")
)

set.seed(20260421)
rung1_rs <- vfold_cv(stage2_data, v = 5, repeats = 10, strata = status_rung1)
set.seed(20260421)
rung2_rs <- vfold_cv(stage2_data, v = 5, repeats = 10, strata = status_rung2)

rung1_repeated_folds <- build_repeated_cv_fold_file(
  rs_obj = rung1_rs,
  outcome_col = "status_rung1",
  rung_name = "rung1",
  layer_name = "repeated_stratified_cv"
)

rung2_repeated_folds <- build_repeated_cv_fold_file(
  rs_obj = rung2_rs,
  outcome_col = "status_rung2",
  rung_name = "rung2",
  layer_name = "repeated_stratified_cv"
)

write_csv(
  rung1_repeated_folds,
  file.path(derived_dir, "model_stage2_rung1_repeated_cv_folds.csv")
)

write_csv(
  rung2_repeated_folds,
  file.path(derived_dir, "model_stage2_rung2_repeated_cv_folds.csv")
)

coord_fold_map <- make_balanced_coord_folds(
  data = stage2_data,
  outcome_col = "status_rung2",
  group_col = "coord_group",
  v = 5L
)

blocked_fold_file <- stage2_data %>%
  left_join(coord_fold_map, by = c("coord_group" = "group_id")) %>%
  transmute(
    validation_layer = "blocked_coordinate_cv",
    split_id = blocked_fold,
    ring,
    coord_group,
    status_rung1 = as.character(status_rung1),
    status_rung2 = as.character(status_rung2)
  )

write_csv(
  blocked_fold_file,
  file.path(derived_dir, "model_stage2_coordinate_block_folds.csv")
)

rung1_blocked_folds <- blocked_fold_file %>%
  transmute(
    rung = "rung1",
    validation_layer,
    split_id,
    repeat_id = NA_character_,
    fold_id = split_id,
    ring,
    coord_group,
    outcome = status_rung1
  )

rung2_blocked_folds <- blocked_fold_file %>%
  transmute(
    rung = "rung2",
    validation_layer,
    split_id,
    repeat_id = NA_character_,
    fold_id = split_id,
    ring,
    coord_group,
    outcome = status_rung2
  )

rung1_eval_primary <- evaluate_candidates(
  data = stage2_data,
  candidates = rung1_candidates,
  fold_definitions = rung1_repeated_folds,
  outcome_col = "status_rung1",
  family = "binomial",
  rung_name = "rung1"
)

rung1_eval_blocked <- evaluate_candidates(
  data = stage2_data,
  candidates = rung1_candidates,
  fold_definitions = rung1_blocked_folds,
  outcome_col = "status_rung1",
  family = "binomial",
  rung_name = "rung1"
)

rung2_eval_primary <- evaluate_candidates(
  data = stage2_data,
  candidates = rung2_candidates,
  fold_definitions = rung2_repeated_folds,
  outcome_col = "status_rung2",
  family = "multinomial",
  rung_name = "rung2"
)

rung2_eval_blocked <- evaluate_candidates(
  data = stage2_data,
  candidates = rung2_candidates,
  fold_definitions = rung2_blocked_folds,
  outcome_col = "status_rung2",
  family = "multinomial",
  rung_name = "rung2"
)

rung1_metrics <- bind_rows(rung1_eval_primary$metrics, rung1_eval_blocked$metrics)
rung1_predictions <- bind_rows(rung1_eval_primary$predictions, rung1_eval_blocked$predictions)
rung2_metrics <- bind_rows(rung2_eval_primary$metrics, rung2_eval_blocked$metrics)
rung2_predictions <- bind_rows(rung2_eval_primary$predictions, rung2_eval_blocked$predictions)

write_csv(
  rung1_metrics,
  file.path(derived_dir, "model_stage2_rung1_resample_metrics.csv")
)

write_csv(
  rung2_metrics,
  file.path(derived_dir, "model_stage2_rung2_resample_metrics.csv")
)

write_csv(
  rung1_predictions,
  file.path(derived_dir, "model_stage2_rung1_assessment_predictions.csv")
)

write_csv(
  rung2_predictions,
  file.path(derived_dir, "model_stage2_rung2_assessment_predictions.csv")
)

rung1_summary <- summarise_model_comparison(rung1_metrics)
rung2_summary <- summarise_model_comparison(rung2_metrics)
rung1_stability <- summarise_ranking_stability(rung1_metrics)
rung2_stability <- summarise_ranking_stability(rung2_metrics)

write_csv(
  rung1_summary,
  file.path(derived_dir, "model_stage2_rung1_model_comparison.csv")
)

write_csv(
  rung2_summary,
  file.path(derived_dir, "model_stage2_rung2_model_comparison.csv")
)

write_csv(
  rung1_stability,
  file.path(derived_dir, "model_stage2_rung1_ranking_stability.csv")
)

write_csv(
  rung2_stability,
  file.path(derived_dir, "model_stage2_rung2_ranking_stability.csv")
)

rung1_calibration <- summarise_calibration(
  prediction_tbl = rung1_predictions,
  rung_name = "rung1",
  class_levels = c("resident", "migrant")
)

rung2_calibration <- summarise_calibration(
  prediction_tbl = rung2_predictions,
  rung_name = "rung2",
  class_levels = c("resident", "NZ migrant", "AU migrant")
)

write_csv(
  rung1_calibration,
  file.path(derived_dir, "model_stage2_rung1_calibration_summary.csv")
)

write_csv(
  rung2_calibration,
  file.path(derived_dir, "model_stage2_rung2_calibration_summary.csv")
)

pick_recommended_models <- function(summary_tbl, stability_tbl, rung_name) {
  primary <- summary_tbl %>%
    filter(validation_layer == "repeated_stratified_cv") %>%
    arrange(mean_log_loss, desc(mean_bal_accuracy))

  blocked <- summary_tbl %>%
    filter(validation_layer == "blocked_coordinate_cv") %>%
    arrange(mean_log_loss, desc(mean_bal_accuracy))

  merged <- primary %>%
    select(candidate_id, candidate_label, caution_level, preferred_default, uses_breast_h, uses_primary_h, primary_mean_log_loss = mean_log_loss, primary_log_loss_rank = log_loss_rank, primary_mean_bal_accuracy = mean_bal_accuracy) %>%
    left_join(
      blocked %>%
        select(candidate_id, blocked_mean_log_loss = mean_log_loss, blocked_log_loss_rank = log_loss_rank, blocked_mean_bal_accuracy = mean_bal_accuracy),
      by = "candidate_id"
    ) %>%
    left_join(
      stability_tbl %>%
        select(candidate_id, log_loss_rank_shift, bal_accuracy_rank_shift),
      by = "candidate_id"
    ) %>%
    mutate(
      uses_any_h = uses_breast_h | uses_primary_h
    )

  best_non_h_primary <- merged %>%
    filter(!uses_any_h) %>%
    summarise(best = min(primary_mean_log_loss, na.rm = TRUE)) %>%
    pull(best)

  best_non_h_blocked <- merged %>%
    filter(!uses_any_h) %>%
    summarise(best = min(blocked_mean_log_loss, na.rm = TRUE)) %>%
    pull(best)

  merged <- merged %>%
    group_by(rung_name = rung_name) %>%
    mutate(
      baseline_primary_rank = rank(ifelse(uses_any_h, Inf, primary_mean_log_loss), ties.method = "min"),
      h_primary_rank = rank(ifelse(!uses_any_h, Inf, primary_mean_log_loss), ties.method = "min")
    ) %>%
    ungroup() %>%
    mutate(
      recommended_for_stage3 = dplyr::case_when(
        !uses_any_h &
          baseline_primary_rank == 1 &
          blocked_log_loss_rank <= 3 ~ TRUE,
        !uses_any_h &
          baseline_primary_rank == 2 &
          primary_log_loss_rank <= 2 &
          blocked_log_loss_rank <= 3 ~ TRUE,
        uses_any_h &
          h_primary_rank == 1 &
          primary_log_loss_rank <= 3 &
          blocked_log_loss_rank <= 2 &
          primary_mean_log_loss + 0.02 < best_non_h_primary &
          blocked_mean_log_loss + 0.02 < best_non_h_blocked ~ TRUE,
        TRUE ~ FALSE
      ),
      recommendation_note = dplyr::case_when(
        recommended_for_stage3 & uses_any_h & uses_primary_h ~
          "Advance only conditionally, alongside the strongest non-H benchmark, because H evidence is provisional and primary H remains high-caution.",
        recommended_for_stage3 & uses_any_h ~
          "Advance only conditionally, alongside the strongest non-H benchmark, because H evidence is provisional pending the museum H triplicate check.",
        recommended_for_stage3 & baseline_primary_rank == 2 ~
          "Advance to Stage 3 as a secondary non-H paired benchmark because it remained near-leading across both validation layers.",
        recommended_for_stage3 ~
          "Advance to Stage 3 as the main non-H benchmark.",
        TRUE ~
          "Do not prioritize for Stage 3."
      )
    ) %>%
    mutate(rung = rung_name)

  merged
}

rung1_recommendations <- pick_recommended_models(rung1_summary, rung1_stability, "rung1")
rung2_recommendations <- pick_recommended_models(rung2_summary, rung2_stability, "rung2")

stage3_recommendations <- bind_rows(rung1_recommendations, rung2_recommendations)

write_csv(
  stage3_recommendations,
  file.path(derived_dir, "model_stage2_stage3_recommendations.csv")
)

summary_lines <- c(
  "# Stage 2 Summary",
  "",
  "## Validation design",
  "",
  "- Primary comparison: repeated stratified 5-fold CV with 10 repeats, fit separately for Rung 1 and Rung 2.",
  "- Sensitivity layer: exact-coordinate blocked 5-fold assignment using available capture coordinates, balanced against a global class-and-size objective on Rung 2 classes.",
  "- Model family: ridge-penalized logistic / multinomial regression with lambda selected inside each training split via `cv.glmnet`.",
  "- Preprocessing: deterministic breast-primary contrasts were constructed observation-wise, while all scaling and paired-distance standardization were estimated inside each training split only.",
  "",
  "## Rung 1 shortlist",
  ""
)

rung1_lines <- rung1_summary %>%
  filter(validation_layer == "repeated_stratified_cv") %>%
  arrange(log_loss_rank) %>%
  transmute(
    line = paste0(
      "- ", candidate_label, ": mean log loss ",
      format_num(mean_log_loss),
      ", mean balanced accuracy ",
      format_num(mean_bal_accuracy),
      ", caution `", caution_level, "`."
    )
  ) %>%
  pull(line)

summary_lines <- c(summary_lines, rung1_lines, "", "## Rung 2 shortlist", "")

rung2_lines <- rung2_summary %>%
  filter(validation_layer == "repeated_stratified_cv") %>%
  arrange(log_loss_rank) %>%
  transmute(
    line = paste0(
      "- ", candidate_label, ": mean log loss ",
      format_num(mean_log_loss),
      ", mean balanced accuracy ",
      format_num(mean_bal_accuracy),
      ", caution `", caution_level, "`."
    )
  ) %>%
  pull(line)

summary_lines <- c(summary_lines, rung2_lines, "", "## Stage 3 recommendations", "")

rec_lines <- stage3_recommendations %>%
  filter(recommended_for_stage3) %>%
  arrange(rung, primary_log_loss_rank) %>%
  transmute(
    line = paste0(
      "- ", rung, " ", candidate_label, ": ",
      recommendation_note,
      " Primary CV rank ",
      primary_log_loss_rank,
      ", blocked rank ",
      blocked_log_loss_rank,
      "."
    )
  ) %>%
  pull(line)

if (length(rec_lines) == 0) {
  rec_lines <- "- No model met the automatic Stage 3 recommendation rule."
}

summary_lines <- c(
  summary_lines,
  rec_lines,
  "",
  "## Constraint notes",
  "",
  "- Breast H entered candidate models by default, but any gains should still be interpreted with technical caution and remain provisional pending confirmation that `data/Master_TCEA_H_SI_Batch 1_triplicate breast feathers_LEH_Nov25_SB.xls` is fully incorporated into the H evidence base.",
  "- Any model using paired H contrasts or primary H was retained only as a higher-caution comparator because Stage 0 could not directly estimate retained primary-H technical repeatability.",
  "- No Stage 3 uncertainty propagation or final operational classifier was implemented in this script."
)

writeLines(
  summary_lines,
  file.path(derived_dir, "model_stage2_summary.md")
)
