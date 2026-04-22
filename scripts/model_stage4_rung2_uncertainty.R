#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
})

source(file.path("R", "model_stage_repeatability_helpers.R"))
source(file.path("R", "model_stage2_helpers.R"))
source(file.path("R", "model_stage3_helpers.R"))
source(file.path("R", "model_stage4_helpers.R"))

derived_dir <- file.path("data", "derived")
ensure_dir(derived_dir)

format_num <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits))
}

n_mc_draws <- 500L
support_threshold <- 0.8
top_prob_threshold <- 0.5
class_levels <- c("resident", "NZ migrant", "AU migrant")

stage4_data <- read_csv(
  file.path(derived_dir, "live_screening_ready_paired_labelled.csv"),
  show_col_types = FALSE
) %>%
  derive_stage2_features()

stage0_variance <- read_csv(
  file.path(derived_dir, "model_stage0_pooled_assay_variance.csv"),
  show_col_types = FALSE
)

stage1_adequacy <- read_csv(
  file.path(derived_dir, "model_stage1_tissue_summary_adequacy.csv"),
  show_col_types = FALSE
)

measurement_spec <- build_stage3_measurement_spec(stage0_variance, stage1_adequacy)

write_csv(
  measurement_spec,
  file.path(derived_dir, "model_stage4_rung2_measurement_spec.csv")
)

rung2_candidates <- list(
  list(
    candidate_id = "r2_a_breast_cn",
    candidate_label = "Rung 2 A: Breast-only C/N",
    caution_level = "baseline_cn_benchmark",
    predictors = c("normalised_d13c_Breast", "normalised_d15n_Breast"),
    distance_vars = character(0),
    distance_name = NA_character_,
    measurement_cols = c("normalised_d13c_Breast", "normalised_d15n_Breast")
  ),
  list(
    candidate_id = "r2_c_breast_cnh_plus_paired_contrast",
    candidate_label = "Rung 2 C: Breast C/N/H + paired contrasts",
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
    measurement_cols = c(
      "normalised_d13c_Breast",
      "normalised_d13c_Primary",
      "normalised_d15n_Breast",
      "normalised_d15n_Primary",
      "normalised_d2h_Breast",
      "normalised_d2h_Primary"
    )
  )
)

candidate_definitions <- bind_rows(
  lapply(
    rung2_candidates,
    function(candidate) {
      tibble(
        candidate_id = candidate$candidate_id,
        candidate_label = candidate$candidate_label,
        caution_level = candidate$caution_level,
        predictors = paste(candidate$predictors, collapse = "|"),
        distance_vars = paste(candidate$distance_vars, collapse = "|"),
        distance_name = ifelse(is.na(candidate$distance_name), "", candidate$distance_name),
        measurement_cols = paste(candidate$measurement_cols, collapse = "|")
      )
    }
  )
)

write_csv(
  candidate_definitions,
  file.path(derived_dir, "model_stage4_rung2_candidate_definitions.csv")
)

required_measurements <- unique(unlist(lapply(rung2_candidates, `[[`, "measurement_cols")))
missing_spec <- measurement_spec %>%
  filter(measurement_col %in% required_measurements, is.na(working_sd))

if (nrow(missing_spec) > 0) {
  stop("Stage 4 could not find working assay SDs for all required measurement columns.")
}

r2a_constraints <- measurement_spec %>%
  filter(measurement_col %in% c("normalised_d13c_Breast", "normalised_d15n_Breast")) %>%
  mutate(ok = adequacy_for_stage2 == "adequate")

if (!all(r2a_constraints$ok)) {
  stop("R2 A violated the Stage 1 adequacy constraints.")
}

r2c_h_constraints <- measurement_spec %>%
  filter(measurement_col %in% c("normalised_d2h_Breast", "normalised_d2h_Primary"))

if (!identical(
  r2c_h_constraints %>%
    filter(measurement_col == "normalised_d2h_Breast") %>%
    pull(adequacy_for_stage2),
  "adequate_with_technical_caution"
)) {
  stop("Stage 4 expected breast H to remain adequate_with_technical_caution.")
}

if (!identical(
  r2c_h_constraints %>%
    filter(measurement_col == "normalised_d2h_Primary") %>%
    pull(adequacy_for_stage2),
  "provisionally_adequate_with_high_caution"
)) {
  stop("Stage 4 expected primary H to remain provisionally_adequate_with_high_caution.")
}

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
  tibble(
    parameter = c(
      "n_mc_draws",
      "support_threshold",
      "top_prob_threshold"
    ),
    value = c(
      as.character(n_mc_draws),
      as.character(support_threshold),
      as.character(top_prob_threshold)
    )
  ),
  file.path(derived_dir, "model_stage4_rung2_uncertainty_config.csv")
)

split_keys <- all_fold_definitions %>%
  distinct(validation_layer, split_id, repeat_id, fold_id) %>%
  arrange(validation_layer, split_id)

assessment_rows <- vector("list", nrow(split_keys) * length(rung2_candidates))
row_idx <- 1L

for (split_idx in seq_len(nrow(split_keys))) {
  split_meta <- split_keys[split_idx, ]

  test_ids <- all_fold_definitions %>%
    filter(
      validation_layer == split_meta$validation_layer[[1]],
      split_id == split_meta$split_id[[1]]
    ) %>%
    pull(ring)

  train_data <- stage4_data %>% filter(!(ring %in% test_ids))
  test_data <- stage4_data %>% filter(ring %in% test_ids)

  for (candidate_idx in seq_along(rung2_candidates)) {
    candidate <- rung2_candidates[[candidate_idx]]

    fit <- fit_stage4_ridge_multinomial(
      train_data = train_data,
      candidate = candidate,
      outcome_col = "status_rung2",
      seed = 9000L + split_idx * 100L + candidate_idx
    )

    if (!identical(fit$status, "ok")) {
      stop(paste("Stage 4 fit failed for", candidate$candidate_id, "in", split_meta$split_id[[1]], "with status", fit$status))
    }

    base_pred <- predict_stage4_ridge_multinomial(fit, test_data)
    base_pred <- bind_cols(
      test_data %>%
        select(ring, coord_group, longitude_capture, latitude_capture, status_rung2),
      base_pred
    ) %>%
      rename(truth = status_rung2)

    draw_prob_array <- array(
      NA_real_,
      dim = c(nrow(test_data), length(class_levels), n_mc_draws),
      dimnames = list(NULL, class_levels, NULL)
    )

    for (draw_idx in seq_len(n_mc_draws)) {
      perturbed_test <- perturb_stage3_measurements(
        data = test_data,
        measurement_spec = measurement_spec,
        measurement_cols = candidate$measurement_cols,
        seed = 1200000L + split_idx * 10000L + candidate_idx * 1000L + draw_idx
      )

      draw_pred <- predict_stage4_ridge_multinomial(fit, perturbed_test)
      draw_prob_array[, , draw_idx] <- as.matrix(draw_pred[, class_levels, drop = FALSE])
    }

    uncertainty_summary <- summarise_stage4_uncertainty(
      base_prediction_tbl = base_pred,
      draw_prob_array = draw_prob_array,
      class_levels = class_levels,
      support_threshold = support_threshold,
      top_prob_threshold = top_prob_threshold
    )

    assessment_rows[[row_idx]] <- bind_cols(
      tibble(
        candidate_id = candidate$candidate_id,
        candidate_label = candidate$candidate_label,
        caution_level = candidate$caution_level,
        validation_layer = split_meta$validation_layer[[1]],
        split_id = split_meta$split_id[[1]],
        repeat_id = split_meta$repeat_id[[1]],
        fold_id = split_meta$fold_id[[1]]
      ),
      base_pred %>%
        select(ring, coord_group, longitude_capture, latitude_capture, truth, all_of(class_levels), base_class),
      uncertainty_summary %>%
        select(-truth, -base_class)
    )

    row_idx <- row_idx + 1L
  }
}

assessment_summary <- bind_rows(assessment_rows)

write_csv(
  assessment_summary,
  file.path(derived_dir, "model_stage4_rung2_uncertainty_assessment_summary.csv")
)

split_summary <- assessment_summary %>%
  group_by(
    candidate_id,
    candidate_label,
    caution_level,
    validation_layer,
    split_id,
    repeat_id,
    fold_id
  ) %>%
  group_modify(~ summarise_stage4_split_metrics(.x, class_levels = class_levels)) %>%
  ungroup()

model_summary <- summarise_stage4_model_metrics(split_summary)
calibration_summary <- summarise_stage4_calibration(assessment_summary, class_levels = class_levels)

confusion_summary <- assessment_summary %>%
  count(candidate_id, candidate_label, validation_layer, truth, operational_call, name = "n")

write_csv(
  split_summary,
  file.path(derived_dir, "model_stage4_rung2_uncertainty_split_summary.csv")
)

write_csv(
  model_summary,
  file.path(derived_dir, "model_stage4_rung2_uncertainty_model_summary.csv")
)

write_csv(
  calibration_summary,
  file.path(derived_dir, "model_stage4_rung2_uncertainty_calibration_summary.csv")
)

write_csv(
  confusion_summary,
  file.path(derived_dir, "model_stage4_rung2_uncertainty_confusion.csv")
)

calibration_overall <- calibration_summary %>%
  filter(class == "overall") %>%
  select(candidate_id, candidate_label, validation_layer, metric, value) %>%
  tidyr::pivot_wider(names_from = metric, values_from = value)

summary_lines <- c(
  "# Stage 4 Rung 2 Uncertainty Summary",
  "",
  "## Configuration",
  "",
  paste0("- Monte Carlo draws per assessment set: ", n_mc_draws, "."),
  paste0("- Indeterminate rule: winning-class support < ", format_num(support_threshold, 2), " or mean top-class probability < ", format_num(top_prob_threshold, 2), "."),
  "- Training/preprocessing remained inside each outer resample split; assay working SDs were read from the Stage 0 pooled-variance table.",
  "- R2 C includes primary H under proxy variance only and should remain high-caution even if it performs better.",
  "",
  "## Model summary",
  ""
)

summary_lines <- c(
  summary_lines,
  model_summary %>%
    left_join(
      calibration_overall,
      by = c("candidate_id", "candidate_label", "validation_layer")
    ) %>%
    transmute(
      line = paste0(
        "- ", candidate_label, " [", validation_layer, "]: uncertainty log loss ",
        format_num(mean_uncertainty_log_loss),
        ", uncertainty balanced accuracy ",
        format_num(mean_uncertainty_bal_accuracy),
        ", macro Brier ",
        format_num(macro_brier),
        ", confidence ECE ",
        format_num(confidence_ece),
        ", indeterminate rate ",
        format_num(mean_indeterminate_rate),
        ", mean switch rate ",
        format_num(mean_switch_rate),
        "."
      )
    ) %>%
    pull(line),
  "",
  "## Interpretation notes",
  "",
  "- R2 A remains the non-H benchmark for winter-signal classification.",
  "- R2 C is the only higher-information candidate carried forward here, but any apparent gain must survive uncertainty propagation under explicit primary-H caution.",
  "- Rung 3 remains out of scope."
)

writeLines(
  summary_lines,
  file.path(derived_dir, "model_stage4_rung2_uncertainty_summary.md")
)
