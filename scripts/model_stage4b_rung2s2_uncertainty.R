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
source(file.path("R", "model_stage2b_rung2_staged_helpers.R"))
source(file.path("R", "model_stage4b_rung2_staged_helpers.R"))

derived_dir <- file.path("data", "derived")
ensure_dir(derived_dir)

format_num <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits))
}

n_mc_draws <- 500L
support_threshold <- 0.8
top_prob_threshold <- 0.5
class_levels <- c("resident", "NZ migrant", "AU migrant")

stage4b_data <- read_csv(
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
  file.path(derived_dir, "model_stage4b_rung2s2_measurement_spec.csv")
)

r1_candidate <- list(
  candidate_id = "r1_c_paired_contrast_cn",
  candidate_label = "Rung 1 C: Paired-contrast C/N",
  predictors = c("abs_delta_d13c", "abs_delta_d15n", "paired_distance_cn"),
  distance_vars = c("abs_delta_d13c", "abs_delta_d15n"),
  distance_name = "paired_distance_cn"
)

migrant_submodel_candidate <- list(
  candidate_id = "r2s2_migrant_submodel_breast_cnh",
  predictors = c(
    "normalised_d13c_Breast",
    "normalised_d15n_Breast",
    "normalised_d2h_Breast"
  ),
  distance_vars = character(0),
  distance_name = NA_character_
)

candidate_meta <- tibble(
  candidate_id = "r2s_b_soft_hierarchical_r1c_then_breast_cnh",
  candidate_label = "Rung 2 S2: Soft hierarchical R1 C then breast C/N/H NZ-vs-AU",
  caution_level = "staged_r1c_then_breast_cnh_no_primary_h",
  uses_r1_upstream = TRUE,
  uses_breast_h = TRUE,
  uses_primary_h = FALSE,
  uses_hard_or_soft_staging = "soft_hierarchical",
  measurement_cols = "normalised_d13c_Breast|normalised_d15n_Breast|normalised_d2h_Breast|normalised_d13c_Primary|normalised_d15n_Primary"
)

write_csv(
  candidate_meta,
  file.path(derived_dir, "model_stage4b_rung2s2_candidate_definition.csv")
)

required_measurements <- c(
  "normalised_d13c_Breast",
  "normalised_d15n_Breast",
  "normalised_d2h_Breast",
  "normalised_d13c_Primary",
  "normalised_d15n_Primary"
)

missing_spec <- measurement_spec %>%
  filter(measurement_col %in% required_measurements, is.na(working_sd))

if (nrow(missing_spec) > 0) {
  stop("Stage 4b could not find working assay SDs for all required measurement columns.")
}

h_constraint <- measurement_spec %>%
  filter(measurement_col == "normalised_d2h_Breast") %>%
  pull(adequacy_for_stage2)

if (!identical(h_constraint, "adequate_with_technical_caution")) {
  stop("Stage 4b expected breast H to remain adequate_with_technical_caution.")
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
  all_fold_definitions,
  file.path(derived_dir, "model_stage4b_rung2s2_outer_folds.csv")
)

write_csv(
  tibble(
    parameter = c("n_mc_draws", "support_threshold", "top_prob_threshold"),
    value = c(
      as.character(n_mc_draws),
      as.character(support_threshold),
      as.character(top_prob_threshold)
    )
  ),
  file.path(derived_dir, "model_stage4b_rung2s2_uncertainty_config.csv")
)

split_keys <- all_fold_definitions %>%
  distinct(validation_layer, split_id, repeat_id, fold_id) %>%
  arrange(validation_layer, split_id)

assessment_rows <- vector("list", nrow(split_keys))

for (split_idx in seq_len(nrow(split_keys))) {
  split_meta <- split_keys[split_idx, ]

  test_ids <- all_fold_definitions %>%
    filter(
      validation_layer == split_meta$validation_layer[[1]],
      split_id == split_meta$split_id[[1]]
    ) %>%
    pull(ring)

  train_data <- stage4b_data %>% filter(!(ring %in% test_ids))
  test_data <- stage4b_data %>% filter(ring %in% test_ids)

  fit <- fit_stage4b_r2s2(
    train_data = train_data,
    r1_candidate = r1_candidate,
    migrant_candidate = migrant_submodel_candidate,
    seed_r1 = 810000L + split_idx,
    seed_submodel = 820000L + split_idx
  )

  if (!identical(fit$status, "ok")) {
    stop(
      paste(
        "Stage 4b fit failed for",
        candidate_meta$candidate_id[[1]],
        "in",
        split_meta$split_id[[1]],
        "with status",
        fit$status
      )
    )
  }

  base_pred <- predict_stage4b_r2s2(fit, test_data) %>%
    bind_cols(
      test_data %>%
        select(ring, coord_group, longitude_capture, latitude_capture, status_rung2)
    ) %>%
    rename(truth = status_rung2) %>%
    relocate(ring, coord_group, longitude_capture, latitude_capture, truth)

  draw_prob_array <- array(
    NA_real_,
    dim = c(nrow(test_data), length(class_levels), n_mc_draws),
    dimnames = list(NULL, class_levels, NULL)
  )

  for (draw_idx in seq_len(n_mc_draws)) {
    perturbed_test <- perturb_stage3_measurements(
      data = test_data,
      measurement_spec = measurement_spec,
      measurement_cols = required_measurements,
      seed = 1300000L + split_idx * 10000L + draw_idx
    )

    draw_pred <- predict_stage4b_r2s2(fit, perturbed_test)
    draw_prob_array[, , draw_idx] <- as.matrix(draw_pred[, class_levels, drop = FALSE])
  }

  uncertainty_summary <- summarise_stage4_uncertainty(
    base_prediction_tbl = base_pred %>%
      select(truth, all_of(class_levels), base_class),
    draw_prob_array = draw_prob_array,
    class_levels = class_levels,
    support_threshold = support_threshold,
    top_prob_threshold = top_prob_threshold
  )

  assessment_rows[[split_idx]] <- bind_cols(
    candidate_meta %>%
      select(
        candidate_id,
        candidate_label,
        caution_level,
        uses_r1_upstream,
        uses_breast_h,
        uses_primary_h,
        uses_hard_or_soft_staging
      ),
    tibble(
      validation_layer = split_meta$validation_layer[[1]],
      split_id = split_meta$split_id[[1]],
      repeat_id = split_meta$repeat_id[[1]],
      fold_id = split_meta$fold_id[[1]]
    ),
    base_pred %>%
      select(
        ring,
        coord_group,
        longitude_capture,
        latitude_capture,
        truth,
        all_of(class_levels),
        base_class,
        r1_prob_migrant,
        p_au_given_migrant
      ),
    uncertainty_summary %>%
      select(
        starts_with("uncertainty_prob_"),
        uncertainty_class,
        class_support,
        top_class_mean_prob,
        top_class_prob_q10,
        top_class_prob_q90,
        any_draw_class_flip,
        switch_rate_from_base,
        top_prob_base,
        top_prob_mean,
        top_prob_abs_shift,
        indeterminate,
        operational_call,
        final_class_changed
      )
  )
}

assessment_summary <- bind_rows(assessment_rows)

write_csv(
  assessment_summary,
  file.path(derived_dir, "model_stage4b_rung2s2_uncertainty_assessment_summary.csv")
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
  file.path(derived_dir, "model_stage4b_rung2s2_uncertainty_split_summary.csv")
)

write_csv(
  model_summary,
  file.path(derived_dir, "model_stage4b_rung2s2_uncertainty_model_summary.csv")
)

write_csv(
  calibration_summary,
  file.path(derived_dir, "model_stage4b_rung2s2_uncertainty_calibration_summary.csv")
)

write_csv(
  confusion_summary,
  file.path(derived_dir, "model_stage4b_rung2s2_uncertainty_confusion.csv")
)

r2c_model_summary <- read_csv(
  file.path(derived_dir, "model_stage4_rung2_uncertainty_model_summary.csv"),
  show_col_types = FALSE
) %>%
  filter(candidate_id == "r2_c_breast_cnh_plus_paired_contrast")

r2c_calibration <- read_csv(
  file.path(derived_dir, "model_stage4_rung2_uncertainty_calibration_summary.csv"),
  show_col_types = FALSE
) %>%
  filter(
    candidate_id == "r2_c_breast_cnh_plus_paired_contrast",
    class == "overall"
  ) %>%
  select(candidate_id, validation_layer, metric, value) %>%
  tidyr::pivot_wider(names_from = metric, values_from = value)

r2s2_calibration <- calibration_summary %>%
  filter(class == "overall") %>%
  select(candidate_id, validation_layer, metric, value) %>%
  tidyr::pivot_wider(names_from = metric, values_from = value)

comparison_tbl <- bind_rows(r2c_model_summary, model_summary) %>%
  left_join(
    bind_rows(r2c_calibration, r2s2_calibration),
    by = c("candidate_id", "validation_layer")
  ) %>%
  mutate(
    comparison_group = case_when(
      candidate_id == "r2_c_breast_cnh_plus_paired_contrast" ~ "accepted_direct",
      TRUE ~ "staged_soft_hierarchical"
    )
  ) %>%
  group_by(validation_layer) %>%
  mutate(
    uncertainty_log_loss_rank = rank(mean_uncertainty_log_loss, ties.method = "min"),
    uncertainty_bal_accuracy_rank = rank(-mean_uncertainty_bal_accuracy, ties.method = "min")
  ) %>%
  ungroup() %>%
  arrange(validation_layer, uncertainty_log_loss_rank, uncertainty_bal_accuracy_rank)

write_csv(
  comparison_tbl,
  file.path(derived_dir, "model_stage4b_rung2s2_vs_r2c_comparison.csv")
)

summary_lines <- c(
  "# Stage 4b Rung 2 S2 Uncertainty Summary",
  "",
  "## Configuration",
  "",
  paste0("- Monte Carlo draws per assessment set: ", n_mc_draws, "."),
  paste0(
    "- Indeterminate rule: winning-class support < ",
    format_num(support_threshold, 2),
    " or mean top-class probability < ",
    format_num(top_prob_threshold, 2),
    "."
  ),
  "- Outer validation splits matched the accepted Rung 2 repeated and blocked fold structure.",
  "- The staged workflow used the accepted R1 C model upstream and a breast C/N/H NZ-vs-AU migrant submodel downstream.",
  "- Assessment rows did not contribute to either the upstream R1 fit or the migrant submodel fit.",
  "",
  "## Comparison against accepted uncertainty-aware R2 C",
  ""
)

summary_lines <- c(
  summary_lines,
  comparison_tbl %>%
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
  "- R2 S2 avoids primary H entirely, so any retained advantage over R2 C should be interpreted as biologically aligned rather than H-driven.",
  "- This run does not replace the accepted direct R2 C result; it evaluates whether the staged soft-hierarchical option remains attractive after uncertainty propagation.",
  "- Rung 3 remains out of scope."
)

writeLines(
  summary_lines,
  file.path(derived_dir, "model_stage4b_rung2s2_uncertainty_summary.md")
)
