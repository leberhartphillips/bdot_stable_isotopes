#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(stringr)
})

source(file.path("R", "model_stage_repeatability_helpers.R"))
source(file.path("R", "model_stage2_helpers.R"))
source(file.path("R", "model_stage3_helpers.R"))

derived_dir <- file.path("data", "derived")
ensure_dir(derived_dir)

format_num <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits))
}

n_mc_draws <- 500L
support_threshold <- 0.8
boundary_interval <- c(0.1, 0.9)

stage3_data <- read_csv(
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
  file.path(derived_dir, "model_stage3_rung1_measurement_spec.csv")
)

rung1_candidates <- list(
  list(
    candidate_id = "r1_c_paired_contrast_cn",
    candidate_label = "Rung 1 C: Paired-contrast C/N",
    caution_level = "paired_cn_benchmark",
    sensitivity_only = FALSE,
    predictors = c("abs_delta_d13c", "abs_delta_d15n", "paired_distance_cn"),
    distance_vars = c("abs_delta_d13c", "abs_delta_d15n"),
    distance_name = "paired_distance_cn",
    measurement_cols = c(
      "normalised_d13c_Breast",
      "normalised_d13c_Primary",
      "normalised_d15n_Breast",
      "normalised_d15n_Primary"
    )
  ),
  list(
    candidate_id = "r1_d_structured_paired_cn",
    candidate_label = "Rung 1 D: Structured paired C/N",
    caution_level = "paired_cn_benchmark",
    sensitivity_only = FALSE,
    predictors = c(
      "normalised_d13c_Breast",
      "normalised_d15n_Breast",
      "abs_delta_d13c",
      "abs_delta_d15n"
    ),
    distance_vars = character(0),
    distance_name = NA_character_,
    measurement_cols = c(
      "normalised_d13c_Breast",
      "normalised_d13c_Primary",
      "normalised_d15n_Breast",
      "normalised_d15n_Primary"
    )
  ),
  list(
    candidate_id = "r1_g_structured_paired_cnh",
    candidate_label = "Rung 1 G: Structured paired C/N/H",
    caution_level = "primary_h_high_caution_tertiary_sensitivity",
    sensitivity_only = TRUE,
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
    rung1_candidates,
    function(candidate) {
      tibble(
        candidate_id = candidate$candidate_id,
        candidate_label = candidate$candidate_label,
        caution_level = candidate$caution_level,
        sensitivity_only = candidate$sensitivity_only,
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
  file.path(derived_dir, "model_stage3_rung1_candidate_definitions.csv")
)

required_measurements <- unique(unlist(lapply(rung1_candidates, `[[`, "measurement_cols")))
missing_spec <- measurement_spec %>%
  filter(measurement_col %in% required_measurements, is.na(working_sd))

if (nrow(missing_spec) > 0) {
  stop("Stage 3 could not find working assay SDs for all required measurement columns.")
}

cn_measurements <- c(
  "normalised_d13c_Breast",
  "normalised_d13c_Primary",
  "normalised_d15n_Breast",
  "normalised_d15n_Primary"
)

cn_constraints <- measurement_spec %>%
  filter(measurement_col %in% cn_measurements) %>%
  mutate(ok = adequacy_for_stage2 %in% c("adequate", "adequate_with_limited_biological_evidence"))

if (!all(cn_constraints$ok)) {
  stop("Stage 3 C/N candidates violated the Stage 1 adequacy constraints.")
}

h_constraints <- measurement_spec %>%
  filter(measurement_col %in% c("normalised_d2h_Breast", "normalised_d2h_Primary"))

if (!identical(
  h_constraints %>%
    filter(measurement_col == "normalised_d2h_Breast") %>%
    pull(adequacy_for_stage2),
  "adequate_with_technical_caution"
)) {
  stop("Stage 3 expected breast H to remain adequate_with_technical_caution.")
}

if (!identical(
  h_constraints %>%
    filter(measurement_col == "normalised_d2h_Primary") %>%
    pull(adequacy_for_stage2),
  "provisionally_adequate_with_high_caution"
)) {
  stop("Stage 3 expected primary H to remain provisionally_adequate_with_high_caution.")
}

repeated_folds <- read_csv(
  file.path(derived_dir, "model_stage2_rung1_repeated_cv_folds.csv"),
  show_col_types = FALSE
)

blocked_fold_file <- read_csv(
  file.path(derived_dir, "model_stage2_coordinate_block_folds.csv"),
  show_col_types = FALSE
)

blocked_folds <- blocked_fold_file %>%
  transmute(
    rung = "rung1",
    validation_layer = "blocked_coordinate_cv",
    split_id,
    repeat_id = NA_character_,
    fold_id = split_id,
    ring,
    coord_group,
    outcome = status_rung1
  )

all_fold_definitions <- bind_rows(repeated_folds, blocked_folds)

write_csv(
  tibble(
    parameter = c(
      "n_mc_draws",
      "support_threshold",
      "boundary_interval_low",
      "boundary_interval_high"
    ),
    value = c(
      as.character(n_mc_draws),
      as.character(support_threshold),
      as.character(boundary_interval[1]),
      as.character(boundary_interval[2])
    )
  ),
  file.path(derived_dir, "model_stage3_rung1_uncertainty_config.csv")
)

split_keys <- all_fold_definitions %>%
  distinct(validation_layer, split_id, repeat_id, fold_id) %>%
  arrange(validation_layer, split_id)

assessment_rows <- vector("list", nrow(split_keys) * length(rung1_candidates))
split_rows <- vector("list", nrow(split_keys) * length(rung1_candidates))
row_idx <- 1L

for (split_idx in seq_len(nrow(split_keys))) {
  split_meta <- split_keys[split_idx, ]

  test_ids <- all_fold_definitions %>%
    filter(
      validation_layer == split_meta$validation_layer[[1]],
      split_id == split_meta$split_id[[1]]
    ) %>%
    pull(ring)

  train_data <- stage3_data %>% filter(!(ring %in% test_ids))
  test_data <- stage3_data %>% filter(ring %in% test_ids)

  for (candidate_idx in seq_along(rung1_candidates)) {
    candidate <- rung1_candidates[[candidate_idx]]

    fit <- fit_stage3_ridge_binomial(
      train_data = train_data,
      candidate = candidate,
      outcome_col = "status_rung1",
      seed = 5000L + split_idx * 100L + candidate_idx
    )

    if (!identical(fit$status, "ok")) {
      stop(paste("Stage 3 fit failed for", candidate$candidate_id, "in", split_meta$split_id[[1]], "with status", fit$status))
    }

    base_pred <- predict_stage3_ridge_binomial(fit, test_data)
    base_pred <- bind_cols(
      test_data %>%
        select(ring, coord_group, longitude_capture, latitude_capture, status_rung1),
      base_pred
    ) %>%
      rename(truth = status_rung1)

    draw_prob_matrix <- matrix(NA_real_, nrow = nrow(test_data), ncol = n_mc_draws)

    for (draw_idx in seq_len(n_mc_draws)) {
      perturbed_test <- perturb_stage3_measurements(
        data = test_data,
        measurement_spec = measurement_spec,
        measurement_cols = candidate$measurement_cols,
        seed = 700000L + split_idx * 10000L + candidate_idx * 1000L + draw_idx
      )

      draw_pred <- predict_stage3_ridge_binomial(fit, perturbed_test)
      draw_prob_matrix[, draw_idx] <- draw_pred$.pred_migrant
    }

    uncertainty_summary <- summarise_stage3_uncertainty(
      base_prediction_tbl = base_pred,
      draw_prob_matrix = draw_prob_matrix,
      support_threshold = support_threshold,
      boundary_interval = boundary_interval
    )

    assessment_rows[[row_idx]] <- bind_cols(
      tibble(
        candidate_id = candidate$candidate_id,
        candidate_label = candidate$candidate_label,
        caution_level = candidate$caution_level,
        sensitivity_only = candidate$sensitivity_only,
        validation_layer = split_meta$validation_layer[[1]],
        split_id = split_meta$split_id[[1]],
        repeat_id = split_meta$repeat_id[[1]],
        fold_id = split_meta$fold_id[[1]]
      ),
      base_pred %>%
        select(ring, coord_group, longitude_capture, latitude_capture),
      uncertainty_summary
    )

    split_rows[[row_idx]] <- bind_cols(
      tibble(
        candidate_id = candidate$candidate_id,
        candidate_label = candidate$candidate_label,
        caution_level = candidate$caution_level,
        sensitivity_only = candidate$sensitivity_only,
        validation_layer = split_meta$validation_layer[[1]],
        split_id = split_meta$split_id[[1]],
        repeat_id = split_meta$repeat_id[[1]],
        fold_id = split_meta$fold_id[[1]]
      ),
      summarise_stage3_split_metrics(assessment_rows[[row_idx]])
    )

    row_idx <- row_idx + 1L
  }
}

assessment_summary <- bind_rows(assessment_rows)
split_summary <- assessment_summary %>%
  group_by(
    candidate_id,
    candidate_label,
    caution_level,
    sensitivity_only,
    validation_layer,
    split_id,
    repeat_id,
    fold_id
  ) %>%
  group_modify(~ summarise_stage3_split_metrics(.x)) %>%
  ungroup()

model_summary <- summarise_stage3_model_metrics(split_summary)

confusion_summary <- assessment_summary %>%
  count(candidate_id, candidate_label, validation_layer, truth, operational_call, name = "n")

write_csv(
  assessment_summary,
  file.path(derived_dir, "model_stage3_rung1_uncertainty_assessment_summary.csv")
)

write_csv(
  split_summary,
  file.path(derived_dir, "model_stage3_rung1_uncertainty_split_summary.csv")
)

write_csv(
  model_summary,
  file.path(derived_dir, "model_stage3_rung1_uncertainty_model_summary.csv")
)

write_csv(
  confusion_summary,
  file.path(derived_dir, "model_stage3_rung1_uncertainty_confusion.csv")
)

summary_lines <- c(
  "# Stage 3 Rung 1 Uncertainty Summary",
  "",
  "## Configuration",
  "",
  paste0("- Monte Carlo draws per assessment set: ", n_mc_draws, "."),
  paste0("- Indeterminate rule: class support < ", format_num(support_threshold, 2), " or central ",
         as.integer((boundary_interval[2] - boundary_interval[1]) * 100),
         "% migrant-probability interval crosses 0.5."),
  "- Training/preprocessing remained inside each outer resample split; assay working SDs were read from the Stage 0 pooled-variance table.",
  "",
  "## Model summary",
  ""
)

summary_lines <- c(
  summary_lines,
  model_summary %>%
    transmute(
      line = paste0(
        "- ", candidate_label, " [", validation_layer, "]: uncertainty log loss ",
        format_num(mean_uncertainty_log_loss),
        ", uncertainty balanced accuracy ",
        format_num(mean_uncertainty_bal_accuracy),
        ", indeterminate rate ",
        format_num(mean_indeterminate_rate),
        ", final-class change rate ",
        format_num(mean_final_class_changed_rate),
        ", mean switch rate ",
        format_num(mean_switch_rate),
        "."
      )
    ) %>%
    pull(line),
  "",
  "## Interpretation notes",
  "",
  "- Breast H support is now treated as explicitly documented by Stage 1 museum triplicate evidence, but primary H remains a proxy-variance, high-caution component.",
  "- The optional C/N/H model is retained only as a tertiary sensitivity analysis and not as a preferred operational target.",
  "- No Rung 2 uncertainty propagation was run here."
)

writeLines(
  summary_lines,
  file.path(derived_dir, "model_stage3_rung1_uncertainty_summary.md")
)
