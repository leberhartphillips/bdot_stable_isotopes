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
source(file.path("R", "model_stage2b_rung2_staged_helpers.R"))
source(file.path("R", "model_stage4_helpers.R"))
source(file.path("R", "model_stage4b_rung2_staged_helpers.R"))

derived_dir <- file.path("data", "derived")
ensure_dir(derived_dir)

format_num <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits))
}

complete_measurements <- function(data, measurement_cols) {
  stats::complete.cases(data[, measurement_cols, drop = FALSE])
}

build_unavailable_rows <- function(data, availability_note, extra_cols = list()) {
  base_tbl <- data %>%
    transmute(
      source_order,
      ring,
      longitude_capture,
      latitude_capture,
      coord_group,
      prediction_available = FALSE,
      availability_note = availability_note
    )

  extra_tbl <- as_tibble(extra_cols)
  if (nrow(extra_tbl) == 1 && nrow(base_tbl) > 1) {
    extra_tbl <- extra_tbl[rep(1, nrow(base_tbl)), , drop = FALSE]
  }

  bind_cols(base_tbl, extra_tbl)
}

unknown_predict_stage3 <- function(
  unknown_data,
  model_obj,
  measurement_spec,
  measurement_cols,
  n_mc_draws,
  support_threshold,
  boundary_interval
) {
  available_mask <- complete_measurements(unknown_data, measurement_cols)
  available_data <- unknown_data[available_mask, , drop = FALSE]
  unavailable_data <- unknown_data[!available_mask, , drop = FALSE]

  available_out <- tibble()
  if (nrow(available_data) > 0) {
    base_pred <- predict_stage3_ridge_binomial(model_obj, available_data)
    base_tbl <- bind_cols(
      tibble(
        truth = factor(rep(NA_character_, nrow(available_data)), levels = c("resident", "migrant"))
      ),
      base_pred
    )

    draw_prob_matrix <- matrix(
      NA_real_,
      nrow = nrow(available_data),
      ncol = n_mc_draws
    )

    for (draw_idx in seq_len(n_mc_draws)) {
      perturbed_data <- perturb_stage3_measurements(
        data = available_data,
        measurement_spec = measurement_spec,
        measurement_cols = measurement_cols,
        seed = 2100000L + draw_idx
      )

      draw_pred <- predict_stage3_ridge_binomial(model_obj, perturbed_data)
      draw_prob_matrix[, draw_idx] <- draw_pred$.pred_migrant
    }

    uncertainty_tbl <- summarise_stage3_uncertainty(
      base_prediction_tbl = base_tbl,
      draw_prob_matrix = draw_prob_matrix,
      support_threshold = support_threshold,
      boundary_interval = boundary_interval
    )

    available_out <- bind_cols(
      available_data %>%
        transmute(
          source_order,
          ring,
          longitude_capture,
          latitude_capture,
          coord_group,
          prediction_available = TRUE,
          availability_note = NA_character_
        ),
      tibble(
        prob_resident_base = base_pred$.pred_resident,
        prob_migrant_base = base_pred$.pred_migrant,
        prob_resident = 1 - uncertainty_tbl$uncertainty_prob_migrant_mean,
        prob_migrant = uncertainty_tbl$uncertainty_prob_migrant_mean,
        prob_migrant_q10 = uncertainty_tbl$uncertainty_prob_migrant_q10,
        prob_migrant_q90 = uncertainty_tbl$uncertainty_prob_migrant_q90,
        top_class = as.character(uncertainty_tbl$uncertainty_class),
        class_support = uncertainty_tbl$class_support,
        indeterminate = uncertainty_tbl$indeterminate,
        operational_call = ifelse(
          uncertainty_tbl$indeterminate,
          "indeterminate",
          as.character(uncertainty_tbl$uncertainty_class)
        ),
        boundary_crosses = uncertainty_tbl$boundary_crosses,
        switch_rate = uncertainty_tbl$switch_rate_from_base,
        top_prob_mean = uncertainty_tbl$top_prob_mean
      )
    )
  }

  unavailable_out <- tibble()
  if (nrow(unavailable_data) > 0) {
    unavailable_out <- build_unavailable_rows(
      data = unavailable_data,
      availability_note = paste(
        "Missing required inputs for accepted R1 C paired-contrast C/N model:",
        paste(measurement_cols, collapse = ", ")
      ),
      extra_cols = list(
        prob_resident_base = NA_real_,
        prob_migrant_base = NA_real_,
        prob_resident = NA_real_,
        prob_migrant = NA_real_,
        prob_migrant_q10 = NA_real_,
        prob_migrant_q90 = NA_real_,
        top_class = NA_character_,
        class_support = NA_real_,
        indeterminate = NA,
        operational_call = NA_character_,
        boundary_crosses = NA,
        switch_rate = NA_real_,
        top_prob_mean = NA_real_
      )
    )
  }

  bind_rows(available_out, unavailable_out) %>%
    arrange(source_order)
}

unknown_predict_stage4_multinomial <- function(
  unknown_data,
  model_obj,
  measurement_spec,
  measurement_cols,
  class_levels,
  n_mc_draws,
  support_threshold,
  top_prob_threshold,
  availability_note_prefix
) {
  available_mask <- complete_measurements(unknown_data, measurement_cols)
  available_data <- unknown_data[available_mask, , drop = FALSE]
  unavailable_data <- unknown_data[!available_mask, , drop = FALSE]

  available_out <- tibble()
  if (nrow(available_data) > 0) {
    base_pred <- predict_stage4_ridge_multinomial(model_obj, available_data)
    base_tbl <- bind_cols(
      tibble(
        truth = factor(rep(NA_character_, nrow(available_data)), levels = class_levels)
      ),
      base_pred[, c(class_levels, "base_class"), drop = FALSE]
    )

    draw_prob_array <- array(
      NA_real_,
      dim = c(nrow(available_data), length(class_levels), n_mc_draws),
      dimnames = list(NULL, class_levels, NULL)
    )

    for (draw_idx in seq_len(n_mc_draws)) {
      perturbed_data <- perturb_stage3_measurements(
        data = available_data,
        measurement_spec = measurement_spec,
        measurement_cols = measurement_cols,
        seed = 2200000L + draw_idx
      )

      draw_pred <- predict_stage4_ridge_multinomial(model_obj, perturbed_data)
      draw_prob_array[, , draw_idx] <- as.matrix(draw_pred[, class_levels, drop = FALSE])
    }

    uncertainty_tbl <- summarise_stage4_uncertainty(
      base_prediction_tbl = base_tbl,
      draw_prob_array = draw_prob_array,
      class_levels = class_levels,
      support_threshold = support_threshold,
      top_prob_threshold = top_prob_threshold
    )

    available_out <- bind_cols(
      available_data %>%
        transmute(
          source_order,
          ring,
          longitude_capture,
          latitude_capture,
          coord_group,
          prediction_available = TRUE,
          availability_note = NA_character_
        ),
      tibble(
        prob_resident_base = base_pred$resident,
        prob_nz_migrant_base = base_pred$`NZ migrant`,
        prob_au_migrant_base = base_pred$`AU migrant`,
        prob_resident = uncertainty_tbl$uncertainty_prob_resident,
        prob_nz_migrant = uncertainty_tbl$uncertainty_prob_NZ_migrant,
        prob_au_migrant = uncertainty_tbl$uncertainty_prob_AU_migrant,
        top_class = as.character(uncertainty_tbl$uncertainty_class),
        class_support = uncertainty_tbl$class_support,
        indeterminate = uncertainty_tbl$indeterminate,
        operational_call = as.character(uncertainty_tbl$operational_call),
        switch_rate = uncertainty_tbl$switch_rate_from_base,
        top_prob_mean = uncertainty_tbl$top_prob_mean
      )
    )
  }

  unavailable_out <- tibble()
  if (nrow(unavailable_data) > 0) {
    unavailable_out <- build_unavailable_rows(
      data = unavailable_data,
      availability_note = paste(
        availability_note_prefix,
        paste(measurement_cols, collapse = ", ")
      ),
      extra_cols = list(
        prob_resident_base = NA_real_,
        prob_nz_migrant_base = NA_real_,
        prob_au_migrant_base = NA_real_,
        prob_resident = NA_real_,
        prob_nz_migrant = NA_real_,
        prob_au_migrant = NA_real_,
        top_class = NA_character_,
        class_support = NA_real_,
        indeterminate = NA,
        operational_call = NA_character_,
        switch_rate = NA_real_,
        top_prob_mean = NA_real_
      )
    )
  }

  bind_rows(available_out, unavailable_out) %>%
    arrange(source_order)
}

unknown_predict_stage4b_r2s2 <- function(
  unknown_data,
  model_obj,
  measurement_spec,
  measurement_cols,
  class_levels,
  n_mc_draws,
  support_threshold,
  top_prob_threshold
) {
  available_mask <- complete_measurements(unknown_data, measurement_cols)
  available_data <- unknown_data[available_mask, , drop = FALSE]
  unavailable_data <- unknown_data[!available_mask, , drop = FALSE]

  available_out <- tibble()
  if (nrow(available_data) > 0) {
    base_pred <- predict_stage4b_r2s2(model_obj, available_data)
    base_tbl <- bind_cols(
      tibble(
        truth = factor(rep(NA_character_, nrow(available_data)), levels = class_levels)
      ),
      base_pred[, c(class_levels, "base_class"), drop = FALSE]
    )

    draw_prob_array <- array(
      NA_real_,
      dim = c(nrow(available_data), length(class_levels), n_mc_draws),
      dimnames = list(NULL, class_levels, NULL)
    )

    for (draw_idx in seq_len(n_mc_draws)) {
      perturbed_data <- perturb_stage3_measurements(
        data = available_data,
        measurement_spec = measurement_spec,
        measurement_cols = measurement_cols,
        seed = 2300000L + draw_idx
      )

      draw_pred <- predict_stage4b_r2s2(model_obj, perturbed_data)
      draw_prob_array[, , draw_idx] <- as.matrix(draw_pred[, class_levels, drop = FALSE])
    }

    uncertainty_tbl <- summarise_stage4_uncertainty(
      base_prediction_tbl = base_tbl,
      draw_prob_array = draw_prob_array,
      class_levels = class_levels,
      support_threshold = support_threshold,
      top_prob_threshold = top_prob_threshold
    )

    available_out <- bind_cols(
      available_data %>%
        transmute(
          source_order,
          ring,
          longitude_capture,
          latitude_capture,
          coord_group,
          prediction_available = TRUE,
          availability_note = NA_character_
        ),
      tibble(
        prob_resident_base = base_pred$resident,
        prob_nz_migrant_base = base_pred$`NZ migrant`,
        prob_au_migrant_base = base_pred$`AU migrant`,
        prob_resident = uncertainty_tbl$uncertainty_prob_resident,
        prob_nz_migrant = uncertainty_tbl$uncertainty_prob_NZ_migrant,
        prob_au_migrant = uncertainty_tbl$uncertainty_prob_AU_migrant,
        top_class = as.character(uncertainty_tbl$uncertainty_class),
        class_support = uncertainty_tbl$class_support,
        indeterminate = uncertainty_tbl$indeterminate,
        operational_call = as.character(uncertainty_tbl$operational_call),
        switch_rate = uncertainty_tbl$switch_rate_from_base,
        top_prob_mean = uncertainty_tbl$top_prob_mean,
        r1_prob_migrant_base = base_pred$r1_prob_migrant,
        p_au_given_migrant_base = base_pred$p_au_given_migrant
      )
    )
  }

  unavailable_out <- tibble()
  if (nrow(unavailable_data) > 0) {
    unavailable_out <- build_unavailable_rows(
      data = unavailable_data,
      availability_note = paste(
        "Missing required inputs for staged R2 S2 workflow:",
        paste(measurement_cols, collapse = ", ")
      ),
      extra_cols = list(
        prob_resident_base = NA_real_,
        prob_nz_migrant_base = NA_real_,
        prob_au_migrant_base = NA_real_,
        prob_resident = NA_real_,
        prob_nz_migrant = NA_real_,
        prob_au_migrant = NA_real_,
        top_class = NA_character_,
        class_support = NA_real_,
        indeterminate = NA,
        operational_call = NA_character_,
        switch_rate = NA_real_,
        top_prob_mean = NA_real_,
        r1_prob_migrant_base = NA_real_,
        p_au_given_migrant_base = NA_real_
      )
    )
  }

  bind_rows(available_out, unavailable_out) %>%
    arrange(source_order)
}

labelled_data <- read_csv(
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

unknown_data <- read_csv(
  file.path(derived_dir, "live_cn_paired_by_ring.csv"),
  show_col_types = FALSE
) %>%
  derive_stage2_features() %>%
  filter(status_known == FALSE) %>%
  mutate(source_order = row_number())

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
  file.path(derived_dir, "model_application_unknown_measurement_spec.csv")
)

n_mc_draws <- 500L
r1_support_threshold <- 0.8
r1_boundary_interval <- c(0.1, 0.9)
r2_support_threshold <- 0.8
r2_top_prob_threshold <- 0.5
r2_class_levels <- c("resident", "NZ migrant", "AU migrant")

write_csv(
  tibble(
    parameter = c(
      "n_mc_draws",
      "r1_support_threshold",
      "r1_boundary_interval_low",
      "r1_boundary_interval_high",
      "r2_support_threshold",
      "r2_top_prob_threshold"
    ),
    value = c(
      as.character(n_mc_draws),
      as.character(r1_support_threshold),
      as.character(r1_boundary_interval[[1]]),
      as.character(r1_boundary_interval[[2]]),
      as.character(r2_support_threshold),
      as.character(r2_top_prob_threshold)
    )
  ),
  file.path(derived_dir, "model_application_unknown_config.csv")
)

r1_candidate <- list(
  candidate_id = "r1_c_paired_contrast_cn",
  candidate_label = "Rung 1 C: Paired-contrast C/N",
  predictors = c("abs_delta_d13c", "abs_delta_d15n", "paired_distance_cn"),
  distance_vars = c("abs_delta_d13c", "abs_delta_d15n"),
  distance_name = "paired_distance_cn"
)

r2a_candidate <- list(
  candidate_id = "r2_a_breast_cn",
  candidate_label = "Rung 2 A: Breast-only C/N",
  predictors = c("normalised_d13c_Breast", "normalised_d15n_Breast"),
  distance_vars = character(0),
  distance_name = NA_character_
)

r2c_candidate <- list(
  candidate_id = "r2_c_breast_cnh_plus_paired_contrast",
  candidate_label = "Rung 2 C: Breast C/N/H + paired contrasts",
  predictors = c(
    "normalised_d13c_Breast",
    "normalised_d15n_Breast",
    "normalised_d2h_Breast",
    "abs_delta_d13c",
    "abs_delta_d15n",
    "abs_delta_d2h"
  ),
  distance_vars = character(0),
  distance_name = NA_character_
)

r2s2_migrant_candidate <- list(
  candidate_id = "r2s2_migrant_submodel_breast_cnh",
  predictors = c(
    "normalised_d13c_Breast",
    "normalised_d15n_Breast",
    "normalised_d2h_Breast"
  ),
  distance_vars = character(0),
  distance_name = NA_character_
)

r1_fit <- fit_stage3_ridge_binomial(
  train_data = labelled_data,
  candidate = r1_candidate,
  outcome_col = "status_rung1",
  seed = 910001L
)

if (!identical(r1_fit$status, "ok")) {
  stop("Full-data fit failed for R1 C.")
}

r2a_fit <- fit_stage4_ridge_multinomial(
  train_data = labelled_data,
  candidate = r2a_candidate,
  outcome_col = "status_rung2",
  seed = 910002L
)

if (!identical(r2a_fit$status, "ok")) {
  stop("Full-data fit failed for R2 A.")
}

r2c_fit <- fit_stage4_ridge_multinomial(
  train_data = labelled_data,
  candidate = r2c_candidate,
  outcome_col = "status_rung2",
  seed = 910003L
)

if (!identical(r2c_fit$status, "ok")) {
  stop("Full-data fit failed for R2 C.")
}

r2s2_fit <- fit_stage4b_r2s2(
  train_data = labelled_data,
  r1_candidate = r1_candidate,
  migrant_candidate = r2s2_migrant_candidate,
  seed_r1 = 910004L,
  seed_submodel = 910005L
)

if (!identical(r2s2_fit$status, "ok")) {
  stop("Full-data fit failed for R2 S2.")
}

r1_predictions <- unknown_predict_stage3(
  unknown_data = unknown_data,
  model_obj = r1_fit,
  measurement_spec = measurement_spec,
  measurement_cols = c(
    "normalised_d13c_Breast",
    "normalised_d13c_Primary",
    "normalised_d15n_Breast",
    "normalised_d15n_Primary"
  ),
  n_mc_draws = n_mc_draws,
  support_threshold = r1_support_threshold,
  boundary_interval = r1_boundary_interval
) %>%
  mutate(
    candidate_id = "r1_c_paired_contrast_cn",
    candidate_label = "Rung 1 C: Paired-contrast C/N"
  ) %>%
  relocate(candidate_id, candidate_label, .before = source_order)

write_csv(
  r1_predictions,
  file.path(derived_dir, "model_application_unknown_rung1_predictions.csv")
)

r2a_predictions <- unknown_predict_stage4_multinomial(
  unknown_data = unknown_data,
  model_obj = r2a_fit,
  measurement_spec = measurement_spec,
  measurement_cols = c("normalised_d13c_Breast", "normalised_d15n_Breast"),
  class_levels = r2_class_levels,
  n_mc_draws = n_mc_draws,
  support_threshold = r2_support_threshold,
  top_prob_threshold = r2_top_prob_threshold,
  availability_note_prefix = "Missing required inputs for R2 A breast-only C/N model:"
) %>%
  mutate(
    candidate_id = "r2_a_breast_cn",
    candidate_label = "Rung 2 A: Breast-only C/N"
  ) %>%
  relocate(candidate_id, candidate_label, .before = source_order)

write_csv(
  r2a_predictions,
  file.path(derived_dir, "model_application_unknown_r2a_predictions.csv")
)

r2c_predictions <- unknown_predict_stage4_multinomial(
  unknown_data = unknown_data,
  model_obj = r2c_fit,
  measurement_spec = measurement_spec,
  measurement_cols = c(
    "normalised_d13c_Breast",
    "normalised_d13c_Primary",
    "normalised_d15n_Breast",
    "normalised_d15n_Primary",
    "normalised_d2h_Breast",
    "normalised_d2h_Primary"
  ),
  class_levels = r2_class_levels,
  n_mc_draws = n_mc_draws,
  support_threshold = r2_support_threshold,
  top_prob_threshold = r2_top_prob_threshold,
  availability_note_prefix = "Missing required inputs for R2 C breast C/N/H plus paired contrasts model:"
) %>%
  mutate(
    candidate_id = "r2_c_breast_cnh_plus_paired_contrast",
    candidate_label = "Rung 2 C: Breast C/N/H + paired contrasts"
  ) %>%
  relocate(candidate_id, candidate_label, .before = source_order)

write_csv(
  r2c_predictions,
  file.path(derived_dir, "model_application_unknown_r2c_predictions.csv")
)

r2s2_predictions <- unknown_predict_stage4b_r2s2(
  unknown_data = unknown_data,
  model_obj = r2s2_fit,
  measurement_spec = measurement_spec,
  measurement_cols = c(
    "normalised_d13c_Breast",
    "normalised_d13c_Primary",
    "normalised_d15n_Breast",
    "normalised_d15n_Primary",
    "normalised_d2h_Breast"
  ),
  class_levels = r2_class_levels,
  n_mc_draws = n_mc_draws,
  support_threshold = r2_support_threshold,
  top_prob_threshold = r2_top_prob_threshold
) %>%
  mutate(
    candidate_id = "r2s_b_soft_hierarchical_r1c_then_breast_cnh",
    candidate_label = "Rung 2 S2: Soft hierarchical R1 C then breast C/N/H NZ-vs-AU"
  ) %>%
  relocate(candidate_id, candidate_label, .before = source_order)

write_csv(
  r2s2_predictions,
  file.path(derived_dir, "model_application_unknown_r2s2_predictions.csv")
)

unknown_meta <- unknown_data %>%
  transmute(
    source_order,
    ring,
    longitude_capture,
    latitude_capture,
    coord_group,
    status_known,
    has_complete_cn_pair,
    has_complete_cnh_pair
  )

combined_predictions <- unknown_meta %>%
  left_join(
    r1_predictions %>%
      select(
        source_order,
        prediction_available,
        availability_note,
        prob_resident,
        prob_migrant,
        top_class,
        class_support,
        indeterminate,
        operational_call,
        switch_rate
      ) %>%
      rename_with(~ paste0("r1_", .x), -source_order),
    by = "source_order"
  ) %>%
  left_join(
    r2a_predictions %>%
      select(
        source_order,
        prediction_available,
        availability_note,
        prob_resident,
        prob_nz_migrant,
        prob_au_migrant,
        top_class,
        class_support,
        indeterminate,
        operational_call,
        switch_rate
      ) %>%
      rename_with(~ paste0("r2a_", .x), -source_order),
    by = "source_order"
  ) %>%
  left_join(
    r2c_predictions %>%
      select(
        source_order,
        prediction_available,
        availability_note,
        prob_resident,
        prob_nz_migrant,
        prob_au_migrant,
        top_class,
        class_support,
        indeterminate,
        operational_call,
        switch_rate
      ) %>%
      rename_with(~ paste0("r2c_", .x), -source_order),
    by = "source_order"
  ) %>%
  left_join(
    r2s2_predictions %>%
      select(
        source_order,
        prediction_available,
        availability_note,
        prob_resident,
        prob_nz_migrant,
        prob_au_migrant,
        top_class,
        class_support,
        indeterminate,
        operational_call,
        switch_rate
      ) %>%
      rename_with(~ paste0("r2s2_", .x), -source_order),
    by = "source_order"
  ) %>%
  rowwise() %>%
  mutate(
    informal_rung2_consensus = {
      calls <- c(r2a_operational_call, r2c_operational_call, r2s2_operational_call)
      calls <- calls[!is.na(calls) & calls != "indeterminate"]
      if (length(calls) == 0) {
        "no determinate Rung 2 call"
      } else {
        tab <- sort(table(calls), decreasing = TRUE)
        if (tab[[1]] >= 2) {
          paste("informal consensus:", names(tab)[1])
        } else if (length(unique(calls)) == 1) {
          paste("single-model lead:", unique(calls))
        } else {
          "models disagree"
        }
      }
    }
  ) %>%
  ungroup() %>%
  arrange(source_order)

write_csv(
  combined_predictions,
  file.path(derived_dir, "model_application_unknown_predictions_wide.csv")
)

probability_long <- bind_rows(
  r1_predictions %>%
    transmute(
      source_order,
      ring,
      model_id = "r1_c",
      model_label = "Rung 1: R1 C paired-contrast C/N",
      prediction_available,
      class = "resident",
      probability = prob_resident,
      operational_call
    ),
  r1_predictions %>%
    transmute(
      source_order,
      ring,
      model_id = "r1_c",
      model_label = "Rung 1: R1 C paired-contrast C/N",
      prediction_available,
      class = "migrant",
      probability = prob_migrant,
      operational_call
    ),
  r2a_predictions %>%
    transmute(
      source_order,
      ring,
      model_id = "r2_a",
      model_label = "Rung 2: R2 A breast-only C/N",
      prediction_available,
      class = "resident",
      probability = prob_resident,
      operational_call
    ),
  r2a_predictions %>%
    transmute(
      source_order,
      ring,
      model_id = "r2_a",
      model_label = "Rung 2: R2 A breast-only C/N",
      prediction_available,
      class = "NZ migrant",
      probability = prob_nz_migrant,
      operational_call
    ),
  r2a_predictions %>%
    transmute(
      source_order,
      ring,
      model_id = "r2_a",
      model_label = "Rung 2: R2 A breast-only C/N",
      prediction_available,
      class = "AU migrant",
      probability = prob_au_migrant,
      operational_call
    ),
  r2c_predictions %>%
    transmute(
      source_order,
      ring,
      model_id = "r2_c",
      model_label = "Rung 2: R2 C direct breast C/N/H + paired contrasts",
      prediction_available,
      class = "resident",
      probability = prob_resident,
      operational_call
    ),
  r2c_predictions %>%
    transmute(
      source_order,
      ring,
      model_id = "r2_c",
      model_label = "Rung 2: R2 C direct breast C/N/H + paired contrasts",
      prediction_available,
      class = "NZ migrant",
      probability = prob_nz_migrant,
      operational_call
    ),
  r2c_predictions %>%
    transmute(
      source_order,
      ring,
      model_id = "r2_c",
      model_label = "Rung 2: R2 C direct breast C/N/H + paired contrasts",
      prediction_available,
      class = "AU migrant",
      probability = prob_au_migrant,
      operational_call
    ),
  r2s2_predictions %>%
    transmute(
      source_order,
      ring,
      model_id = "r2_s2",
      model_label = "Rung 2: R2 S2 staged hierarchical approach",
      prediction_available,
      class = "resident",
      probability = prob_resident,
      operational_call
    ),
  r2s2_predictions %>%
    transmute(
      source_order,
      ring,
      model_id = "r2_s2",
      model_label = "Rung 2: R2 S2 staged hierarchical approach",
      prediction_available,
      class = "NZ migrant",
      probability = prob_nz_migrant,
      operational_call
    ),
  r2s2_predictions %>%
    transmute(
      source_order,
      ring,
      model_id = "r2_s2",
      model_label = "Rung 2: R2 S2 staged hierarchical approach",
      prediction_available,
      class = "AU migrant",
      probability = prob_au_migrant,
      operational_call
    )
) %>%
  arrange(source_order, model_id, class)

write_csv(
  probability_long,
  file.path(derived_dir, "model_application_unknown_probability_long.csv")
)

summary_lines <- c(
  "# Unknown-status Live Bird Prediction Summary",
  "",
  "These outputs are model-based predictions for the 18 live birds with unknown migratory phenotype.",
  "They are application-layer predictions only and were not used in the labelled evaluation set.",
  "",
  "## Availability",
  "",
  paste0(
    "- R1 C predictions available for ",
    sum(r1_predictions$prediction_available, na.rm = TRUE),
    " of ",
    nrow(r1_predictions),
    " birds."
  ),
  paste0(
    "- R2 A predictions available for ",
    sum(r2a_predictions$prediction_available, na.rm = TRUE),
    " of ",
    nrow(r2a_predictions),
    " birds."
  ),
  paste0(
    "- R2 C predictions available for ",
    sum(r2c_predictions$prediction_available, na.rm = TRUE),
    " of ",
    nrow(r2c_predictions),
    " birds."
  ),
  paste0(
    "- R2 S2 predictions available for ",
    sum(r2s2_predictions$prediction_available, na.rm = TRUE),
    " of ",
    nrow(r2s2_predictions),
    " birds."
  ),
  "",
  "## Notes",
  "",
  "- R1 and Rung 2 calls use the same accepted indeterminate logic as the labelled-bird analyses.",
  "- CP19938 lacks the paired inputs needed for R1 C, R2 C, and R2 S2, but still carries enough breast C/N data for R2 A.",
  "- These outputs are useful as biological leads, not as validation evidence."
)

writeLines(
  summary_lines,
  file.path(derived_dir, "model_application_unknown_summary.md")
)
