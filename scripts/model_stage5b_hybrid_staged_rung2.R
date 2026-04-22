#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(tidyr)
})

source(file.path("R", "model_stage_repeatability_helpers.R"))
source(file.path("R", "model_stage2_helpers.R"))
source(file.path("R", "model_stage2b_rung2_staged_helpers.R"))
source(file.path("R", "model_stage3_helpers.R"))
source(file.path("R", "model_stage4b_rung2_staged_helpers.R"))
source(file.path("R", "model_stage5_siteaware_helpers.R"))

derived_dir <- file.path("data", "derived")
ensure_dir(derived_dir)

format_num <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits))
}

class_levels_r2 <- c("resident", "NZ migrant", "AU migrant")

evaluate_hybrid_staged_split <- function(train_data, test_data, upstream_candidate, migrant_candidate, seed_r1, seed_sub) {
  fit <- fit_stage4b_r2s2(
    train_data = train_data,
    r1_candidate = upstream_candidate,
    migrant_candidate = migrant_candidate,
    seed_r1 = seed_r1,
    seed_submodel = seed_sub
  )

  if (!identical(fit$status, "ok")) {
    return(list(status = fit$status, message = fit$detail$message %||% NA_character_))
  }

  pred <- predict_stage4b_r2s2(fit, test_data) %>%
    bind_cols(
      test_data %>%
        select(ring, coord_group, longitude_capture, latitude_capture, truth = status_rung2)
    ) %>%
    relocate(ring, coord_group, longitude_capture, latitude_capture, truth) %>%
    rename(.pred_class = base_class) %>%
    mutate(
      top_prob = pmax(resident, `NZ migrant`, `AU migrant`),
      weak_support = top_prob < 0.5
    )

  metrics <- stage2b_prediction_metrics(pred, class_levels_r2) %>%
    rename(
      log_loss = log_loss,
      balanced_accuracy = balanced_accuracy,
      accuracy = accuracy,
      brier = brier,
      mean_top_prob = mean_top_prob,
      weak_support_rate = weak_support_rate
    )

  cal <- multiclass_calibration_summary(
    predictions = pred %>%
      select(truth, all_of(class_levels_r2)),
    class_levels = class_levels_r2
  )

  metrics <- metrics %>%
    mutate(
      ece = cal %>% filter(metric == "confidence_ece", class == "overall") %>% pull(value)
    )

  list(
    status = "ok",
    predictions = pred,
    metrics = metrics,
    predictors_used = paste(c(upstream_candidate$predictors, migrant_candidate$predictors), collapse = "|")
  )
}

multiclass_pooled_summary <- function(prediction_tbl) {
  prob_mat <- as.matrix(prediction_tbl[, class_levels_r2, drop = FALSE])
  cal <- multiclass_calibration_summary(
    predictions = prediction_tbl %>%
      transmute(
        truth = factor(truth, levels = class_levels_r2),
        resident = resident,
        `NZ migrant` = `NZ migrant`,
        `AU migrant` = `AU migrant`
      ),
    class_levels = class_levels_r2
  )

  tibble(
    pooled_log_loss = multiclass_log_loss(prediction_tbl$truth, prob_mat),
    pooled_bal_accuracy = macro_bal_accuracy(
      factor(prediction_tbl$truth, levels = class_levels_r2),
      factor(prediction_tbl$.pred_class, levels = class_levels_r2)
    ),
    pooled_accuracy = mean(prediction_tbl$.pred_class == prediction_tbl$truth),
    pooled_brier = multiclass_brier(prediction_tbl$truth, prob_mat),
    pooled_ece = cal %>% filter(metric == "confidence_ece", class == "overall") %>% pull(value),
    pooled_mean_top_prob = mean(prediction_tbl$top_prob, na.rm = TRUE),
    pooled_weak_support_rate = mean(prediction_tbl$weak_support, na.rm = TRUE)
  )
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

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

stage5_baselines <- read_csv(
  file.path(derived_dir, "model_stage5_siteaware_model_summary.csv"),
  show_col_types = FALSE
)

fold_definitions <- bind_rows(site_block_folds, site_loso_folds) %>%
  arrange(validation_layer, split_id, ring)

candidate_def <- tibble(
  rung = "Rung 2",
  candidate_id = "r2h_b_hybrid_siteaware_r1_then_breast_cnh",
  candidate_label = "Rung 2 H: Site-aware R1 then breast C/N/H NZ-vs-AU",
  model_type = "staged_soft_hybrid",
  candidate_group = "exploratory_extension",
  caution_level = "siteaware_r1_primary_reference_then_non_siteaware_breast_cnh",
  upstream_stage_label = "Site-aware upstream R1: paired-contrast C/N + primary-site reference distance",
  downstream_stage_label = "Non-site-aware downstream migrant NZ-vs-AU model: breast C/N/H",
  uses_primary_site_reference_upstream = TRUE,
  uses_winter_reference_downstream = FALSE,
  uses_r1_upstream = TRUE,
  uses_breast_h = TRUE,
  uses_primary_h = FALSE
)

write_csv(
  candidate_def,
  file.path(derived_dir, "model_stage5b_hybrid_candidate_definition.csv")
)

upstream_candidate <- list(
  candidate_id = "r1_x_paired_contrast_cn_plus_primary_site_reference",
  candidate_label = "Site-aware upstream R1: paired-contrast C/N + primary-site reference distance",
  predictors = c("abs_delta_d13c", "abs_delta_d15n", "paired_distance_cn", "primary_ref_min_site_dist_cn"),
  distance_vars = c("abs_delta_d13c", "abs_delta_d15n"),
  distance_name = "paired_distance_cn"
)

downstream_candidate <- list(
  candidate_id = "r2h_migrant_submodel_breast_cnh",
  predictors = c(
    "normalised_d13c_Breast",
    "normalised_d15n_Breast",
    "normalised_d2h_Breast"
  ),
  distance_vars = character(0),
  distance_name = NA_character_
)

split_keys <- fold_definitions %>%
  distinct(validation_layer, split_id) %>%
  arrange(validation_layer, split_id)

split_metric_rows <- vector("list", nrow(split_keys))
prediction_rows <- vector("list", nrow(split_keys))

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
    add_primary_site_reference_features(primary_reference_train)

  test_aug <- test_data %>%
    add_primary_site_reference_features(primary_reference_train)

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
      rung2_all_classes_present
    )

  res <- evaluate_hybrid_staged_split(
    train_data = train_aug,
    test_data = test_aug,
    upstream_candidate = upstream_candidate,
    migrant_candidate = downstream_candidate,
    seed_r1 = 950000L + split_idx,
    seed_sub = 960000L + split_idx
  )

  split_metric_rows[[split_idx]] <- bind_cols(
    candidate_def %>%
      select(
        rung,
        candidate_id,
        candidate_label,
        model_type,
        candidate_group,
        caution_level,
        upstream_stage_label,
        downstream_stage_label,
        uses_primary_site_reference_upstream,
        uses_winter_reference_downstream,
        uses_r1_upstream,
        uses_breast_h,
        uses_primary_h
      ),
    split_context,
    tibble(
      n_analysis = nrow(train_data),
      primary_reference_rows = nrow(primary_reference_train),
      primary_reference_sites = n_distinct(primary_reference_train$sampling_site_label),
      fit_status = res$status
    ),
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
        mean_top_prob = NA_real_,
        weak_support_rate = NA_real_,
        ece = NA_real_,
        predictors_used = NA_character_
      )
    }
  )

  if (identical(res$status, "ok")) {
    prediction_rows[[split_idx]] <- bind_cols(
      candidate_def %>%
        select(rung, candidate_id, candidate_label, model_type, candidate_group),
      split_context %>% select(validation_layer, split_id),
      res$predictions
    )
  }
}

split_metrics <- bind_rows(split_metric_rows) %>%
  arrange(validation_layer, split_id)

assessment_predictions <- bind_rows(prediction_rows) %>%
  arrange(validation_layer, split_id, ring)

write_csv(
  split_metrics,
  file.path(derived_dir, "model_stage5b_hybrid_split_metrics.csv")
)

write_csv(
  assessment_predictions,
  file.path(derived_dir, "model_stage5b_hybrid_assessment_predictions.csv")
)

model_summary <- split_metrics %>%
  group_by(
    rung,
    candidate_id,
    candidate_label,
    model_type,
    candidate_group,
    caution_level,
    upstream_stage_label,
    downstream_stage_label,
    uses_primary_site_reference_upstream,
    uses_winter_reference_downstream,
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
    assessment_predictions %>%
      group_by(rung, candidate_id, candidate_label, model_type, candidate_group, validation_layer) %>%
      group_modify(~ multiclass_pooled_summary(.x)) %>%
      ungroup(),
    by = c("rung", "candidate_id", "candidate_label", "model_type", "candidate_group", "validation_layer")
  )

write_csv(
  model_summary,
  file.path(derived_dir, "model_stage5b_hybrid_model_summary.csv")
)

calibration_summary <- assessment_predictions %>%
  group_by(rung, candidate_id, candidate_label, validation_layer) %>%
  group_modify(~ multiclass_calibration_summary(
    predictions = .x %>%
      transmute(
        truth = factor(truth, levels = class_levels_r2),
        resident = resident,
        `NZ migrant` = `NZ migrant`,
        `AU migrant` = `AU migrant`
      ),
    class_levels = class_levels_r2
  )) %>%
  ungroup()

write_csv(
  calibration_summary,
  file.path(derived_dir, "model_stage5b_hybrid_calibration_summary.csv")
)

comparison_vs_siteaware <- bind_rows(
  stage5_baselines %>%
    filter(candidate_id %in% c(
      "r2_c_breast_cnh_plus_paired_contrast",
      "r2s_b_soft_hierarchical_r1c_then_breast_cnh"
    )) %>%
    transmute(
      candidate_id,
      candidate_label,
      validation_layer,
      pooled_log_loss,
      pooled_bal_accuracy,
      pooled_brier,
      pooled_ece,
      pooled_mean_top_prob,
      pooled_weak_support_rate
    ),
  model_summary %>%
    transmute(
      candidate_id,
      candidate_label,
      validation_layer,
      pooled_log_loss,
      pooled_bal_accuracy,
      pooled_brier,
      pooled_ece,
      pooled_mean_top_prob,
      pooled_weak_support_rate
    )
) %>%
  group_by(validation_layer) %>%
  mutate(
    pooled_log_loss_rank = min_rank(pooled_log_loss),
    pooled_bal_accuracy_rank = min_rank(desc(pooled_bal_accuracy))
  ) %>%
  ungroup() %>%
  arrange(validation_layer, pooled_log_loss_rank, pooled_bal_accuracy_rank)

write_csv(
  comparison_vs_siteaware,
  file.path(derived_dir, "model_stage5b_hybrid_vs_siteaware_coleaders.csv")
)

summary_lines <- c(
  "# Stage 5b Hybrid Staged Rung 2 Summary",
  "",
  "This exploratory candidate keeps all accepted and Stage 5 results intact.",
  "",
  "## Hybrid structure",
  "",
  "- Upstream stage: a site-aware Rung 1 model using paired-contrast C/N plus a training-fold-only primary-site reference distance.",
  "- Downstream stage: a non-site-aware migrant NZ-vs-AU breast C/N/H model.",
  "- Held-out assessment sites never contribute to the upstream primary-site reference cloud.",
  "",
  "## Pooled comparison",
  "",
  paste0(
    "- Grouped site-block pooled log loss: ",
    format_num(model_summary %>% filter(validation_layer == 'blocked_site_cv') %>% pull(pooled_log_loss))
  ),
  paste0(
    "- Grouped site-block pooled balanced accuracy: ",
    format_num(model_summary %>% filter(validation_layer == 'blocked_site_cv') %>% pull(pooled_bal_accuracy))
  ),
  paste0(
    "- LOSO pooled log loss: ",
    format_num(model_summary %>% filter(validation_layer == 'leave_one_site_out_cv') %>% pull(pooled_log_loss))
  ),
  paste0(
    "- LOSO pooled balanced accuracy: ",
    format_num(model_summary %>% filter(validation_layer == 'leave_one_site_out_cv') %>% pull(pooled_bal_accuracy))
  ),
  "",
  "See `model_stage5b_hybrid_vs_siteaware_coleaders.csv` for the direct comparison against `R2 S2` and `R2 C`."
)

writeLines(
  summary_lines,
  file.path(derived_dir, "model_stage5b_hybrid_summary.md")
)
