suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(glmnet)
})

stage3_measurement_column_map <- function() {
  tibble(
    measurement_col = c(
      "normalised_d13c_Breast",
      "normalised_d15n_Breast",
      "normalised_d2h_Breast",
      "normalised_d13c_Primary",
      "normalised_d15n_Primary",
      "normalised_d2h_Primary"
    ),
    feather_type = c("Breast", "Breast", "Breast", "Primary", "Primary", "Primary"),
    isotope = c(
      "normalised_d13c",
      "normalised_d15n",
      "normalised_d2h",
      "normalised_d13c",
      "normalised_d15n",
      "normalised_d2h"
    )
  )
}

build_stage3_measurement_spec <- function(stage0_variance, stage1_adequacy) {
  stage3_measurement_column_map() %>%
    left_join(
      stage0_variance %>%
        select(
          feather_type,
          isotope,
          direct_estimate_available,
          working_variance,
          working_sd,
          working_source,
          working_notes
        ),
      by = c("feather_type", "isotope")
    ) %>%
    left_join(
      stage1_adequacy %>%
        select(
          feather_type,
          isotope,
          adequacy_for_stage2,
          adequacy_notes
        ),
      by = c("feather_type", "isotope")
    ) %>%
    mutate(
      measurement_support = dplyr::case_when(
        feather_type == "Primary" & isotope == "normalised_d2h" ~
          "primary_h_proxy_variance_high_caution",
        feather_type == "Breast" & isotope == "normalised_d2h" ~
          "breast_h_stage1_repeatability_plus_weak_stage0_assay_support",
        TRUE ~ "direct_stage0_working_variance"
      )
    )
}

prep_stage3_candidate <- function(train_data, candidate) {
  train <- train_data

  if (length(candidate$distance_vars) > 0) {
    distance_means <- vapply(
      candidate$distance_vars,
      function(var) mean(train[[var]], na.rm = TRUE),
      numeric(1)
    )
    distance_sds <- vapply(
      candidate$distance_vars,
      function(var) stats::sd(train[[var]], na.rm = TRUE),
      numeric(1)
    )
    distance_sds[is.na(distance_sds) | distance_sds == 0] <- 1

    scaled_distance <- sweep(
      as.matrix(train[, candidate$distance_vars, drop = FALSE]),
      2,
      distance_means,
      "-"
    )
    scaled_distance <- sweep(scaled_distance, 2, distance_sds, "/")
    train[[candidate$distance_name]] <- sqrt(rowSums(scaled_distance^2))
  } else {
    distance_means <- numeric(0)
    distance_sds <- numeric(0)
  }

  predictor_names <- candidate$predictors
  train_x_raw <- as.matrix(train[, predictor_names, drop = FALSE])

  predictor_means <- colMeans(train_x_raw, na.rm = TRUE)
  predictor_sds <- apply(train_x_raw, 2, stats::sd, na.rm = TRUE)
  predictor_sds[is.na(predictor_sds) | predictor_sds == 0] <- 1

  train_x <- sweep(train_x_raw, 2, predictor_means, "-")
  train_x <- sweep(train_x, 2, predictor_sds, "/")

  keep_mask <- apply(train_x, 2, function(x) stats::sd(x, na.rm = TRUE) > 0)
  if (!any(keep_mask)) {
    keep_mask <- rep(TRUE, length(predictor_names))
  }

  list(
    predictor_names = predictor_names,
    predictor_means = predictor_means,
    predictor_sds = predictor_sds,
    kept_predictors = predictor_names[keep_mask],
    keep_mask = keep_mask,
    distance_vars = candidate$distance_vars,
    distance_name = candidate$distance_name,
    distance_means = distance_means,
    distance_sds = distance_sds
  )
}

transform_stage3_candidate <- function(data, prep) {
  transformed <- data

  if (length(prep$distance_vars) > 0) {
    scaled_distance <- sweep(
      as.matrix(transformed[, prep$distance_vars, drop = FALSE]),
      2,
      prep$distance_means,
      "-"
    )
    scaled_distance <- sweep(scaled_distance, 2, prep$distance_sds, "/")
    transformed[[prep$distance_name]] <- sqrt(rowSums(scaled_distance^2))
  }

  x_raw <- as.matrix(transformed[, prep$predictor_names, drop = FALSE])
  x <- sweep(x_raw, 2, prep$predictor_means, "-")
  x <- sweep(x, 2, prep$predictor_sds, "/")
  x <- x[, prep$keep_mask, drop = FALSE]

  list(
    data = transformed,
    x = x
  )
}

fit_stage3_ridge_binomial <- function(train_data, candidate, outcome_col, seed) {
  prep <- prep_stage3_candidate(train_data = train_data, candidate = candidate)
  train_x <- transform_stage3_candidate(train_data, prep)$x

  if (ncol(train_x) == 0) {
    return(list(status = "no_predictors"))
  }

  y_train <- train_data[[outcome_col]]
  if (length(unique(y_train)) < length(levels(y_train))) {
    return(list(status = "missing_training_class"))
  }

  y_train_num <- ifelse(y_train == levels(y_train)[2], 1, 0)

  set.seed(seed)
  fit <- tryCatch(
    cv.glmnet(
      x = train_x,
      y = y_train_num,
      family = "binomial",
      alpha = 0,
      nfolds = choose_inner_folds(y_train),
      type.measure = "deviance",
      standardize = FALSE
    ),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    return(list(status = "fit_error", message = fit$message))
  }

  list(
    status = "ok",
    fit = fit,
    prep = prep,
    class_levels = levels(y_train)
  )
}

predict_stage3_ridge_binomial <- function(model_obj, new_data) {
  transformed <- transform_stage3_candidate(new_data, model_obj$prep)
  prob_pos <- as.numeric(
    predict(
      model_obj$fit,
      newx = transformed$x,
      s = "lambda.1se",
      type = "response"
    )
  )

  class_levels <- model_obj$class_levels
  pred_class <- ifelse(prob_pos >= 0.5, class_levels[2], class_levels[1])

  tibble(
    .pred_resident = 1 - prob_pos,
    .pred_migrant = prob_pos,
    base_class = factor(pred_class, levels = class_levels)
  )
}

perturb_stage3_measurements <- function(data, measurement_spec, measurement_cols, seed) {
  perturbed <- data

  for (col_name in measurement_cols) {
    col_spec <- measurement_spec %>%
      filter(measurement_col == !!col_name)

    if (nrow(col_spec) != 1 || is.na(col_spec$working_sd[[1]])) {
      stop(paste("No usable working SD found for", col_name))
    }

    set.seed(seed + match(col_name, measurement_cols) * 1000L)
    perturbed[[col_name]] <- perturbed[[col_name]] +
      stats::rnorm(nrow(perturbed), mean = 0, sd = col_spec$working_sd[[1]])
  }

  derive_stage2_features(perturbed)
}

summarise_stage3_uncertainty <- function(
  base_prediction_tbl,
  draw_prob_matrix,
  support_threshold = 0.8,
  boundary_interval = c(0.1, 0.9)
) {
  prob_mean <- rowMeans(draw_prob_matrix)
  prob_sd <- apply(draw_prob_matrix, 1, stats::sd)
  prob_q_low <- apply(draw_prob_matrix, 1, stats::quantile, probs = boundary_interval[1], names = FALSE)
  prob_q_high <- apply(draw_prob_matrix, 1, stats::quantile, probs = boundary_interval[2], names = FALSE)

  prop_migrant_draws <- rowMeans(draw_prob_matrix >= 0.5)
  uncertainty_class <- ifelse(prob_mean >= 0.5, "migrant", "resident")
  class_support <- ifelse(
    uncertainty_class == "migrant",
    prop_migrant_draws,
    1 - prop_migrant_draws
  )

  boundary_crosses <- prob_q_low < 0.5 & prob_q_high > 0.5
  indeterminate <- class_support < support_threshold | boundary_crosses

  base_class_chr <- as.character(base_prediction_tbl$base_class)
  top_prob_base <- pmax(base_prediction_tbl$.pred_migrant, base_prediction_tbl$.pred_resident)
  top_prob_mean <- rowMeans(pmax(draw_prob_matrix, 1 - draw_prob_matrix))

  tibble(
    truth = factor(base_prediction_tbl$truth, levels = c("resident", "migrant")),
    base_prob_migrant = base_prediction_tbl$.pred_migrant,
    base_class = factor(base_class_chr, levels = c("resident", "migrant")),
    uncertainty_prob_migrant_mean = prob_mean,
    uncertainty_prob_migrant_sd = prob_sd,
    uncertainty_prob_migrant_q10 = prob_q_low,
    uncertainty_prob_migrant_q90 = prob_q_high,
    uncertainty_class = factor(uncertainty_class, levels = c("resident", "migrant")),
    class_support = class_support,
    prop_migrant_draws = prop_migrant_draws,
    any_draw_class_flip = ifelse(
      base_class_chr == "migrant",
      prop_migrant_draws < 1,
      prop_migrant_draws > 0
    ),
    switch_rate_from_base = ifelse(
      base_class_chr == "migrant",
      1 - prop_migrant_draws,
      prop_migrant_draws
    ),
    top_prob_base = top_prob_base,
    top_prob_mean = top_prob_mean,
    top_prob_abs_shift = abs(top_prob_mean - top_prob_base),
    boundary_crosses = boundary_crosses,
    indeterminate = indeterminate,
    operational_call = factor(
      ifelse(indeterminate, "indeterminate", uncertainty_class),
      levels = c("resident", "migrant", "indeterminate")
    ),
    final_class_changed = !indeterminate & uncertainty_class != base_class_chr
  )
}

summarise_stage3_split_metrics <- function(assessment_summary) {
  truth <- factor(assessment_summary$truth, levels = c("resident", "migrant"))
  base_class <- factor(assessment_summary$base_class, levels = c("resident", "migrant"))
  uncertainty_class <- factor(assessment_summary$uncertainty_class, levels = c("resident", "migrant"))

  determinate_idx <- !assessment_summary$indeterminate

  tibble(
    n_assessment = nrow(assessment_summary),
    base_log_loss = binary_log_loss(truth, assessment_summary$base_prob_migrant, positive_class = "migrant"),
    uncertainty_log_loss = binary_log_loss(truth, assessment_summary$uncertainty_prob_migrant_mean, positive_class = "migrant"),
    base_bal_accuracy = macro_bal_accuracy(truth, base_class),
    uncertainty_bal_accuracy = macro_bal_accuracy(truth, uncertainty_class),
    determinate_rate = mean(!assessment_summary$indeterminate),
    indeterminate_rate = mean(assessment_summary$indeterminate),
    final_class_changed_rate = mean(assessment_summary$final_class_changed),
    any_flip_rate = mean(assessment_summary$any_draw_class_flip),
    base_to_indeterminate_rate = mean(assessment_summary$indeterminate),
    boundary_cross_rate = mean(assessment_summary$boundary_crosses),
    mean_switch_rate = mean(assessment_summary$switch_rate_from_base),
    median_switch_rate = stats::median(assessment_summary$switch_rate_from_base),
    mean_abs_top_prob_shift = mean(assessment_summary$top_prob_abs_shift),
    top_prob_shift_gt_0_05_rate = mean(assessment_summary$top_prob_abs_shift > 0.05)
  )
}

summarise_stage3_model_metrics <- function(split_summary_tbl) {
  split_summary_tbl %>%
    group_by(candidate_id, candidate_label, validation_layer, caution_level, sensitivity_only) %>%
    summarise(
      n_splits = n(),
      mean_base_log_loss = mean(base_log_loss, na.rm = TRUE),
      sd_base_log_loss = stats::sd(base_log_loss, na.rm = TRUE),
      mean_uncertainty_log_loss = mean(uncertainty_log_loss, na.rm = TRUE),
      sd_uncertainty_log_loss = stats::sd(uncertainty_log_loss, na.rm = TRUE),
      mean_base_bal_accuracy = mean(base_bal_accuracy, na.rm = TRUE),
      mean_uncertainty_bal_accuracy = mean(uncertainty_bal_accuracy, na.rm = TRUE),
      mean_determinate_rate = mean(determinate_rate, na.rm = TRUE),
      mean_indeterminate_rate = mean(indeterminate_rate, na.rm = TRUE),
      mean_final_class_changed_rate = mean(final_class_changed_rate, na.rm = TRUE),
      mean_any_flip_rate = mean(any_flip_rate, na.rm = TRUE),
      mean_base_to_indeterminate_rate = mean(base_to_indeterminate_rate, na.rm = TRUE),
      mean_boundary_cross_rate = mean(boundary_cross_rate, na.rm = TRUE),
      mean_switch_rate = mean(mean_switch_rate, na.rm = TRUE),
      mean_abs_top_prob_shift = mean(mean_abs_top_prob_shift, na.rm = TRUE),
      mean_top_prob_shift_gt_0_05_rate = mean(top_prob_shift_gt_0_05_rate, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(validation_layer, mean_uncertainty_log_loss, desc(mean_uncertainty_bal_accuracy))
}
