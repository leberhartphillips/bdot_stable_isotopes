suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(glmnet)
})

stage4_safe_class_name <- function(x) {
  gsub("_+$", "", gsub("^_+", "", gsub("[^A-Za-z0-9]+", "_", x)))
}

stage4_prob_col_map <- function(prefix, class_levels) {
  stats::setNames(
    paste0(prefix, vapply(class_levels, stage4_safe_class_name, character(1))),
    class_levels
  )
}

fit_stage4_ridge_multinomial <- function(train_data, candidate, outcome_col, seed) {
  prep <- prep_stage3_candidate(train_data = train_data, candidate = candidate)
  train_x <- transform_stage3_candidate(train_data, prep)$x

  if (ncol(train_x) == 0) {
    return(list(status = "no_predictors"))
  }

  y_train <- train_data[[outcome_col]]
  if (length(unique(y_train)) < length(levels(y_train))) {
    return(list(status = "missing_training_class"))
  }

  set.seed(seed)
  fit <- tryCatch(
    cv.glmnet(
      x = train_x,
      y = y_train,
      family = "multinomial",
      alpha = 0,
      nfolds = choose_inner_folds(y_train),
      type.measure = "deviance",
      standardize = FALSE,
      type.multinomial = "grouped"
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

predict_stage4_ridge_multinomial <- function(model_obj, new_data) {
  transformed <- transform_stage3_candidate(new_data, model_obj$prep)
  pred_arr <- predict(
    model_obj$fit,
    newx = transformed$x,
    s = "lambda.1se",
    type = "response"
  )
  pred_arr <- drop(pred_arr)

  if (is.null(dim(pred_arr))) {
    pred_arr <- matrix(pred_arr, nrow = 1)
  }

  if (nrow(pred_arr) != nrow(new_data)) {
    pred_arr <- t(pred_arr)
  }

  class_levels <- model_obj$class_levels
  colnames(pred_arr) <- class_levels

  pred_class <- factor(
    class_levels[max.col(pred_arr, ties.method = "first")],
    levels = class_levels
  )

  prob_tbl <- as_tibble(pred_arr, .name_repair = "minimal")
  prob_tbl$base_class <- pred_class
  prob_tbl
}

summarise_stage4_uncertainty <- function(
  base_prediction_tbl,
  draw_prob_array,
  class_levels,
  support_threshold = 0.8,
  top_prob_threshold = 0.5
) {
  n_assessment <- dim(draw_prob_array)[1]
  n_draws <- dim(draw_prob_array)[3]

  mean_prob_mat <- apply(draw_prob_array, c(1, 2), mean)
  if (is.null(dim(mean_prob_mat))) {
    mean_prob_mat <- matrix(mean_prob_mat, nrow = n_assessment)
  }
  colnames(mean_prob_mat) <- class_levels

  base_prob_mat <- as.matrix(base_prediction_tbl[, class_levels, drop = FALSE])
  base_class_chr <- as.character(base_prediction_tbl$base_class)
  top_prob_base <- apply(base_prob_mat, 1, max)

  uncertainty_rows <- lapply(
    seq_len(n_assessment),
    function(i) {
      draw_prob_mat <- t(draw_prob_array[i, , , drop = FALSE][1, , ])
      if (is.null(dim(draw_prob_mat))) {
        draw_prob_mat <- matrix(draw_prob_mat, nrow = 1)
      }
      colnames(draw_prob_mat) <- class_levels

      draw_winners <- class_levels[max.col(draw_prob_mat, ties.method = "first")]
      winner_counts <- table(factor(draw_winners, levels = class_levels))
      modal_class <- names(winner_counts)[which.max(winner_counts)]
      class_support <- max(winner_counts) / n_draws

      modal_prob_draws <- draw_prob_mat[, modal_class]
      top_prob_mean <- mean(modal_prob_draws)
      top_prob_q10 <- stats::quantile(modal_prob_draws, probs = 0.1, names = FALSE)
      top_prob_q90 <- stats::quantile(modal_prob_draws, probs = 0.9, names = FALSE)

      indeterminate <- class_support < support_threshold || top_prob_mean < top_prob_threshold

      tibble(
        uncertainty_class = factor(modal_class, levels = class_levels),
        class_support = class_support,
        top_class_mean_prob = top_prob_mean,
        top_class_prob_q10 = top_prob_q10,
        top_class_prob_q90 = top_prob_q90,
        any_draw_class_flip = any(draw_winners != base_class_chr[i]),
        switch_rate_from_base = mean(draw_winners != base_class_chr[i]),
        top_prob_base = top_prob_base[i],
        top_prob_mean = top_prob_mean,
        top_prob_abs_shift = abs(top_prob_mean - top_prob_base[i]),
        indeterminate = indeterminate,
        operational_call = factor(
          ifelse(indeterminate, "indeterminate", modal_class),
          levels = c(class_levels, "indeterminate")
        ),
        final_class_changed = !indeterminate && modal_class != base_class_chr[i]
      )
    }
  )

  mean_prob_tbl <- as_tibble(mean_prob_mat, .name_repair = "minimal")
  mean_prob_cols <- stage4_prob_col_map("uncertainty_prob_", class_levels)
  colnames(mean_prob_tbl) <- unname(mean_prob_cols[class_levels])

  bind_cols(
    tibble(
      truth = factor(base_prediction_tbl$truth, levels = class_levels),
      base_class = factor(base_class_chr, levels = class_levels)
    ),
    mean_prob_tbl,
    bind_rows(uncertainty_rows)
  )
}

summarise_stage4_split_metrics <- function(assessment_summary, class_levels) {
  truth <- factor(assessment_summary$truth, levels = class_levels)
  base_class <- factor(assessment_summary$base_class, levels = class_levels)
  uncertainty_class <- factor(assessment_summary$uncertainty_class, levels = class_levels)

  base_prob_mat <- as.matrix(assessment_summary[, class_levels, drop = FALSE])

  uncertainty_prob_cols <- unname(stage4_prob_col_map("uncertainty_prob_", class_levels)[class_levels])
  uncertainty_prob_mat <- as.matrix(assessment_summary[, uncertainty_prob_cols, drop = FALSE])
  colnames(uncertainty_prob_mat) <- class_levels

  tibble(
    n_assessment = nrow(assessment_summary),
    base_log_loss = multiclass_log_loss(truth, base_prob_mat),
    uncertainty_log_loss = multiclass_log_loss(truth, uncertainty_prob_mat),
    base_bal_accuracy = macro_bal_accuracy(truth, base_class),
    uncertainty_bal_accuracy = macro_bal_accuracy(truth, uncertainty_class),
    determinate_rate = mean(!assessment_summary$indeterminate),
    indeterminate_rate = mean(assessment_summary$indeterminate),
    final_class_changed_rate = mean(assessment_summary$final_class_changed),
    any_flip_rate = mean(assessment_summary$any_draw_class_flip),
    base_to_indeterminate_rate = mean(assessment_summary$indeterminate),
    mean_switch_rate = mean(assessment_summary$switch_rate_from_base),
    median_switch_rate = stats::median(assessment_summary$switch_rate_from_base),
    mean_class_support = mean(assessment_summary$class_support),
    mean_top_prob_shift = mean(assessment_summary$top_prob_abs_shift, na.rm = TRUE),
    top_prob_shift_gt_0_05_rate = mean(assessment_summary$top_prob_abs_shift > 0.05, na.rm = TRUE)
  )
}

summarise_stage4_model_metrics <- function(split_summary_tbl) {
  split_summary_tbl %>%
    group_by(candidate_id, candidate_label, validation_layer, caution_level) %>%
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
      mean_switch_rate = mean(mean_switch_rate, na.rm = TRUE),
      mean_class_support = mean(mean_class_support, na.rm = TRUE),
      mean_top_prob_shift = mean(mean_top_prob_shift, na.rm = TRUE),
      mean_top_prob_shift_gt_0_05_rate = mean(top_prob_shift_gt_0_05_rate, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(validation_layer, mean_uncertainty_log_loss, desc(mean_uncertainty_bal_accuracy))
}

summarise_stage4_calibration <- function(assessment_summary, class_levels) {
  grouped <- split(
    assessment_summary,
    interaction(assessment_summary$candidate_id, assessment_summary$validation_layer, drop = TRUE)
  )

  bind_rows(
    lapply(
      grouped,
      function(df) {
        cal_input <- tibble(truth = factor(df$truth, levels = class_levels))
        for (class_name in class_levels) {
          cal_input[[class_name]] <- df[[stage4_prob_col_map("uncertainty_prob_", class_levels)[[class_name]]]]
        }

        cal <- multiclass_calibration_summary(cal_input, class_levels = class_levels)
        base_info <- tibble(
          candidate_id = df$candidate_id[[1]],
          candidate_label = df$candidate_label[[1]],
          validation_layer = df$validation_layer[[1]],
          caution_level = df$caution_level[[1]]
        )

        bind_cols(
          base_info[rep(1, nrow(cal)), ],
          cal
        )
      }
    )
  )
}
