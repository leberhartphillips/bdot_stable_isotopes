suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

safe_divide <- function(numerator, denominator) {
  ifelse(is.na(denominator) | denominator == 0, NA_real_, numerator / denominator)
}

binary_au_confusion_metrics <- function(truth, pred_class, positive_class = "AU") {
  truth_chr <- as.character(truth)
  pred_chr <- as.character(pred_class)

  tp <- sum(truth_chr == positive_class & pred_chr == positive_class, na.rm = TRUE)
  fn <- sum(truth_chr == positive_class & pred_chr != positive_class, na.rm = TRUE)
  tn <- sum(truth_chr != positive_class & pred_chr != positive_class, na.rm = TRUE)
  fp <- sum(truth_chr != positive_class & pred_chr == positive_class, na.rm = TRUE)

  tibble(
    au_sensitivity = safe_divide(tp, tp + fn),
    au_specificity = safe_divide(tn, tn + fp),
    au_ppv = safe_divide(tp, tp + fp),
    au_npv = safe_divide(tn, tn + fn)
  )
}

augment_binary_predictions <- function(predictions, positive_class = "AU", indeterminate_threshold = 0.8) {
  out <- predictions

  if (!(".pred_positive" %in% colnames(out)) || !(".pred_negative" %in% colnames(out))) {
    stop("Binary prediction table must contain .pred_positive and .pred_negative.")
  }

  if (!(".pred_class" %in% colnames(out))) {
    out <- out %>%
      mutate(
        .pred_class = factor(
          ifelse(.pred_positive >= 0.5, positive_class, setdiff(levels(truth), positive_class)[1]),
          levels = levels(truth)
        )
      )
  }

  out %>%
    mutate(
      top_prob = pmax(.pred_positive, .pred_negative),
      weak_support = top_prob < 0.5,
      indeterminate = top_prob < indeterminate_threshold
    )
}

binary_au_prediction_metrics <- function(prediction_tbl, positive_class = "AU", indeterminate_threshold = 0.8) {
  pred <- augment_binary_predictions(
    predictions = prediction_tbl,
    positive_class = positive_class,
    indeterminate_threshold = indeterminate_threshold
  )

  conf <- binary_au_confusion_metrics(pred$truth, pred$.pred_class, positive_class = positive_class)

  tibble(
    log_loss = binary_log_loss(pred$truth, pred$.pred_positive, positive_class = positive_class),
    balanced_accuracy = macro_bal_accuracy(pred$truth, pred$.pred_class),
    accuracy = mean(pred$.pred_class == pred$truth),
    brier = binary_brier(pred$truth, pred$.pred_positive, positive_class = positive_class),
    mean_top_prob = mean(pred$top_prob, na.rm = TRUE),
    weak_support_rate = mean(pred$weak_support, na.rm = TRUE),
    indeterminate_rate = mean(pred$indeterminate, na.rm = TRUE)
  ) %>%
    bind_cols(conf)
}

binary_au_pooled_summary <- function(prediction_tbl, positive_class = "AU", indeterminate_threshold = 0.8) {
  pred <- augment_binary_predictions(
    predictions = prediction_tbl,
    positive_class = positive_class,
    indeterminate_threshold = indeterminate_threshold
  )

  cal <- binary_calibration_summary(pred, positive_class = positive_class)

  binary_au_prediction_metrics(
    prediction_tbl = pred,
    positive_class = positive_class,
    indeterminate_threshold = indeterminate_threshold
  ) %>%
    transmute(
      pooled_log_loss = log_loss,
      pooled_bal_accuracy = balanced_accuracy,
      pooled_accuracy = accuracy,
      pooled_brier = brier,
      pooled_mean_top_prob = mean_top_prob,
      pooled_weak_support_rate = weak_support_rate,
      pooled_indeterminate_rate = indeterminate_rate,
      pooled_au_sensitivity = au_sensitivity,
      pooled_au_specificity = au_specificity,
      pooled_au_ppv = au_ppv,
      pooled_au_npv = au_npv,
      pooled_ece = cal %>% filter(metric == "ece") %>% pull(value),
      pooled_calibration_intercept = cal %>% filter(metric == "calibration_intercept") %>% pull(value),
      pooled_calibration_slope = cal %>% filter(metric == "calibration_slope") %>% pull(value)
    )
}

summarise_binary_au_models <- function(split_summary_tbl) {
  split_summary_tbl %>%
    group_by(candidate_id, candidate_label, model_type, candidate_group, caution_level, validation_layer) %>%
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
      mean_split_indeterminate_rate = mean(indeterminate_rate, na.rm = TRUE),
      mean_split_au_sensitivity = mean(au_sensitivity, na.rm = TRUE),
      mean_split_au_specificity = mean(au_specificity, na.rm = TRUE),
      mean_split_au_ppv = mean(au_ppv, na.rm = TRUE),
      mean_split_au_npv = mean(au_npv, na.rm = TRUE),
      .groups = "drop"
    )
}

collapse_multiclass_predictions_to_au_binary <- function(prediction_tbl, indeterminate_threshold = 0.8) {
  prediction_tbl %>%
    transmute(
      ring,
      coord_group,
      longitude_capture,
      latitude_capture,
      truth = factor(
        ifelse(truth == "AU migrant", "AU", "non_AU"),
        levels = c("non_AU", "AU")
      ),
      .pred_class = factor(
        ifelse(`AU migrant` >= (resident + `NZ migrant`), "AU", "non_AU"),
        levels = c("non_AU", "AU")
      ),
      .pred_negative = resident + `NZ migrant`,
      .pred_positive = `AU migrant`
    ) %>%
    augment_binary_predictions(
      positive_class = "AU",
      indeterminate_threshold = indeterminate_threshold
    )
}
