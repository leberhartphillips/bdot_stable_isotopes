suppressPackageStartupMessages({
  library(dplyr)
  library(glmnet)
  library(rsample)
  library(tibble)
})

fit_stage2b_ridge_binomial_generic <- function(train_data, candidate, outcome_col, seed) {
  prep <- prep_stage3_candidate(train_data = train_data, candidate = candidate)
  train_x <- transform_stage3_candidate(train_data, prep)$x

  if (ncol(train_x) == 0) {
    return(list(status = "no_predictors"))
  }

  y_train <- factor(train_data[[outcome_col]])
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

predict_stage2b_ridge_binomial_generic <- function(model_obj, new_data) {
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
    prob_negative = 1 - prob_pos,
    prob_positive = prob_pos,
    pred_class = factor(pred_class, levels = class_levels)
  )
}

crossfit_stage2b_r1_probabilities <- function(train_data, r1_candidate, seed) {
  set.seed(seed)
  inner_rs <- vfold_cv(
    train_data,
    v = choose_inner_folds(train_data$status_rung1),
    strata = status_rung1
  )

  pred_rows <- lapply(
    seq_len(nrow(inner_rs)),
    function(i) {
      analysis_dat <- analysis(inner_rs$splits[[i]])
      assessment_dat <- assessment(inner_rs$splits[[i]])

      fit <- fit_stage3_ridge_binomial(
        train_data = analysis_dat,
        candidate = r1_candidate,
        outcome_col = "status_rung1",
        seed = seed + i
      )

      if (!identical(fit$status, "ok")) {
        stop(
          paste(
            "Cross-fitted R1 model failed with status",
            fit$status,
            "during staged Rung 2 evaluation."
          )
        )
      }

      pred <- predict_stage3_ridge_binomial(fit, assessment_dat)

      tibble(
        ring = assessment_dat$ring,
        r1_prob_migrant = pred$.pred_migrant
      )
    }
  )

  bind_rows(pred_rows) %>%
    group_by(ring) %>%
    summarise(r1_prob_migrant = mean(r1_prob_migrant), .groups = "drop")
}

build_stage2b_hierarchical_probabilities <- function(p_migrant, p_au_given_migrant) {
  prob_mat <- cbind(
    resident = 1 - p_migrant,
    `NZ migrant` = p_migrant * (1 - p_au_given_migrant),
    `AU migrant` = p_migrant * p_au_given_migrant
  )

  rownames(prob_mat) <- NULL
  prob_mat
}

stage2b_prediction_table <- function(test_data, prob_mat) {
  class_levels <- colnames(prob_mat)
  pred_class <- factor(class_levels[max.col(prob_mat, ties.method = "first")], levels = class_levels)

  bind_cols(
    test_data %>% select(ring, coord_group, longitude_capture, latitude_capture),
    tibble(
      truth = test_data$status_rung2,
      .pred_class = pred_class
    ),
    as_tibble(prob_mat, .name_repair = "minimal")
  ) %>%
    mutate(
      top_prob = apply(prob_mat, 1, max),
      weak_support = top_prob < 0.5
    )
}

stage2b_prediction_metrics <- function(prediction_tbl, class_levels) {
  prob_mat <- as.matrix(prediction_tbl[, class_levels, drop = FALSE])

  tibble(
    log_loss = multiclass_log_loss(prediction_tbl$truth, prob_mat),
    balanced_accuracy = macro_bal_accuracy(prediction_tbl$truth, prediction_tbl$.pred_class),
    accuracy = mean(prediction_tbl$.pred_class == prediction_tbl$truth),
    brier = multiclass_brier(prediction_tbl$truth, prob_mat),
    mean_top_prob = mean(prediction_tbl$top_prob, na.rm = TRUE),
    weak_support_rate = mean(prediction_tbl$weak_support, na.rm = TRUE)
  )
}

stage2b_calibration_summary <- function(prediction_tbl, class_levels) {
  grouped <- split(
    prediction_tbl,
    interaction(prediction_tbl$candidate_id, prediction_tbl$validation_layer, drop = TRUE)
  )

  bind_rows(
    lapply(
      grouped,
      function(df) {
        base_info <- df %>%
          summarise(
            candidate_id = first(candidate_id),
            candidate_label = first(candidate_label),
            validation_layer = first(validation_layer),
            strategy_family = first(strategy_family),
            caution_level = first(caution_level)
          )

        cal <- multiclass_calibration_summary(
          predictions = df %>%
            select(truth, all_of(class_levels)),
          class_levels = class_levels
        )

        bind_cols(base_info[rep(1, nrow(cal)), ], cal)
      }
    )
  )
}
