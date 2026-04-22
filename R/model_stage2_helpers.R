suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(rsample)
  library(glmnet)
})

derive_stage2_features <- function(data) {
  data %>%
    mutate(
      delta_d13c_breast_minus_primary = normalised_d13c_Breast - normalised_d13c_Primary,
      delta_d15n_breast_minus_primary = normalised_d15n_Breast - normalised_d15n_Primary,
      delta_d2h_breast_minus_primary = normalised_d2h_Breast - normalised_d2h_Primary,
      abs_delta_d13c = abs(delta_d13c_breast_minus_primary),
      abs_delta_d15n = abs(delta_d15n_breast_minus_primary),
      abs_delta_d2h = abs(delta_d2h_breast_minus_primary),
      status_rung1 = factor(status_bin, levels = c("resident", "migrant")),
      status_rung2 = factor(
        dplyr::case_when(
          status == "resident" ~ "resident",
          status == "AU migrant" ~ "AU migrant",
          status %in% c("SI migrant", "NI migrant") ~ "NZ migrant",
          TRUE ~ NA_character_
        ),
        levels = c("resident", "NZ migrant", "AU migrant")
      ),
      coord_group = paste(longitude_capture, latitude_capture, sep = ",")
    )
}

choose_inner_folds <- function(y) {
  class_counts <- table(y)
  min_class_n <- min(class_counts)
  max(3L, min(5L, as.integer(min_class_n)))
}

binary_log_loss <- function(y, p, positive_class = "migrant") {
  y_num <- ifelse(y == positive_class, 1, 0)
  p <- pmax(pmin(p, 1 - 1e-6), 1e-6)
  -mean(y_num * log(p) + (1 - y_num) * log(1 - p))
}

binary_brier <- function(y, p, positive_class = "migrant") {
  y_num <- ifelse(y == positive_class, 1, 0)
  mean((y_num - p)^2)
}

multiclass_log_loss <- function(y, prob_mat) {
  class_index <- match(as.character(y), colnames(prob_mat))
  prob_true <- prob_mat[cbind(seq_len(nrow(prob_mat)), class_index)]
  prob_true <- pmax(pmin(prob_true, 1 - 1e-6), 1e-6)
  -mean(log(prob_true))
}

multiclass_brier <- function(y, prob_mat) {
  class_names <- colnames(prob_mat)
  truth_mat <- sapply(class_names, function(class_name) as.numeric(y == class_name))
  mean(rowMeans((truth_mat - prob_mat)^2))
}

macro_bal_accuracy <- function(y_true, y_pred) {
  class_levels <- levels(y_true)
  mean(
    vapply(
      class_levels,
      function(class_name) {
        mean(y_pred[y_true == class_name] == class_name)
      },
      numeric(1)
    ),
    na.rm = TRUE
  )
}

confidence_ece <- function(confidence, correct, n_bins = 10L) {
  if (length(confidence) == 0) {
    return(NA_real_)
  }

  confidence <- pmax(pmin(confidence, 1), 0)
  breaks <- seq(0, 1, length.out = n_bins + 1L)
  bin_id <- cut(confidence, breaks = breaks, include.lowest = TRUE, labels = FALSE)

  bin_summary <- tibble(
    bin_id = bin_id,
    confidence = confidence,
    correct = correct
  ) %>%
    group_by(bin_id) %>%
    summarise(
      n = n(),
      mean_confidence = mean(confidence),
      mean_correct = mean(correct),
      .groups = "drop"
    )

  sum((bin_summary$n / length(confidence)) * abs(bin_summary$mean_confidence - bin_summary$mean_correct))
}

binary_calibration_summary <- function(predictions, positive_class = "migrant") {
  p <- pmax(pmin(predictions$.pred_positive, 1 - 1e-6), 1e-6)
  y_num <- ifelse(predictions$truth == positive_class, 1, 0)
  logit_p <- qlogis(p)

  intercept_fit <- tryCatch(
    stats::glm(y_num ~ offset(logit_p), family = stats::binomial()),
    error = function(e) NULL
  )
  slope_fit <- tryCatch(
    stats::glm(y_num ~ logit_p, family = stats::binomial()),
    error = function(e) NULL
  )

  tibble(
    metric = c("brier", "ece", "calibration_intercept", "calibration_slope"),
    value = c(
      binary_brier(predictions$truth, p, positive_class = positive_class),
      confidence_ece(p, y_num),
      if (is.null(intercept_fit)) NA_real_ else stats::coef(intercept_fit)[1],
      if (is.null(slope_fit)) NA_real_ else stats::coef(slope_fit)[2]
    ),
    class = positive_class
  )
}

multiclass_calibration_summary <- function(predictions, class_levels) {
  prob_mat <- as.matrix(predictions[, class_levels, drop = FALSE])
  max_prob <- apply(prob_mat, 1, max)
  pred_class <- factor(class_levels[max.col(prob_mat)], levels = class_levels)
  correct <- as.numeric(pred_class == predictions$truth)

  overall <- tibble(
    metric = c("macro_brier", "confidence_ece"),
    value = c(
      multiclass_brier(predictions$truth, prob_mat),
      confidence_ece(max_prob, correct)
    ),
    class = "overall"
  )

  classwise <- bind_rows(
    lapply(
      class_levels,
      function(class_name) {
        class_prob <- prob_mat[, class_name]
        tibble(
          metric = c("one_vs_rest_brier", "one_vs_rest_ece"),
          value = c(
            mean((as.numeric(predictions$truth == class_name) - class_prob)^2),
            confidence_ece(class_prob, as.numeric(predictions$truth == class_name))
          ),
          class = class_name
        )
      }
    )
  )

  bind_rows(overall, classwise)
}

prepare_candidate_matrices <- function(train_data, test_data, candidate) {
  train <- train_data
  test <- test_data

  if (!is.null(candidate$distance_vars) && length(candidate$distance_vars) > 0) {
    distance_name <- candidate$distance_name
    distance_vars <- candidate$distance_vars

    means <- vapply(distance_vars, function(var) mean(train[[var]], na.rm = TRUE), numeric(1))
    sds <- vapply(distance_vars, function(var) stats::sd(train[[var]], na.rm = TRUE), numeric(1))
    sds[is.na(sds) | sds == 0] <- 1

    scaled_train <- sweep(as.matrix(train[, distance_vars, drop = FALSE]), 2, means, "-")
    scaled_train <- sweep(scaled_train, 2, sds, "/")
    scaled_test <- sweep(as.matrix(test[, distance_vars, drop = FALSE]), 2, means, "-")
    scaled_test <- sweep(scaled_test, 2, sds, "/")

    train[[distance_name]] <- sqrt(rowSums(scaled_train^2))
    test[[distance_name]] <- sqrt(rowSums(scaled_test^2))
  }

  predictors <- candidate$predictors
  train_x <- as.matrix(train[, predictors, drop = FALSE])
  test_x <- as.matrix(test[, predictors, drop = FALSE])

  means <- colMeans(train_x, na.rm = TRUE)
  sds <- apply(train_x, 2, stats::sd, na.rm = TRUE)
  sds[is.na(sds) | sds == 0] <- 1

  train_x <- sweep(train_x, 2, means, "-")
  train_x <- sweep(train_x, 2, sds, "/")
  test_x <- sweep(test_x, 2, means, "-")
  test_x <- sweep(test_x, 2, sds, "/")

  zero_var <- apply(train_x, 2, function(x) stats::sd(x, na.rm = TRUE) == 0)
  if (any(zero_var)) {
    keep <- !zero_var
    train_x <- train_x[, keep, drop = FALSE]
    test_x <- test_x[, keep, drop = FALSE]
  }

  list(
    train_x = train_x,
    test_x = test_x,
    predictors_used = colnames(train_x)
  )
}

fit_predict_ridge <- function(train_data, test_data, candidate, outcome_col, family, seed) {
  prep <- prepare_candidate_matrices(train_data, test_data, candidate)
  x_train <- prep$train_x
  x_test <- prep$test_x

  if (ncol(x_train) == 0) {
    return(list(status = "no_predictors"))
  }

  y_train <- train_data[[outcome_col]]
  y_test <- test_data[[outcome_col]]

  if (length(unique(y_train)) < length(levels(y_train))) {
    return(list(status = "missing_training_class"))
  }

  set.seed(seed)
  inner_folds <- choose_inner_folds(y_train)

  if (family == "binomial") {
    y_train_num <- ifelse(y_train == levels(y_train)[2], 1, 0)

    fit <- tryCatch(
      cv.glmnet(
        x = x_train,
        y = y_train_num,
        family = "binomial",
        alpha = 0,
        nfolds = inner_folds,
        type.measure = "deviance",
        standardize = FALSE
      ),
      error = function(e) e
    )

    if (inherits(fit, "error")) {
      return(list(status = "fit_error", message = fit$message))
    }

    prob_pos <- as.numeric(predict(fit, newx = x_test, s = "lambda.1se", type = "response"))
    pred_class <- ifelse(prob_pos >= 0.5, levels(y_train)[2], levels(y_train)[1])
    pred_class <- factor(pred_class, levels = levels(y_train))

    prob_tbl <- tibble(
      .pred_negative = 1 - prob_pos,
      .pred_positive = prob_pos
    )

    metrics <- tibble(
      log_loss = binary_log_loss(y_test, prob_pos, positive_class = levels(y_train)[2]),
      balanced_accuracy = macro_bal_accuracy(y_test, pred_class),
      accuracy = mean(pred_class == y_test),
      brier = binary_brier(y_test, prob_pos, positive_class = levels(y_train)[2]),
      selected_lambda = fit$lambda.1se
    )

    predictions <- bind_cols(
      test_data %>% select(ring, coord_group, longitude_capture, latitude_capture),
      tibble(
        truth = y_test,
        .pred_class = pred_class
      ),
      prob_tbl
    )

    return(list(
      status = "ok",
      predictions = predictions,
      metrics = metrics,
      predictors_used = prep$predictors_used
    ))
  }

  fit <- tryCatch(
    cv.glmnet(
      x = x_train,
      y = y_train,
      family = "multinomial",
      alpha = 0,
      nfolds = inner_folds,
      type.measure = "deviance",
      standardize = FALSE,
      type.multinomial = "grouped"
    ),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    return(list(status = "fit_error", message = fit$message))
  }

  pred_arr <- predict(fit, newx = x_test, s = "lambda.1se", type = "response")
  pred_arr <- drop(pred_arr)

  if (is.null(dim(pred_arr))) {
    pred_arr <- matrix(pred_arr, nrow = 1)
  }

  if (nrow(pred_arr) != nrow(test_data)) {
    pred_arr <- t(pred_arr)
  }

  class_levels <- levels(y_train)
  colnames(pred_arr) <- class_levels
  pred_class <- factor(class_levels[max.col(pred_arr, ties.method = "first")], levels = class_levels)

  metrics <- tibble(
    log_loss = multiclass_log_loss(y_test, pred_arr),
    balanced_accuracy = macro_bal_accuracy(y_test, pred_class),
    accuracy = mean(pred_class == y_test),
    brier = multiclass_brier(y_test, pred_arr),
    selected_lambda = fit$lambda.1se
  )

  predictions <- bind_cols(
    test_data %>% select(ring, coord_group, longitude_capture, latitude_capture),
    tibble(
      truth = y_test,
      .pred_class = pred_class
    ),
    as_tibble(pred_arr, .name_repair = "minimal")
  )

  list(
    status = "ok",
    predictions = predictions,
    metrics = metrics,
    predictors_used = prep$predictors_used
  )
}

make_balanced_coord_folds <- function(data, outcome_col, group_col = "coord_group", v = 5L) {
  outcome_sym <- rlang::sym(outcome_col)
  group_sym <- rlang::sym(group_col)

  group_counts <- data %>%
    count(!!group_sym, !!outcome_sym, name = "n") %>%
    tidyr::pivot_wider(
      names_from = !!outcome_sym,
      values_from = n,
      values_fill = 0
    ) %>%
    rename(group_id = !!rlang::sym(group_col))

  class_cols <- setdiff(colnames(group_counts), "group_id")
  group_counts$total_n <- rowSums(group_counts[, class_cols, drop = FALSE])
  group_counts <- group_counts %>%
    arrange(desc(total_n))

  target_class <- colSums(group_counts[, class_cols, drop = FALSE]) / v
  target_total <- sum(group_counts$total_n) / v
  target_mat <- matrix(rep(target_class, each = v), nrow = v)

  fold_class <- matrix(0, nrow = v, ncol = length(class_cols))
  colnames(fold_class) <- class_cols
  fold_total <- rep(0, v)
  assigned_fold <- integer(nrow(group_counts))

  for (i in seq_len(nrow(group_counts))) {
    row_counts <- as.numeric(group_counts[i, class_cols, drop = FALSE])
    row_total <- group_counts$total_n[i]

    scores <- vapply(
      seq_len(v),
      function(fold_id) {
        new_fold_class <- fold_class
        new_fold_class[fold_id, ] <- new_fold_class[fold_id, ] + row_counts

        new_fold_total <- fold_total
        new_fold_total[fold_id] <- new_fold_total[fold_id] + row_total

        class_score <- sum(((new_fold_class - target_mat) / pmax(target_mat, 1))^2)
        total_score <- sum(((new_fold_total - target_total) / pmax(target_total, 1))^2)
        empty_penalty <- sum(new_fold_total == 0) * 100

        class_score + total_score + empty_penalty
      },
      numeric(1)
    )

    best_folds <- which(scores == min(scores))
    if (length(best_folds) > 1) {
      best_fold <- best_folds[which.min(fold_total[best_folds])]
    } else {
      best_fold <- best_folds
    }

    assigned_fold[i] <- best_fold
    fold_class[best_fold, ] <- fold_class[best_fold, ] + row_counts
    fold_total[best_fold] <- fold_total[best_fold] + row_total
  }

  tibble(
    group_id = group_counts$group_id,
    blocked_fold = paste0("Fold", assigned_fold)
  )
}
