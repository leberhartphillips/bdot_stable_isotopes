suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

fit_stage4b_r2s2 <- function(train_data, r1_candidate, migrant_candidate, seed_r1, seed_submodel) {
  r1_fit <- fit_stage3_ridge_binomial(
    train_data = train_data,
    candidate = r1_candidate,
    outcome_col = "status_rung1",
    seed = seed_r1
  )

  if (!identical(r1_fit$status, "ok")) {
    return(list(status = paste0("r1_", r1_fit$status), detail = r1_fit))
  }

  migrant_train <- train_data %>%
    filter(status_rung2 != "resident")

  sub_fit <- fit_stage2b_ridge_binomial_generic(
    train_data = migrant_train,
    candidate = migrant_candidate,
    outcome_col = "status_nzau",
    seed = seed_submodel
  )

  if (!identical(sub_fit$status, "ok")) {
    return(list(status = paste0("submodel_", sub_fit$status), detail = sub_fit))
  }

  list(
    status = "ok",
    r1_fit = r1_fit,
    sub_fit = sub_fit,
    class_levels = c("resident", "NZ migrant", "AU migrant")
  )
}

predict_stage4b_r2s2 <- function(model_obj, new_data) {
  r1_pred <- predict_stage3_ridge_binomial(model_obj$r1_fit, new_data)
  sub_pred <- predict_stage2b_ridge_binomial_generic(model_obj$sub_fit, new_data)

  prob_mat <- build_stage2b_hierarchical_probabilities(
    p_migrant = r1_pred$.pred_migrant,
    p_au_given_migrant = sub_pred$prob_positive
  )

  class_levels <- model_obj$class_levels
  pred_class <- factor(
    class_levels[max.col(prob_mat, ties.method = "first")],
    levels = class_levels
  )

  bind_cols(
    as_tibble(prob_mat, .name_repair = "minimal"),
    tibble(
      base_class = pred_class,
      r1_prob_migrant = r1_pred$.pred_migrant,
      p_au_given_migrant = sub_pred$prob_positive
    )
  )
}
