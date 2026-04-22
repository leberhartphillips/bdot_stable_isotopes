suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

safe_sd <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) <= 1) {
    return(NA_real_)
  }
  stats::sd(x)
}

safe_range <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(NA_real_)
  }
  max(x) - min(x)
}

format_num <- function(x, digits = 3) {
  ifelse(
    is.na(x),
    "NA",
    formatC(x, format = "f", digits = digits)
  )
}

classify_evidence_strength <- function(n_groups, total_within_df) {
  if (is.na(n_groups) || n_groups == 0 || is.na(total_within_df) || total_within_df <= 0) {
    return("none")
  }
  if (n_groups <= 2 || total_within_df <= 2) {
    return("weak")
  }
  if (n_groups <= 4 || total_within_df <= 6) {
    return("moderate")
  }
  "stronger"
}

pooled_within_summary <- function(data, group_col, value_col) {
  value_sym <- rlang::sym(value_col)
  group_sym <- rlang::sym(group_col)

  group_summary <- data %>%
    filter(!is.na(!!group_sym), !is.na(!!value_sym)) %>%
    group_by(!!group_sym) %>%
    summarise(
      n_runs = n(),
      mean_value = mean(!!value_sym),
      within_ss = sum((!!value_sym - mean(!!value_sym))^2),
      within_sd = safe_sd(!!value_sym),
      within_range = safe_range(!!value_sym),
      mean_abs_deviation = mean(abs(!!value_sym - mean(!!value_sym))),
      .groups = "drop"
    ) %>%
    filter(n_runs >= 2)

  if (nrow(group_summary) == 0) {
    return(
      tibble(
        n_groups = 0L,
        n_rows = 0L,
        total_within_df = 0L,
        min_runs = NA_integer_,
        median_runs = NA_real_,
        max_runs = NA_integer_,
        pooled_within_variance = NA_real_,
        pooled_within_sd = NA_real_,
        mean_within_sd = NA_real_,
        median_within_sd = NA_real_,
        mean_within_range = NA_real_,
        median_within_range = NA_real_,
        mean_abs_deviation = NA_real_
      )
    )
  }

  total_within_df <- sum(group_summary$n_runs - 1L)
  pooled_var <- sum(group_summary$within_ss) / total_within_df

  tibble(
    n_groups = nrow(group_summary),
    n_rows = sum(group_summary$n_runs),
    total_within_df = total_within_df,
    min_runs = min(group_summary$n_runs),
    median_runs = stats::median(group_summary$n_runs),
    max_runs = max(group_summary$n_runs),
    pooled_within_variance = pooled_var,
    pooled_within_sd = sqrt(pooled_var),
    mean_within_sd = mean(group_summary$within_sd, na.rm = TRUE),
    median_within_sd = stats::median(group_summary$within_sd, na.rm = TRUE),
    mean_within_range = mean(group_summary$within_range, na.rm = TRUE),
    median_within_range = stats::median(group_summary$within_range, na.rm = TRUE),
    mean_abs_deviation = mean(group_summary$mean_abs_deviation, na.rm = TRUE)
  )
}

one_way_random_effects_metrics <- function(data, group_col, value_col) {
  value_sym <- rlang::sym(value_col)
  group_sym <- rlang::sym(group_col)

  group_summary <- data %>%
    filter(!is.na(!!group_sym), !is.na(!!value_sym)) %>%
    group_by(!!group_sym) %>%
    summarise(
      n_runs = n(),
      mean_value = mean(!!value_sym),
      within_ss = sum((!!value_sym - mean(!!value_sym))^2),
      .groups = "drop"
    ) %>%
    filter(n_runs >= 2)

  if (nrow(group_summary) < 2) {
    return(
      tibble(
        n_groups_for_icc = nrow(group_summary),
        total_rows_for_icc = sum(group_summary$n_runs),
        n0 = NA_real_,
        ms_between = NA_real_,
        ms_within = NA_real_,
        between_group_variance = NA_real_,
        within_group_variance = NA_real_,
        repeatability_icc = NA_real_
      )
    )
  }

  total_rows <- sum(group_summary$n_runs)
  total_within_df <- total_rows - nrow(group_summary)
  if (total_within_df <= 0) {
    return(
      tibble(
        n_groups_for_icc = nrow(group_summary),
        total_rows_for_icc = total_rows,
        n0 = NA_real_,
        ms_between = NA_real_,
        ms_within = NA_real_,
        between_group_variance = NA_real_,
        within_group_variance = NA_real_,
        repeatability_icc = NA_real_
      )
    )
  }

  grand_mean <- sum(group_summary$n_runs * group_summary$mean_value) / total_rows
  ss_between <- sum(group_summary$n_runs * (group_summary$mean_value - grand_mean)^2)
  ss_within <- sum(group_summary$within_ss)

  ms_between <- ss_between / (nrow(group_summary) - 1L)
  ms_within <- ss_within / total_within_df
  n0 <- (total_rows - sum(group_summary$n_runs^2) / total_rows) / (nrow(group_summary) - 1L)

  between_var <- max((ms_between - ms_within) / n0, 0)
  within_var <- ms_within
  icc <- if ((between_var + within_var) > 0) {
    between_var / (between_var + within_var)
  } else {
    NA_real_
  }

  tibble(
    n_groups_for_icc = nrow(group_summary),
    total_rows_for_icc = total_rows,
    n0 = n0,
    ms_between = ms_between,
    ms_within = ms_within,
    between_group_variance = between_var,
    within_group_variance = within_var,
    repeatability_icc = icc
  )
}
