#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
  library(readr)
  library(stringr)
  library(tidyr)
})

derived_dir <- file.path("data", "derived")

extract_sample_date <- function(x) {
  as.Date(str_extract(x, "20[0-9]{2}-[0-9]{2}-[0-9]{2}"))
}

safe_first_date <- function(x) {
  x <- sort(unique(x[!is.na(x)]))
  if (length(x) == 0) {
    return(as.Date(NA))
  }
  x[[1]]
}

raw_live <- read_csv(
  file.path(derived_dir, "live_cn_measurements_raw_qc_long.csv"),
  show_col_types = FALSE
)

adjudicated_long <- read_csv(
  file.path(derived_dir, "live_cn_measurements_adjudicated_long.csv"),
  show_col_types = FALSE
)

adjudicated_by_sample <- read_csv(
  file.path(derived_dir, "live_cn_measurements_adjudicated_by_sample.csv"),
  show_col_types = FALSE
)

live_tissue <- read_csv(
  file.path(derived_dir, "live_cn_tissue_summary.csv"),
  show_col_types = FALSE
)

live_paired <- read_csv(
  file.path(derived_dir, "live_cn_paired_by_ring.csv"),
  show_col_types = FALSE
)

live_labelled <- read_csv(
  file.path(derived_dir, "live_screening_ready_paired_labelled.csv"),
  show_col_types = FALSE
)

museum_winter <- read_csv(
  file.path(derived_dir, "museum_screening_ready_breast_homogenate_winter.csv"),
  show_col_types = FALSE
)

qc_issues <- read_csv(
  file.path(derived_dir, "qc_issues.csv"),
  show_col_types = FALSE
)

model_final_leading <- read_csv(
  file.path(derived_dir, "model_final_leading_models.csv"),
  show_col_types = FALSE
)

model_final_rung_status <- read_csv(
  file.path(derived_dir, "model_final_rung_status.csv"),
  show_col_types = FALSE
)

stage3_rung1_assessment <- read_csv(
  file.path(derived_dir, "model_stage3_rung1_uncertainty_assessment_summary.csv"),
  show_col_types = FALSE
)

stage3_rung1_model_summary <- read_csv(
  file.path(derived_dir, "model_stage3_rung1_uncertainty_model_summary.csv"),
  show_col_types = FALSE
)

stage4_rung2_assessment <- read_csv(
  file.path(derived_dir, "model_stage4_rung2_uncertainty_assessment_summary.csv"),
  show_col_types = FALSE
)

stage4_rung2_model_summary <- read_csv(
  file.path(derived_dir, "model_stage4_rung2_uncertainty_model_summary.csv"),
  show_col_types = FALSE
)

stage4b_rung2s2_assessment <- read_csv(
  file.path(derived_dir, "model_stage4b_rung2s2_uncertainty_assessment_summary.csv"),
  show_col_types = FALSE
)

stage4b_rung2s2_model_summary <- read_csv(
  file.path(derived_dir, "model_stage4b_rung2s2_uncertainty_model_summary.csv"),
  show_col_types = FALSE
)

stage4_rung2_calibration <- read_csv(
  file.path(derived_dir, "model_stage4_rung2_uncertainty_calibration_summary.csv"),
  show_col_types = FALSE
)

stage4b_rung2s2_calibration <- read_csv(
  file.path(derived_dir, "model_stage4b_rung2s2_uncertainty_calibration_summary.csv"),
  show_col_types = FALSE
)

unknown_predictions_wide <- read_csv(
  file.path(derived_dir, "model_application_unknown_predictions_wide.csv"),
  show_col_types = FALSE
)

unknown_probability_long <- read_csv(
  file.path(derived_dir, "model_application_unknown_probability_long.csv"),
  show_col_types = FALSE
)

live_ring_dates <- adjudicated_by_sample %>%
  mutate(sample_date = extract_sample_date(sample_id_base)) %>%
  group_by(ring) %>%
  summarise(
    sample_date = safe_first_date(sample_date),
    n_unique_capture_dates = n_distinct(sample_date, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    capture_year = year(sample_date),
    capture_month = month(sample_date),
    capture_month_label = if_else(
      is.na(capture_month),
      NA_character_,
      as.character(month(sample_date, label = TRUE, abbr = TRUE))
    )
  )

report_sample_flow <- tibble(
  stage_order = 1:10,
  stage = c(
    "Raw live assay-isotope rows",
    "Audit-retained live assay-isotope rows",
    "Model-contributing live assay-isotope rows",
    "Adjudicated live sample-isotope values",
    "Adjudicated live sample rows",
    "Live tissue summaries",
    "Paired rings",
    "Complete paired C/N rings",
    "Complete paired C/N/H rings",
    "Labelled paired birds"
  ),
  n = c(
    nrow(raw_live),
    sum(raw_live$assay_qc_audit_retained, na.rm = TRUE),
    sum(raw_live$assay_qc_contributes_to_model_value, na.rm = TRUE),
    nrow(adjudicated_long),
    nrow(adjudicated_by_sample),
    nrow(live_tissue),
    nrow(live_paired),
    sum(live_paired$has_complete_cn_pair, na.rm = TRUE),
    sum(live_paired$has_complete_cnh_pair, na.rm = TRUE),
    nrow(live_labelled)
  ),
  unit = c(
    rep("assay-isotope rows", 3),
    "sample-isotope values",
    "samples",
    "ring x feather summaries",
    rep("rings", 3),
    "birds"
  ),
  note = c(
    "One row per live assay x isotope value before adjudication.",
    "Retained for audit after explicit workbook-backed QC exclusions only.",
    "Retained rows that contributed to the adjudicated modelling values.",
    "One adjudicated modelling value per sample x isotope.",
    "Sample-level rows after recombining isotopes by adjudicated sample.",
    "One row per ring x feather type summary used for pairing.",
    "One ring-level row with breast and primary tissues side by side.",
    "Paired rings with complete breast and primary C/N.",
    "Paired rings with complete breast and primary C/N/H.",
    "Known-status birds retained for accepted screening-model comparison."
  )
)

report_live_pairing_status <- live_paired %>%
  transmute(
    status_known = if_else(status_known, "Known status", "Unknown status"),
    paired_cn = if_else(has_complete_cn_pair, "Complete C/N pair", "Incomplete C/N pair"),
    paired_cnh = if_else(has_complete_cnh_pair, "Complete C/N/H pair", "Incomplete C/N/H pair")
  ) %>%
  count(status_known, paired_cn, paired_cnh, name = "n_rings")

report_isotope_completeness <- bind_rows(
  live_tissue %>%
    transmute(
      source = "live",
      subset = "tissue_summary",
      feather_type = feather_type,
      n_total = n(),
      d13c = !is.na(normalised_d13c),
      d15n = !is.na(normalised_d15n),
      d2h = !is.na(normalised_d2h)
    ) %>%
    pivot_longer(
      cols = c(d13c, d15n, d2h),
      names_to = "isotope",
      values_to = "present"
    ) %>%
    group_by(source, subset, feather_type, isotope) %>%
    summarise(
      n_present = sum(present, na.rm = TRUE),
      n_total = first(n_total),
      prop_present = n_present / n_total,
      .groups = "drop"
    ),
  museum_winter %>%
    transmute(
      source = "museum",
      subset = "winter_breast_homogenate",
      feather_type = "Breast",
      n_total = n(),
      d13c = has_d13c,
      d15n = has_d15n,
      d18o = has_d18o,
      d2h = has_d2h
    ) %>%
    pivot_longer(
      cols = c(d13c, d15n, d18o, d2h),
      names_to = "isotope",
      values_to = "present"
    ) %>%
    group_by(source, subset, feather_type, isotope) %>%
    summarise(
      n_present = sum(present, na.rm = TRUE),
      n_total = first(n_total),
      prop_present = n_present / n_total,
      .groups = "drop"
    )
)

report_rung_class_counts <- bind_rows(
  live_labelled %>%
    count(status_bin, name = "n") %>%
    transmute(rung = "Rung 1", class = status_bin, n),
  live_labelled %>%
    mutate(
      rung2_class = case_when(
        status == "resident" ~ "resident",
        status %in% c("SI migrant", "NI migrant") ~ "NZ migrant",
        status == "AU migrant" ~ "AU migrant",
        TRUE ~ NA_character_
      )
    ) %>%
    count(rung2_class, name = "n") %>%
    transmute(rung = "Rung 2", class = rung2_class, n),
  live_labelled %>%
    count(status, name = "n") %>%
    transmute(rung = "Rung 3", class = status, n)
)

report_live_status_counts <- live_labelled %>%
  count(status, status_bin, name = "n_birds") %>%
  arrange(desc(n_birds), status)

report_live_sampling_points <- live_paired %>%
  left_join(live_ring_dates, by = "ring") %>%
  transmute(
    ring,
    status,
    status_bin,
    status_known,
    longitude_capture,
    latitude_capture,
    sample_date,
    capture_year,
    capture_month,
    has_complete_cn_pair,
    has_complete_cnh_pair
  )

report_museum_sampling_points <- museum_winter %>%
  transmute(
    specimen_id,
    geo_region_bin,
    collection_year,
    collection_month,
    lon,
    lat,
    has_d13c,
    has_d15n,
    has_d18o,
    has_d2h,
    has_all_4_isotopes,
    has_core_3_isotopes
  )

report_live_temporal_summary <- live_ring_dates %>%
  mutate(
    capture_year = year(sample_date),
    capture_month = month(sample_date),
    capture_month_label = if_else(
      is.na(capture_month),
      NA_character_,
      as.character(month(sample_date, label = TRUE, abbr = TRUE))
    )
  ) %>%
  count(capture_year, capture_month, capture_month_label, name = "n_rings") %>%
  arrange(capture_year, capture_month)

report_museum_temporal_summary <- museum_winter %>%
  count(collection_year, collection_month, name = "n_specimens") %>%
  arrange(collection_year, collection_month)

report_qc_issue_summary <- qc_issues %>%
  count(dataset, severity, issue_type, wt = n_records, name = "n_records") %>%
  arrange(desc(n_records), dataset, severity, issue_type)

report_leading_models_long <- model_final_leading %>%
  select(
    rung,
    role,
    candidate_id,
    candidate_label,
    repeat_uncertainty_log_loss,
    blocked_uncertainty_log_loss,
    repeat_uncertainty_bal_accuracy,
    blocked_uncertainty_bal_accuracy,
    repeat_indeterminate_rate,
    blocked_indeterminate_rate,
    repeat_mean_switch_rate,
    blocked_mean_switch_rate,
    repeat_brier_or_macro_brier,
    blocked_brier_or_macro_brier,
    repeat_ece,
    blocked_ece
  ) %>%
  pivot_longer(
    cols = -c(rung, role, candidate_id, candidate_label),
    names_to = c("validation_layer", "metric"),
    names_pattern = "^(repeat|blocked)_(.*)$",
    values_to = "value"
  ) %>%
  mutate(
    validation_layer = recode(
      validation_layer,
      `repeat` = "Repeated stratified CV",
      blocked = "Blocked coordinate CV"
    )
  ) %>%
  pivot_wider(
    names_from = metric,
    values_from = value
  )

report_rung_status <- model_final_rung_status %>%
  rename(status = recommended_status)

report_rung1_features <- live_labelled %>%
  transmute(
    ring,
    rung1_class = status_bin,
    status,
    abs_delta_d13c = abs(delta_d13c_breast_minus_primary),
    abs_delta_d15n = abs(delta_d15n_breast_minus_primary),
    normalised_d13c_Breast,
    normalised_d13c_Primary,
    normalised_d15n_Breast,
    normalised_d15n_Primary
  )

report_rung1_breast_primary_long <- report_rung1_features %>%
  transmute(
    ring,
    rung1_class,
    d13c_breast = normalised_d13c_Breast,
    d13c_primary = normalised_d13c_Primary,
    d15n_breast = normalised_d15n_Breast,
    d15n_primary = normalised_d15n_Primary
  ) %>%
  pivot_longer(
    cols = c(d13c_breast, d13c_primary, d15n_breast, d15n_primary),
    names_to = c("isotope", "tissue"),
    names_sep = "_",
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from = tissue,
    values_from = value
  )

report_rung1_uncertainty_cases <- stage3_rung1_assessment %>%
  filter(candidate_id == "r1_c_paired_contrast_cn") %>%
  transmute(
    validation_layer,
    ring,
    truth,
    base_prob_migrant,
    uncertainty_prob_migrant_mean,
    uncertainty_prob_migrant_q10,
    uncertainty_prob_migrant_q90,
    class_support,
    boundary_crosses,
    indeterminate,
    operational_call,
    switch_rate_from_base,
    top_prob_abs_shift
  )

report_rung2_features <- live_labelled %>%
  mutate(
    rung2_class = case_when(
      status == "resident" ~ "resident",
      status %in% c("SI migrant", "NI migrant") ~ "NZ migrant",
      status == "AU migrant" ~ "AU migrant",
      TRUE ~ NA_character_
    ),
    abs_delta_d13c = abs(delta_d13c_breast_minus_primary),
    abs_delta_d15n = abs(delta_d15n_breast_minus_primary),
    abs_delta_d2h = abs(normalised_d2h_Breast - normalised_d2h_Primary)
  ) %>%
  transmute(
    ring,
    rung2_class,
    status,
    normalised_d13c_Breast,
    normalised_d15n_Breast,
    normalised_d2h_Breast,
    abs_delta_d13c,
    abs_delta_d15n,
    abs_delta_d2h
  )

report_rung2_breast_long <- report_rung2_features %>%
  select(ring, rung2_class, starts_with("normalised_")) %>%
  pivot_longer(
    cols = starts_with("normalised_"),
    names_to = "isotope",
    values_to = "value"
  )

report_rung2_uncertainty_cases <- stage4_rung2_assessment %>%
  bind_rows(stage4b_rung2s2_assessment) %>%
  filter(candidate_id %in% c(
    "r2_a_breast_cn",
    "r2_c_breast_cnh_plus_paired_contrast",
    "r2s_b_soft_hierarchical_r1c_then_breast_cnh"
  )) %>%
  transmute(
    candidate_id,
    candidate_label,
    validation_layer,
    ring,
    truth,
    top_class_mean_prob,
    class_support,
    indeterminate,
    operational_call,
    switch_rate_from_base,
    top_prob_abs_shift
  )

report_metric_interpretation <- tibble(
  metric = c(
    "Balanced accuracy",
    "Log loss",
    "Macro Brier score",
    "Confidence ECE",
    "Indeterminate rate",
    "Class-switch rate"
  ),
  plain_language_meaning = c(
    "How well the model separates the classes overall, without letting the biggest class dominate the score.",
    "How much probability the model places on the right answer, with a strong penalty for confident wrong answers.",
    "How close the predicted probabilities are to what actually happened across the classes.",
    "How believable the model's stated confidence is. If it says 80% confidence often, about 80% of those cases should be right.",
    "How often the model declines to make a confident call after uncertainty is taken into account.",
    "How often technical uncertainty changes the predicted class across Monte Carlo draws."
  ),
  better_direction = c(
    "Higher is better",
    "Lower is better",
    "Lower is better",
    "Lower is better",
    "Not automatically lower",
    "Lower is generally better"
  ),
  good_pattern_here = c(
    "Higher than the benchmark while still holding up in blocked validation.",
    "Low values for both repeated and blocked validation, especially for the leading candidates.",
    "Low values alongside sensible calibration and not just sharper but brittle probabilities.",
    "Low enough that the model's confidence looks believable rather than exaggerated.",
    "Moderate rather than extreme: enough indeterminate cases to avoid false certainty, but not so many that the model becomes unusable.",
    "Low enough that assay noise usually softens borderline calls instead of flipping birds into different classes."
  ),
  tradeoff_captured = c(
    "Discrimination between classes.",
    "Probability quality and overconfidence penalties.",
    "Overall probability accuracy.",
    "Calibration of stated confidence.",
    "Operational coverage versus caution.",
    "Stability under measurement uncertainty."
  )
)

report_rung1_candidate_metrics <- stage3_rung1_model_summary %>%
  filter(candidate_id %in% c(
    "r1_c_paired_contrast_cn",
    "r1_d_structured_paired_cn",
    "r1_g_structured_paired_cnh"
  )) %>%
  transmute(
    rung = "Rung 1",
    candidate_id,
    candidate_label,
    validation_layer = if_else(
      validation_layer == "repeated_stratified_cv",
      "Repeated stratified CV",
      "Blocked coordinate CV"
    ),
    uncertainty_log_loss = mean_uncertainty_log_loss,
    uncertainty_bal_accuracy = mean_uncertainty_bal_accuracy,
    indeterminate_rate = mean_indeterminate_rate,
    mean_switch_rate = mean_switch_rate,
    sensitivity_only = sensitivity_only,
    caution_level
  )

report_rung2_candidate_metrics <- bind_rows(
  stage4_rung2_model_summary %>%
    filter(candidate_id %in% c("r2_a_breast_cn", "r2_c_breast_cnh_plus_paired_contrast")),
  stage4b_rung2s2_model_summary
) %>%
  left_join(
    bind_rows(stage4_rung2_calibration, stage4b_rung2s2_calibration) %>%
      filter(class == "overall", metric %in% c("macro_brier", "confidence_ece")) %>%
      select(candidate_id, validation_layer, metric, value) %>%
      pivot_wider(names_from = metric, values_from = value),
    by = c("candidate_id", "validation_layer")
  ) %>%
  transmute(
    rung = "Rung 2",
    candidate_id,
    candidate_label,
    validation_layer = if_else(
      validation_layer == "repeated_stratified_cv",
      "Repeated stratified CV",
      "Blocked coordinate CV"
    ),
    uncertainty_log_loss = mean_uncertainty_log_loss,
    uncertainty_bal_accuracy = mean_uncertainty_bal_accuracy,
    macro_brier = macro_brier,
    confidence_ece = confidence_ece,
    indeterminate_rate = mean_indeterminate_rate,
    mean_switch_rate = mean_switch_rate,
    caution_level
  )

report_model_outcome_composition <- bind_rows(
  stage3_rung1_assessment %>%
    filter(candidate_id %in% c(
      "r1_c_paired_contrast_cn",
      "r1_d_structured_paired_cn",
      "r1_g_structured_paired_cnh"
    )) %>%
    transmute(
      rung = "Rung 1",
      candidate_id,
      candidate_label,
      validation_layer = if_else(
        validation_layer == "repeated_stratified_cv",
        "Repeated stratified CV",
        "Blocked coordinate CV"
      ),
      truth,
      operational_call = if_else(
        indeterminate,
        "indeterminate",
        as.character(uncertainty_class)
      ),
      indeterminate
    ),
  bind_rows(stage4_rung2_assessment, stage4b_rung2s2_assessment) %>%
    filter(candidate_id %in% c(
      "r2_a_breast_cn",
      "r2_c_breast_cnh_plus_paired_contrast",
      "r2s_b_soft_hierarchical_r1c_then_breast_cnh"
    )) %>%
    transmute(
      rung = "Rung 2",
      candidate_id,
      candidate_label,
      validation_layer = if_else(
        validation_layer == "repeated_stratified_cv",
        "Repeated stratified CV",
        "Blocked coordinate CV"
      ),
      truth,
      operational_call,
      indeterminate
    )
) %>%
  mutate(
    outcome_bucket = case_when(
      indeterminate ~ "Indeterminate",
      operational_call == truth ~ "Determinate correct",
      TRUE ~ "Determinate incorrect"
    )
  ) %>%
  count(rung, candidate_id, candidate_label, validation_layer, outcome_bucket, name = "n") %>%
  group_by(rung, candidate_id, candidate_label, validation_layer) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

report_unknown_predictions <- unknown_predictions_wide %>%
  mutate(
    r1_summary = case_when(
      !r1_prediction_available ~ "Unavailable",
      TRUE ~ r1_operational_call
    ),
    r2a_summary = case_when(
      !r2a_prediction_available ~ "Unavailable",
      TRUE ~ r2a_operational_call
    ),
    r2c_summary = case_when(
      !r2c_prediction_available ~ "Unavailable",
      TRUE ~ r2c_operational_call
    ),
    r2s2_summary = case_when(
      !r2s2_prediction_available ~ "Unavailable",
      TRUE ~ r2s2_operational_call
    )
  ) %>%
  arrange(source_order)

report_unknown_probability_long <- unknown_probability_long %>%
  left_join(
    unknown_predictions_wide %>%
      select(
        source_order,
        ring,
        longitude_capture,
        latitude_capture,
        informal_rung2_consensus
      ),
    by = c("source_order", "ring")
  ) %>%
  mutate(
    model_short = recode(
      model_id,
      r1_c = "R1 C",
      r2_a = "R2 A",
      r2_c = "R2 C",
      r2_s2 = "R2 S2"
    )
  ) %>%
  arrange(source_order, model_id, class)

report_unknown_map <- report_unknown_predictions %>%
  transmute(
    source_order,
    ring,
    longitude_capture,
    latitude_capture,
    has_complete_cn_pair,
    has_complete_cnh_pair,
    r1_map_status = case_when(
      !r1_prediction_available ~ "R1 unavailable",
      TRUE ~ paste("R1:", r1_operational_call)
    ),
    r2a_map_status = case_when(
      !r2a_prediction_available ~ "R2 A unavailable",
      TRUE ~ paste("R2 A:", r2a_operational_call)
    ),
    r2c_map_status = case_when(
      !r2c_prediction_available ~ "R2 C unavailable",
      TRUE ~ paste("R2 C:", r2c_operational_call)
    ),
    r2s2_map_status = case_when(
      !r2s2_prediction_available ~ "R2 S2 unavailable",
      TRUE ~ paste("R2 S2:", r2s2_operational_call)
    ),
    informal_rung2_consensus
  ) %>%
  arrange(source_order)

report_unknown_prediction_counts <- bind_rows(
  report_unknown_predictions %>%
    transmute(
      model = "R1 C",
      prediction_available = r1_prediction_available,
      operational_call = r1_operational_call
    ),
  report_unknown_predictions %>%
    transmute(
      model = "R2 A",
      prediction_available = r2a_prediction_available,
      operational_call = r2a_operational_call
    ),
  report_unknown_predictions %>%
    transmute(
      model = "R2 C",
      prediction_available = r2c_prediction_available,
      operational_call = r2c_operational_call
    ),
  report_unknown_predictions %>%
    transmute(
      model = "R2 S2",
      prediction_available = r2s2_prediction_available,
      operational_call = r2s2_operational_call
    )
) %>%
  mutate(
    call_bucket = case_when(
      !prediction_available ~ "Unavailable",
      operational_call == "indeterminate" ~ "Indeterminate",
      TRUE ~ operational_call
    )
  ) %>%
  count(model, call_bucket, name = "n_birds") %>%
  group_by(model) %>%
  mutate(prop = n_birds / sum(n_birds)) %>%
  ungroup()

write_csv(report_sample_flow, file.path(derived_dir, "report_sample_flow.csv"))
write_csv(report_live_pairing_status, file.path(derived_dir, "report_live_pairing_status.csv"))
write_csv(report_isotope_completeness, file.path(derived_dir, "report_isotope_completeness.csv"))
write_csv(report_rung_class_counts, file.path(derived_dir, "report_rung_class_counts.csv"))
write_csv(report_live_status_counts, file.path(derived_dir, "report_live_status_counts.csv"))
write_csv(report_live_sampling_points, file.path(derived_dir, "report_live_sampling_points.csv"))
write_csv(report_museum_sampling_points, file.path(derived_dir, "report_museum_sampling_points.csv"))
write_csv(report_live_temporal_summary, file.path(derived_dir, "report_live_temporal_summary.csv"))
write_csv(report_museum_temporal_summary, file.path(derived_dir, "report_museum_temporal_summary.csv"))
write_csv(report_qc_issue_summary, file.path(derived_dir, "report_qc_issue_summary.csv"))
write_csv(report_leading_models_long, file.path(derived_dir, "report_leading_models_long.csv"))
write_csv(report_rung_status, file.path(derived_dir, "report_rung_status.csv"))
write_csv(report_rung1_features, file.path(derived_dir, "report_rung1_features.csv"))
write_csv(report_rung1_breast_primary_long, file.path(derived_dir, "report_rung1_breast_primary_long.csv"))
write_csv(report_rung1_uncertainty_cases, file.path(derived_dir, "report_rung1_uncertainty_cases.csv"))
write_csv(report_rung2_features, file.path(derived_dir, "report_rung2_features.csv"))
write_csv(report_rung2_breast_long, file.path(derived_dir, "report_rung2_breast_long.csv"))
write_csv(report_rung2_uncertainty_cases, file.path(derived_dir, "report_rung2_uncertainty_cases.csv"))
write_csv(report_metric_interpretation, file.path(derived_dir, "report_metric_interpretation.csv"))
write_csv(report_rung1_candidate_metrics, file.path(derived_dir, "report_rung1_candidate_metrics.csv"))
write_csv(report_rung2_candidate_metrics, file.path(derived_dir, "report_rung2_candidate_metrics.csv"))
write_csv(report_model_outcome_composition, file.path(derived_dir, "report_model_outcome_composition.csv"))
write_csv(report_unknown_predictions, file.path(derived_dir, "report_unknown_predictions.csv"))
write_csv(report_unknown_probability_long, file.path(derived_dir, "report_unknown_probability_long.csv"))
write_csv(report_unknown_map, file.path(derived_dir, "report_unknown_map.csv"))
write_csv(report_unknown_prediction_counts, file.path(derived_dir, "report_unknown_prediction_counts.csv"))
