#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(stringr)
})

derived_dir <- file.path("data", "derived")

format_num <- function(x, digits = 3) {
  ifelse(is.na(x), "", formatC(x, format = "f", digits = digits))
}

md_escape <- function(x) {
  x <- ifelse(is.na(x), "", x)
  str_replace_all(x, "\\|", "\\\\|")
}

df_to_markdown <- function(df) {
  header <- paste0("| ", paste(names(df), collapse = " | "), " |")
  separator <- paste0("| ", paste(rep("---", ncol(df)), collapse = " | "), " |")
  rows <- apply(df, 1, function(row) {
    paste0("| ", paste(md_escape(as.character(row)), collapse = " | "), " |")
  })
  paste(c(header, separator, rows), collapse = "\n")
}

stage3_r1 <- read_csv(
  file.path(derived_dir, "model_stage3_rung1_uncertainty_model_summary.csv"),
  show_col_types = FALSE
)

stage4_r2 <- read_csv(
  file.path(derived_dir, "model_stage4_rung2_uncertainty_model_summary.csv"),
  show_col_types = FALSE
)

stage4b_r2s2 <- read_csv(
  file.path(derived_dir, "model_stage4b_rung2s2_uncertainty_model_summary.csv"),
  show_col_types = FALSE
)

stage2_r1_cal <- read_csv(
  file.path(derived_dir, "model_stage2_rung1_calibration_summary.csv"),
  show_col_types = FALSE
)

stage4_r2_cal <- read_csv(
  file.path(derived_dir, "model_stage4_rung2_uncertainty_calibration_summary.csv"),
  show_col_types = FALSE
)

stage4b_r2s2_cal <- read_csv(
  file.path(derived_dir, "model_stage4b_rung2s2_uncertainty_calibration_summary.csv"),
  show_col_types = FALSE
)

rung_status_tbl <- tibble(
  rung = c("Rung 1", "Rung 2", "Rung 3"),
  recommended_status = c("operational", "promising but not operational", "exploratory only"),
  lead_model = c(
    "R1 C: Paired-contrast C/N",
    "Co-leading: R2 C direct breast C/N/H + paired contrasts; R2 S2 staged hierarchical R1 C then breast C/N/H NZ-vs-AU",
    "None"
  ),
  key_reasons = c(
    "Paired breast-primary C/N contrast outperformed baselines and remained the strongest validated option after uncertainty propagation with an explicit indeterminate rule.",
    "Two Rung 2 candidates now lead the field. R2 C remains the strongest direct high-information model, while R2 S2 stays competitive after uncertainty propagation, avoids direct dependence on primary H, and is easier to explain biologically; both still have high indeterminate rates and neither is operational.",
    "Rung 3 still lacks enough support for operational use because class counts are small, deployment-relevant blocking is weak, and no validated uncertainty-aware win has been shown."
  )
)

r1_final <- stage3_r1 %>%
  filter(candidate_id %in% c("r1_c_paired_contrast_cn", "r1_d_structured_paired_cn")) %>%
  select(
    candidate_id,
    candidate_label,
    validation_layer,
    mean_uncertainty_log_loss,
    mean_uncertainty_bal_accuracy,
    mean_indeterminate_rate,
    mean_switch_rate
  ) %>%
  mutate(
    layer_prefix = ifelse(validation_layer == "repeated_stratified_cv", "repeat", "blocked")
  ) %>%
  select(-validation_layer) %>%
  pivot_wider(
    names_from = layer_prefix,
    names_glue = "{layer_prefix}_{.value}",
    values_from = c(
      mean_uncertainty_log_loss,
      mean_uncertainty_bal_accuracy,
      mean_indeterminate_rate,
      mean_switch_rate
    )
  )

r1_cal <- stage2_r1_cal %>%
  filter(
    candidate_id %in% c("r1_c_paired_contrast_cn", "r1_d_structured_paired_cn"),
    metric %in% c("brier", "ece")
  ) %>%
  mutate(
    metric_name = ifelse(metric == "brier", "brier_or_macro_brier", "ece"),
    layer_prefix = ifelse(validation_layer == "repeated_stratified_cv", "repeat", "blocked")
  ) %>%
  select(candidate_id, metric_name, layer_prefix, value) %>%
  pivot_wider(
    names_from = c(layer_prefix, metric_name),
    values_from = value
  )

r2_final <- bind_rows(stage4_r2, stage4b_r2s2) %>%
  filter(candidate_id %in% c(
    "r2_a_breast_cn",
    "r2_c_breast_cnh_plus_paired_contrast",
    "r2s_b_soft_hierarchical_r1c_then_breast_cnh"
  )) %>%
  select(
    candidate_id,
    candidate_label,
    validation_layer,
    mean_uncertainty_log_loss,
    mean_uncertainty_bal_accuracy,
    mean_indeterminate_rate,
    mean_switch_rate
  ) %>%
  mutate(
    layer_prefix = ifelse(validation_layer == "repeated_stratified_cv", "repeat", "blocked")
  ) %>%
  select(-validation_layer) %>%
  pivot_wider(
    names_from = layer_prefix,
    names_glue = "{layer_prefix}_{.value}",
    values_from = c(
      mean_uncertainty_log_loss,
      mean_uncertainty_bal_accuracy,
      mean_indeterminate_rate,
      mean_switch_rate
    )
  )

r2_cal <- bind_rows(stage4_r2_cal, stage4b_r2s2_cal) %>%
  filter(
    candidate_id %in% c(
      "r2_a_breast_cn",
      "r2_c_breast_cnh_plus_paired_contrast",
      "r2s_b_soft_hierarchical_r1c_then_breast_cnh"
    ),
    class == "overall",
    metric %in% c("macro_brier", "confidence_ece")
  ) %>%
  mutate(
    metric_name = ifelse(metric == "confidence_ece", "ece", "brier_or_macro_brier"),
    layer_prefix = ifelse(validation_layer == "repeated_stratified_cv", "repeat", "blocked")
  ) %>%
  select(candidate_id, metric_name, layer_prefix, value) %>%
  pivot_wider(
    names_from = c(layer_prefix, metric_name),
    values_from = value
  )

leading_models_tbl <- bind_rows(
  r1_final %>%
    left_join(r1_cal, by = "candidate_id") %>%
    mutate(
      rung = "Rung 1",
      role = case_when(
        candidate_id == "r1_c_paired_contrast_cn" ~ "best-supported operational model",
        candidate_id == "r1_d_structured_paired_cn" ~ "secondary non-H comparator",
        TRUE ~ "candidate"
      ),
      calibration_scale = "binary_brier"
    ),
  r2_final %>%
    left_join(r2_cal, by = "candidate_id") %>%
    mutate(
      rung = "Rung 2",
      role = case_when(
        candidate_id == "r2_a_breast_cn" ~ "non-H benchmark",
        candidate_id == "r2_c_breast_cnh_plus_paired_contrast" ~ "co-leading direct candidate",
        candidate_id == "r2s_b_soft_hierarchical_r1c_then_breast_cnh" ~ "co-leading staged candidate (avoids primary H)",
        TRUE ~ "candidate"
      ),
      calibration_scale = "multiclass_macro_brier"
    )
) %>%
  arrange(
    rung,
    match(
      candidate_id,
      c(
        "r1_c_paired_contrast_cn",
        "r1_d_structured_paired_cn",
        "r2_a_breast_cn",
        "r2_c_breast_cnh_plus_paired_contrast",
        "r2s_b_soft_hierarchical_r1c_then_breast_cnh"
      )
    )
  ) %>%
  select(
    rung,
    role,
    candidate_id,
    candidate_label,
    repeat_mean_uncertainty_log_loss,
    blocked_mean_uncertainty_log_loss,
    repeat_mean_uncertainty_bal_accuracy,
    blocked_mean_uncertainty_bal_accuracy,
    repeat_mean_indeterminate_rate,
    blocked_mean_indeterminate_rate,
    repeat_mean_switch_rate,
    blocked_mean_switch_rate,
    calibration_scale,
    repeat_brier_or_macro_brier,
    blocked_brier_or_macro_brier,
    repeat_ece,
    blocked_ece
  ) %>%
  rename(
    repeat_uncertainty_log_loss = repeat_mean_uncertainty_log_loss,
    blocked_uncertainty_log_loss = blocked_mean_uncertainty_log_loss,
    repeat_uncertainty_bal_accuracy = repeat_mean_uncertainty_bal_accuracy,
    blocked_uncertainty_bal_accuracy = blocked_mean_uncertainty_bal_accuracy,
    repeat_indeterminate_rate = repeat_mean_indeterminate_rate,
    blocked_indeterminate_rate = blocked_mean_indeterminate_rate,
    repeat_mean_switch_rate = repeat_mean_switch_rate,
    blocked_mean_switch_rate = blocked_mean_switch_rate
  )

indeterminate_rules_tbl <- tibble(
  rung = c("Rung 1", "Rung 1", "Rung 2", "Rung 2"),
  stage = c("Stage 3", "Stage 3", "Stage 4", "Stage 4"),
  rule_component = c(
    "winning-class support",
    "probability-boundary crossing",
    "winning-class support",
    "mean top-class probability"
  ),
  threshold = c("< 0.80", "central 80% interval crosses 0.5", "< 0.80", "< 0.50"),
  why_it_matters = c(
    "Prevents a hard resident/migrant call when Monte Carlo draws do not strongly support one class.",
    "Flags binary calls that sit too close to the resident/migrant boundary under assay perturbation.",
    "Prevents a hard 3-class call when the winning class is not stable across uncertainty draws.",
    "Prevents a hard 3-class call when even the average winning class is weakly supported."
  )
)

limitations_tbl <- tibble(
  limitation = c(
    "Primary H remains the main H-related weak point",
    "No canonical site field for stronger blocking",
    "Multinomial sample-size limits remain important",
    "Conservative indeterminate rules reduce immediate coverage"
  ),
  scope = c(
    "Rung 2 and any future H-heavy extensions",
    "All rungs, especially deployment-relevant validation",
    "Rung 2 and especially Rung 3",
    "Rung 1 and Rung 2 operational use"
  ),
  why_it_matters = c(
    "Primary H still relies on a breast-H proxy variance and lacks direct retained repeat support, so H-driven gains need caution.",
    "Current blocked validation uses coordinate groups as a proxy rather than a canonical site variable, which weakens the strength of deployment-relevant blocking.",
    "Some multinomial training splits still have very small class counts, which increases instability and limits how strongly Rung 2 and Rung 3 can be claimed.",
    "The current rules are intentionally cautious, so promising models can still yield high indeterminate rates and limited immediate operational coverage."
  )
)

next_improvements_tbl <- tibble(
  rank = 1:5,
  improvement = c(
    "Collect direct primary-H technical repeatability and repeat-primary biological evidence",
    "Add a canonical site field to live modelling tables",
    "Increase labelled live sample size for Rung 2 and Rung 3, especially NZ subgroups",
    "Assemble an independent external or temporal holdout set",
    "Evaluate additional isotopes only after repeatability support is in place"
  ),
  why_priority = c(
    "This directly targets the main weakness behind the current Rung 2 caution and would make H-based gains much more defensible.",
    "This would allow stronger blocked validation that better matches intended deployment conditions.",
    "This would reduce multinomial instability and is essential before any serious Rung 3 operational attempt.",
    "This would test whether current results survive outside the present labelled development set and would support threshold tuning.",
    "Additional isotope expansion could help, but only after measurement support is demonstrated; sulfur is the clearest next candidate if validated."
  )
)

write_csv(
  rung_status_tbl,
  file.path(derived_dir, "model_final_rung_status.csv")
)

write_csv(
  leading_models_tbl,
  file.path(derived_dir, "model_final_leading_models.csv")
)

write_csv(
  indeterminate_rules_tbl,
  file.path(derived_dir, "model_final_indeterminate_rules.csv")
)

write_csv(
  limitations_tbl,
  file.path(derived_dir, "model_final_limitations.csv")
)

write_csv(
  next_improvements_tbl,
  file.path(derived_dir, "model_final_next_improvements.csv")
)

md_rung_status <- rung_status_tbl %>%
  transmute(
    rung = rung,
    status = recommended_status,
    lead_model = lead_model
  )

md_leading_models <- leading_models_tbl %>%
  transmute(
    rung = rung,
    role = role,
    model = candidate_label,
    repeat_log_loss = format_num(repeat_uncertainty_log_loss),
    blocked_log_loss = format_num(blocked_uncertainty_log_loss),
    repeat_bal_acc = format_num(repeat_uncertainty_bal_accuracy),
    blocked_bal_acc = format_num(blocked_uncertainty_bal_accuracy),
    repeat_indeterminate = format_num(repeat_indeterminate_rate),
    blocked_indeterminate = format_num(blocked_indeterminate_rate)
  )

md_rules <- indeterminate_rules_tbl %>%
  transmute(
    rung = rung,
    component = rule_component,
    threshold = threshold
  )

md_next <- next_improvements_tbl %>%
  transmute(
    rank = rank,
    improvement = improvement
  )

summary_lines <- c(
  "# Final Modelling Summary",
  "",
  "## Current conclusion",
  "",
  "- Best-supported Rung 1 model: `R1 C` paired-contrast C/N. It stayed strongest across the paired C/N candidates and remained defensible after assay-uncertainty propagation with an explicit indeterminate rule.",
  "- Current status of Rung 2: promising but not operational. Two candidates now lead this rung: `R2 C`, the direct breast C/N/H plus paired-contrast model, and `R2 S2`, the staged hierarchical model that first uses the accepted Rung 1 migrant signal and then uses breast C/N/H to separate NZ from AU migrants. `R2 S2` is especially attractive because it avoids direct dependence on primary H and is easier to explain biologically, but neither candidate is operational yet.",
  "- Current status of Rung 3: exploratory only. Class sizes, blocking strength, and uncertainty support are still too weak for operational screening.",
  "",
  "## Final Rung Status",
  "",
  df_to_markdown(md_rung_status),
  "",
  "## Plain-language Rung 2 logic",
  "",
  "The direct `R2 C` model tries to solve the full three-class problem in one step by combining breast C/N/H with paired seasonal contrasts. The staged `R2 S2` model breaks the question into two biologically simpler steps: first, use the accepted Rung 1 paired-contrast C/N model to ask whether a bird looks migrant-like at all; second, for birds with migrant-like signal, use breast C/N/H to separate NZ from AU. This staged logic is easier to explain and reduces reliance on primary H, but it still needs the same conservative uncertainty treatment and remains provisional.",
  "",
  "## Leading Models",
  "",
  df_to_markdown(md_leading_models),
  "",
  "## Indeterminate Rules",
  "",
  df_to_markdown(md_rules),
  "",
  "Indeterminate calls matter operationally because they prevent the pipeline from turning weakly supported or instability-prone predictions into false certainty.",
  "",
  "## Key Evidence Limitations",
  "",
  paste0("- ", limitations_tbl$limitation, ": ", limitations_tbl$why_it_matters),
  "",
  "## Ranked Next Improvements",
  "",
  df_to_markdown(md_next)
)

writeLines(
  summary_lines,
  file.path(derived_dir, "model_final_summary.md")
)
