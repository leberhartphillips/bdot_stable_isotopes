# Canonical Derived Datasets

This directory contains non-destructive, analysis-ready tables built from in-repo inputs used by `bdot_stable_isotopes2.qmd`.

## Museum data

- `museum_measurements_raw_qc.csv`
  - one row per raw museum isotope assay before metadata joins
  - preserves raw value, adjudicated usable value, assay-QC action, and explicit inclusion/exclusion flags
- `museum_measurements_long.csv`
  - one row per isotope assay
  - keeps only adjudicated usable museum assays after QC and metadata joins
  - preserves source file, sheet, row, specimen label, assay branch, replicate label, and joined metadata
- `museum_metadata_canonical.csv`
  - one row per `specimen_id + feather_type`
  - preserves source inventory row ids and joined checkpoint-derived coordinates
  - includes flags showing whether checkpoint location/date agree with current inventory metadata
- `museum_specimen_region_isotopes_long.csv`
  - one row per `specimen_id + feather_type + region + isotope`
  - averages repeated assays within the same tissue-region-isotope unit
- `museum_specimen_region_isotopes_wide.csv`
  - one row per `specimen_id + feather_type + region`
  - keeps isotope values in wide format with completeness flags
- `museum_screening_ready_breast_homogenate_winter.csv`
  - screening-ready candidate subset
  - restricted to breast feathers, homogenate measurements, and May-August collection months
  - no final model has been fit here

## Live bird data

- `live_phenotypes_canonical.csv`
  - one row per ring from the in-repo LIMS workbook
- `live_cn_measurements_by_sample.csv`
  - one row per live bird assay row from batch 5 C/N results
  - includes raw, audit-retained, and adjudicated C/N values, repeat-assay flags, phenotype joins, and assay-QC fields
- `live_cn_measurements_raw_qc_long.csv`
  - one row per raw live assay-isotope value
  - preserves assay-level QC actions, repeat-pair status, audit-retention flags, and whether each retained row contributed to the adjudicated modelling value
- `live_cn_measurements_adjudicated_long.csv`
  - one row per adjudicated live sample-isotope modelling value
  - derived from retained assay rows after applying explicit workbook exclusions only
  - uses explicit selected rows, explicit green averages, or the mean of retained valid repeats where no explicit workbook instruction chooses one repeat
- `live_cn_measurements_adjudicated_long_strict.csv`
  - strict adjudicated live sample-isotope values
  - excludes unresolved repeat pairs from model-input use
- `live_cn_measurements_adjudicated_long_sensitivity.csv`
  - sensitivity adjudicated live sample-isotope values
  - retains a deterministic non-`_rpt` default where an unresolved repeat pair has exactly one original row
- `live_cn_measurements_adjudicated_by_sample.csv`
  - one row per adjudicated live assay/sample after isotope values are recombined
- `live_cn_measurements_adjudicated_by_sample_strict.csv`
  - one row per strict adjudicated live assay/sample after isotope values are recombined
- `live_cn_measurements_adjudicated_by_sample_sensitivity.csv`
  - one row per sensitivity adjudicated live assay/sample after isotope values are recombined
- `live_cn_tissue_summary.csv`
  - one row per `ring + feather_type`
  - averages adjudicated live assays within tissue
- `live_cn_tissue_summary_strict.csv`
  - strict tissue-level C/N summaries for model-input use
- `live_cn_tissue_summary_sensitivity.csv`
  - sensitivity tissue-level C/N summaries for model-input use
- `live_cn_paired_by_ring.csv`
  - one row per ring with breast and primary C/N values side by side
  - includes completeness flags and tissue deltas
- `live_cn_paired_by_ring_strict.csv`
  - strict ring-level paired C/N table
  - unresolved repeat-pair tissues are removed before pairing
- `live_cn_paired_by_ring_sensitivity.csv`
  - sensitivity ring-level paired C/N table
  - unresolved repeat-pair tissues use the deterministic non-`_rpt` default when available
- `live_screening_ready_paired_labelled.csv`
  - labelled live-bird candidate subset
  - restricted to known migratory status and complete primary + breast C/N pairs
- `live_screening_ready_paired_labelled_strict.csv`
  - strict labelled live-bird candidate subset for conservative model-input use
- `live_screening_ready_paired_labelled_sensitivity.csv`
  - sensitivity labelled live-bird candidate subset for deterministic repeat-pair sensitivity checks

## QC outputs

- `assay_qc_adjudication.csv`
  - one row per raw assay value that received an explicit QC disposition
  - records whether the raw value was retained for audit, excluded by explicit lab instruction, or contributed to the adjudicated modelling value
  - also records repeat-pair detection, resolution status, and the modelling-value method
- `assay_qc_affected_samples.csv`
  - one row per affected sample/tissue/isotope showing the final adjudicated outcome
- `assay_qc_manual_adjudication.csv`
  - sample-isotope cases that still require manual review because no reproducible workbook rule could be parsed
- `assay_qc_repeat_pairs.csv`
  - one row per detected repeat-analysis pair in the live batch 5 C/N workbook
  - records whether the pair was detected by green formatting, `_rpt` naming, comments, or combinations of those signals
- `assay_qc_repeat_pair_summary.csv`
  - counts of repeat-analysis pairs by detection mechanism and green-formatting pattern
- `assay_qc_sheet_rules.csv`
  - sheet-level notes and whether they were applied as rule-based adjudications or only recorded as context
- `assay_qc_summary.csv`
  - counts of excluded values, average replacements, explicit keeps, strict/sensitivity inclusions, and pending manual-review rows
- `qc_step_log.csv`
  - row-count log for each filter, join, reshape, and subset step
- `qc_issues.csv`
  - explicit warnings and information on exclusions, harmonization, and unresolved ambiguity
- `qc_missingness_summary.csv`
  - variable-level missingness summaries for key derived tables
- `qc_completeness_summary.csv`
  - isotope/tissue completeness summaries for the main modelling-ready subsets

## Important current limitations

- Museum coordinates are inherited from `data/CNOH_museum_samples_201125.rds`, which reflects the current notebook's historical geocoding state.
- Live telemetry-derived winter-distance fields and database-backed site fields are not included here because they depend on external resources outside this repository.
- `AV21986_A feather` oxygen and `AV705_D feather` hydrogen are now excluded by explicit workbook-backed lab QC, not by unexplained carry-forward notebook logic.
- `CP19311_2024-10-11_primary feather` is now retained in the assay-level audit tables as a valid repeat pair and contributes to the modelling tables via the mean of retained valid repeats, because no workbook instruction explicitly excludes either row or selects one row over the other.
