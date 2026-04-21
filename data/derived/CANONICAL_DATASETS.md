# Canonical Derived Datasets

This directory contains non-destructive, analysis-ready tables built from in-repo inputs used by `bdot_stable_isotopes2.qmd`.

## Museum data

- `museum_measurements_long.csv`
  - one row per isotope assay
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
  - includes repeat-assay flags and phenotype joins
- `live_cn_tissue_summary.csv`
  - one row per `ring + feather_type`
  - averages repeated assays within tissue
- `live_cn_paired_by_ring.csv`
  - one row per ring with breast and primary C/N values side by side
  - includes completeness flags and tissue deltas
- `live_screening_ready_paired_labelled.csv`
  - labelled live-bird candidate subset
  - restricted to known migratory status and complete primary + breast C/N pairs

## QC outputs

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
- The explicit exclusion of `AV21986_A feather` in triplicate oxygen data is preserved from the current QMD but remains scientifically unresolved.
