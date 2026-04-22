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

Authoritative live assay sources currently folded into these tables:
- `data/Master_EA_C_N_SI_Batch 5 Breast and Primary Feathers_LEH_Repeats added_April25_v1_SB_v2.xls`
- `data/Master_TCEA_H_SI_Batch 5 breast and primary feathers_LEH_v1_SB.xls`

Note:
- The `live_cn_*` filenames are now legacy names. These tables now carry live `d13C`, `d15N`, and `d2H` where available.

- `live_phenotypes_canonical.csv`
  - one row per ring from the in-repo LIMS workbook
  - now preserves `sampling_location_general`, raw `sampling_location_specific`, and a plotting-friendly `sampling_site_label`
  - `sampling_site_label` is derived from `Sampling Location (Specific)` by taking the first comma-delimited site name, with `Kaitorete Spit` harmonized to `Kaitorete`
- `live_cn_measurements_by_sample.csv`
  - one row per live bird assay row from the current batch 5 C/N and H workbooks
  - includes raw, audit-retained, and adjudicated `d13C`, `d15N`, and `d2H` values, repeat-assay flags, phenotype joins, assay-QC fields, and site columns from the submission form
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
  - includes `d13C`, `d15N`, and `d2H` columns plus provenance showing which retained assay rows contributed to the modelling value
- `live_cn_measurements_adjudicated_by_sample_strict.csv`
  - one row per strict adjudicated live assay/sample after isotope values are recombined
- `live_cn_measurements_adjudicated_by_sample_sensitivity.csv`
  - one row per sensitivity adjudicated live assay/sample after isotope values are recombined
- `live_cn_tissue_summary.csv`
  - one row per `ring + feather_type`
  - averages adjudicated live assays within tissue
  - now includes tissue-level `d2H` plus `sampling_location_general`, `sampling_location_specific`, and `sampling_site_label`
- `live_cn_tissue_summary_strict.csv`
  - strict tissue-level live-isotope summaries for model-input use
- `live_cn_tissue_summary_sensitivity.csv`
  - sensitivity tissue-level live-isotope summaries for model-input use
- `live_cn_paired_by_ring.csv`
  - one row per ring with breast and primary `d13C`, `d15N`, and `d2H` values side by side
  - includes completeness flags for paired C/N, paired H, and paired C/N/H
  - now preserves ring-level site identity for population/site summaries and exploratory plots
- `live_cn_paired_by_ring_strict.csv`
  - strict ring-level paired live-isotope table
  - unresolved repeat-pair tissues are removed before pairing
- `live_cn_paired_by_ring_sensitivity.csv`
  - sensitivity ring-level paired live-isotope table
  - unresolved repeat-pair tissues use the deterministic non-`_rpt` default when available
- `live_screening_ready_paired_labelled.csv`
  - labelled live-bird candidate subset
  - restricted to known migratory status and complete primary + breast C/N pairs
  - now carries `d2H` columns, site identity, and C/N/H completeness flags for downstream screening decisions
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
  - one row per detected repeat-analysis pair across the current live batch 5 C/N and H workbooks
  - records whether the pair was detected by green formatting, `_rpt` naming, comments, average columns, or combinations of those signals
- `assay_qc_repeat_pair_summary.csv`
  - counts of repeat-analysis pairs by detection mechanism and green-formatting pattern
- `assay_qc_sheet_rules.csv`
  - sheet-level notes and whether they were applied as rule-based adjudications or only recorded as context
- `assay_qc_summary.csv`
  - counts of excluded values, retained repeat rows, model-contributing rows, and pending manual-review rows by isotope
- `qc_step_log.csv`
  - row-count log for each filter, join, reshape, and subset step
- `qc_issues.csv`
  - explicit warnings and information on exclusions, harmonization, and unresolved ambiguity
- `qc_missingness_summary.csv`
  - variable-level missingness summaries for key derived tables
- `qc_completeness_summary.csv`
  - isotope/tissue completeness summaries for the main modelling-ready subsets

## Stage evidence outputs

- `model_stage0_assay_repeatability_repeat_groups.csv`
  - Stage 0 technical repeat-run groups from the retained live assay ledger
  - now records source files and labels these groups explicitly as technical repeats of the same sample:tissue
- `model_stage0_assay_repeatability_summary.csv`
  - Stage 0 pooled assay-repeatability summaries by tissue and isotope
  - preserves source-file provenance for each direct technical repeat estimate
- `model_stage1_breast_feather_repeatability_detail.csv`
  - Stage 1 repeated-breast-feather detail table from museum data
  - now records source files and labels this evidence as distinct breast feathers from the same bird, not technical reruns
- `model_stage1_breast_feather_repeatability_summary.csv`
  - Stage 1 biological/tissue repeatability summary for repeated museum breast feathers
- `model_stage1_tissue_summary_adequacy.csv`
  - joined Stage 0 and Stage 1 adequacy summary for later screening stages
  - now distinguishes Stage 0 technical-repeat evidence sources from Stage 1 repeated-feather evidence sources
- `model_stage_evidence_sources.csv`
  - compact provenance table mapping each raw evidence workbook to the Stage and evidence stream it informs

## Site Validation Outputs

- `model_stage2_site_lookup.csv`
  - one row per unique live submission-form site string retained in the canonical phenotype table
  - records the raw specific location, the harmonized `sampling_site_label`, and whether any grouping ambiguity was detected
- `model_stage2_site_unit_summary.csv`
  - one row per live site/population validation unit
  - summarises phenotype coverage, paired-isotope coverage, labelled Rung 1 and Rung 2 counts, class composition, primary-reference coverage, and LOSO caveats
- `model_stage2_site_validation_completeness.csv`
  - site/coordinate completeness audit for the live phenotype, paired, labelled, and primary-reference inputs
- `model_stage2_site_block_folds.csv`
  - recommended grouped site-block outer-fold definition for later modelling work
  - built on `sampling_site_label` and balanced across the current Rung 2 class structure while keeping sites intact
- `model_stage2_site_loso_folds.csv`
  - leave-one-site-out fold definition for site stress tests
  - useful for sensitivity analysis, but many held-out sites are too small or too class-pure for stable performance estimation
- `model_stage2_site_fold_summary.csv`
  - per-fold assessment counts, site membership, and class composition for both grouped site-block and LOSO designs

## Reference Summary Outputs

- `reference_live_primary_site_summary.csv`
  - full-data descriptive primary-feather breeding-site reference summary by `sampling_site_label`
  - includes labelled and unknown-status live birds when a primary tissue row exists
  - provides counts, centroids, variances, coordinate centroids, and contributor provenance
- `reference_live_primary_site_covariance_long.csv`
  - pairwise primary-feather covariance and correlation summaries within each live site reference
- `reference_museum_winter_breast_region_summary.csv`
  - broad winter breast-feather reference summary by museum `geo_region`
  - New Zealand winter museum birds remain only New Zealand winter references here; they are not relabelled as residents or NZ migrants
- `reference_museum_winter_breast_region_covariance_long.csv`
  - pairwise covariance and correlation summaries for the broad museum winter breast reference groups
- `reference_object_manifest.csv`
  - documents each reference-summary table, its source contributors, intended interpretation, and how to use it later without leakage

## Parallel Target Extensions

- `extension_au_non_au_matched_paired_labelled.csv`
  - frozen binary-extension benchmark subset for exploratory AU-vs-non_AU screening
  - preserves the original 3-class winter target alongside the derived binary target on the same paired complete-C/N/H labelled birds
- `extension_au_non_au_breast_features.csv`
  - breast-only modelling-ready feature table for the AU-vs-non_AU extension
  - uses the same frozen ring set as the paired-contrast binary table for matched downstream comparison
- `extension_au_non_au_paired_contrast_features.csv`
  - paired-contrast modelling-ready feature table for the AU-vs-non_AU extension
  - keeps breast and primary values plus paired contrast features on the same frozen ring set
- `extension_au_non_au_candidate_definitions.csv`
  - compact description of the exploratory binary candidate feature sets and their source tables
- `extension_au_non_au_class_counts.csv`
  - overall AU-vs-non_AU class counts for the frozen binary extension subset
  - also carries the matched 3-class benchmark counts for direct reporting
- `extension_au_non_au_site_class_counts.csv`
  - site/population-level class counts for the frozen binary extension subset
  - useful for checking site imbalance before later blocked validation
- `extension_au_non_au_site_block_folds.csv`
  - binary-target copy of the accepted grouped site-block fold definitions
  - reuses the existing fold assignments without changing them
- `extension_au_non_au_site_loso_folds.csv`
  - binary-target copy of the accepted LOSO site fold definitions
  - intended for later stress-test reporting, not as a replacement for the accepted fold objects
- `extension_au_non_au_subset_check.csv`
  - ring-level audit confirming that the frozen benchmark table, the binary feature tables, and the copied fold files all use the same comparison subset

### South-Island-Only AU-vs-non_AU Extension

- `extension_au_non_au_south_island_site_eligibility.csv`
  - explicit site-level inclusion table for the South-Island-only extension branch
  - defines South Island breeder status from `sampling_site_label` plus raw `sampling_location_specific`, because `sampling_location_general` is not discriminative in the current live dataset
- `extension_au_non_au_south_island_matched_paired_labelled.csv`
  - frozen South-Island-only benchmark subset for exploratory AU-vs-non_AU screening
  - preserves the original `status`, the current 3-class winter target, and the South-Island-only binary target on the same paired complete-C/N/H labelled birds
- `extension_au_non_au_south_island_breast_features.csv`
  - breast-only modelling-ready feature table for the South-Island-only AU-vs-non_AU extension
- `extension_au_non_au_south_island_paired_contrast_features.csv`
  - paired-contrast modelling-ready feature table for the South-Island-only AU-vs-non_AU extension
- `extension_au_non_au_south_island_candidate_definitions.csv`
  - compact description of the South-Island-only binary candidate feature sets and their source tables
- `extension_au_non_au_south_island_class_counts.csv`
  - overall class counts for the South-Island-only frozen subset
  - also carries the matched 3-class benchmark counts and the full binary-extension reference counts
- `extension_au_non_au_south_island_site_class_counts.csv`
  - site/population-level class counts for the South-Island-only frozen subset
- `extension_au_non_au_south_island_site_block_folds.csv`
  - South-Island-only copy of the accepted grouped site-block fold definitions
  - preserves the original fold assignments after filtering to South Island breeding sites
- `extension_au_non_au_south_island_site_loso_folds.csv`
  - South-Island-only copy of the accepted LOSO site fold definitions
- `extension_au_non_au_south_island_fold_summary.csv`
  - per-fold assessment-set counts and binary class composition after South Island filtering
- `extension_au_non_au_south_island_subset_check.csv`
  - ring-level audit confirming that the frozen South-Island-only benchmark table, the South-Island-only binary feature tables, and the copied fold files all use the same comparison subset
- `extension_au_non_au_south_island_imbalance_comparison.csv`
  - compact comparison of class balance and site imbalance between the full AU-vs-non_AU extension and the South-Island-only branch

## Important current limitations

- Museum coordinates are inherited from `data/CNOH_museum_samples_201125.rds`, which reflects the current notebook's historical geocoding state.
- Live telemetry-derived winter-distance fields and database-backed site fields are not included here because they depend on external resources outside this repository.
- `AV21986_A feather` oxygen and `AV705_D feather` hydrogen are now excluded by explicit workbook-backed lab QC, not by unexplained carry-forward notebook logic.
- `CP19311_2024-10-11_primary feather` is now retained in the assay-level audit tables as a valid repeat pair and contributes to the modelling tables via the mean of retained valid repeats, because no workbook instruction explicitly excludes either row or selects one row over the other.
- The live H workbook differs from the C/N repeat-analysis workbook: it contains two `_rpt` repeat pairs and worksheet average columns, but no explicit sheet note instructing the analyst to prefer those averages. The build therefore retains the H repeat rows for audit and uses the mean of retained valid repeats for modelling; the displayed `Avg d2H` values match that mean for both H repeat pairs.
