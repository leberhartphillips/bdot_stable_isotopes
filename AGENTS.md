# AGENTS.md

## Project scope
This repository develops a reproducible stable-isotope screening pipeline for migratory phenotype in banded dotterels, with emphasis on defensible validation, uncertainty, and biological interpretability.

## Working directory rules
- Work only inside this repository.
- Do not access or modify files outside this repository unless I explicitly approve it.
- Prefer small, reviewable diffs over broad rewrites.
- Do not delete files unless I explicitly ask.

## Data safety rules
- Never modify raw data files.
- Treat all raw lab spreadsheets and original source files as immutable source data, even if they are not stored in `data/raw`.
- Save cleaned or derived outputs to `data/derived` or another clearly named derived-data folder.
- Preserve traceability from every derived table back to the original source files.
- Never silently drop rows.
- Every row loss must be counted and explained.
- Flag empty joins, duplicated IDs, conflicting labels, and unexpected missingness.

## Scientific and statistical rules
- Do not invent citations, methods, or package behavior.
- Do not make scientific or statistical decisions silently.
- When uncertain, state the uncertainty explicitly.
- Prioritize correctness and robustness over optimism.
- Optimize for defensible screening performance, not maximal apparent classification accuracy.
- Be critical about leakage, overfitting, circularity, and misuse of imputation.
- Keep accepted results separate from exploratory extension branches unless a new approach clearly justifies updating conclusions.

## Repository-specific context
- Exploratory analysis notebook: `bdot_stable_isotopes2.qmd`
- Report-focused integration document: `bdot_stable_isotope_screening_report.qmd`
- The repository began messy and exploratory; keep non-destructive habits.
- Raw stable isotope data originate from Excel files and are wrangled into modelling-ready formats in R.
- The project includes museum and live/field datasets.
- Breast feathers reflect wintering-ground signal; primary/flight feathers reflect breeding-ground signal.

## Expected workflow for substantial tasks
For any substantial task, return:
1. Plan
2. Actions taken
3. Verification
4. Unresolved issues

## Modelling guardrails
- Compare a small set of justified candidate models rather than jumping to one approach.
- Report assumptions clearly.
- Use validation appropriate for the sample size, class structure, and site structure.
- Report class-wise performance, confusion matrices, calibration, indeterminate burden, and uncertainty where relevant.
- Explicitly state where inference is reliable and where it is not.
- Keep preprocessing, scaling, feature construction, and tuning inside resampling when required to avoid leakage.
- Site/population-aware blocking is now an important validation layer and should be preferred over weaker spatial shortcuts when the question is transfer across sites.
- Reference or site-informed features must be built inside training folds only.

## Current accepted project conclusions
- Rung 1 is currently the strongest and most defensible part of the framework.
- The accepted Rung 1 operational model is the paired-contrast C/N approach (`R1 C`), with an indeterminate option.
- Rung 2 remains promising but not operational.
- Under stricter site-aware validation, the staged Rung 2 approach (`R2 S2`) is the strongest current Rung 2 direction.
- Direct Rung 2 alternatives and site/reference-informed extensions may still be useful as exploratory comparisons, but they should not replace the accepted interpretation without clear evidence.
- Rung 3 remains exploratory only.

## Project-specific modelling rules
- For Rung 1, paired breast-primary contrast models must be explicitly compared against breast-only and primary-only baselines.
- For Rung 2 and Rung 3, breast-feather winter signal is the default basis of classification unless paired features show validated incremental value.
- Museum samples are support/reference data unless a task explicitly tests a new modelling role for them.
- Unknown-status live birds must not be given pseudo-labels; they may be used for prediction-only outputs or as unlabeled/reference contributors where methodologically justified.
- Direct AU-vs-non-AU simplifications have been tested and are not currently operationally useful; treat them as secondary exploratory results unless explicitly revisited.
- Site-aware and reference-informed extension branches should remain clearly separated from accepted mainline results unless they clearly outperform under appropriate validation.
- Continuous isoscapes/geospatial surfaces are exploratory only for now and are not the primary operational model for this repository.

## Writing rules
- Keep prose direct and non-hyped.
- Do not overclaim model performance.
- Separate evidence from interpretation.
- For unlabelled birds, use language like “predicted status of unlabelled birds” and avoid wording that implies validation.