# AGENTS.md

## Project scope
This repository develops a reproducible stable-isotope screening pipeline for migratory phenotype in banded dotterels.

## Working directory rules
- Work only inside this repository.
- Do not access or modify files outside this repository unless I explicitly approve it.
- Prefer small, reviewable diffs over broad rewrites.
- Do not delete files unless I explicitly ask.

## Data safety rules
- Never modify raw data files.
- Treat all files in `data/raw` as immutable source data.
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
- Be critical about leakage, overfitting, and misuse of imputation.

## Repository-specific context
- Main working Quarto file so far: `bdot_stable_isotopes2.qmd`
- The repository is currently messy and exploratory.
- Before large edits, inspect the existing structure and propose a non-destructive plan.
- Raw stable isotope data originate from Excel files and are wrangled into modelling-ready formats in R.
- The project includes museum and field datasets.
- Breast feathers reflect wintering-ground signal; flight feathers reflect breeding-ground signal.

## Expected workflow for substantial tasks
For any substantial task, return:
1. Plan
2. Actions taken
3. Verification
4. Unresolved issues

## Modelling guardrails
- Compare a small set of justified candidate models rather than jumping to one approach.
- Report assumptions clearly.
- Use validation that is appropriate for the sample size and class structure.
- Report class-wise performance, confusion matrices, and uncertainty.
- Explicitly state where inference is reliable and where it is not.

## Writing rules
- Keep prose direct and non-hyped.
- Do not overclaim model performance.
- Separate evidence from interpretation.