# bdot_stable_isotopes

Stable-isotope screening pipeline work for migratory phenotype in banded dotterels.

## Current analysis status

- Canonical active analysis file: `bdot_stable_isotopes2.qmd`
- Legacy analysis file retained for reference only: `bdot_stable_isotopes.qmd`
- This cleanup pass is documentation-only. It does not refactor scientific logic, move existing files, or modify raw inputs.

## Working conventions for this repository

- Treat existing source spreadsheets and other raw inputs in `data/` as immutable.
- Put new derived outputs in `data/derived/` going forward unless there is a clearly documented reason to write elsewhere.
- Keep traceability from each derived object back to its source file(s) and analysis step(s).
- Count and explain row loss, empty joins, duplicated IDs, conflicting labels, and unexpected missingness when the workflow is refactored further.

## Current workflow snapshot

- The active notebook is `bdot_stable_isotopes2.qmd`.
- `bdot_stable_isotopes2.qmd` currently mixes imports, wrangling, live geocoding, modeling, and reporting in one document.
- The notebook currently reads and writes a derived checkpoint at `data/CNOH_museum_samples_201125.rds`.
- A bootstrap ROC object is currently read from `output/roc_all.rds`.
- A conservative inventory of in-repo data files is documented in `DATA_MANIFEST.md`.

## External dependencies currently used by the canonical QMD

These are documented only in this pass. They are not being replaced yet.

- Google Sheet used in sample management via `googlesheets4::read_sheet()`
- Sibling-repo Excel file at `../../stable_isotopes/data/for_Sarah/master_museum_sample_banded_dotterels.xlsx`
- External telemetry RDS at `/Users/luketheduke2/ownCloud/kemp_projects/bdot/R_projects/bdot_db/output/bdot_argos_proc_2024_lonE.rds`
- Database-backed fields from `FIELD_2024_BADOatNZ` via `dbo::dbcon()` and `dbq(con, "SELECT * FROM CAPTURES")`
- Live geocoding calls inside `bdot_stable_isotopes2.qmd`; these should eventually be replaced by a committed lookup table

## Documentation added in this pass

- `README.md` updated with current repository status and workflow boundaries
- `DATA_MANIFEST.md` added with conservative file classification and dependency notes
- `data/derived/README.md` added to reserve a home for future derived outputs
