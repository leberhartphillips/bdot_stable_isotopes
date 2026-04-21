# Data Manifest

This is a conservative, documentation-only inventory of current repository data as of 2026-04-21.

- Canonical active analysis file: `bdot_stable_isotopes2.qmd`
- Legacy analysis file: `bdot_stable_isotopes.qmd`
- No files were moved, deleted, or reclassified on disk in this pass.
- Labels below are provisional and meant to support cleanup planning, not to make irreversible provenance claims.

## Status labels used here

- `likely raw`: likely source data or direct lab/export input that should be treated as immutable
- `likely derived`: likely output created from one or more upstream inputs
- `likely reference`: likely administrative, template, planning, or support material
- `likely superseded`: likely older version, historical copy, or interim file retained for traceability

## External dependencies currently used by `bdot_stable_isotopes2.qmd`

| Dependency | Current use in canonical QMD | Notes |
| --- | --- | --- |
| Google Sheets auth | `gs4_auth()` at lines 89-92 | Render depends on interactive/local auth state |
| Google Sheet | `read_sheet()` at line 1979 | Used for sample management input |
| Sibling-repo Excel | `../../stable_isotopes/data/for_Sarah/master_museum_sample_banded_dotterels.xlsx` at lines 1830 and 2151 | Outside this repository |
| External telemetry RDS | `/Users/luketheduke2/ownCloud/kemp_projects/bdot/R_projects/bdot_db/output/bdot_argos_proc_2024_lonE.rds` at line 4135 | Outside this repository |
| Database-backed fields | `dbo::dbcon()` / `dbq(con, "SELECT * FROM CAPTURES")` at lines 4424-4442 | Depends on local DB access and schema |

## Live geocoding currently performed inside `bdot_stable_isotopes2.qmd`

These calls should eventually be replaced by a committed lookup table, but are only documented here in this pass.

- `geocode(unique_locations)` around lines 359, 1854, 2165, and 2371
- `geocode(museum_sample_loc)` around line 1998

## Likely derived files

| Path | Used by canonical QMD | Notes |
| --- | --- | --- |
| `data/CNOH_museum_samples_111125.rds` | no | Likely historical derived checkpoint |
| `data/CNOH_museum_samples_181125.rds` | no | Likely historical derived checkpoint |
| `data/CNOH_museum_samples_201125.rds` | yes | Likely current derived checkpoint; read at lines 462, 2172, and 2449, and written at line 2435 |
| `data/museum_samples_to_analyze_161224.csv` | yes | Likely derived sample-management export written at line 2078 |

## Likely raw files

| Path | Used by canonical QMD | Notes |
| --- | --- | --- |
| `data/LIMS_251105_MaxPlanck_Eberhart-Hertel.xlsx` | yes | Read at line 4052 |
| `data/MASTER_museum_ALL SAMPLES INVENTORY_banded_dotterels_v4.xlsx` | yes | Read at lines 333 and 2343 |
| `data/Master_EA_C_N_SI_Batch 3 Breast feathers_shaft and barb comparison_LEH_Sept25_v1_SB.xls` | yes | Read at line 2117 |
| `data/Master_EA_C_N_SI_Batch 3 vane_barb and Batch 4 Breast Feathers_LEH_Sept25_v1_SB.xls` | yes | Read at line 2234 |
| `data/Master_EA_C_N_SI_Batch 4 Breast Feathers_LEH_Sept25_v1_SB_v2-1.xls` | no | Likely source/interim lab file; not currently read |
| `data/Master_EA_C_N_SI_Batch 4 Breast Feathers_LEH_Sept25_v1_SB_v2.xls` | no | Likely source/interim lab file; not currently read |
| `data/Master_EA_C_N_SI_Batch 5 Breast and Primary Feathers_LEH_Dec25_SB.xls` | no | Likely source/interim lab file; not currently read |
| `data/Master_EA_C_N_SI_Batch 5 Breast and Primary Feathers_LEH_Dec25_SB_v2.xls` | no | Likely source/interim lab file; not currently read |
| `data/Master_EA_C_N_SI_Batch 5 Breast and Primary Feathers_LEH_Repeats added_April25_v1_SB_v2.xls` | yes | Read at line 4081 |
| `data/Master_SI_Replicate test Feathers_N_9_LEH_Jan25_v1_SB.xls` | yes | Read at line 1772 |
| `data/Master_TCEA_H_SI_BATCH 4 BREAST Feathers_LE_Nov25_v1.xls` | yes | Read at line 2303 |
| `data/Master_TCEA_H_SI_Batch 1_triplicate breast feathers_LEH_Nov25_SB.xls` | yes | Read at line 265 |
| `data/Master_TCEA_O_SI_Batch 4 Breast Feathers_Batch 1 Replication_LEH_Sept25_SB.xls` | no | Likely source/interim lab file; not currently read |
| `data/Master_TCEA_O_SI_Breast Feathers_Replication_Batch 4_LEH_Sept25_SB.xls` | yes | Read at lines 293 and 2274 |
| `data/RAW DATA_SI_Feathers_LEH_Nov25.xls` | no | Likely raw source; referenced in commented code at lines 4078-4080 |
| `data/old file versions/Master_SI_Feathers_LEH_Sept24_v1_SB_v3_sample wts added.xls` | yes | Read at line 234 from an `old file versions` path; treated as immutable historical raw input in this pass |

## Likely reference files

| Path | Used by canonical QMD | Notes |
| --- | --- | --- |
| `data/Charges for first contract.xlsx` | no | Administrative/reference |
| `data/LIMS Sample Submission Template_updated Nov24_FINAL.xlsx` | no | Template/reference |
| `data/LIMS_analysis_BATCH 4 Breast and Primary_LEH_Aug25.xlsx` | no | Analysis support/reference |
| `data/Luke Eberhart Hertel raw data_UrFLMA_250520_shaft_barb feather tests_v2.xls` | no | Reference/supporting work file; not currently read |
| `data/Luke Eberhart Hertel_LIMS Sample Submission_Otago_Canterbury museum_Banded dotterel feathers_v6_weight.xlsx` | no | Submission/reference |
| `data/Luke Eberhart_Hertel_Summary analyses FY2024_25 and FY2025_26_V3_12 Nov 2025.xlsx` | no | Summary/reference |
| `data/copy_museum_samples_to_analyze_161224_feather nos_weights_v3.xlsx` | no | Derived support/reference workbook based on sample-management outputs |

## Likely superseded files

| Path | Used by canonical QMD | Notes |
| --- | --- | --- |
| `data/Luke Eberhart Hertel raw data_UrFLMA_250520_shaft_barb feather tests.xls` | no | Older version beside a `v2` file |
| `data/Luke Eberhart_Hertel_Summary analyses FY2024_25 and FY2025_26_V2_Sept 2025_v2.xlsx` | no | Older version beside a `V3` file |
| `data/Master_SI_Feathers_LEH_Sept24_v1_SB_v3_sample wts added._13 Jan re-analysis strategy_v2.xls` | no | Older version beside a `v3` file |
| `data/Master_SI_Feathers_LEH_Sept24_v1_SB_v3_sample wts added._13 Jan re-analysis strategy_v3.xls` | no | Historical re-analysis strategy file retained for traceability |
| `data/museum_samples_to_analyze_161224_feather nos_weights_v2.xlsx` | no | Older version beside a `v3` file |
| `data/old file versions/MASTER_museum_ALL SAMPLES INVENTORY_banded_dotterels_v3.xlsx` | no | Historical older inventory version |
| `data/old file versions/Master_SI_Feathers_LEH_Sept24_v1_SB_v3_sample wts added._13 Jan re-analysis strategy.xls` | no | Historical older version |
| `data/old file versions/museum_samples_to_analyze_161224.xlsx` | no | Historical older version |

## System artifacts

| Path | Notes |
| --- | --- |
| `data/.DS_Store` | macOS filesystem metadata; not project data |

## Notes on current checkpoint authority

- The active QMD clearly reads and writes `data/CNOH_museum_samples_201125.rds`, so this is the best current candidate for the operative derived checkpoint.
- The earlier dated `CNOH_museum_samples_111125.rds` and `CNOH_museum_samples_181125.rds` are retained here as likely historical checkpoints.
- This pass does not rename, move, or delete any of these files.
