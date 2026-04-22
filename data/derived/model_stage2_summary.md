# Stage 2 Summary

## Validation design

- Primary comparison: repeated stratified 5-fold CV with 10 repeats, fit separately for Rung 1 and Rung 2.
- Sensitivity layer: exact-coordinate blocked 5-fold assignment using available capture coordinates, balanced against a global class-and-size objective on Rung 2 classes.
- Model family: ridge-penalized logistic / multinomial regression with lambda selected inside each training split via `cv.glmnet`.
- Preprocessing: deterministic breast-primary contrasts were constructed observation-wise, while all scaling and paired-distance standardization were estimated inside each training split only.

## Rung 1 shortlist

- Rung 1 C: Paired-contrast C/N: mean log loss 0.517, mean balanced accuracy 0.721, caution `paired_cn_benchmark`.
- Rung 1 D: Structured paired C/N: mean log loss 0.530, mean balanced accuracy 0.682, caution `paired_cn_benchmark`.
- Rung 1 G: Structured paired C/N/H: mean log loss 0.542, mean balanced accuracy 0.690, caution `primary_h_provisional_pending_triplicate_check`.
- Rung 1 F: Paired-contrast C/N/H: mean log loss 0.553, mean balanced accuracy 0.664, caution `primary_h_provisional_pending_triplicate_check`.
- Rung 1 B: Primary-only C/N: mean log loss 0.660, mean balanced accuracy 0.502, caution `primary_cn_limited_biological_support`.
- Rung 1 H: Primary-only C/N/H: mean log loss 0.662, mean balanced accuracy 0.504, caution `primary_h_provisional_pending_triplicate_check`.
- Rung 1 A: Breast-only C/N: mean log loss 0.670, mean balanced accuracy 0.500, caution `baseline_cn_benchmark`.
- Rung 1 E: Breast-only C/N/H: mean log loss 0.670, mean balanced accuracy 0.500, caution `breast_h_provisional_pending_triplicate_check`.

## Rung 2 shortlist

- Rung 2 C: Breast C/N/H + paired contrasts: mean log loss 0.942, mean balanced accuracy 0.540, caution `primary_h_provisional_pending_triplicate_check`.
- Rung 2 D: Breast C/N/H + paired distance: mean log loss 1.092, mean balanced accuracy 0.329, caution `primary_h_provisional_pending_triplicate_check`.
- Rung 2 A: Breast-only C/N: mean log loss 1.093, mean balanced accuracy 0.334, caution `baseline_cn_benchmark`.
- Rung 2 B: Breast-only C/N/H: mean log loss 1.095, mean balanced accuracy 0.326, caution `breast_h_provisional_pending_triplicate_check`.

## Stage 3 recommendations

- rung1 Rung 1 C: Paired-contrast C/N: Advance to Stage 3 as the main non-H benchmark. Primary CV rank 1, blocked rank 3.
- rung1 Rung 1 D: Structured paired C/N: Advance to Stage 3 as a secondary non-H paired benchmark because it remained near-leading across both validation layers. Primary CV rank 2, blocked rank 1.
- rung2 Rung 2 C: Breast C/N/H + paired contrasts: Advance only conditionally, alongside the strongest non-H benchmark, because H evidence is provisional and primary H remains high-caution. Primary CV rank 1, blocked rank 1.
- rung2 Rung 2 A: Breast-only C/N: Advance to Stage 3 as the main non-H benchmark. Primary CV rank 3, blocked rank 3.

## Constraint notes

- Breast H entered candidate models by default, but any gains should still be interpreted with technical caution and remain provisional pending confirmation that `data/Master_TCEA_H_SI_Batch 1_triplicate breast feathers_LEH_Nov25_SB.xls` is fully incorporated into the H evidence base.
- Any model using paired H contrasts or primary H was retained only as a higher-caution comparator because Stage 0 could not directly estimate retained primary-H technical repeatability.
- No Stage 3 uncertainty propagation or final operational classifier was implemented in this script.
