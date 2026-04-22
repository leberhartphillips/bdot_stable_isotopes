# Stage 4 Rung 2 Uncertainty Summary

## Configuration

- Monte Carlo draws per assessment set: 500.
- Indeterminate rule: winning-class support < 0.80 or mean top-class probability < 0.50.
- Training/preprocessing remained inside each outer resample split; assay working SDs were read from the Stage 0 pooled-variance table.
- R2 C includes primary H under proxy variance only and should remain high-caution even if it performs better.

## Model summary

- Rung 2 C: Breast C/N/H + paired contrasts [blocked_coordinate_cv]: uncertainty log loss 0.921, uncertainty balanced accuracy 0.589, macro Brier 0.183, confidence ECE 0.137, indeterminate rate 0.694, mean switch rate 0.061.
- Rung 2 A: Breast-only C/N [blocked_coordinate_cv]: uncertainty log loss 1.095, uncertainty balanced accuracy 0.333, macro Brier 0.221, confidence ECE 0.032, indeterminate rate 1.000, mean switch rate 0.000.
- Rung 2 C: Breast C/N/H + paired contrasts [repeated_stratified_cv]: uncertainty log loss 0.961, uncertainty balanced accuracy 0.542, macro Brier 0.193, confidence ECE 0.124, indeterminate rate 0.711, mean switch rate 0.078.
- Rung 2 A: Breast-only C/N [repeated_stratified_cv]: uncertainty log loss 1.094, uncertainty balanced accuracy 0.334, macro Brier 0.221, confidence ECE 0.020, indeterminate rate 1.000, mean switch rate 0.002.

## Interpretation notes

- R2 A remains the non-H benchmark for winter-signal classification.
- R2 C is the only higher-information candidate carried forward here, but any apparent gain must survive uncertainty propagation under explicit primary-H caution.
- Rung 3 remains out of scope.
