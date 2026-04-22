# Stage 4b Rung 2 S2 Uncertainty Summary

## Configuration

- Monte Carlo draws per assessment set: 500.
- Indeterminate rule: winning-class support < 0.80 or mean top-class probability < 0.50.
- Outer validation splits matched the accepted Rung 2 repeated and blocked fold structure.
- The staged workflow used the accepted R1 C model upstream and a breast C/N/H NZ-vs-AU migrant submodel downstream.
- Assessment rows did not contribute to either the upstream R1 fit or the migrant submodel fit.

## Comparison against accepted uncertainty-aware R2 C

- Rung 2 S2: Soft hierarchical R1 C then breast C/N/H NZ-vs-AU [blocked_coordinate_cv]: uncertainty log loss 0.899, uncertainty balanced accuracy 0.642, macro Brier 0.180, confidence ECE 0.197, indeterminate rate 0.605, mean switch rate 0.083.
- Rung 2 C: Breast C/N/H + paired contrasts [blocked_coordinate_cv]: uncertainty log loss 0.921, uncertainty balanced accuracy 0.589, macro Brier 0.183, confidence ECE 0.137, indeterminate rate 0.694, mean switch rate 0.061.
- Rung 2 S2: Soft hierarchical R1 C then breast C/N/H NZ-vs-AU [repeated_stratified_cv]: uncertainty log loss 0.924, uncertainty balanced accuracy 0.531, macro Brier 0.185, confidence ECE 0.111, indeterminate rate 0.627, mean switch rate 0.064.
- Rung 2 C: Breast C/N/H + paired contrasts [repeated_stratified_cv]: uncertainty log loss 0.961, uncertainty balanced accuracy 0.542, macro Brier 0.193, confidence ECE 0.124, indeterminate rate 0.711, mean switch rate 0.078.

## Interpretation notes

- R2 S2 avoids primary H entirely, so any retained advantage over R2 C should be interpreted as biologically aligned rather than H-driven.
- This run does not replace the accepted direct R2 C result; it evaluates whether the staged soft-hierarchical option remains attractive after uncertainty propagation.
- Rung 3 remains out of scope.
