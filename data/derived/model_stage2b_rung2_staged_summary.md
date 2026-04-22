# Staged Rung 2 Summary

## Leakage control

- The probability-augmented direct model used cross-fitted Rung 1 migrant probabilities inside each outer training split and a separate outer-trained Rung 1 fit for each assessment split.
- The soft hierarchical model fit the accepted Rung 1 model and the NZ-vs-AU migrant submodel on the outer analysis set only; no assessment rows contributed to either fit.

## Staged candidates

- Rung 2 S1: Breast C/N/H + out-of-fold R1 migrant probability: strategy family `probability_augmented_direct`, caution `staged_r1c_prob_plus_breast_cnh_no_primary_h`.
- Rung 2 S2: Soft hierarchical R1 C then breast C/N/H NZ-vs-AU: strategy family `soft_hierarchical_two_step`, caution `staged_r1c_then_breast_cnh_no_primary_h`.

## Headline comparison versus accepted direct Rung 2 models

- Best repeated-CV log-loss rank: Rung 2 S2: Soft hierarchical R1 C then breast C/N/H NZ-vs-AU (log loss 0.924, balanced accuracy 0.532).
- Best blocked-CV log-loss rank: Rung 2 S2: Soft hierarchical R1 C then breast C/N/H NZ-vs-AU (log loss 0.903, balanced accuracy 0.617).

## Interpretation notes

- These staged approaches were evaluated as additional Stage 2-style candidates only. They do not replace the accepted Stage 4 uncertainty-aware Rung 2 results.
- Deterministic weak-support rates are reported here from top-class probabilities < 0.5, but Stage 4-style uncertainty-aware indeterminate behavior was not re-estimated for these staged variants.
- The soft hierarchical model avoids primary H entirely; any gain from that variant should therefore be interpreted as improved biological alignment rather than as proof that the accepted direct R2 C model should be displaced.
