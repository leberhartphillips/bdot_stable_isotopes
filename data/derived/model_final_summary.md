# Final Modelling Summary

## Current conclusion

- Best-supported Rung 1 model: `R1 C` paired-contrast C/N. It stayed strongest across the paired C/N candidates and remained defensible after assay-uncertainty propagation with an explicit indeterminate rule.
- Current status of Rung 2: promising but not operational. Two candidates now lead this rung: `R2 C`, the direct breast C/N/H plus paired-contrast model, and `R2 S2`, the staged hierarchical model that first uses the accepted Rung 1 migrant signal and then uses breast C/N/H to separate NZ from AU migrants. `R2 S2` is especially attractive because it avoids direct dependence on primary H and is easier to explain biologically, but neither candidate is operational yet.
- Current status of Rung 3: exploratory only. Class sizes, blocking strength, and uncertainty support are still too weak for operational screening.

## Final Rung Status

| rung | status | lead_model |
| --- | --- | --- |
| Rung 1 | operational | R1 C: Paired-contrast C/N |
| Rung 2 | promising but not operational | Co-leading: R2 C direct breast C/N/H + paired contrasts; R2 S2 staged hierarchical R1 C then breast C/N/H NZ-vs-AU |
| Rung 3 | exploratory only | None |

## Plain-language Rung 2 logic

The direct `R2 C` model tries to solve the full three-class problem in one step by combining breast C/N/H with paired seasonal contrasts. The staged `R2 S2` model breaks the question into two biologically simpler steps: first, use the accepted Rung 1 paired-contrast C/N model to ask whether a bird looks migrant-like at all; second, for birds with migrant-like signal, use breast C/N/H to separate NZ from AU. This staged logic is easier to explain and reduces reliance on primary H, but it still needs the same conservative uncertainty treatment and remains provisional.

## Leading Models

| rung | role | model | repeat_log_loss | blocked_log_loss | repeat_bal_acc | blocked_bal_acc | repeat_indeterminate | blocked_indeterminate |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Rung 1 | best-supported operational model | Rung 1 C: Paired-contrast C/N | 0.512 | 0.494 | 0.735 | 0.764 | 0.281 | 0.283 |
| Rung 1 | secondary non-H comparator | Rung 1 D: Structured paired C/N | 0.533 | 0.509 | 0.674 | 0.719 | 0.278 | 0.373 |
| Rung 2 | non-H benchmark | Rung 2 A: Breast-only C/N | 1.094 | 1.095 | 0.334 | 0.333 | 1.000 | 1.000 |
| Rung 2 | co-leading direct candidate | Rung 2 C: Breast C/N/H + paired contrasts | 0.961 | 0.921 | 0.542 | 0.589 | 0.711 | 0.694 |
| Rung 2 | co-leading staged candidate (avoids primary H) | Rung 2 S2: Soft hierarchical R1 C then breast C/N/H NZ-vs-AU | 0.924 | 0.899 | 0.531 | 0.642 | 0.627 | 0.605 |

## Indeterminate Rules

| rung | component | threshold |
| --- | --- | --- |
| Rung 1 | winning-class support | < 0.80 |
| Rung 1 | probability-boundary crossing | central 80% interval crosses 0.5 |
| Rung 2 | winning-class support | < 0.80 |
| Rung 2 | mean top-class probability | < 0.50 |

Indeterminate calls matter operationally because they prevent the pipeline from turning weakly supported or instability-prone predictions into false certainty.

## Key Evidence Limitations

- Primary H remains the main H-related weak point: Primary H still relies on a breast-H proxy variance and lacks direct retained repeat support, so H-driven gains need caution.
- No canonical site field for stronger blocking: Current blocked validation uses coordinate groups as a proxy rather than a canonical site variable, which weakens the strength of deployment-relevant blocking.
- Multinomial sample-size limits remain important: Some multinomial training splits still have very small class counts, which increases instability and limits how strongly Rung 2 and Rung 3 can be claimed.
- Conservative indeterminate rules reduce immediate coverage: The current rules are intentionally cautious, so promising models can still yield high indeterminate rates and limited immediate operational coverage.

## Ranked Next Improvements

| rank | improvement |
| --- | --- |
| 1 | Collect direct primary-H technical repeatability and repeat-primary biological evidence |
| 2 | Add a canonical site field to live modelling tables |
| 3 | Increase labelled live sample size for Rung 2 and Rung 3, especially NZ subgroups |
| 4 | Assemble an independent external or temporal holdout set |
| 5 | Evaluate additional isotopes only after repeatability support is in place |
