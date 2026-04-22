# Stage 3 Rung 1 Uncertainty Summary

## Configuration

- Monte Carlo draws per assessment set: 500.
- Indeterminate rule: class support < 0.80 or central 80% migrant-probability interval crosses 0.5.
- Training/preprocessing remained inside each outer resample split; assay working SDs were read from the Stage 0 pooled-variance table.

## Model summary

- Rung 1 C: Paired-contrast C/N [blocked_coordinate_cv]: uncertainty log loss 0.494, uncertainty balanced accuracy 0.764, indeterminate rate 0.283, final-class change rate 0.000, mean switch rate 0.111.
- Rung 1 D: Structured paired C/N [blocked_coordinate_cv]: uncertainty log loss 0.509, uncertainty balanced accuracy 0.719, indeterminate rate 0.373, final-class change rate 0.000, mean switch rate 0.102.
- Rung 1 G: Structured paired C/N/H [blocked_coordinate_cv]: uncertainty log loss 0.553, uncertainty balanced accuracy 0.647, indeterminate rate 0.217, final-class change rate 0.000, mean switch rate 0.074.
- Rung 1 C: Paired-contrast C/N [repeated_stratified_cv]: uncertainty log loss 0.512, uncertainty balanced accuracy 0.735, indeterminate rate 0.281, final-class change rate 0.000, mean switch rate 0.087.
- Rung 1 D: Structured paired C/N [repeated_stratified_cv]: uncertainty log loss 0.533, uncertainty balanced accuracy 0.674, indeterminate rate 0.278, final-class change rate 0.000, mean switch rate 0.094.
- Rung 1 G: Structured paired C/N/H [repeated_stratified_cv]: uncertainty log loss 0.534, uncertainty balanced accuracy 0.676, indeterminate rate 0.276, final-class change rate 0.000, mean switch rate 0.089.

## Interpretation notes

- Breast H support is now treated as explicitly documented by Stage 1 museum triplicate evidence, but primary H remains a proxy-variance, high-caution component.
- The optional C/N/H model is retained only as a tertiary sensitivity analysis and not as a preferred operational target.
- No Rung 2 uncertainty propagation was run here.
