# Stage 5 Site-Aware Validation and Reference Extension Summary

This extension branch keeps all previously accepted results intact and adds a separate site-aware validation pass.

## Validation layers

- Primary new evaluation layer: grouped 4-fold site-block validation (`sampling_site_label`).
- Secondary stress-test layer: leave-one-site-out validation.
- Formal site-level partial-pooling models were not fitted here because the held-out-site validation target and the small number of labelled birds per site make that extension unstable and difficult to interpret.

## Grouped site-block highlights

- Best grouped-site Rung 1 pooled log loss: r1_c_paired_contrast_cn = 0.533.
- Best grouped-site Rung 2 pooled log loss: r2s_b_soft_hierarchical_r1c_then_breast_cnh = 1.066.
- Best grouped-site Rung 2 pooled balanced accuracy: r2x_b_staged_r1c_then_breast_cnh_plus_winter_reference = 0.406.

## LOSO stress-test highlights

- LOSO should be read cautiously because most held-out sites do not contain all Rung 2 classes.
- Best LOSO Rung 2 pooled log loss: r2s_b_soft_hierarchical_r1c_then_breast_cnh = 1.091.

## Extension scope

- `R1 X` adds a training-fold-only primary-site reference distance built from live primary feathers. This is where the unknown-status live birds contribute: they help characterize breeding-site primary structure without being given phenotype labels.
- `R2 X` adds museum winter breast reference distances to the direct Rung 2 model.
- `R2 X2` adds the same museum winter reference distances to the staged `R2 S2` migrant submodel.

See the CSV outputs in `data/derived/` for pooled summaries, split metrics, confusion tables, and the comparison against previous validation layers.
