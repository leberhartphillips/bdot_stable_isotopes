# Stage 5b Hybrid Staged Rung 2 Summary

This exploratory candidate keeps all accepted and Stage 5 results intact.

## Hybrid structure

- Upstream stage: a site-aware Rung 1 model using paired-contrast C/N plus a training-fold-only primary-site reference distance.
- Downstream stage: a non-site-aware migrant NZ-vs-AU breast C/N/H model.
- Held-out assessment sites never contribute to the upstream primary-site reference cloud.

## Pooled comparison

- Grouped site-block pooled log loss: 1.197
- Grouped site-block pooled balanced accuracy: 0.299
- LOSO pooled log loss: 1.182
- LOSO pooled balanced accuracy: 0.292

See `model_stage5b_hybrid_vs_siteaware_coleaders.csv` for the direct comparison against `R2 S2` and `R2 C`.
