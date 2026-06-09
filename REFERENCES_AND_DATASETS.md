# ContamGuard-ST — references (with code) & datasets

Consolidated reference + dataset index. Paper DOIs verified via Crossref and code
repositories via the GitHub API on 2026-06-09. See `BASELINE_REFERENCES.md` for the
full provenance and audit boundary.

## Reference papers & method baselines (with public code)

| Role | Method | Venue / year | DOI | Code |
|------|--------|--------------|-----|------|
| Primary | SPLIT — Resolving sensitivity, specificity and signal contamination in Xenium ST | Nature Methods 2026 | `10.1038/s41592-026-03089-8` | https://github.com/bdsc-tds/SPLIT |
| Companion | Xenium analysis pipeline | — | — | https://github.com/bdsc-tds/xenium_analysis_pipeline |
| Companion | Bilous2026 reproducibility | — | — | https://github.com/bdsc-tds/Bilous2026 |
| Baseline | proseg — probabilistic cell segmentation | — | — | https://github.com/dcjones/proseg |
| Baseline | Baysor — Bayesian molecule segmentation | — | — | https://github.com/kharchenkolab/Baysor |
| Baseline | FastReseg — segmentation-error correction | — | — | https://github.com/Nanostring-Biostats/FastReseg |

## Datasets

This repo runs on user-supplied local imaging-ST inputs (no shipped dataset catalog).
Reference data comes from the baselines above:
- Xenium control / contamination cohorts (e.g. `xenium_cervix_controls`).
- SPLIT vignette Xenium panels (183-gene panel; ATERA comparison from ~1,600 shared genes).

> Verification: SPLIT DOI confirmed in Crossref; SPLIT / proseg / Baysor / FastReseg /
> xenium_analysis_pipeline repos confirmed live via GitHub API (2026-06-09).
