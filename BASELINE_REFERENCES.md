# ContamGuard-ST baseline references

Verification date: 2026-06-07

## Baseline decision summary

| Role | Baseline | Decision |
|---|---|---|
| Primary | SPLIT / Xenium contamination analysis | Use as the first open-code/public-artifact starting point because it directly matches this track's input-output problem. |
| Secondary | See table below | Use only for comparison, adapter design, and ablation inspiration; do not copy implementation. |

## Primary baseline

- Paper title: Resolving sensitivity, specificity and signal contamination in spatial transcriptomics
- Venue/date: Nature Methods, 2026
- Article URL: https://www.nature.com/articles/s41592-026-03089-8
- Code/artifact URL: https://github.com/bdsc-tds/SPLIT
- Verification date: 2026-06-07
- Default branch/artifact: main
- Observed HEAD SHA or DOI: f226987f0fecf2847646bda50a64886e0ea432a9
- Local audit checkout/artifact: `baselines/SPLIT-original`
- License note: GPL-3.0 license file observed locally
- Local use: Primary public-code reference for estimating signal contamination, sensitivity/specificity limits, and blank/negative-control-aware molecule reliability.
- Gate-2 status (2026-06-08): `REFERENCE_REPORTED` only for this track round. The frozen vignette reports a Xenium purification example with a 183-gene Xenium panel, an ATERA comparison/reference derived from ~1,600 shared genes, and `chunk_size=10000`; it describes lower but detectable transcript spillover and spatial diffusion-aware purification, but it does not provide a same-card numeric parity metric for `data/processed/contamguard_xenium_cervix_controls/anndata.h5ad`.
- Same-card run status: not run; needs META to build an isolated R/Bioconductor/RCTD/Seurat environment for the GPL package if a `RAN` external baseline is required.
- Fallback: If this public code/artifact becomes unavailable, mark this track `deferred-unverified` until a comparable open-code baseline is found.
- Verification command/evidence:
  - `git ls-remote --symref <repo> HEAD` for GitHub/Git/Bioconductor repositories.
  - `zenodo.org/api/records/<record>` plus local file checks for Zenodo artifacts.

## Secondary verified references

| Baseline | Role | Code/artifact URL | Local audit checkout | Branch/artifact | Observed SHA/DOI | License note |
|---|---|---|---|---|---:|---|
| Xenium analysis pipeline | Best-practice Xenium workflow companion | https://github.com/bdsc-tds/xenium_analysis_pipeline | `baselines/xenium_analysis_pipeline-original` | `main` | `5e61b1c5b1a1` | GPL-3.0 license file observed locally |
| Bilous2026 reproducibility | Paper-linked reproducibility repository | https://github.com/bdsc-tds/Bilous2026 | `not cloned yet` | `main` | `adc9563bed26` | BSD-3-Clause reported by GitHub API |
| proseg | Probabilistic cell segmentation for in situ ST | https://github.com/dcjones/proseg | `not cloned yet` | `main` | `bfecb6fe33d2` | license requires re-check |
| Baysor | Bayesian molecule-based segmentation | https://github.com/kharchenkolab/Baysor | `not cloned yet` | `cpp` | `b6e289565d74` | MIT reported by GitHub API |
| FastReseg | Transcript-profile-based segmentation error correction | https://github.com/Nanostring-Biostats/FastReseg | `not cloned yet` | `main` | `7e0fe3650332` | license requires re-check |

## Brand independence note

Reference names in this file are provenance labels only. Local package names, CLI commands, figure labels, and manuscript novelty claims must use `ContamGuard-ST` terminology and the independent refinements in `README.md`, not upstream branding.
