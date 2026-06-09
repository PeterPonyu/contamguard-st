# ContamGuard-ST

ContamGuard-ST scores control-channel contamination and tests spatial spillover recovery against planted-alpha benchmarks.

This repository is a conservative public code surface: method implementation, command-line entry points, tests, and the byte-locked results schema. Background citations are listed in `BASELINE_REFERENCES.md`.

## Install

```bash
python -m pip install -e .
```

The lightweight unit tests run without bundled datasets. Real-data commands expect local spatial-omics inputs and expose `--help` for path overrides.

## Command-line usage

```bash
python -m contamguard_st.cli smoke-synthetic
python -m contamguard_st.cli smoke-real --help
python -m contamguard_st.cli gate2-parity --help
python -m contamguard_st.cli gate2-robustness --help
python -m contamguard_st.cli gate3-analysis --help
python -m contamguard_st.cli claim-status
```

Commands emit JSON to stdout. Gate commands also write uniform contract outputs under `results/<project>/` via the vendored `results_contract.py`.

## Evidence-derived status

`python -m contamguard_st.cli claim-status` reads `evidence/summary.json`, derives the public claim label from the visible evidence and stated validation bar, and currently prints `preliminary` plus the missing-evidence list. It does not rely on private governance documents or a hardcoded validation constant.

The current public-safe scope is preliminary: the multi-seed planted-spillover effect is robust but synthetic/modest, and the external SPLIT-style same-card baseline has not been run.

## Citations and references

See `BASELINE_REFERENCES.md` for papers, code references, and citation context.
