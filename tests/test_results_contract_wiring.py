from __future__ import annotations

import hashlib
import json
import tempfile
import unittest
from pathlib import Path

from contamguard_st.contracts import ClaimGateEvidence
from contamguard_st.parity import Gate2RobustnessResult
from contamguard_st.real_smoke import RealSmokeConfig
from contamguard_st.result_wiring import emit_gate2_robustness_results


class ResultsContractWiringTests(unittest.TestCase):
    def test_vendored_contract_matches_recorded_sha256(self):
        package_dir = Path(__file__).resolve().parents[1] / "src" / "contamguard_st"
        contract_path = package_dir / "results_contract.py"
        expected = (package_dir / "results_contract.sha256").read_text(encoding="utf-8").strip().split()[0]
        observed = hashlib.sha256(contract_path.read_bytes()).hexdigest()
        self.assertEqual(observed, expected)

    def test_gate2_robustness_emits_contract_files(self):
        report = Gate2RobustnessResult(
            rows=[
                {
                    "plant_seed": 20260608,
                    "provenance": "RAN",
                    "pearson_gain_over_naive": 0.036804,
                    "spearman_gain_over_naive": 0.057327,
                    "mae_reduction_vs_naive": 0.04794,
                }
            ],
            aggregate={
                "verdict": "robust_real_differentiator",
                "seed_count": 20.0,
                "metric_summaries": {
                    "pearson_gain_over_naive": {
                        "mean": 0.036804,
                        "std": 0.010898,
                        "ci_low": 0.031703,
                        "ci_high": 0.041904,
                        "positive_seed_count": 20.0,
                    }
                },
            },
            reference_rows=[],
            evidence=ClaimGateEvidence(public_data_smoke=True, baseline_comparison=True, ablation=True, failure_modes=True),
        )
        config = RealSmokeConfig(adata_path=Path("data/processed/example/anndata.h5ad"))
        with tempfile.TemporaryDirectory() as tmp:
            paths = emit_gate2_robustness_results(report, config, results_dir=Path(tmp))
            metrics = json.loads(paths["metrics"].read_text(encoding="utf-8"))
            self.assertEqual(metrics["project"], "contamguard-st")
            self.assertEqual(metrics["metrics"]["gate2.robustness.seed_count"], 20.0)
            self.assertEqual(metrics["metrics"]["gate2.robustness.metric_summaries.pearson_gain_over_naive.mean"], 0.036804)
            self.assertEqual(metrics["metrics"]["gate2.robustness.seed_20260608.pearson_gain_over_naive"], 0.036804)
            self.assertTrue(paths["run_metadata"].is_file())


if __name__ == "__main__":
    unittest.main()
