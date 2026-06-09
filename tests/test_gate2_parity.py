from __future__ import annotations

import unittest
from pathlib import Path

import numpy as np

from contamguard_st.contracts import ClaimGateEvidence
from contamguard_st.parity import (
    Gate2ParityConfig,
    Gate2ParityResult,
    Gate2RobustnessConfig,
    Gate2RobustnessResult,
    _planted_alpha,
    _summary_stats,
    spatial_control_score,
    spatial_neighbor_control_fraction,
)
from contamguard_st.real_smoke import RealSmokeConfig


class Gate2ParityUnitTests(unittest.TestCase):
    def test_spatial_neighbor_fraction_uses_nearby_cells(self):
        control_fraction = np.array([0.0, 1.0, 0.0, 0.0])
        xy = np.array([[0.0, 0.0], [0.2, 0.0], [10.0, 10.0], [11.0, 10.0]])
        neighbor_fraction, distance = spatial_neighbor_control_fraction(control_fraction, xy, neighbor_count=2)
        self.assertGreater(neighbor_fraction[0], neighbor_fraction[3])
        self.assertTrue(np.all(np.isfinite(distance)))

    def test_spatial_score_blends_local_and_neighbor_anchor(self):
        score = spatial_control_score(np.array([0.2]), np.array([0.8]), spatial_weight=0.5)
        self.assertAlmostEqual(float(score[0]), 0.5)

    def test_planted_alpha_is_seeded_and_bounded(self):
        config = Gate2ParityConfig(
            smoke=RealSmokeConfig(adata_path=Path("fixture.h5ad")),
            plant_seed=7,
            max_alpha=0.8,
        )
        alpha = _planted_alpha(config, 32)
        self.assertEqual(len(alpha), 32)
        self.assertGreater(float(np.std(alpha)), 0.0)
        self.assertGreaterEqual(float(np.min(alpha)), 0.0)
        self.assertLessEqual(float(np.max(alpha)), 0.8)

    def test_full_gate_evidence_stays_locked_without_license(self):
        result = Gate2ParityResult(
            rows=[{"method_id": "local", "provenance": "RAN"}],
            differentiator={"verdict": "real_differentiator"},
            evidence=ClaimGateEvidence(public_data_smoke=True, baseline_comparison=True, ablation=True, failure_modes=True),
        )
        payload = result.to_jsonable()
        self.assertEqual(payload["claim_status"], "locked")
        self.assertEqual(payload["missing_claim_evidence"], ["license_review"])

    def test_robustness_config_requires_at_least_fifteen_seeds(self):
        with self.assertRaises(ValueError):
            Gate2RobustnessConfig(
                smoke=RealSmokeConfig(adata_path=Path("fixture.h5ad")),
                seed_count=14,
            ).validate()
        self.assertEqual(
            len(
                Gate2RobustnessConfig(
                    smoke=RealSmokeConfig(adata_path=Path("fixture.h5ad")),
                    seed_count=15,
                    seed_start=10,
                ).seeds()
            ),
            15,
        )

    def test_summary_stats_reports_positive_seed_fraction_and_ci(self):
        summary = _summary_stats([0.1, 0.2, 0.3, -0.1], ci_level=0.95)
        self.assertEqual(summary["positive_seed_count"], 3.0)
        self.assertEqual(summary["positive_seed_fraction"], 0.75)
        self.assertLess(summary["ci_low"], summary["ci_high"])

    def test_robustness_result_stays_locked_without_license(self):
        result = Gate2RobustnessResult(
            rows=[{"plant_seed": 1, "provenance": "RAN"}],
            aggregate={"verdict": "robust_real_differentiator"},
            reference_rows=[],
            evidence=ClaimGateEvidence(public_data_smoke=True, baseline_comparison=True, ablation=True, failure_modes=True),
        )
        payload = result.to_jsonable()
        self.assertEqual(payload["claim_status"], "locked")
        self.assertEqual(payload["missing_claim_evidence"], ["license_review"])


if __name__ == "__main__":
    unittest.main()
