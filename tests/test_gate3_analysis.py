from __future__ import annotations

import unittest
from pathlib import Path

from contamguard_st.contracts import ClaimGateEvidence
from contamguard_st.gate3 import Gate3Config, Gate3Result
from contamguard_st.real_smoke import RealSmokeConfig


class Gate3AnalysisUnitTests(unittest.TestCase):
    def test_config_builds_matching_gate2_parameters(self):
        config = Gate3Config(
            smoke=RealSmokeConfig(adata_path=Path("fixture.h5ad")),
            neighbor_count=4,
            seed=99,
            max_alpha=0.5,
            anchor_bins=8,
        )
        gate2 = config.to_gate2_config()
        self.assertEqual(gate2.neighbor_count, 4)
        self.assertEqual(gate2.plant_seed, 99)
        self.assertEqual(gate2.max_alpha, 0.5)
        self.assertEqual(gate2.anchor_bins, 8)

    def test_gate3_result_reports_license_as_only_missing_item(self):
        result = Gate3Result(
            ablation={"removed_component": "spatial_neighbor_lag"},
            failure_mode={"floor_statement": "low-count floor"},
            evidence=ClaimGateEvidence(public_data_smoke=True, baseline_comparison=True, ablation=True, failure_modes=True),
        )
        self.assertEqual(result.claim_status.value, "locked")
        self.assertEqual(result.to_jsonable()["missing_claim_evidence"], ["license_review"])


if __name__ == "__main__":
    unittest.main()
