import tempfile
import unittest
from pathlib import Path

from contamguard_st.contracts import ClaimGateEvidence
from contamguard_st.data_paths import find_repo_root, processed_data_path
from contamguard_st.real_smoke import RealSmokeConfig, RealSmokeResult, _robust_flag_threshold


class RealDataSmokeUnitTests(unittest.TestCase):
    def test_result_keeps_claim_locked(self):
        result = RealSmokeResult({"control_vs_gene_log10_mean_gap": 1.0}, ClaimGateEvidence(public_data_smoke=True))
        self.assertEqual(result.claim_status.value, "locked")
        self.assertEqual(
            result.to_jsonable()["missing_claim_evidence"],
            ["baseline_comparison", "ablation", "failure_modes", "license_review"],
        )

    def test_config_rejects_tiny_sample(self):
        with self.assertRaises(ValueError):
            RealSmokeConfig(adata_path=Path("fixture.h5ad"), max_cells=4).validate()

    def test_repo_root_path_resolution(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "data" / "processed").mkdir(parents=True)
            anchor = root / "src" / "package" / "module.py"
            anchor.parent.mkdir(parents=True)
            anchor.touch()
            self.assertEqual(find_repo_root(anchor), root)
            self.assertEqual(
                processed_data_path("fixture_card", anchor=anchor),
                root / "data" / "processed" / "fixture_card",
            )

    def test_threshold_has_floor(self):
        import numpy as np

        self.assertGreaterEqual(_robust_flag_threshold(np.array([0.0, 0.001, 0.002]), 0.01), 0.01)


if __name__ == "__main__":
    unittest.main()
