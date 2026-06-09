import unittest

from contamguard_st.smoke import build_synthetic_molecules, estimate_control_spillover, run_synthetic_smoke


class SmokeTests(unittest.TestCase):
    def test_spillover_smoke_flags_control_heavy_cells(self):
        report = run_synthetic_smoke()
        self.assertEqual(report.claim_status.value, "locked")
        self.assertEqual(report.metrics["control_channels_retained"], 1.0)
        self.assertGreaterEqual(report.metrics["flagged_cells"], 2.0)

    def test_dropping_controls_fails_loudly(self):
        records = [r for r in build_synthetic_molecules() if r.channel_type == "gene"]
        with self.assertRaises(ValueError):
            estimate_control_spillover(records)

if __name__ == "__main__":
    unittest.main()
