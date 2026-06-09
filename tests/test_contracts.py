import unittest

from contamguard_st import ClaimGateEvidence, ClaimStatus, ControlChannelSpec, evaluate_claim_gate


class ContractTests(unittest.TestCase):
    def test_claim_gate_requires_all_evidence(self):
        self.assertEqual(evaluate_claim_gate(ClaimGateEvidence(public_data_smoke=True)), ClaimStatus.LOCKED)

    def test_human_signed_full_evidence_validates_claim(self):
        evidence = ClaimGateEvidence(
            public_data_smoke=True,
            baseline_comparison=True,
            ablation=True,
            failure_modes=True,
            license_review=True,
        )
        self.assertEqual(evaluate_claim_gate(evidence, human_signed=True), ClaimStatus.VALIDATED)

    def test_control_channel_spec_keeps_anchors(self):
        spec = ControlChannelSpec()
        self.assertTrue(spec.is_control("negative_control"))
        self.assertFalse(spec.is_control("gene"))

if __name__ == "__main__":
    unittest.main()
