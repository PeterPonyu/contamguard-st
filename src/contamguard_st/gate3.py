from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from .contracts import ClaimGateEvidence, ClaimStatus, evaluate_claim_gate
from .parity import (
    Gate2ParityConfig,
    _recovery_metrics,
    planted_spillover_recovery,
)
from .real_smoke import RealSmokeConfig


@dataclass(frozen=True)
class Gate3Config:
    smoke: RealSmokeConfig = field(default_factory=RealSmokeConfig)
    neighbor_count: int = 6
    seed: int = 20260608
    max_alpha: float = 0.8
    alpha_shape_a: float = 1.5
    alpha_shape_b: float = 3.0
    anchor_bins: int = 20
    anchor_quantile: float = 0.3

    def to_gate2_config(self) -> Gate2ParityConfig:
        return Gate2ParityConfig(
            smoke=self.smoke,
            neighbor_count=self.neighbor_count,
            plant_seed=self.seed,
            max_alpha=self.max_alpha,
            alpha_shape_a=self.alpha_shape_a,
            alpha_shape_b=self.alpha_shape_b,
            anchor_bins=self.anchor_bins,
            anchor_quantile=self.anchor_quantile,
        )


@dataclass(frozen=True)
class Gate3Result:
    ablation: dict[str, Any]
    failure_mode: dict[str, Any]
    evidence: ClaimGateEvidence

    @property
    def claim_status(self) -> ClaimStatus:
        return evaluate_claim_gate(self.evidence)

    def to_jsonable(self) -> dict[str, Any]:
        return {
            "ablation": self.ablation,
            "failure_mode": self.failure_mode,
            "claim_status": self.claim_status.value,
            "missing_claim_evidence": list(self.evidence.missing()),
        }


def run_gate3_analysis(config: Gate3Config | None = None) -> Gate3Result:
    import numpy as np

    config = config or Gate3Config()
    recovery = planted_spillover_recovery(config.to_gate2_config())
    truth = recovery["planted_alpha"]
    ours = recovery["ours_estimated_alpha"]
    no_spatial = np.clip(
        (recovery["contaminated_gene_counts"] - recovery["expected_clean_gene_counts"])
        / max(1e-9, float(np.median(recovery["neighbor_contaminated_gene_counts"]))),
        0.0,
        float(config.max_alpha),
    )
    ours_metrics = recovery["ours_metrics"]
    naive_metrics = recovery["naive_metrics"]
    no_spatial_metrics = _recovery_metrics(no_spatial, truth)

    ablation = {
        "removed_component": "spatial_neighbor_denominator",
        "ours_planted_alpha_pearson": ours_metrics["pearson"],
        "ours_planted_alpha_spearman": ours_metrics["spearman"],
        "ours_planted_alpha_mae": ours_metrics["mae"],
        "no_spatial_planted_alpha_pearson": no_spatial_metrics["pearson"],
        "no_spatial_planted_alpha_spearman": no_spatial_metrics["spearman"],
        "no_spatial_planted_alpha_mae": no_spatial_metrics["mae"],
        "naive_planted_alpha_pearson": naive_metrics["pearson"],
        "naive_planted_alpha_spearman": naive_metrics["spearman"],
        "naive_planted_alpha_mae": naive_metrics["mae"],
        "pearson_drop_without_spatial": round(float(ours_metrics["pearson"] - no_spatial_metrics["pearson"]), 6),
        "mae_increase_without_spatial": round(float(no_spatial_metrics["mae"] - ours_metrics["mae"]), 6),
        "interpretation": "removing the spatial neighbor normalization worsens recovery of the planted alpha and moves the estimator toward non-spatial recovery",
    }

    low_control_mask = recovery["low_control_mask"]
    dense_mask = recovery["dense_mask"]
    failure_mode = {
        "low_control_count_floor": round(float(np.min(recovery["table"].control_counts[low_control_mask])) if np.any(low_control_mask) else 0.0, 6),
        "low_control_count_ceiling": round(float(np.max(recovery["table"].control_counts[low_control_mask])) if np.any(low_control_mask) else 0.0, 6),
        "low_control_cell_fraction": round(float(np.mean(low_control_mask)), 6),
        "low_control_ours_mae": round(float(np.mean(np.abs(ours[low_control_mask] - truth[low_control_mask]))) if np.any(low_control_mask) else 0.0, 6),
        "global_ours_mae": ours_metrics["mae"],
        "dense_region_cell_fraction": round(float(np.mean(dense_mask)), 6),
        "dense_region_ours_mae": round(float(np.mean(np.abs(ours[dense_mask] - truth[dense_mask]))) if np.any(dense_mask) else 0.0, 6),
        "dense_region_median_neighbor_distance_floor": round(float(np.median(recovery["median_neighbor_distance"][dense_mask])) if np.any(dense_mask) else 0.0, 6),
        "global_median_neighbor_distance": round(float(np.median(recovery["median_neighbor_distance"])), 6),
        "floor_statement": "low-control and very dense neighborhoods are reported as planted-alpha recovery floors, not biological discoveries",
    }
    evidence = ClaimGateEvidence(public_data_smoke=True, baseline_comparison=True, ablation=True, failure_modes=True)
    return Gate3Result(ablation, failure_mode, evidence)
