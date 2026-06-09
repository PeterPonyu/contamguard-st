from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

from .contracts import ClaimGateEvidence
from .gate3 import Gate3Result
from .parity import Gate2ParityResult, Gate2RobustnessResult
from .real_smoke import RealSmokeConfig, RealSmokeResult
from .results_contract import dataset_card_id, write_results

PROJECT_ID = "contamguard-st"


def _numeric_items(payload: Mapping[str, Any], prefix: str = "") -> dict[str, float | None]:
    metrics: dict[str, float | None] = {}
    for key, value in payload.items():
        name = f"{prefix}.{key}" if prefix else str(key)
        if value is None:
            metrics[name] = None
        elif isinstance(value, bool):
            continue
        elif isinstance(value, (int, float)):
            metrics[name] = float(value)
        elif isinstance(value, Mapping):
            metrics.update(_numeric_items(value, name))
    return metrics


def _method_id(value: Any) -> str:
    return str(value or "method").replace(".", "_").replace(" ", "_").replace("/", "_")


def _config_paths(config: RealSmokeConfig) -> list[str]:
    return [str(config.adata_path)]


def _metadata(config: RealSmokeConfig, report: RealSmokeResult | None, *, notes: str) -> dict[str, Any]:
    metrics = report.metrics if report is not None else {}
    return {
        "dataset_paths": _config_paths(config),
        "n_obs": metrics.get("sampled_cell_count"),
        "n_vars": None,
        "seed": None,
        "deterministic": True,
        "num_threads": 1,
        "reproducibility_level": "seeded",
        "normalization": {"applied": False, "method": "raw control/gene counts"},
        "interpretability": {
            "claim_scope": "robust_modest_planted_spillover_recovery_not_real_contamination_gt",
            "control_anchor": True,
            "spatial_neighbor_estimator": True,
        },
        "notes": notes,
        "provenance": {
            "max_cells": config.max_cells,
            "flag_floor": config.flag_floor,
        },
    }


def emit_real_smoke_results(
    report: RealSmokeResult,
    config: RealSmokeConfig,
    *,
    results_dir: Path | None = None,
    outputs: Mapping[str, Any] | None = None,
) -> dict[str, Path]:
    return write_results(
        PROJECT_ID,
        dataset_card_id(_config_paths(config)),
        report.metrics,
        outputs=outputs,
        run_metadata=_metadata(config, report, notes="real-data control-channel smoke emitted through the vendored results contract"),
        results_dir=results_dir,
    )


def emit_gate2_results(
    report: Gate2ParityResult,
    config: RealSmokeConfig,
    *,
    results_dir: Path | None = None,
    outputs: Mapping[str, Any] | None = None,
) -> dict[str, Path]:
    metrics: dict[str, float | None] = {}
    metrics.update(_numeric_items(report.differentiator, "gate2"))
    for row in report.rows:
        method = _method_id(row.get("method_id"))
        metrics.update(_numeric_items(row, f"gate2.{method}"))
    synthetic_report = RealSmokeResult(dict(metrics), evidence=ClaimGateEvidence(public_data_smoke=True))
    return write_results(
        PROJECT_ID,
        dataset_card_id(_config_paths(config)),
        metrics,
        outputs=outputs,
        run_metadata=_metadata(
            config,
            synthetic_report,
            notes="single-seed planted-spillover parity emitted through the vendored results contract",
        ),
        results_dir=results_dir,
    )


def emit_gate2_robustness_results(
    report: Gate2RobustnessResult,
    config: RealSmokeConfig,
    *,
    results_dir: Path | None = None,
    outputs: Mapping[str, Any] | None = None,
) -> dict[str, Path]:
    metrics: dict[str, float | None] = {}
    metrics.update(_numeric_items(report.aggregate, "gate2.robustness"))
    for row in report.rows:
        seed_label = f"seed_{int(row.get('plant_seed', 0))}"
        metrics.update(_numeric_items(row, f"gate2.robustness.{seed_label}"))
    synthetic_report = RealSmokeResult(dict(metrics), evidence=ClaimGateEvidence(public_data_smoke=True))
    return write_results(
        PROJECT_ID,
        dataset_card_id(_config_paths(config)),
        metrics,
        outputs=outputs,
        run_metadata=_metadata(
            config,
            synthetic_report,
            notes="multi-seed planted-spillover robustness emitted through the vendored results contract",
        ),
        results_dir=results_dir,
    )


def emit_gate3_results(
    report: Gate3Result,
    config: RealSmokeConfig,
    *,
    results_dir: Path | None = None,
    outputs: Mapping[str, Any] | None = None,
) -> dict[str, Path]:
    metrics: dict[str, float | None] = {}
    metrics.update(_numeric_items(report.ablation, "gate3.ablation"))
    metrics.update(_numeric_items(report.failure_mode, "gate3.failure_mode"))
    synthetic_report = RealSmokeResult(dict(metrics), evidence=ClaimGateEvidence(public_data_smoke=True))
    return write_results(
        PROJECT_ID,
        dataset_card_id(_config_paths(config)),
        metrics,
        outputs=outputs,
        run_metadata=_metadata(
            config,
            synthetic_report,
            notes="planted-spillover ablation and failure-mode floor emitted through the vendored results contract",
        ),
        results_dir=results_dir,
    )
