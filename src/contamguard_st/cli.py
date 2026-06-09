from __future__ import annotations

import argparse
import json
from pathlib import Path

from .contracts import ClaimGateEvidence, evaluate_claim_gate
from .claim_status import graduation_claim_status
from .real_smoke import default_adata_path
from .smoke import run_synthetic_smoke


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="contamguard-st")
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("smoke-synthetic")
    real = sub.add_parser("smoke-real")
    real.add_argument("--adata-path", "--st-path", dest="adata_path", type=Path, default=None)
    real.add_argument("--max-cells", type=int, default=5000)
    real.add_argument("--flag-floor", type=float, default=0.01)
    real.add_argument("--results-dir", type=Path, default=None)
    gate2 = sub.add_parser("gate2-parity")
    gate2.add_argument("--adata-path", "--st-path", dest="adata_path", type=Path, default=None)
    gate2.add_argument("--max-cells", type=int, default=5000)
    gate2.add_argument("--flag-floor", type=float, default=0.01)
    gate2.add_argument("--neighbor-count", type=int, default=6)
    gate2.add_argument("--plant-seed", type=int, default=20260608)
    gate2.add_argument("--max-alpha", type=float, default=0.8)
    gate2.add_argument("--alpha-shape-a", type=float, default=1.5)
    gate2.add_argument("--alpha-shape-b", type=float, default=3.0)
    gate2.add_argument("--anchor-bins", type=int, default=20)
    gate2.add_argument("--anchor-quantile", type=float, default=0.3)
    gate2.add_argument("--reference-path", type=Path, default=None)
    gate2.add_argument("--out-path", type=Path, default=None)
    gate2.add_argument("--results-dir", type=Path, default=None)
    robust = sub.add_parser("gate2-robustness")
    robust.add_argument("--adata-path", "--st-path", dest="adata_path", type=Path, default=None)
    robust.add_argument("--max-cells", type=int, default=5000)
    robust.add_argument("--flag-floor", type=float, default=0.01)
    robust.add_argument("--neighbor-count", type=int, default=6)
    robust.add_argument("--seed-start", type=int, default=20260608)
    robust.add_argument("--seed-count", type=int, default=20)
    robust.add_argument("--seed-step", type=int, default=1)
    robust.add_argument("--max-alpha", type=float, default=0.8)
    robust.add_argument("--alpha-shape-a", type=float, default=1.5)
    robust.add_argument("--alpha-shape-b", type=float, default=3.0)
    robust.add_argument("--anchor-bins", type=int, default=20)
    robust.add_argument("--anchor-quantile", type=float, default=0.3)
    robust.add_argument("--ci-level", type=float, default=0.95)
    robust.add_argument("--consistency-threshold", type=float, default=0.8)
    robust.add_argument("--reference-path", type=Path, default=None)
    robust.add_argument("--out-path", type=Path, default=None)
    robust.add_argument("--results-dir", type=Path, default=None)
    gate3 = sub.add_parser("gate3-analysis")
    gate3.add_argument("--adata-path", "--st-path", dest="adata_path", type=Path, default=None)
    gate3.add_argument("--max-cells", type=int, default=5000)
    gate3.add_argument("--flag-floor", type=float, default=0.01)
    gate3.add_argument("--neighbor-count", type=int, default=6)
    gate3.add_argument("--plant-seed", type=int, default=20260608)
    gate3.add_argument("--max-alpha", type=float, default=0.8)
    gate3.add_argument("--alpha-shape-a", type=float, default=1.5)
    gate3.add_argument("--alpha-shape-b", type=float, default=3.0)
    gate3.add_argument("--anchor-bins", type=int, default=20)
    gate3.add_argument("--anchor-quantile", type=float, default=0.3)
    gate3.add_argument("--out-path", type=Path, default=None)
    gate3.add_argument("--results-dir", type=Path, default=None)
    sub.add_parser("claim-status")
    args = parser.parse_args(argv)
    if args.command == "smoke-synthetic":
        report = run_synthetic_smoke()
        print(json.dumps({"metrics": report.metrics, "claim_status": report.claim_status.value}, sort_keys=True))
        return 0
    if args.command == "smoke-real":
        from .real_smoke import RealSmokeConfig, run_real_data_smoke
        from .result_wiring import emit_real_smoke_results

        config = RealSmokeConfig(
            adata_path=args.adata_path or default_adata_path(),
            max_cells=args.max_cells,
            flag_floor=args.flag_floor,
        )
        report = run_real_data_smoke(config)
        payload = report.to_jsonable()
        payload["contract_results"] = _paths_payload(emit_real_smoke_results(report, config, results_dir=args.results_dir))
        print(json.dumps(payload, sort_keys=True))
        return 0
    if args.command == "gate2-parity":
        from .parity import Gate2ParityConfig, run_gate2_parity
        from .result_wiring import emit_gate2_results

        config = Gate2ParityConfig(
            smoke=_real_config_from_args(args),
            reference_path=args.reference_path,
            neighbor_count=args.neighbor_count,
            plant_seed=args.plant_seed,
            max_alpha=args.max_alpha,
            alpha_shape_a=args.alpha_shape_a,
            alpha_shape_b=args.alpha_shape_b,
            anchor_bins=args.anchor_bins,
            anchor_quantile=args.anchor_quantile,
        )
        report = run_gate2_parity(config)
        payload = report.to_jsonable()
        outputs = {}
        if args.out_path is not None:
            outputs["parity_table"] = args.out_path
        payload["contract_results"] = _paths_payload(
            emit_gate2_results(report, config.smoke, results_dir=args.results_dir, outputs=outputs)
        )
        _write_json_if_requested(args.out_path, payload)
        print(json.dumps(payload, sort_keys=True))
        return 0
    if args.command == "gate2-robustness":
        from .parity import Gate2RobustnessConfig, run_gate2_robustness
        from .result_wiring import emit_gate2_robustness_results

        config = Gate2RobustnessConfig(
            smoke=_real_config_from_args(args),
            reference_path=args.reference_path,
            neighbor_count=args.neighbor_count,
            seed_start=args.seed_start,
            seed_count=args.seed_count,
            seed_step=args.seed_step,
            max_alpha=args.max_alpha,
            alpha_shape_a=args.alpha_shape_a,
            alpha_shape_b=args.alpha_shape_b,
            anchor_bins=args.anchor_bins,
            anchor_quantile=args.anchor_quantile,
            ci_level=args.ci_level,
            consistency_threshold=args.consistency_threshold,
        )
        report = run_gate2_robustness(config)
        payload = report.to_jsonable()
        outputs = {}
        if args.out_path is not None:
            outputs["robustness_table"] = args.out_path
        payload["contract_results"] = _paths_payload(
            emit_gate2_robustness_results(report, config.smoke, results_dir=args.results_dir, outputs=outputs)
        )
        _write_json_if_requested(args.out_path, payload)
        print(json.dumps(payload, sort_keys=True))
        return 0
    if args.command == "gate3-analysis":
        from .gate3 import Gate3Config, run_gate3_analysis
        from .result_wiring import emit_gate3_results

        config = Gate3Config(
            smoke=_real_config_from_args(args),
            neighbor_count=args.neighbor_count,
            seed=args.plant_seed,
            max_alpha=args.max_alpha,
            alpha_shape_a=args.alpha_shape_a,
            alpha_shape_b=args.alpha_shape_b,
            anchor_bins=args.anchor_bins,
            anchor_quantile=args.anchor_quantile,
        )
        report = run_gate3_analysis(config)
        payload = report.to_jsonable()
        outputs = {}
        if args.out_path is not None:
            outputs["gate3_table"] = args.out_path
        payload["contract_results"] = _paths_payload(emit_gate3_results(report, config.smoke, results_dir=args.results_dir, outputs=outputs))
        _write_json_if_requested(args.out_path, payload)
        print(json.dumps(payload, sort_keys=True))
        return 0
    if args.command == "claim-status":
        print(graduation_claim_status().value)
        return 0
    print(evaluate_claim_gate(ClaimGateEvidence()).value)
    return 0


def _real_config_from_args(args: argparse.Namespace):
    from .real_smoke import RealSmokeConfig

    return RealSmokeConfig(
        adata_path=args.adata_path or default_adata_path(),
        max_cells=args.max_cells,
        flag_floor=args.flag_floor,
    )


def _write_json_if_requested(path: Path | None, payload: dict[str, object]) -> None:
    if path is None:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _paths_payload(paths: dict[str, Path]) -> dict[str, str]:
    return {key: str(path) for key, path in paths.items()}


if __name__ == "__main__":
    raise SystemExit(main())
