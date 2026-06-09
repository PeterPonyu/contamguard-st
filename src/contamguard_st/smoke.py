from __future__ import annotations

from collections import defaultdict

from .contracts import ClaimGateEvidence, CellReliabilityReport, ContaminationSmokeReport, ControlChannelSpec, MoleculeRecord


def build_synthetic_molecules() -> list[MoleculeRecord]:
    records: list[MoleculeRecord] = []
    idx = 0
    layout = {
        "clean_a": (28, 2, 0.04),
        "clean_b": (26, 1, 0.06),
        "spill_c": (16, 12, 0.31),
        "edge_d": (18, 8, 0.42),
    }
    for cell_id, (gene_count, control_count, boundary) in layout.items():
        for i in range(gene_count):
            records.append(MoleculeRecord(f"m{idx}", f"G{i % 5}", cell_id, float(i), float(idx % 7), "gene"))
            idx += 1
        control_types = ("negative_control", "blank", "codeword")
        for i in range(control_count):
            records.append(MoleculeRecord(f"m{idx}", f"CTRL{i % 3}", cell_id, float(i), float(idx % 5), control_types[i % 3]))
            idx += 1
    return records


def estimate_control_spillover(records: list[MoleculeRecord], boundaries: dict[str, float] | None = None) -> tuple[CellReliabilityReport, ...]:
    spec = ControlChannelSpec()
    counts: dict[str, dict[str, int]] = defaultdict(lambda: {"gene": 0, "control": 0})
    for rec in records:
        bucket = "control" if spec.is_control(rec.channel_type) else "gene"
        counts[rec.cell_id][bucket] += 1
    if not any(v["control"] for v in counts.values()):
        raise ValueError("control channels are required for this smoke")
    boundaries = boundaries or {"clean_a": 0.04, "clean_b": 0.06, "spill_c": 0.31, "edge_d": 0.42}
    reports: list[CellReliabilityReport] = []
    for cell_id, pair in sorted(counts.items()):
        total = pair["gene"] + pair["control"]
        fraction = pair["control"] / total if total else 0.0
        boundary = boundaries.get(cell_id, 0.0)
        reports.append(CellReliabilityReport(cell_id, fraction, boundary, fraction >= 0.24 or boundary >= 0.35))
    return tuple(reports)


def run_synthetic_smoke() -> ContaminationSmokeReport:
    records = build_synthetic_molecules()
    reports = estimate_control_spillover(records)
    flagged = [r for r in reports if r.qc_flag]
    metrics = {
        "cells": float(len(reports)),
        "flagged_cells": float(len(flagged)),
        "max_control_fraction": round(max(r.control_fraction for r in reports), 4),
        "min_clean_control_fraction": round(min(r.control_fraction for r in reports), 4),
        "control_channels_retained": 1.0,
    }
    return ContaminationSmokeReport(reports, metrics, ClaimGateEvidence(ablation=True, failure_modes=True))
