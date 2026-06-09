from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .contracts import ClaimGateEvidence, ClaimStatus, evaluate_claim_gate
from .data_paths import processed_data_path


def default_adata_path() -> Path:
    card = "contam" + "guard_xenium_cervix_controls"
    return processed_data_path(card, "anndata.h5ad")


@dataclass(frozen=True)
class RealSmokeConfig:
    adata_path: Path = field(default_factory=default_adata_path)
    max_cells: int = 5000
    flag_floor: float = 0.01

    def validate(self) -> None:
        if self.max_cells < 16:
            raise ValueError("max_cells must be at least 16")
        if self.flag_floor < 0.0:
            raise ValueError("flag_floor must be non-negative")


@dataclass(frozen=True)
class RealSmokeResult:
    metrics: dict[str, float]
    evidence: ClaimGateEvidence

    @property
    def claim_status(self) -> ClaimStatus:
        return evaluate_claim_gate(self.evidence)

    def to_jsonable(self) -> dict[str, object]:
        return {
            "metrics": self.metrics,
            "claim_status": self.claim_status.value,
            "missing_claim_evidence": list(self.evidence.missing()),
        }


def _resolve_h5ad(path: str | Path) -> Path:
    resolved = Path(path).expanduser()
    if not resolved.exists():
        raise FileNotFoundError(f"AnnData file does not exist: {resolved}")
    if resolved.suffix != ".h5ad":
        raise ValueError(f"expected a .h5ad file, got: {resolved}")
    return resolved


def _read_h5ad(path: str | Path) -> Any:
    try:
        import anndata as ad  # type: ignore[import-untyped]
    except ImportError as exc:  # pragma: no cover - depends on optional real-data env
        raise RuntimeError("real-data smoke requires the optional anndata package") from exc
    return ad.read_h5ad(_resolve_h5ad(path), backed="r")


def _close(adata: Any) -> None:
    close = getattr(getattr(adata, "file", None), "close", None)
    if close is not None:
        close()


def _sample_indices(total: int, limit: int) -> Any:
    import numpy as np

    count = min(total, limit)
    if count == total:
        return np.arange(total, dtype=int)
    return np.unique(np.linspace(0, total - 1, count, dtype=int))


def _as_1d(values: Any) -> Any:
    import numpy as np

    return np.asarray(values, dtype=np.float64).reshape(-1)


def _axis_sum(matrix: Any, axis: int) -> Any:
    if hasattr(matrix, "sum"):
        return _as_1d(matrix.sum(axis=axis))
    import numpy as np

    return _as_1d(np.asarray(matrix).sum(axis=axis))


def _safe_fraction(numer: Any, denom: Any) -> Any:
    import numpy as np

    numer = np.asarray(numer, dtype=np.float64)
    denom = np.asarray(denom, dtype=np.float64)
    return np.divide(numer, denom, out=np.zeros_like(numer, dtype=np.float64), where=denom > 0.0)


def _robust_flag_threshold(values: Any, floor: float) -> float:
    import numpy as np

    median = float(np.median(values))
    mad = float(np.median(np.abs(values - median)))
    return max(floor, median + 3.0 * mad)


def run_real_data_smoke(config: RealSmokeConfig) -> RealSmokeResult:
    import numpy as np

    config.validate()
    adata = _read_h5ad(config.adata_path)
    try:
        if "is_control" not in adata.var:
            raise KeyError("real-data smoke requires var['is_control']")
        control_mask = np.asarray(adata.var["is_control"].astype(bool).to_numpy())
        control_idx = np.flatnonzero(control_mask)
        gene_idx = np.flatnonzero(~control_mask)
        if len(control_idx) == 0 or len(gene_idx) == 0:
            raise ValueError("real-data smoke requires both control and gene-expression features")

        rows = _sample_indices(int(adata.n_obs), config.max_cells)
        control_matrix = adata[rows, control_idx].X
        gene_matrix = adata[rows, gene_idx].X
        control_cell_counts = _axis_sum(control_matrix, axis=1)
        gene_cell_counts = _axis_sum(gene_matrix, axis=1)
        total_cell_counts = control_cell_counts + gene_cell_counts
        cell_fraction = _safe_fraction(control_cell_counts, total_cell_counts)
        threshold = _robust_flag_threshold(cell_fraction, config.flag_floor)
        control_feature_means = _axis_sum(control_matrix, axis=0) / float(len(rows))
        gene_feature_means = _axis_sum(gene_matrix, axis=0) / float(len(rows))
        eps = 1e-9
        mean_gap = float(np.log10(gene_feature_means.mean() + eps) - np.log10(control_feature_means.mean() + eps))
        median_gap = float(np.log10(np.median(gene_feature_means) + eps) - np.log10(np.median(control_feature_means) + eps))
        flagged = cell_fraction >= threshold
        metrics = {
            "cell_count_total": float(adata.n_obs),
            "sampled_cell_count": float(len(rows)),
            "gene_feature_count": float(len(gene_idx)),
            "control_feature_count": float(len(control_idx)),
            "median_cell_control_fraction": round(float(np.median(cell_fraction)), 6),
            "p95_cell_control_fraction": round(float(np.quantile(cell_fraction, 0.95)), 6),
            "flag_threshold_control_fraction": round(float(threshold), 6),
            "flagged_cell_fraction": round(float(np.mean(flagged)), 6),
            "control_vs_gene_log10_mean_gap": round(mean_gap, 6),
            "control_vs_gene_log10_median_gap": round(median_gap, 6),
            "mean_gene_count_per_sampled_cell": round(float(np.mean(gene_cell_counts)), 6),
            "mean_control_count_per_sampled_cell": round(float(np.mean(control_cell_counts)), 6),
        }
        return RealSmokeResult(metrics, ClaimGateEvidence(public_data_smoke=True))
    finally:
        _close(adata)
