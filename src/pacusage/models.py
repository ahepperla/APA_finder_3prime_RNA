"""Shared typed records used by the scientific modules."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any


@dataclass(frozen=True)
class Sample:
    sample_id: str
    alignment: str
    condition: str
    control: str
    control_condition: str
    replicate: str = ""
    batch: str = ""
    donor: str = ""
    layout: str = "auto"
    strandedness: str = "auto"
    library_profile: str = "generic_3prime"
    evidence_source: str = "auto"
    covariates: dict[str, str] | None = None

    def as_dict(self) -> dict[str, Any]:
        row = asdict(self)
        covariates = row.pop("covariates") or {}
        row.update(covariates)
        return row


@dataclass(frozen=True)
class EvidenceObservation:
    sample_id: str
    contig: str
    strand: str
    coordinate: int
    count: int
    poly_a_clip_count: int = 0
    evidence_source: str = "read_3p"


@dataclass(frozen=True)
class PacCandidate:
    contig: str
    strand: str
    coordinate: int
    total_count: int
    supporting_samples: int
    capped_support: int
    member_coordinates: tuple[int, ...]
    status: str = "primary"
    rejection_reason: str = ""
    fraction_within_2nt: float = 0.0
    width_90: int = 0
    poly_a_clip_fraction: float = 0.0
    local_enrichment: float = 0.0
    coordinate_interval_low: int | None = None
    coordinate_interval_high: int | None = None
    coordinate_bootstrap_successes: int = 0
