"""Protocol strategies kept separate from PAC discovery and quantification."""

from __future__ import annotations

from dataclasses import dataclass

from .errors import PacusageError


@dataclass(frozen=True)
class ProtocolProfile:
    name: str
    default_layout: str
    default_strandedness: str
    default_evidence_source: str
    default_endpoint_model: str
    read_selector: str = "read1"
    expected_terminal_span: int | None = None

    def candidate_sources(self, layout: str) -> tuple[str, ...]:
        if self.name == "plasmidsaurus_3prime":
            return ("read_3p",)
        if self.name == "exact_boundary":
            return (
                ("polyA_junction", "read_3p")
                if layout == "SE"
                else (
                    "polyA_junction",
                    "fragment_3p",
                    "read_3p",
                )
            )
        if layout == "PE":
            return ("read_3p", "read_5p", "fragment_3p", "polyA_junction")
        return ("read_3p", "read_5p", "polyA_junction")


PROFILES = {
    "plasmidsaurus_3prime": ProtocolProfile(
        name="plasmidsaurus_3prime",
        default_layout="SE",
        default_strandedness="forward",
        default_evidence_source="read_3p",
        default_endpoint_model="proximal_tag",
        expected_terminal_span=400,
    ),
    "exact_boundary": ProtocolProfile(
        name="exact_boundary",
        default_layout="auto",
        default_strandedness="auto",
        default_evidence_source="auto",
        default_endpoint_model="exact_boundary",
    ),
    "generic_3prime": ProtocolProfile(
        name="generic_3prime",
        default_layout="auto",
        default_strandedness="auto",
        default_evidence_source="auto",
        default_endpoint_model="auto",
    ),
}


def get_profile(name: str) -> ProtocolProfile:
    try:
        return PROFILES[name]
    except KeyError as error:
        raise PacusageError(
            f"Unsupported library_profile={name!r}; expected one of: " + ", ".join(sorted(PROFILES))
        ) from error


def resolve_profile_defaults(
    profile_name: str,
    layout: str,
    strandedness: str,
    evidence_source: str,
    endpoint_model: str,
) -> dict[str, str]:
    profile = get_profile(profile_name)
    return {
        "library_profile": profile.name,
        "layout": profile.default_layout if layout == "auto" else layout,
        "strandedness": (profile.default_strandedness if strandedness == "auto" else strandedness),
        "evidence_source": (
            profile.default_evidence_source if evidence_source == "auto" else evidence_source
        ),
        "endpoint_model": (
            profile.default_endpoint_model if endpoint_model == "auto" else endpoint_model
        ),
    }
