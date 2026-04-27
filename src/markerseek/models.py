"""Core data models for MarkerSeek."""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class RegionSegment:
    """A named region on a circular genome."""

    name: str
    start: int
    end: int

    @property
    def wraps_origin(self) -> bool:
        return self.end <= self.start

    def spans(self, genome_length: int) -> list[tuple[int, int]]:
        if not self.wraps_origin:
            return [(self.start, self.end)]
        return [(self.start, genome_length), (0, self.end)]

    def length(self, genome_length: int) -> int:
        if not self.wraps_origin:
            return self.end - self.start
        return (genome_length - self.start) + self.end

    def contains(self, position: int, genome_length: int) -> bool:
        for start, end in self.spans(genome_length):
            if start <= position < end:
                return True
        return False


@dataclass(frozen=True)
class AnnotatedInterval:
    """A gene, tRNA, rRNA, or intergenic spacer projected onto reference coordinates."""

    feature_id: str
    label: str
    feature_type: str
    parent_gene: str
    start: int
    end: int
    strand: int = 0
    region: str = ""

    @property
    def wraps_origin(self) -> bool:
        return self.end <= self.start

    def spans(self, genome_length: int) -> list[tuple[int, int]]:
        if not self.wraps_origin:
            return [(self.start, self.end)]
        return [(self.start, genome_length), (0, self.end)]

    def length(self, genome_length: int) -> int:
        if not self.wraps_origin:
            return self.end - self.start
        return (genome_length - self.start) + self.end

    def midpoint(self, genome_length: int) -> int:
        return (self.start + (self.length(genome_length) // 2)) % genome_length


@dataclass
class WindowResult:
    window_id: str
    start: int
    end: int
    midpoint: int
    pi: float | None
    valid_sites: int
    region: str
    label_name: str
    is_hotspot: bool = False
    hotspot_rank: int | None = None


@dataclass
class FeatureResult:
    feature_id: str
    feature_type: str
    parent_gene: str
    label_name: str
    start: int
    end: int
    strand: int
    length_bp: int
    region: str
    pi: float | None
    valid_sites: int
    spans_origin: bool
    is_hotspot: bool = False
    hotspot_rank: int | None = None


@dataclass
class AnalysisResult:
    reference_name: str
    genome_length: int
    sample_count: int
    regions: list[RegionSegment]
    position_pi: list[float | None]
    windows: list[WindowResult]
    features: list[FeatureResult]
    aligned_sequences: dict[str, str] = field(default_factory=dict)
    exon_intervals: list[AnnotatedInterval] = field(default_factory=list)
    sample_order: list[str] = field(default_factory=list)
