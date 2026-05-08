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
    variable_sites: int = 0
    indel_sites: int = 0
    conserved_left_bp: int = 0
    conserved_right_bp: int = 0
    primer_available: str = "NA"
    species_resolution: float | None = None
    unique_haplotype_count: int = 0
    species_specific_haplotype_ratio: float | None = None
    interspecific_divergence: float | None = None
    intraspecific_divergence: float | None = None
    nearest_neighbor_discrimination: float | None = None
    barcoding_gap: float | None = None
    misclassification_risk: float | None = None
    alignment_reliability: float | None = None
    missing_ambig_ratio: float | None = None
    markerseek_score: float | None = None
    haplotypes: list[str] = field(default_factory=list)

    def spans(self, genome_length: int) -> list[tuple[int, int]]:
        """Return 0-based half-open spans for this feature on a circular genome."""

        start = self.start - 1
        end = self.end
        if not self.spans_origin:
            return [(start, end)]
        return [(start, genome_length), (0, end)]


@dataclass
class PrimerResult:
    label_name: str
    rank: int
    fwd_seq: str
    rev_seq: str
    fwd_len: int
    rev_len: int
    fwd_gc: float
    rev_gc: float
    fwd_tm: float
    rev_tm: float
    fwd_self_any_th: float
    rev_self_any_th: float
    primer3_penalty: float
    target_start: int
    target_end: int
    amplicon_min_len: int
    amplicon_max_len: int
    amplicon_mean_len: float
    cross_species_success_rate: float
    amplicon_variable_sites: int
    amplicon_indel_sites: int
    primer_score: float
    primer_id: str = field(init=False)

    def __post_init__(self) -> None:
        self.primer_id = f"{self.label_name}_p{self.rank}"


@dataclass
class SampleMetadata:
    sample_name: str
    species: str
    source_path: str


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
    sample_metadata: list[SampleMetadata] = field(default_factory=list)
    score_weights: dict[str, float] = field(default_factory=dict)
    primers: list[PrimerResult] = field(default_factory=list)
    primer_amplicon_samples: dict[str, set[str]] = field(default_factory=dict)
