"""Species metadata helpers for MarkerSeek inputs."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Iterable

from .models import SampleMetadata


def infer_species_from_record(record, fallback_stem: str) -> str:
    """Infer a species label from GenBank organism metadata or a filename stem."""

    annotations = getattr(record, "annotations", {}) or {}
    organism = str(annotations.get("organism", "")).strip()
    species = organism or "_".join(fallback_stem.split("_")[:2])
    return re.sub(r"\s+", "_", species.strip()) or fallback_stem


def build_sample_metadata(genome_records: Iterable, sample_order: list[str]) -> list[SampleMetadata]:
    """Build sample metadata keyed by the sanitized sample names used in analysis."""

    records_by_sample = {genome_record.sample_name: genome_record for genome_record in genome_records}
    metadata: list[SampleMetadata] = []
    for sample_name in sample_order:
        genome_record = records_by_sample.get(sample_name)
        if genome_record is None:
            metadata.append(
                SampleMetadata(
                    sample_name=sample_name,
                    species=infer_species_from_record(object(), sample_name),
                    source_path="",
                )
            )
            continue
        path = Path(genome_record.path)
        metadata.append(
            SampleMetadata(
                sample_name=sample_name,
                species=infer_species_from_record(genome_record.record, path.stem),
                source_path=str(path),
            )
        )
    return metadata
