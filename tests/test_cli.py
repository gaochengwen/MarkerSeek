from __future__ import annotations

import csv
import struct
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from markerseek.cli import build_parser, main


def test_cli_smoke_with_fake_mafft(tmp_path: Path):
    inputs_dir = tmp_path / "inputs"
    inputs_dir.mkdir()
    for index in range(3):
        record = build_record(index)
        path = inputs_dir / f"sample_{index + 1}.gb"
        SeqIO.write(record, path, "genbank")

    fake_mafft = tmp_path / "fake_mafft.py"
    fake_mafft.write_text(
        "\n".join(
            [
                "#!/usr/bin/env python3",
                "from Bio import SeqIO",
                "import sys",
                "records = list(SeqIO.parse(sys.argv[-1], 'fasta'))",
                "max_len = max(len(record.seq) for record in records)",
                "for record in records:",
                "    seq = str(record.seq)",
                "    record.seq = type(record.seq)(seq + ('-' * (max_len - len(seq))))",
                "SeqIO.write(records, sys.stdout, 'fasta')",
            ]
        ),
        encoding="utf-8",
    )
    fake_mafft.chmod(0o755)

    outdir = tmp_path / "output"
    exit_code = main(
        [
            "analyze",
            str(inputs_dir),
            "--outdir",
            str(outdir),
            "--mafft-bin",
            str(fake_mafft),
            "--label-mode",
            "peak-only",
            "--label-max",
            "3",
            "--label-min-distance",
            "50",
            "--lsc",
            "751:300",
            "--irb",
            "301:450",
            "--ssc",
            "451:600",
            "--ira",
            "601:750",
        ]
    )

    assert exit_code == 0
    assert (outdir / "pi_plot.pdf").exists()
    assert (outdir / "pi_plot.png").exists()
    assert (outdir / "pi_windows.tsv").exists()
    assert (outdir / "pi_features.tsv").exists()

    with (outdir / "pi_windows.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows
    assert rows[0]["window_id"] == "W0001"
    hotspot_count = sum(1 for row in rows if row["is_hotspot"] == "yes")
    assert hotspot_count == 1

    png_width, png_height = png_dimensions(outdir / "pi_plot.png")
    assert png_width > 4200
    assert png_height > 1800


def test_cli_parses_window_step_threshold_and_label_options():
    args = build_parser().parse_args(
        [
            "analyze",
            "test_data",
            "--window",
            "800",
            "--step",
            "100",
            "--hotspot-mode",
            "threshold",
            "--hotspot-value",
            "0.02",
            "--label-mode",
            "peak-only",
            "--label-max",
            "6",
            "--label-min-distance",
            "2500",
        ]
    )

    assert args.window == 800
    assert args.step == 100
    assert args.hotspot_mode == "threshold"
    assert args.hotspot_value == 0.02
    assert args.label_mode == "peak-only"
    assert args.label_max == 6
    assert args.label_min_distance == 2500


def test_cli_default_hotspot_and_label_values():
    args = build_parser().parse_args(["analyze", "test_data"])

    assert args.hotspot_mode == "top-percent"
    assert args.hotspot_value == 3.0
    assert args.label_mode == "peak-only"
    assert args.label_max is None


def build_record(index: int) -> SeqRecord:
    irb = "ATGCGTAC" * 18 + "ATGCGA"
    ira = reverse_complement(irb)
    sequence = list(("A" * 300) + irb + ("C" * 150) + ira + ("G" * 150))
    if index > 0:
        sequence[320 + index] = "T"
        sequence[470 + index] = "G"
        sequence[690 + index] = "A"
    record = SeqRecord(Seq("".join(sequence)), id=f"sample_{index + 1}", name=f"sample_{index + 1}")
    record.annotations["molecule_type"] = "DNA"
    record.annotations["topology"] = "circular"
    record.features = [
        SeqFeature(FeatureLocation(20, 120, strand=1), type="gene", qualifiers={"gene": ["psaA"]}),
        SeqFeature(FeatureLocation(330, 380, strand=1), type="gene", qualifiers={"gene": ["rrn16S"]}),
        SeqFeature(FeatureLocation(635, 685, strand=-1), type="gene", qualifiers={"gene": ["rrn16S"]}),
        SeqFeature(FeatureLocation(410, 440, strand=1), type="gene", qualifiers={"gene": ["ndhB"]}),
        SeqFeature(FeatureLocation(715, 745, strand=-1), type="gene", qualifiers={"gene": ["ndhB"]}),
    ]
    return record


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def png_dimensions(path: Path) -> tuple[int, int]:
    with path.open("rb") as handle:
        signature = handle.read(8)
        assert signature == b"\x89PNG\r\n\x1a\n"
        _chunk_length = handle.read(4)
        chunk_type = handle.read(4)
        assert chunk_type == b"IHDR"
        width, height = struct.unpack(">II", handle.read(8))
        return width, height
