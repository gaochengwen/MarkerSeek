"""Genus-level Database module for MarkerSeek.

Indexes per-genus reference GenBank collections and pre-computed MarkerSeek
analysis outputs into a single SQLite catalogue so the web UI can browse,
search, and serve results across thousands of genera.

Directory convention under MARKERSEEK_DB_ROOT:

    <root>/
    ├── taxonomy.csv              # optional; columns: genus,family,order
    ├── <Genus>_ref/              # required for "reference" genera; *.gb files
    └── <Genus>_results/          # optional; MarkerSeek analyse-output drop

A genus is included when at least one of <Genus>_ref/ or <Genus>_results/
exists. Species lists are read on demand from <Genus>_ref/ at request time
to keep the index compact.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sqlite3
import sys
import threading
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Iterable

DEFAULT_DB_FILENAME = "markerseek_database.sqlite3"
DEFAULT_DB_ROOT_ENV = "MARKERSEEK_DB_ROOT"
DEFAULT_DB_PATH_ENV = "MARKERSEEK_DATABASE_DB"

REF_SUFFIX = "_ref"
RESULTS_SUFFIX = "_results"
TAXONOMY_FILENAME = "taxonomy.csv"

GENBANK_SUFFIXES = {".gb", ".gbk", ".genbank"}

DB_LOCK = threading.Lock()

LENGTH_BUCKETS_BP = (
    (0, 120_000),
    (120_000, 130_000),
    (130_000, 140_000),
    (140_000, 150_000),
    (150_000, 160_000),
    (160_000, 170_000),
    (170_000, 200_000),
    (200_000, 250_000),
    (250_000, 350_000),
    (350_000, 10_000_000),
)
GC_BUCKETS_PCT = (
    (0.0, 30.0),
    (30.0, 32.0),
    (32.0, 34.0),
    (34.0, 36.0),
    (36.0, 38.0),
    (38.0, 40.0),
    (40.0, 42.0),
    (42.0, 44.0),
    (44.0, 46.0),
    (46.0, 100.0),
)
SCORE_BUCKETS = (
    (0.0, 10.0),
    (10.0, 20.0),
    (20.0, 30.0),
    (30.0, 40.0),
    (40.0, 50.0),
    (50.0, 60.0),
    (60.0, 70.0),
    (70.0, 80.0),
    (80.0, 90.0),
    (90.0, 100.01),
)

LOCUS_LENGTH_RE = re.compile(r"^LOCUS\s+\S+\s+(\d+)\s+bp", re.IGNORECASE)
ACCESSION_RE = re.compile(r"^ACCESSION\s+(\S+)", re.IGNORECASE)
ORGANISM_RE = re.compile(r"^\s*ORGANISM\s+(.+)$", re.IGNORECASE)


@dataclass(frozen=True)
class SpeciesEntry:
    species: str
    accession: str | None
    length_bp: int | None
    gc_percent: float | None
    file_path: str
    file_name: str


@dataclass(frozen=True)
class GenusRecord:
    genus: str
    family: str | None
    taxonomic_order: str | None
    species_count: int
    ref_dir: str | None
    results_dir: str | None
    has_primers: bool
    genome_length_min: int | None
    genome_length_max: int | None
    top_marker_label: str | None
    top_marker_score: float | None
    marker_count: int
    primer_count: int
    primer_marker_count: int
    summary_json: str
    indexed_at: str


# ----- schema ---------------------------------------------------------------

SCHEMA_SQL = """
CREATE TABLE IF NOT EXISTS genera (
    genus TEXT PRIMARY KEY,
    family TEXT,
    taxonomic_order TEXT,
    species_count INTEGER NOT NULL,
    ref_dir TEXT,
    results_dir TEXT,
    has_primers INTEGER NOT NULL,
    genome_length_min INTEGER,
    genome_length_max INTEGER,
    top_marker_label TEXT,
    top_marker_score REAL,
    marker_count INTEGER NOT NULL,
    primer_count INTEGER NOT NULL,
    primer_marker_count INTEGER NOT NULL DEFAULT 0,
    summary_json TEXT NOT NULL,
    indexed_at TEXT NOT NULL
);
"""

INDEX_SQL = (
    "CREATE INDEX IF NOT EXISTS idx_genera_family ON genera(family);",
    "CREATE INDEX IF NOT EXISTS idx_genera_order ON genera(taxonomic_order);",
    "CREATE INDEX IF NOT EXISTS idx_genera_genus_lower ON genera(genus COLLATE NOCASE);",
)


def init_database_db(db_path: Path) -> None:
    db_path.parent.mkdir(parents=True, exist_ok=True)
    with DB_LOCK, sqlite3.connect(db_path) as conn:
        conn.execute(SCHEMA_SQL)
        for stmt in INDEX_SQL:
            conn.execute(stmt)
        existing_cols = {row[1] for row in conn.execute("PRAGMA table_info(genera)")}
        if "primer_marker_count" not in existing_cols:
            conn.execute("ALTER TABLE genera ADD COLUMN primer_marker_count INTEGER NOT NULL DEFAULT 0")
        conn.commit()


def connect(db_path: Path) -> sqlite3.Connection:
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    return conn


# ----- GenBank lightweight parsing ------------------------------------------


def parse_genbank_summary(path: Path) -> SpeciesEntry | None:
    """Parse only the summary fields needed for the index.

    Reads LOCUS, ACCESSION, ORGANISM, and counts G/C bases in the ORIGIN
    section. Returns None on unreadable files so a single broken record does
    not break a reindex pass.
    """

    length_bp: int | None = None
    accession: str | None = None
    organism: str | None = None
    gc_count = 0
    at_count = 0
    n_count = 0
    in_origin = False

    try:
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if not in_origin:
                    if length_bp is None:
                        match = LOCUS_LENGTH_RE.match(line)
                        if match:
                            length_bp = int(match.group(1))
                            continue
                    if accession is None:
                        match = ACCESSION_RE.match(line)
                        if match:
                            accession = match.group(1)
                            continue
                    if organism is None:
                        match = ORGANISM_RE.match(line)
                        if match:
                            organism = match.group(1).strip()
                            continue
                    if line.startswith("ORIGIN"):
                        in_origin = True
                    continue
                if line.startswith("//"):
                    break
                for char in line:
                    lower = char.lower()
                    if lower in ("g", "c"):
                        gc_count += 1
                    elif lower in ("a", "t", "u"):
                        at_count += 1
                    elif lower.isalpha():
                        n_count += 1
    except OSError:
        return None

    species = _species_from_filename(path) if not organism else organism
    canonical_count = gc_count + at_count
    gc_percent: float | None = None
    if canonical_count > 0:
        gc_percent = round(100.0 * gc_count / canonical_count, 2)

    return SpeciesEntry(
        species=species,
        accession=accession,
        length_bp=length_bp,
        gc_percent=gc_percent,
        file_path=str(path),
        file_name=path.name,
    )


def _species_from_filename(path: Path) -> str:
    stem = path.stem
    parts = stem.split("_")
    if len(parts) >= 2:
        # "Astragalus_membranaceus_var_membranaceus" -> "Astragalus membranaceus var membranaceus"
        return " ".join(parts)
    return stem


def list_species_entries(ref_dir: Path) -> list[SpeciesEntry]:
    if not ref_dir.exists() or not ref_dir.is_dir():
        return []
    entries: list[SpeciesEntry] = []
    for path in sorted(ref_dir.iterdir()):
        if path.suffix.lower() not in GENBANK_SUFFIXES:
            continue
        if not path.is_file():
            continue
        entry = parse_genbank_summary(path)
        if entry is not None:
            entries.append(entry)
    return entries


# ----- TSV ingestion --------------------------------------------------------


def _read_marker_features(tsv_path: Path) -> list[dict[str, str]]:
    if not tsv_path.exists():
        return []
    with tsv_path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _to_float(value: Any) -> float | None:
    if value is None:
        return None
    text = str(value).strip()
    if not text or text.upper() == "NA":
        return None
    try:
        return float(text)
    except ValueError:
        return None


def _count_lines(path: Path) -> int:
    if not path.exists():
        return 0
    count = 0
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for _ in handle:
            count += 1
    return max(count - 1, 0)  # subtract header


def _count_primer_markers(primer_path: Path | None) -> int:
    """Number of distinct ``label_name`` values present in ``primers.tsv``."""

    if primer_path is None or not primer_path.exists():
        return 0
    seen: set[str] = set()
    with primer_path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            label = (row.get("label_name") or "").strip()
            if label:
                seen.add(label)
    return len(seen)


# ----- summary aggregation --------------------------------------------------


def _bucket_distribution(values: Iterable[float], buckets: tuple[tuple[float, float], ...]) -> list[dict[str, Any]]:
    counts = [0] * len(buckets)
    for value in values:
        if value is None:
            continue
        for index, (lo, hi) in enumerate(buckets):
            if lo <= value < hi:
                counts[index] += 1
                break
    return [
        {"label": _bucket_label(lo, hi), "lo": lo, "hi": hi, "count": count}
        for (lo, hi), count in zip(buckets, counts, strict=False)
    ]


def _bucket_label(lo: float, hi: float) -> str:
    if hi >= 1_000_000:
        return f">{int(lo / 1000)} kb"
    if hi > 1000:
        return f"{int(lo / 1000)}–{int(hi / 1000)} kb"
    if isinstance(lo, float) and isinstance(hi, float):
        return f"{lo:g}–{hi:g}"
    return f"{lo}–{hi}"


def aggregate_genus(
    *,
    genus: str,
    ref_dir: Path | None,
    results_dir: Path | None,
    species_entries: list[SpeciesEntry],
    family: str | None,
    taxonomic_order: str | None,
    indexed_at: str,
) -> GenusRecord:
    species_count = len(species_entries)
    lengths = [e.length_bp for e in species_entries if e.length_bp is not None]
    gcs = [e.gc_percent for e in species_entries if e.gc_percent is not None]

    marker_path = (results_dir / "candidate_marker_features.tsv") if results_dir else None
    primer_path = (results_dir / "primers.tsv") if results_dir else None

    markers = _read_marker_features(marker_path) if marker_path else []
    marker_count = len(markers)
    top_label: str | None = None
    top_score: float | None = None
    score_values: list[float] = []
    for row in markers:
        score = _to_float(row.get("markerseek_score"))
        if score is None:
            continue
        score_values.append(score)
        if top_score is None or score > top_score:
            top_score = score
            top_label = row.get("label_name") or row.get("feature_id")

    primer_count = _count_lines(primer_path) if primer_path else 0
    primer_marker_count = _count_primer_markers(primer_path)

    summary = {
        "length_distribution": _bucket_distribution(lengths, LENGTH_BUCKETS_BP),
        "gc_distribution": _bucket_distribution(gcs, GC_BUCKETS_PCT),
        "score_distribution": _bucket_distribution(score_values, SCORE_BUCKETS) if score_values else [],
        "lengths_bp": lengths,
        "gc_percent": gcs,
    }

    return GenusRecord(
        genus=genus,
        family=family,
        taxonomic_order=taxonomic_order,
        species_count=species_count,
        ref_dir=str(ref_dir) if ref_dir else None,
        results_dir=str(results_dir) if results_dir else None,
        has_primers=primer_count > 0,
        genome_length_min=min(lengths) if lengths else None,
        genome_length_max=max(lengths) if lengths else None,
        top_marker_label=top_label,
        top_marker_score=top_score,
        marker_count=marker_count,
        primer_count=primer_count,
        primer_marker_count=primer_marker_count,
        summary_json=json.dumps(summary, sort_keys=True),
        indexed_at=indexed_at,
    )


# ----- discovery ------------------------------------------------------------


def discover_genera(root: Path) -> dict[str, dict[str, Path]]:
    """Return {genus: {"ref": Path|None, "results": Path|None}}."""

    found: dict[str, dict[str, Path]] = {}
    if not root.exists() or not root.is_dir():
        return found
    for entry in sorted(root.iterdir()):
        if not entry.is_dir():
            continue
        name = entry.name
        if name.endswith(REF_SUFFIX):
            genus = name[: -len(REF_SUFFIX)]
            found.setdefault(genus, {"ref": None, "results": None})["ref"] = entry
        elif name.endswith(RESULTS_SUFFIX):
            genus = name[: -len(RESULTS_SUFFIX)]
            found.setdefault(genus, {"ref": None, "results": None})["results"] = entry
    return {genus: dirs for genus, dirs in found.items() if genus}


def load_taxonomy(path: Path) -> dict[str, dict[str, str | None]]:
    if not path.exists():
        return {}
    out: dict[str, dict[str, str | None]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            return {}
        lower_headers = {name.lower(): name for name in reader.fieldnames}
        if "genus" not in lower_headers:
            return {}
        for row in reader:
            genus = (row.get(lower_headers["genus"]) or "").strip()
            if not genus:
                continue
            family = (row.get(lower_headers.get("family", "")) or "").strip() or None
            order = (row.get(lower_headers.get("order", "")) or "").strip() or None
            out[genus] = {"family": family, "order": order}
    return out


# ----- reindex orchestration ------------------------------------------------


UPSERT_SQL = """
INSERT INTO genera (
    genus, family, taxonomic_order, species_count,
    ref_dir, results_dir, has_primers,
    genome_length_min, genome_length_max,
    top_marker_label, top_marker_score,
    marker_count, primer_count, primer_marker_count,
    summary_json, indexed_at
) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
ON CONFLICT(genus) DO UPDATE SET
    family = excluded.family,
    taxonomic_order = excluded.taxonomic_order,
    species_count = excluded.species_count,
    ref_dir = excluded.ref_dir,
    results_dir = excluded.results_dir,
    has_primers = excluded.has_primers,
    genome_length_min = excluded.genome_length_min,
    genome_length_max = excluded.genome_length_max,
    top_marker_label = excluded.top_marker_label,
    top_marker_score = excluded.top_marker_score,
    marker_count = excluded.marker_count,
    primer_count = excluded.primer_count,
    primer_marker_count = excluded.primer_marker_count,
    summary_json = excluded.summary_json,
    indexed_at = excluded.indexed_at;
"""


def reindex(
    *,
    root: Path,
    db_path: Path,
    only_genus: str | None = None,
    log: Any = None,
) -> dict[str, int]:
    init_database_db(db_path)
    taxonomy = load_taxonomy(root / TAXONOMY_FILENAME)
    discovered = discover_genera(root)
    if only_genus is not None:
        discovered = {g: d for g, d in discovered.items() if g == only_genus}
    indexed_at = datetime.now(UTC).replace(microsecond=0).isoformat().replace("+00:00", "Z")
    inserted = 0
    skipped = 0
    with DB_LOCK, sqlite3.connect(db_path) as conn:
        for genus, dirs in sorted(discovered.items()):
            ref_dir = dirs.get("ref")
            results_dir = dirs.get("results")
            species_entries = list_species_entries(ref_dir) if ref_dir else _species_from_results(results_dir)
            tax = taxonomy.get(genus, {})
            record = aggregate_genus(
                genus=genus,
                ref_dir=ref_dir,
                results_dir=results_dir,
                species_entries=species_entries,
                family=tax.get("family"),
                taxonomic_order=tax.get("order"),
                indexed_at=indexed_at,
            )
            if record.species_count == 0 and not record.results_dir:
                skipped += 1
                if log:
                    log(f"skip {genus}: empty ref directory and no results")
                continue
            conn.execute(
                UPSERT_SQL,
                (
                    record.genus,
                    record.family,
                    record.taxonomic_order,
                    record.species_count,
                    record.ref_dir,
                    record.results_dir,
                    1 if record.has_primers else 0,
                    record.genome_length_min,
                    record.genome_length_max,
                    record.top_marker_label,
                    record.top_marker_score,
                    record.marker_count,
                    record.primer_count,
                    record.primer_marker_count,
                    record.summary_json,
                    record.indexed_at,
                ),
            )
            inserted += 1
            if log:
                log(
                    f"indexed {genus}: species={record.species_count} "
                    f"markers={record.marker_count} primers={record.primer_count}"
                )
        conn.commit()
    return {"indexed": inserted, "skipped": skipped, "total_seen": len(discovered)}


def _species_from_results(results_dir: Path | None) -> list[SpeciesEntry]:
    """Fallback: derive a synthetic species list from sample_metadata.tsv when no _ref/ exists."""

    if results_dir is None:
        return []
    meta_path = results_dir / "sample_metadata.tsv"
    if not meta_path.exists():
        return []
    entries: list[SpeciesEntry] = []
    with meta_path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            species = (row.get("species") or row.get("sample_name") or "").strip()
            if not species:
                continue
            entries.append(
                SpeciesEntry(
                    species=species,
                    accession=None,
                    length_bp=None,
                    gc_percent=None,
                    file_path=row.get("source_path", ""),
                    file_name=Path(row.get("source_path", "")).name,
                )
            )
    return entries


# ----- read access used by web layer ----------------------------------------


def list_genera(
    db_path: Path,
    *,
    query: str | None = None,
    family: str | None = None,
    taxonomic_order: str | None = None,
    page: int = 1,
    page_size: int = 50,
) -> tuple[list[sqlite3.Row], int]:
    page = max(1, page)
    page_size = max(1, min(page_size, 200))
    where: list[str] = []
    params: list[Any] = []
    if query:
        where.append("genus LIKE ? COLLATE NOCASE")
        params.append(f"{query}%")
    if family:
        where.append("family = ?")
        params.append(family)
    if taxonomic_order:
        where.append("taxonomic_order = ?")
        params.append(taxonomic_order)
    where_clause = (" WHERE " + " AND ".join(where)) if where else ""
    with connect(db_path) as conn:
        total = conn.execute(f"SELECT COUNT(*) FROM genera{where_clause}", params).fetchone()[0]
        rows = list(
            conn.execute(
                f"SELECT * FROM genera{where_clause} ORDER BY genus LIMIT ? OFFSET ?",
                (*params, page_size, (page - 1) * page_size),
            )
        )
    return rows, total


def get_genus(db_path: Path, genus: str) -> sqlite3.Row | None:
    with connect(db_path) as conn:
        return conn.execute("SELECT * FROM genera WHERE genus = ?", (genus,)).fetchone()


def list_taxonomy_facets(db_path: Path) -> dict[str, list[str]]:
    with connect(db_path) as conn:
        families = [
            row[0]
            for row in conn.execute(
                "SELECT DISTINCT family FROM genera WHERE family IS NOT NULL ORDER BY family"
            )
        ]
        orders = [
            row[0]
            for row in conn.execute(
                "SELECT DISTINCT taxonomic_order FROM genera WHERE taxonomic_order IS NOT NULL ORDER BY taxonomic_order"
            )
        ]
    return {"families": families, "orders": orders}


def database_stats(db_path: Path) -> dict[str, int]:
    if not db_path.exists():
        return {
            "genus_count": 0,
            "with_results": 0,
            "with_primers": 0,
            "species_total": 0,
            "primer_marker_total": 0,
        }
    with connect(db_path) as conn:
        genus_count = conn.execute("SELECT COUNT(*) FROM genera").fetchone()[0]
        with_results = conn.execute(
            "SELECT COUNT(*) FROM genera WHERE marker_count > 0"
        ).fetchone()[0]
        with_primers = conn.execute(
            "SELECT COUNT(*) FROM genera WHERE has_primers = 1"
        ).fetchone()[0]
        species_total = conn.execute(
            "SELECT COALESCE(SUM(species_count), 0) FROM genera"
        ).fetchone()[0]
        primer_marker_total = conn.execute(
            "SELECT COALESCE(SUM(primer_marker_count), 0) FROM genera"
        ).fetchone()[0]
    return {
        "genus_count": int(genus_count),
        "with_results": int(with_results),
        "with_primers": int(with_primers),
        "species_total": int(species_total),
        "primer_marker_total": int(primer_marker_total),
    }


def alphabet_index(db_path: Path) -> list[dict[str, Any]]:
    """Return [{letter, count}, ...] for A-Z plus '#' for non-alpha."""

    if not db_path.exists():
        return []
    counts: dict[str, int] = {}
    with connect(db_path) as conn:
        for (genus,) in conn.execute("SELECT genus FROM genera"):
            initial = (genus[:1] or "#").upper()
            if not initial.isalpha():
                initial = "#"
            counts[initial] = counts.get(initial, 0) + 1
    letters = ["#"] + [chr(c) for c in range(ord("A"), ord("Z") + 1)]
    return [{"letter": ch, "count": counts.get(ch, 0)} for ch in letters]


# ----- CLI ------------------------------------------------------------------


def _resolve_db_path(args: argparse.Namespace, default_root: Path) -> Path:
    if args.db_path:
        return Path(args.db_path).expanduser().resolve()
    env_db = os.environ.get(DEFAULT_DB_PATH_ENV)
    if env_db:
        return Path(env_db).expanduser().resolve()
    return (default_root / DEFAULT_DB_FILENAME).resolve()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="markerseek-db",
        description="Index per-genus MarkerSeek references and results into a SQLite catalogue.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    reindex_p = sub.add_parser("reindex", help="Scan the database root and (re)index every discovered genus.")
    reindex_p.add_argument(
        "--root",
        default=os.environ.get(DEFAULT_DB_ROOT_ENV, "/home/ubuntu/plant_plastid"),
        help=f"Directory holding <Genus>_ref / <Genus>_results subfolders. "
             f"Defaults to ${DEFAULT_DB_ROOT_ENV} or /home/ubuntu/plant_plastid.",
    )
    reindex_p.add_argument("--db-path", help="SQLite catalogue path (default: <root>/markerseek_database.sqlite3).")
    reindex_p.add_argument("--genus", help="Reindex only this genus (incremental update).")
    reindex_p.add_argument("--quiet", action="store_true", help="Suppress per-genus progress lines.")

    stats_p = sub.add_parser("stats", help="Print catalogue counts.")
    stats_p.add_argument("--root", default=os.environ.get(DEFAULT_DB_ROOT_ENV, "/home/ubuntu/plant_plastid"))
    stats_p.add_argument("--db-path")

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    root = Path(args.root).expanduser().resolve()
    db_path = _resolve_db_path(args, root)

    if args.command == "reindex":
        log = (lambda msg: None) if args.quiet else (lambda msg: print(msg, file=sys.stderr))
        if not root.exists():
            print(f"error: root directory does not exist: {root}", file=sys.stderr)
            return 2
        result = reindex(root=root, db_path=db_path, only_genus=args.genus, log=log)
        print(
            f"Indexed {result['indexed']} genera "
            f"(skipped {result['skipped']}, seen {result['total_seen']}). DB: {db_path}"
        )
        return 0

    if args.command == "stats":
        stats = database_stats(db_path)
        print(json.dumps(stats, indent=2))
        return 0

    parser.error(f"Unknown command: {args.command}")
    return 2


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
