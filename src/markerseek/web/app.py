"""FastAPI web application for MarkerSeek."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import queue
import re
import shutil
import sqlite3
import threading
import traceback
from contextlib import asynccontextmanager
from datetime import UTC, datetime, timedelta
from io import StringIO
from pathlib import Path
from secrets import randbelow
from statistics import median
from typing import Any
from urllib.parse import quote

import uvicorn
from fastapi import FastAPI, File, Form, HTTPException, Request, UploadFile
from fastapi.responses import FileResponse, HTMLResponse, JSONResponse, RedirectResponse, Response
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates

from markerseek import database as db_catalog
from markerseek.analysis import MarkerSeekError, run_analysis
from markerseek.models import AnalysisResult, WindowResult
from markerseek.output import (
    MARKER_FEATURE_COLUMNS,
    PRIMER_COLUMNS,
    RESULT_FILENAMES,
    create_results_zip,
    feature_id_safe,
    format_float,
    write_analysis_outputs,
    write_json,
)
from markerseek.plotting import select_label_windows
from markerseek.visualisation import render_alignment_svg_block

PACKAGE_DIR = Path(__file__).resolve().parent
TEMPLATE_DIR = PACKAGE_DIR / "templates"
STATIC_DIR = PACKAGE_DIR / "static"

SUPPORTED_UPLOAD_SUFFIXES = {".gb", ".gbk", ".genbank"}
JOB_STATUS_ORDER = {"queued", "running", "succeeded", "failed", "expired"}
DOWNLOADABLE_FILES = set(RESULT_FILENAMES) | {"results.zip"}
DOWNLOAD_LIST_HIDDEN_SUFFIXES = (".png",)
FEATURE_DATA_KEYS = ("pi_curve", "snp_positions", "indel_positions", "haplotype_network", "species_pca", "primers")
INTRA_DEPENDENT_METRICS = {
    "intraspecific_divergence",
    "nearest_neighbor_discrimination",
    "barcoding_gap",
    "misclassification_risk",
}
INSUFFICIENT_REPLICATES_TEXT = "Insufficient replicates"
NA_DIAGNOSTIC_EXPLANATION = (
    "nearest_neighbor_discrimination and barcoding_gap require intraspecific data. They are reported as NA "
    "when no species has at least two samples, or when the alignment does not provide enough valid "
    "within-species pairwise distances to estimate the statistic."
)
CANDIDATE_MARKER_EXPLANATION = (
    "Candidate markers include all hotspot labels marked in the diversity plot, listed in genomic order, "
    "plus the top 10 remaining regions by MarkerSeek score. Primer design is attempted for every listed region."
)
PRIMER_AVAILABLE_EXPLANATION = (
    "Primer available: yes means at least one valid primer pair was found; no means no valid pair passed "
    "screening; NA means primer design was skipped."
)
DEMO_JOB_ID = "MSK-EXAMPLE-DEMO"
FEATURE_METADATA_KEYS = (
    "label_name",
    "feature_type",
    "region",
    "length_bp",
    "variable_sites",
    "pi",
    "unique_haplotype_count",
    "primer_available",
    "alignment_reliability",
    "markerseek_score",
)
INTEGER_FEATURE_METADATA_KEYS = {"length_bp", "variable_sites", "unique_haplotype_count"}
TWO_DECIMAL_FEATURE_METADATA_KEYS = {"alignment_reliability", "markerseek_score"}

MAX_FILES = 20
MIN_FILES = 2
DEFAULT_RETENTION_DAYS = 7
DEFAULT_MAX_UPLOAD_BYTES = 20 * 1024 * 1024
DEFAULT_MAFFT_THREADS = min(4, os.cpu_count() or 1)
DEFAULT_PRIMER_TM_MIN = 52.0
DEFAULT_PRIMER_TM_OPT = 58.0
DEFAULT_PRIMER_TM_MAX = 70.0
DEFAULT_PRIMER_LEN_MIN = 18
DEFAULT_PRIMER_LEN_MAX = 27
DEFAULT_PRIMER_AMPLICON_MIN = 80
DEFAULT_PRIMER_AMPLICON_MAX = 3000
DEFAULT_PRIMER_MISMATCH = 1
DEFAULT_PRIMER_ANCHOR_BP = 5
DEFAULT_PRIMER_NUM_RETURN = 5
DEFAULT_PRIMER_FLANK_MAX_LEN = 120
SAFE_FILENAME_RE = re.compile(r"[^A-Za-z0-9_.-]+")


class WebSettings:
    def __init__(self) -> None:
        self.data_dir = Path(os.environ.get("MARKERSEEK_WEB_DATA", "markerseek_jobs")).resolve()
        self.retention_days = int(os.environ.get("MARKERSEEK_RETENTION_DAYS", str(DEFAULT_RETENTION_DAYS)))
        self.max_upload_bytes = int(os.environ.get("MARKERSEEK_MAX_UPLOAD_BYTES", str(DEFAULT_MAX_UPLOAD_BYTES)))
        self.mafft_bin = os.environ.get("MARKERSEEK_MAFFT_BIN", "mafft")
        self.database_root = Path(
            os.environ.get(db_catalog.DEFAULT_DB_ROOT_ENV, "/home/ubuntu/plant_plastid")
        ).resolve()

    @property
    def jobs_dir(self) -> Path:
        return self.data_dir / "jobs"

    @property
    def db_path(self) -> Path:
        return self.data_dir / "jobs.sqlite3"

    @property
    def database_db_path(self) -> Path:
        custom = os.environ.get(db_catalog.DEFAULT_DB_PATH_ENV)
        if custom:
            return Path(custom).expanduser().resolve()
        return self.database_root / db_catalog.DEFAULT_DB_FILENAME


settings = WebSettings()
templates = Jinja2Templates(directory=str(TEMPLATE_DIR))
templates.env.filters["fmt1"] = lambda value: format_template_number(value, digits=1)
templates.env.filters["fmt2"] = lambda value: format_template_number(value, digits=2)
templates.env.filters["fmt3"] = lambda value: format_template_number(value, digits=3)
templates.env.filters["fmt4"] = lambda value: format_template_number(value, digits=4)
templates.env.filters["sig3"] = lambda value: format_significant_number(value, digits=3)
templates.env.filters["is_na"] = lambda value: is_na_value(value)

JOB_QUEUE: queue.Queue[str] = queue.Queue()
DB_LOCK = threading.Lock()
WORKER_STARTED = False


@asynccontextmanager
async def lifespan(_app: FastAPI):
    startup_services()
    yield


app = FastAPI(title="MarkerSeek Analyzer", lifespan=lifespan)
app.mount("/static", StaticFiles(directory=str(STATIC_DIR)), name="static")


def utcnow() -> datetime:
    return datetime.now(UTC)


def isoformat(value: datetime) -> str:
    return value.replace(microsecond=0).isoformat().replace("+00:00", "Z")


def parse_iso(value: str) -> datetime:
    return datetime.fromisoformat(value.replace("Z", "+00:00"))


def safe_filename(raw_name: str) -> str:
    name = Path(raw_name or "").name
    name = SAFE_FILENAME_RE.sub("_", name).strip("._")
    return name or "input.gb"


def job_dir(job_id: str) -> Path:
    return settings.jobs_dir / job_id


def inputs_dir(job_id: str) -> Path:
    return job_dir(job_id) / "inputs"


def outputs_dir(job_id: str) -> Path:
    return job_dir(job_id) / "outputs"


def feature_payload_dir(job_id: str) -> Path:
    return outputs_dir(job_id) / "feature_payload"


def init_db() -> None:
    settings.data_dir.mkdir(parents=True, exist_ok=True)
    settings.jobs_dir.mkdir(parents=True, exist_ok=True)
    with DB_LOCK, sqlite3.connect(settings.db_path) as conn:
        conn.execute(
            """
            CREATE TABLE IF NOT EXISTS jobs (
                job_id TEXT PRIMARY KEY,
                status TEXT NOT NULL,
                project_name TEXT NOT NULL,
                created_at TEXT NOT NULL,
                updated_at TEXT NOT NULL,
                started_at TEXT,
                completed_at TEXT,
                expires_at TEXT NOT NULL,
                params_json TEXT NOT NULL,
                input_files_json TEXT NOT NULL,
                reference_file TEXT,
                error TEXT,
                permanent INTEGER NOT NULL DEFAULT 0
            )
            """
        )
        columns = {row[1] for row in conn.execute("PRAGMA table_info(jobs)").fetchall()}
        if "permanent" not in columns:
            conn.execute("ALTER TABLE jobs ADD COLUMN permanent INTEGER NOT NULL DEFAULT 0")
        conn.commit()


def db_execute(sql: str, values: tuple[Any, ...] = ()) -> None:
    with DB_LOCK, sqlite3.connect(settings.db_path) as conn:
        conn.execute(sql, values)
        conn.commit()


def db_fetchone(sql: str, values: tuple[Any, ...] = ()) -> sqlite3.Row | None:
    with DB_LOCK, sqlite3.connect(settings.db_path) as conn:
        conn.row_factory = sqlite3.Row
        return conn.execute(sql, values).fetchone()


def db_fetchall(sql: str, values: tuple[Any, ...] = ()) -> list[sqlite3.Row]:
    with DB_LOCK, sqlite3.connect(settings.db_path) as conn:
        conn.row_factory = sqlite3.Row
        return list(conn.execute(sql, values).fetchall())


def create_job_record(
    *,
    job_id: str,
    project_name: str,
    params: dict[str, Any],
    input_files: list[str],
    reference_file: str | None,
    permanent: bool = False,
) -> None:
    now = utcnow()
    expires_at = now + timedelta(days=settings.retention_days)
    db_execute(
        """
        INSERT INTO jobs (
            job_id, status, project_name, created_at, updated_at, expires_at,
            params_json, input_files_json, reference_file, permanent
        )
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            job_id,
            "queued",
            project_name,
            isoformat(now),
            isoformat(now),
            "permanent" if permanent else isoformat(expires_at),
            json.dumps(params, sort_keys=True),
            json.dumps(input_files),
            reference_file,
            1 if permanent else 0,
        ),
    )


def update_job(job_id: str, **fields: Any) -> None:
    if not fields:
        return
    fields["updated_at"] = isoformat(utcnow())
    assignments = ", ".join(f"{name} = ?" for name in fields)
    values = tuple(fields.values()) + (job_id,)
    db_execute(f"UPDATE jobs SET {assignments} WHERE job_id = ?", values)


def get_job(job_id: str) -> sqlite3.Row | None:
    return db_fetchone("SELECT * FROM jobs WHERE job_id = ?", (job_id,))


def generate_job_id() -> str:
    today = utcnow().strftime("%Y%m%d")
    while True:
        job_id = f"MSK-{today}{randbelow(10_000_000):07d}"
        if get_job(job_id) is None:
            return job_id


def parse_int(value: str, field_name: str, *, minimum: int = 1) -> int:
    try:
        parsed = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{field_name} must be an integer.") from exc
    if parsed < minimum:
        raise ValueError(f"{field_name} must be at least {minimum}.")
    return parsed


def parse_optional_int(value: str, field_name: str, *, minimum: int = 1) -> int | None:
    if value is None or value.strip() == "":
        return None
    return parse_int(value, field_name, minimum=minimum)


def parse_float(value: str, field_name: str, *, minimum: float | None = None, maximum: float | None = None) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{field_name} must be a number.") from exc
    if minimum is not None and parsed < minimum:
        raise ValueError(f"{field_name} must be at least {minimum}.")
    if maximum is not None and parsed > maximum:
        raise ValueError(f"{field_name} must be at most {maximum}.")
    return parsed


def parse_job_params(form: dict[str, str]) -> dict[str, Any]:
    hotspot_mode = form["hotspot_mode"]
    label_mode = form["label_mode"]
    if hotspot_mode not in {"top-percent", "top-n", "threshold"}:
        raise ValueError("Unsupported hotspot mode.")
    if label_mode not in {"peak-only", "all", "none"}:
        raise ValueError("Unsupported label mode.")
    hotspot_value = parse_float(form["hotspot_value"], "Hotspot value", minimum=0.0)
    if hotspot_mode == "top-percent" and hotspot_value <= 0:
        raise ValueError("Top-percent hotspot value must be greater than 0.")
    if hotspot_mode == "top-n" and hotspot_value < 1:
        raise ValueError("Top-N hotspot value must be at least 1.")

    max_threads = os.cpu_count() or 1
    mafft_threads = parse_int(form["mafft_threads"], "MAFFT threads", minimum=1)
    if mafft_threads > max_threads:
        raise ValueError(f"MAFFT threads cannot exceed available CPU cores ({max_threads}).")

    primer_tm_min = parse_float(form["primer_tm_min"], "Primer minimum Tm", minimum=0.0)
    primer_tm_opt = parse_float(form["primer_tm_opt"], "Primer optimal Tm", minimum=0.0)
    primer_tm_max = parse_float(form["primer_tm_max"], "Primer maximum Tm", minimum=0.0)
    if primer_tm_min > primer_tm_max:
        raise ValueError("Primer minimum Tm cannot exceed primer maximum Tm.")
    if not primer_tm_min <= primer_tm_opt <= primer_tm_max:
        raise ValueError("Primer optimal Tm must be between primer minimum and maximum Tm.")

    primer_len_min = parse_int(form["primer_len_min"], "Primer minimum length", minimum=1)
    primer_len_opt = parse_int(form["primer_len_opt"], "Primer optimal length", minimum=1)
    primer_len_max = parse_int(form["primer_len_max"], "Primer maximum length", minimum=1)
    if primer_len_min > primer_len_max:
        raise ValueError("Primer minimum length cannot exceed primer maximum length.")
    if not primer_len_min <= primer_len_opt <= primer_len_max:
        raise ValueError("Primer optimal length must be between primer minimum and maximum length.")

    primer_amplicon_min = parse_int(form["primer_amplicon_min"], "Primer minimum amplicon", minimum=1)
    primer_amplicon_max = parse_int(form["primer_amplicon_max"], "Primer maximum amplicon", minimum=1)
    if primer_amplicon_min > primer_amplicon_max:
        raise ValueError("Primer minimum amplicon cannot exceed primer maximum amplicon.")

    return {
        "window": parse_int(form["window"], "Window size", minimum=1),
        "step": parse_int(form["step"], "Step size", minimum=1),
        "hotspot_mode": hotspot_mode,
        "hotspot_value": hotspot_value,
        "label_mode": label_mode,
        "label_max": parse_optional_int(form.get("label_max", ""), "Label maximum", minimum=1),
        "label_min_distance": parse_int(form["label_min_distance"], "Label minimum distance", minimum=0),
        "similarity_window": parse_int(form["similarity_window"], "Similarity window", minimum=1),
        "similarity_step": parse_int(form["similarity_step"], "Similarity step", minimum=1),
        "similarity_floor": parse_float(form["similarity_floor"], "Similarity floor", minimum=0.0, maximum=1.0),
        "no_similarity_plot": form.get("no_similarity_plot") == "on",
        "primer_design": form.get("skip_primer_design") != "on",
        "primer_tm_min": primer_tm_min,
        "primer_tm_opt": primer_tm_opt,
        "primer_tm_max": primer_tm_max,
        "primer_len_min": primer_len_min,
        "primer_len_opt": primer_len_opt,
        "primer_len_max": primer_len_max,
        "primer_amplicon_min": primer_amplicon_min,
        "primer_amplicon_max": primer_amplicon_max,
        "primer_mismatch": parse_int(form["primer_mismatch"], "Primer mismatch", minimum=0),
        "primer_anchor_bp": parse_int(form["primer_anchor_bp"], "Primer 3-prime anchor", minimum=1),
        "primer_num_return": parse_int(form["primer_num_return"], "Primer pairs to evaluate", minimum=1),
        "primer_flank_max_len": parse_int(form["primer_flank_max_len"], "Primer flank maximum length", minimum=1),
        "mafft_threads": mafft_threads,
    }


def default_form_values() -> dict[str, str]:
    return {
        "project_name": "MarkerSeek analysis",
        "window": "600",
        "step": "200",
        "hotspot_mode": "top-percent",
        "hotspot_value": "3.0",
        "label_mode": "peak-only",
        "label_max": "",
        "label_min_distance": "0",
        "similarity_window": "200",
        "similarity_step": "60",
        "similarity_floor": "0.5",
        "mafft_threads": str(DEFAULT_MAFFT_THREADS),
        "no_similarity_plot": "",
        "skip_primer_design": "",
        "primer_tm_min": str(DEFAULT_PRIMER_TM_MIN),
        "primer_tm_opt": str(DEFAULT_PRIMER_TM_OPT),
        "primer_tm_max": str(DEFAULT_PRIMER_TM_MAX),
        "primer_len_min": str(DEFAULT_PRIMER_LEN_MIN),
        "primer_len_opt": "20",
        "primer_len_max": str(DEFAULT_PRIMER_LEN_MAX),
        "primer_amplicon_min": str(DEFAULT_PRIMER_AMPLICON_MIN),
        "primer_amplicon_max": str(DEFAULT_PRIMER_AMPLICON_MAX),
        "primer_mismatch": str(DEFAULT_PRIMER_MISMATCH),
        "primer_anchor_bp": str(DEFAULT_PRIMER_ANCHOR_BP),
        "primer_num_return": str(DEFAULT_PRIMER_NUM_RETURN),
        "primer_flank_max_len": str(DEFAULT_PRIMER_FLANK_MAX_LEN),
    }


def render_home(request: Request, *, error: str | None = None, form_values: dict[str, str] | None = None) -> HTMLResponse:
    values = default_form_values()
    if form_values:
        values.update(form_values)
    return templates.TemplateResponse(
        request,
        "home.html",
        {
            "error": error,
            "values": values,
            "max_files": MAX_FILES,
            "min_files": MIN_FILES,
            "max_upload_mb": settings.max_upload_bytes // (1024 * 1024),
            "retention_days": settings.retention_days,
            "max_threads": os.cpu_count() or 1,
        },
    )


def display_job_params(params: dict[str, Any]) -> list[tuple[str, Any]]:
    hidden = {
        "label_max",
        "label_min_distance",
        "label_mode",
        "mafft_threads",
        "similarity_floor",
        "primer_tm_min",
        "primer_tm_opt",
        "primer_tm_max",
        "primer_len_opt",
        "primer_anchor_bp",
        "primer_flank_max_len",
    }
    labels = {
        "window": "hotspot_window",
        "step": "hotspot_step",
        "primer_len_min": "primer_length_min",
        "primer_len_max": "primer_length_max",
        "primer_amplicon_min": "primer_amplicon_min",
        "primer_amplicon_max": "primer_amplicon_max",
        "primer_mismatch": "primer_mismatch",
        "primer_num_return": "primer_pairs_returned",
    }
    priority = [
        "hotspot_mode",
        "hotspot_value",
        "window",
        "step",
        "similarity_window",
        "similarity_step",
        "primer_design",
        "primer_len_min",
        "primer_len_max",
        "primer_amplicon_min",
        "primer_amplicon_max",
        "primer_mismatch",
        "primer_num_return",
        "no_similarity_plot",
    ]
    ordered_keys = [key for key in priority if key in params]
    ordered_keys.extend(key for key in params if key not in hidden and key not in ordered_keys)
    return [(labels.get(key, key), params[key]) for key in ordered_keys if key not in hidden]


def primer_settings_from_params(params: dict[str, Any]) -> dict[str, Any]:
    return {
        "PRIMER_MIN_TM": params.get("primer_tm_min", DEFAULT_PRIMER_TM_MIN),
        "PRIMER_OPT_TM": params.get("primer_tm_opt", DEFAULT_PRIMER_TM_OPT),
        "PRIMER_MAX_TM": params.get("primer_tm_max", DEFAULT_PRIMER_TM_MAX),
        "PRIMER_MIN_SIZE": params.get("primer_len_min", DEFAULT_PRIMER_LEN_MIN),
        "PRIMER_OPT_SIZE": params.get("primer_len_opt", 20),
        "PRIMER_MAX_SIZE": params.get("primer_len_max", DEFAULT_PRIMER_LEN_MAX),
        "PRIMER_PRODUCT_SIZE_RANGE": [[
            params.get("primer_amplicon_min", DEFAULT_PRIMER_AMPLICON_MIN),
            params.get("primer_amplicon_max", DEFAULT_PRIMER_AMPLICON_MAX),
        ]],
        "PRIMER_NUM_RETURN": params.get("primer_num_return", DEFAULT_PRIMER_NUM_RETURN),
    }


def in_silico_settings_from_params(params: dict[str, Any]) -> dict[str, Any]:
    return {
        "min_amplicon": params.get("primer_amplicon_min", DEFAULT_PRIMER_AMPLICON_MIN),
        "max_amplicon": params.get("primer_amplicon_max", DEFAULT_PRIMER_AMPLICON_MAX),
        "max_mismatch": params.get("primer_mismatch", DEFAULT_PRIMER_MISMATCH),
        "anchor_3prime": params.get("primer_anchor_bp", DEFAULT_PRIMER_ANCHOR_BP),
        "flank_max_len": params.get("primer_flank_max_len", DEFAULT_PRIMER_FLANK_MAX_LEN),
    }


def is_na_value(value: Any) -> bool:
    if value is None:
        return True
    if isinstance(value, str):
        return value.strip().upper() in {"", "NA", "NAN", "NONE", "NULL"}
    if isinstance(value, float):
        return math.isnan(value)
    return False


def format_template_number(value: Any, *, digits: int = 2) -> str:
    if is_na_value(value):
        return "NA"
    if isinstance(value, bool):
        return str(value)
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return str(value)
    return f"{numeric:.{digits}f}"


def format_significant_number(value: Any, *, digits: int = 3) -> str:
    if is_na_value(value):
        return "NA"
    if isinstance(value, bool):
        return str(value)
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return str(value)
    return f"{numeric:.{digits}g}"


def format_integer_value(value: Any) -> str:
    if is_na_value(value):
        return "NA"
    try:
        return str(int(float(value)))
    except (TypeError, ValueError):
        return str(value)


def format_feature_metadata_value(key: str, value: Any) -> str:
    if key in INTEGER_FEATURE_METADATA_KEYS:
        return format_integer_value(value)
    if key == "pi":
        return format_significant_number(value, digits=3)
    if key in TWO_DECIMAL_FEATURE_METADATA_KEYS:
        return format_template_number(value, digits=2)
    if is_na_value(value):
        return "NA"
    return str(value)


def feature_metadata_rows(feature: dict[str, Any]) -> list[tuple[str, str]]:
    return [
        (key, format_feature_metadata_value(key, feature.get(key)))
        for key in FEATURE_METADATA_KEYS
        if key in feature
    ]


def parse_tsv_float(value: Any) -> float | None:
    if is_na_value(value):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def parse_tsv_int(value: Any) -> int:
    try:
        return int(float(str(value)))
    except (TypeError, ValueError):
        return 0


def original_hotspot_ranks(
    rows: list[dict[str, str]],
    params: dict[str, Any] | None,
) -> dict[int, int]:
    mode = str((params or {}).get("hotspot_mode", "top-percent"))
    value = parse_tsv_float((params or {}).get("hotspot_value", 3.0))
    if value is None or mode not in {"top-percent", "top-n", "threshold"}:
        return {}

    valid_rows: list[tuple[int, dict[str, str], float]] = []
    for index, row in enumerate(rows):
        pi = parse_tsv_float(row.get("pi"))
        if pi is not None:
            valid_rows.append((index, row, pi))

    ordered = sorted(
        valid_rows,
        key=lambda item: (
            -item[2],
            parse_tsv_int(item[1].get("start")),
            parse_tsv_int(item[1].get("end")),
            item[1].get("label_name", ""),
        ),
    )
    if mode == "threshold":
        selected_indices = {index for index, _row, pi in ordered if pi >= value}
    elif mode == "top-n":
        count = max(1, int(value))
        selected_indices = {index for index, _row, _pi in ordered[:count]}
    else:
        count = max(1, math.ceil(len(ordered) * (value / 100.0)))
        selected_indices = {index for index, _row, _pi in ordered[:count]}

    ranks: dict[int, int] = {}
    rank = 1
    for index, _row, _pi in ordered:
        if index in selected_indices:
            ranks[index] = rank
            rank += 1
    return ranks


def parse_label_max(value: Any) -> int | None:
    if is_na_value(value):
        return None
    parsed = parse_tsv_int(value)
    return parsed if parsed > 0 else None


def read_pi_windows(job_id: str) -> list[WindowResult]:
    path = outputs_dir(job_id) / "pi_windows.tsv"
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    windows: list[WindowResult] = []
    for index, row in enumerate(rows, start=1):
        start = parse_tsv_int(row.get("start"))
        end = parse_tsv_int(row.get("end"))
        midpoint = parse_tsv_int(row.get("midpoint"))
        valid_sites = parse_tsv_int(row.get("valid_sites"))
        if start <= 0 or end <= 0 or midpoint <= 0:
            continue
        windows.append(
            WindowResult(
                window_id=row.get("window_id") or f"W{index:04d}",
                start=start,
                end=end,
                midpoint=midpoint,
                pi=parse_tsv_float(row.get("pi")),
                valid_sites=valid_sites,
                region=row.get("region", ""),
                label_name=row.get("label_name", ""),
                is_hotspot=row.get("is_hotspot", "").strip().lower() == "yes",
            )
        )
    return windows


def row_spans(row: dict[str, Any], genome_length: int) -> list[tuple[int, int]]:
    start = parse_tsv_int(row.get("start")) - 1
    end = parse_tsv_int(row.get("end"))
    if start < 0 or end <= 0:
        return []
    if end > start:
        return [(start, end)]
    return [(start, genome_length), (0, end)]


def window_spans(window: WindowResult, genome_length: int) -> list[tuple[int, int]]:
    if window.end >= window.start:
        return [(window.start - 1, window.end)]
    return [(window.start - 1, genome_length), (0, window.end)]


def row_window_overlap(row: dict[str, Any], window: WindowResult, genome_length: int) -> int:
    overlap = 0
    for row_start, row_end in row_spans(row, genome_length):
        for window_start, window_end in window_spans(window, genome_length):
            left = max(row_start, window_start)
            right = min(row_end, window_end)
            if right > left:
                overlap += right - left
    return overlap


def row_midpoint(row: dict[str, Any], genome_length: int) -> int:
    spans = row_spans(row, genome_length)
    if not spans:
        return 0
    start = parse_tsv_int(row.get("start"))
    length = sum(end - start_ for start_, end in spans)
    return ((start - 1) + (length // 2)) % genome_length + 1


def best_row_index_for_window(
    rows: list[dict[str, str]],
    window: WindowResult,
    genome_length: int,
    selected_feature_ids: set[str],
) -> int | None:
    exact_matches = [
        (index, row)
        for index, row in enumerate(rows)
        if row.get("label_name") == window.label_name
    ]
    candidates = exact_matches or [
        (index, row)
        for index, row in enumerate(rows)
        if row_window_overlap(row, window, genome_length) > 0
    ]
    if not candidates:
        return None
    index, _row = min(
        candidates,
        key=lambda item: (
            item[1].get("feature_id", "") in selected_feature_ids,
            -row_window_overlap(item[1], window, genome_length),
            abs(row_midpoint(item[1], genome_length) - window.midpoint),
            parse_tsv_int(item[1].get("start")),
            parse_tsv_int(item[1].get("end")),
            item[1].get("feature_id", ""),
        ),
    )
    return index


def plot_label_hotspot_ranks(
    job_id: str,
    rows: list[dict[str, str]],
    params: dict[str, Any] | None,
) -> dict[int, int]:
    windows = read_pi_windows(job_id)
    if not windows:
        return original_hotspot_ranks(rows, params)

    genome_length = max((window.end for window in windows), default=0)
    if genome_length <= 0:
        return original_hotspot_ranks(rows, params)

    label_mode = str((params or {}).get("label_mode", "peak-only"))
    if label_mode not in {"peak-only", "all", "none"}:
        label_mode = "peak-only"
    label_result = AnalysisResult(
        reference_name="",
        genome_length=genome_length,
        sample_count=0,
        regions=[],
        position_pi=[],
        windows=windows,
        features=[],
    )
    labeled_windows = select_label_windows(
        label_result,
        label_mode=label_mode,
        label_max=parse_label_max((params or {}).get("label_max")),
        label_min_distance_bp=parse_tsv_int((params or {}).get("label_min_distance", 0)),
    )

    ranks: dict[int, int] = {}
    selected_feature_ids: set[str] = set()
    rank = 1
    for window in labeled_windows:
        index = best_row_index_for_window(rows, window, genome_length, selected_feature_ids)
        if index is None:
            continue
        feature_id = rows[index].get("feature_id", "")
        if feature_id in selected_feature_ids:
            continue
        ranks[index] = rank
        selected_feature_ids.add(feature_id)
        rank += 1
    return ranks or original_hotspot_ranks(rows, params)


def candidate_marker_sort_key(row: dict[str, str | bool | int]) -> tuple[float, ...]:
    hotspot_rank = row.get("candidate_hotspot_rank")
    if isinstance(hotspot_rank, int):
        return (
            0.0,
            float(parse_tsv_int(row.get("start"))),
            float(parse_tsv_int(row.get("end"))),
        )

    score = parse_tsv_float(row.get("markerseek_score"))
    score_key = -score if score is not None else math.inf
    return (
        1.0,
        score_key,
        float(parse_tsv_int(row.get("start"))),
        float(parse_tsv_int(row.get("end"))),
    )


def read_marker_feature_rows(job_id: str, params: dict[str, Any] | None = None) -> list[dict[str, str | bool | int]]:
    path = outputs_dir(job_id) / "candidate_marker_features.tsv"
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    hotspot_ranks = plot_label_hotspot_ranks(job_id, rows, params)
    visible_rows: list[dict[str, str | bool | int]] = []
    for index, row in enumerate(rows):
        row["has_detail"] = (feature_payload_dir(job_id) / f"{feature_id_safe(row.get('feature_id', ''))}.json").exists()
        row["has_missing_replicate_metrics"] = any(is_na_value(row.get(metric)) for metric in INTRA_DEPENDENT_METRICS)
        if index in hotspot_ranks:
            row["candidate_hotspot_rank"] = hotspot_ranks[index]
        if row["has_detail"]:
            visible_rows.append(row)
    return sorted(visible_rows, key=candidate_marker_sort_key)


def load_feature_payload(job_id: str, feature_id: str) -> tuple[dict[str, Any], dict[str, Any]]:
    row = get_job(job_id)
    if row is None or row["status"] != "succeeded":
        raise HTTPException(status_code=404, detail="Results are not available.")
    path = feature_payload_dir(job_id) / f"{feature_id_safe(feature_id)}.json"
    if not path.exists():
        raise HTTPException(status_code=404, detail="Feature detail is not available.")
    payload = json.loads(path.read_text(encoding="utf-8"))
    return row_to_job(row), payload


def feature_fasta(payload: dict[str, Any]) -> str:
    lines: list[str] = []
    for sample_name, sequence in payload.get("aligned_block", {}).items():
        lines.append(f">{sample_name}")
        lines.extend(sequence[index : index + 80] for index in range(0, len(sequence), 80))
    return "\n".join(lines) + ("\n" if lines else "")


def feature_csv(payload: dict[str, Any]) -> str:
    feature = payload["feature"]
    primers = payload.get("primers", [])
    output = StringIO()
    writer = csv.writer(output)
    writer.writerow(MARKER_FEATURE_COLUMNS + PRIMER_COLUMNS)
    feature_values = [_feature_csv_value(feature.get(column)) for column in MARKER_FEATURE_COLUMNS]
    if primers:
        for primer in primers:
            writer.writerow(feature_values + [_feature_csv_value(primer.get(column)) for column in PRIMER_COLUMNS])
    else:
        writer.writerow(feature_values + [""] * len(PRIMER_COLUMNS))
    return output.getvalue()


def _feature_csv_value(value: Any) -> str | int:
    if isinstance(value, float) or value is None:
        return format_float(value)
    return value


def input_file_metrics(job_id: str, input_files: list[str]) -> dict[str, Any]:
    total_bytes = 0
    found_files = 0
    for name in input_files:
        path = inputs_dir(job_id) / safe_filename(name)
        if not path.exists():
            continue
        found_files += 1
        total_bytes += path.stat().st_size
    return {
        "file_count": len(input_files),
        "found_files": found_files,
        "total_bytes": total_bytes,
        "total_mb": total_bytes / (1024 * 1024),
    }


def format_data_size(byte_count: int) -> str:
    if byte_count >= 1024 * 1024:
        return f"{byte_count / (1024 * 1024):.1f} MB"
    if byte_count >= 1024:
        return f"{byte_count / 1024:.1f} KB"
    return f"{byte_count} B"


def runtime_work_units(file_count: int, total_bytes: int) -> float:
    count = max(file_count, 1)
    measured_mb = total_bytes / (1024 * 1024)
    estimated_mb = max(measured_mb, count * 0.16)
    return max(1.0, estimated_mb * math.sqrt(count) + count * 0.2)


def seconds_between(start: str | None, end: str | None) -> float | None:
    if not start or not end:
        return None
    try:
        seconds = (parse_iso(end) - parse_iso(start)).total_seconds()
    except ValueError:
        return None
    return seconds if seconds > 0 else None


def recent_runtime_rates(*, exclude_job_id: str | None = None, limit: int = 30) -> list[float]:
    rows = db_fetchall(
        """
        SELECT job_id, started_at, completed_at, input_files_json
        FROM jobs
        WHERE status = 'succeeded' AND started_at IS NOT NULL AND completed_at IS NOT NULL
        ORDER BY completed_at DESC
        LIMIT ?
        """,
        (limit,),
    )
    rates: list[float] = []
    for row in rows:
        if row["job_id"] == exclude_job_id:
            continue
        duration = seconds_between(row["started_at"], row["completed_at"])
        if duration is None:
            continue
        try:
            input_files = json.loads(row["input_files_json"])
        except json.JSONDecodeError:
            continue
        metrics = input_file_metrics(row["job_id"], input_files)
        units = runtime_work_units(metrics["file_count"], metrics["total_bytes"])
        rates.append(duration / units)
    return rates


def server_load_info() -> dict[str, Any]:
    cores = os.cpu_count() or 1
    try:
        load_1m = os.getloadavg()[0]
    except (AttributeError, OSError):
        return {"label": "unknown", "factor": 1.0, "text": f"{cores} CPU cores"}
    ratio = load_1m / cores
    if ratio >= 1.25:
        label = "high"
        factor = 1.35
    elif ratio >= 0.85:
        label = "moderate"
        factor = 1.15
    else:
        label = "normal"
        factor = 1.0
    return {"label": label, "factor": factor, "text": f"{cores} CPU cores, {label} load"}


def estimate_compute_seconds(
    metrics: dict[str, Any],
    params: dict[str, Any] | None,
    rates: list[float],
    load_factor: float,
) -> float:
    file_count = metrics["file_count"]
    total_bytes = metrics["total_bytes"]
    total_mb = total_bytes / (1024 * 1024)
    units = runtime_work_units(file_count, total_bytes)
    fallback = 25 + (file_count ** 1.18 * 1.8) + (total_mb * 22) + (units * 9)
    if params and params.get("no_similarity_plot"):
        fallback *= 0.85
    if rates:
        estimate = (20 + units * median(rates)) * 0.75 + fallback * 0.25
    else:
        estimate = fallback
    return max(20.0, estimate * load_factor)


def format_duration(seconds: float | None) -> str:
    if seconds is None:
        return "not available"
    seconds = max(0, seconds)
    if seconds < 45:
        return "under 1 min"
    if seconds < 90:
        return "about 1 min"
    minutes = seconds / 60
    if minutes < 60:
        rounded = max(1, round(minutes))
        return f"about {rounded} min"
    hours = minutes / 60
    if hours < 24:
        return f"about {hours:.1f} h"
    return f"about {hours / 24:.1f} days"


def format_duration_range(seconds: float | None) -> str:
    if seconds is None:
        return "not available"
    seconds = max(0, seconds)
    if seconds < 45:
        return "under 1 min"
    lower = max(1, math.floor(seconds * 0.75 / 60))
    upper = max(lower, math.ceil(seconds * 1.35 / 60))
    if upper < 60:
        return f"about {lower}-{upper} min" if lower != upper else f"about {upper} min"
    lower_h = max(1, math.floor(lower / 60))
    upper_h = max(lower_h, math.ceil(upper / 60))
    return f"about {lower_h}-{upper_h} h" if lower_h != upper_h else f"about {upper_h} h"


def active_job_rows() -> list[sqlite3.Row]:
    return db_fetchall(
        """
        SELECT job_id, status, created_at, started_at, input_files_json, params_json
        FROM jobs
        WHERE status IN ('queued', 'running')
        ORDER BY created_at, job_id
        """
    )


def estimate_job_runtime(job: dict[str, Any]) -> dict[str, Any]:
    metrics = input_file_metrics(job["job_id"], job["input_files"])
    rates = recent_runtime_rates(exclude_job_id=job["job_id"])
    load = server_load_info()
    compute_seconds = estimate_compute_seconds(metrics, job["params"], rates, load["factor"])
    active_rows = active_job_rows()
    running_count = sum(1 for row in active_rows if row["status"] == "running")
    queued_count = sum(1 for row in active_rows if row["status"] == "queued")

    wait_seconds = 0.0
    jobs_ahead = 0
    for row in active_rows:
        if row["job_id"] == job["job_id"]:
            break
        jobs_ahead += 1
        try:
            row_files = json.loads(row["input_files_json"])
            row_params = json.loads(row["params_json"])
        except json.JSONDecodeError:
            row_files = []
            row_params = {}
        row_metrics = input_file_metrics(row["job_id"], row_files)
        row_compute = estimate_compute_seconds(row_metrics, row_params, rates, load["factor"])
        if row["status"] == "running" and row["started_at"]:
            try:
                elapsed = (utcnow() - parse_iso(row["started_at"])).total_seconds()
            except ValueError:
                elapsed = 0.0
            wait_seconds += max(30.0, row_compute - elapsed)
        else:
            wait_seconds += row_compute

    elapsed_seconds = None
    remaining_seconds = None
    primary_label = "Estimated runtime"
    primary_seconds = compute_seconds
    if job["status"] == "succeeded":
        actual = seconds_between(job["started_at"], job["completed_at"])
        if actual is not None:
            primary_label = "Actual runtime"
            primary_seconds = actual
    elif job["status"] == "running" and job["started_at"]:
        try:
            elapsed_seconds = max(0.0, (utcnow() - parse_iso(job["started_at"])).total_seconds())
        except ValueError:
            elapsed_seconds = None
        if elapsed_seconds is None:
            remaining_seconds = compute_seconds
        else:
            remaining_seconds = max(30.0, compute_seconds - elapsed_seconds)
        primary_label = "Estimated remaining time"
        primary_seconds = remaining_seconds
    elif job["status"] == "queued":
        remaining_seconds = wait_seconds + compute_seconds
        primary_label = "Estimated time to completion"
        primary_seconds = remaining_seconds

    if rates:
        basis = f"Calibrated from {len(rates)} recent completed jobs."
    else:
        basis = "Estimated from uploaded data size and current server load."

    return {
        "primary_label": primary_label,
        "primary_text": format_duration(primary_seconds) if primary_label == "Actual runtime" else format_duration_range(primary_seconds),
        "compute_text": format_duration_range(compute_seconds),
        "wait_text": format_duration_range(wait_seconds) if wait_seconds > 0 else "no queue wait",
        "elapsed_text": format_duration(elapsed_seconds) if elapsed_seconds is not None else None,
        "data_text": f"{metrics['file_count']} files, {format_data_size(metrics['total_bytes'])} uploaded",
        "server_text": f"{running_count} running, {queued_count} queued; {load['text']}",
        "queue_text": f"{jobs_ahead} job{'s' if jobs_ahead != 1 else ''} ahead",
        "basis": basis,
    }


def row_to_job(row: sqlite3.Row) -> dict[str, Any]:
    params = json.loads(row["params_json"])
    input_files = json.loads(row["input_files_json"])
    return {
        "job_id": row["job_id"],
        "status": row["status"],
        "project_name": row["project_name"],
        "created_at": row["created_at"],
        "updated_at": row["updated_at"],
        "started_at": row["started_at"],
        "completed_at": row["completed_at"],
        "expires_at": row["expires_at"],
        "params": params,
        "input_files": input_files,
        "reference_file": row["reference_file"],
        "error": row["error"],
        "permanent": bool(row["permanent"]),
    }


def append_log(job_id: str, message: str) -> None:
    path = outputs_dir(job_id) / "run.log"
    path.parent.mkdir(parents=True, exist_ok=True)
    timestamp = isoformat(utcnow())
    with path.open("a", encoding="utf-8") as handle:
        handle.write(f"[{timestamp}] {message}\n")


async def save_uploads(job_id: str, files: list[UploadFile]) -> list[str]:
    destination = inputs_dir(job_id)
    destination.mkdir(parents=True, exist_ok=True)
    saved_names: list[str] = []
    seen_names: set[str] = set()
    total_size = 0

    for upload in files:
        name = safe_filename(upload.filename or "")
        suffix = Path(name).suffix.lower()
        if suffix not in SUPPORTED_UPLOAD_SUFFIXES:
            raise ValueError(f"Unsupported file type for {upload.filename or name}. Use .gb, .gbk, or .genbank.")
        if name in seen_names:
            raise ValueError(f"Duplicate uploaded filename after sanitizing: {name}.")
        seen_names.add(name)
        saved_names.append(name)

        target = destination / name
        with target.open("wb") as handle:
            while True:
                chunk = await upload.read(1024 * 1024)
                if not chunk:
                    break
                total_size += len(chunk)
                if total_size > settings.max_upload_bytes:
                    raise ValueError(f"Total upload size exceeds {settings.max_upload_bytes // (1024 * 1024)} MB.")
                handle.write(chunk)
        await upload.close()

    return sorted(saved_names)


def recover_unfinished_jobs() -> None:
    rows = db_fetchall("SELECT job_id, status FROM jobs WHERE status IN ('queued', 'running')")
    for row in rows:
        if row["status"] == "running":
            update_job(
                row["job_id"],
                status="failed",
                completed_at=isoformat(utcnow()),
                error="Server restarted while this job was running. Please submit the analysis again.",
            )
        else:
            JOB_QUEUE.put(row["job_id"])


def cleanup_expired_jobs() -> None:
    now = utcnow()
    rows = db_fetchall("SELECT job_id, status, expires_at, permanent FROM jobs WHERE status != 'expired'")
    for row in rows:
        if row["permanent"]:
            continue
        try:
            expired = parse_iso(row["expires_at"]) <= now
        except ValueError:
            expired = False
        if not expired:
            continue
        shutil.rmtree(job_dir(row["job_id"]), ignore_errors=True)
        update_job(row["job_id"], status="expired", error="This job has expired and its files were removed.")


def run_job(job_id: str) -> None:
    row = get_job(job_id)
    if row is None:
        return
    job = row_to_job(row)
    params = job["params"]
    output_dir = outputs_dir(job_id)
    output_dir.mkdir(parents=True, exist_ok=True)
    update_job(job_id, status="running", started_at=isoformat(utcnow()), error=None)
    append_log(job_id, "Job started.")

    try:
        input_paths = [inputs_dir(job_id) / name for name in job["input_files"]]
        reference = inputs_dir(job_id) / job["reference_file"] if job["reference_file"] else None
        append_log(job_id, f"Input files: {len(input_paths)}.")
        append_log(job_id, f"Reference: {reference.name if reference else 'first sorted input'}")
        append_log(job_id, f"MAFFT threads: {params['mafft_threads']}.")

        result = run_analysis(
            input_paths,
            reference=reference,
            window=params["window"],
            step=params["step"],
            hotspot_mode=params["hotspot_mode"],
            hotspot_value=params["hotspot_value"],
            mafft_bin=settings.mafft_bin,
            mafft_threads=params["mafft_threads"],
            label_mode=params["label_mode"],
            label_max=params["label_max"],
            label_min_distance_bp=params["label_min_distance"],
            primer_design=params.get("primer_design", False),
            primer_settings=primer_settings_from_params(params),
            in_silico_settings=in_silico_settings_from_params(params),
        )
        write_analysis_outputs(
            result,
            output_dir,
            hotspot_mode=params["hotspot_mode"],
            hotspot_value=params["hotspot_value"],
            label_mode=params["label_mode"],
            label_max=params["label_max"],
            label_min_distance_bp=params["label_min_distance"],
            similarity_window=params["similarity_window"],
            similarity_step=params["similarity_step"],
            similarity_floor=params["similarity_floor"],
            include_similarity_plot=not params["no_similarity_plot"],
            primer_design=params.get("primer_design", False),
            mafft_bin=settings.mafft_bin,
        )
        summary = {
            "job_id": job_id,
            "project_name": job["project_name"],
            "status": "succeeded",
            "created_at": job["created_at"],
            "completed_at": isoformat(utcnow()),
            "expires_at": job["expires_at"],
            "params": params,
            "input_files": job["input_files"],
            "reference_file": job["reference_file"],
            "reference_name": result.reference_name,
            "sample_count": result.sample_count,
            "genome_length": result.genome_length,
            "window_count": len(result.windows),
            "hotspot_window_count": sum(1 for window in result.windows if window.is_hotspot),
            "feature_count": len(result.features),
            "primer_pair_count": len(result.primers),
        }
        write_json(output_dir / "job.json", summary)
        append_log(job_id, "Analysis completed successfully.")
        create_results_zip(output_dir, output_dir / "results.zip")
        update_job(job_id, status="succeeded", completed_at=summary["completed_at"], error=None)
    except (MarkerSeekError, Exception) as exc:
        message = str(exc) or exc.__class__.__name__
        append_log(job_id, f"Job failed: {message}")
        append_log(job_id, traceback.format_exc())
        update_job(job_id, status="failed", completed_at=isoformat(utcnow()), error=message)


def worker_loop() -> None:
    while True:
        try:
            job_id = JOB_QUEUE.get(timeout=60)
        except queue.Empty:
            cleanup_expired_jobs()
            continue
        try:
            run_job(job_id)
        finally:
            JOB_QUEUE.task_done()


def startup_services() -> None:
    global WORKER_STARTED
    init_db()
    init_database_catalog()
    cleanup_expired_jobs()
    recover_unfinished_jobs()
    if not WORKER_STARTED:
        thread = threading.Thread(target=worker_loop, name="markerseek-worker", daemon=True)
        thread.start()
        WORKER_STARTED = True


def init_database_catalog() -> None:
    db_path = settings.database_db_path
    try:
        db_catalog.init_database_db(db_path)
    except OSError:
        pass


@app.get("/", include_in_schema=False)
def root() -> RedirectResponse:
    return RedirectResponse(url="/markerseek", status_code=303)


@app.get("/markerseek", response_class=HTMLResponse)
def home(request: Request) -> HTMLResponse:
    cleanup_expired_jobs()
    return render_home(request)


@app.get("/MarkerSeek", include_in_schema=False)
def legacy_capitalized_home() -> RedirectResponse:
    return RedirectResponse(url="/markerseek", status_code=301)


@app.get("/analyzer/home", include_in_schema=False)
def legacy_home() -> RedirectResponse:
    return RedirectResponse(url="/markerseek", status_code=301)


@app.post("/analyzer/jobs", response_model=None)
async def submit_job(
    request: Request,
    files: list[UploadFile] = File(...),
    project_name: str = Form("MarkerSeek analysis"),
    reference_file: str = Form(""),
    window: str = Form("600"),
    step: str = Form("200"),
    hotspot_mode: str = Form("top-percent"),
    hotspot_value: str = Form("3.0"),
    label_mode: str = Form("peak-only"),
    label_max: str = Form(""),
    label_min_distance: str = Form("0"),
    similarity_window: str = Form("200"),
    similarity_step: str = Form("60"),
    similarity_floor: str = Form("0.5"),
    no_similarity_plot: str | None = Form(None),
    skip_primer_design: str | None = Form(None),
    primer_tm_min: str = Form(str(DEFAULT_PRIMER_TM_MIN)),
    primer_tm_opt: str = Form(str(DEFAULT_PRIMER_TM_OPT)),
    primer_tm_max: str = Form(str(DEFAULT_PRIMER_TM_MAX)),
    primer_len_min: str = Form(str(DEFAULT_PRIMER_LEN_MIN)),
    primer_len_opt: str = Form("20"),
    primer_len_max: str = Form(str(DEFAULT_PRIMER_LEN_MAX)),
    primer_amplicon_min: str = Form(str(DEFAULT_PRIMER_AMPLICON_MIN)),
    primer_amplicon_max: str = Form(str(DEFAULT_PRIMER_AMPLICON_MAX)),
    primer_mismatch: str = Form(str(DEFAULT_PRIMER_MISMATCH)),
    primer_anchor_bp: str = Form(str(DEFAULT_PRIMER_ANCHOR_BP)),
    primer_num_return: str = Form(str(DEFAULT_PRIMER_NUM_RETURN)),
    primer_flank_max_len: str = Form(str(DEFAULT_PRIMER_FLANK_MAX_LEN)),
    mafft_threads: str = Form(str(DEFAULT_MAFFT_THREADS)),
) -> Response:
    form_values = {
        "project_name": project_name,
        "window": window,
        "step": step,
        "hotspot_mode": hotspot_mode,
        "hotspot_value": hotspot_value,
        "label_mode": label_mode,
        "label_max": label_max,
        "label_min_distance": label_min_distance,
        "similarity_window": similarity_window,
        "similarity_step": similarity_step,
        "similarity_floor": similarity_floor,
        "no_similarity_plot": no_similarity_plot or "",
        "skip_primer_design": skip_primer_design or "",
        "primer_tm_min": primer_tm_min,
        "primer_tm_opt": primer_tm_opt,
        "primer_tm_max": primer_tm_max,
        "primer_len_min": primer_len_min,
        "primer_len_opt": primer_len_opt,
        "primer_len_max": primer_len_max,
        "primer_amplicon_min": primer_amplicon_min,
        "primer_amplicon_max": primer_amplicon_max,
        "primer_mismatch": primer_mismatch,
        "primer_anchor_bp": primer_anchor_bp,
        "primer_num_return": primer_num_return,
        "primer_flank_max_len": primer_flank_max_len,
        "mafft_threads": mafft_threads,
    }
    if len(files) < MIN_FILES or len(files) > MAX_FILES:
        return render_home(
            request,
            error=f"Upload between {MIN_FILES} and {MAX_FILES} GenBank files.",
            form_values=form_values,
        )

    try:
        params = parse_job_params(form_values)
        job_id = generate_job_id()
        saved_files = await save_uploads(job_id, files)
        selected_reference = safe_filename(reference_file) if reference_file.strip() else None
        if selected_reference and selected_reference not in saved_files:
            raise ValueError("Selected reference file is not one of the uploaded files.")
        clean_project_name = project_name.strip() or "MarkerSeek analysis"
        create_job_record(
            job_id=job_id,
            project_name=clean_project_name,
            params=params,
            input_files=saved_files,
            reference_file=selected_reference,
        )
        JOB_QUEUE.put(job_id)
    except ValueError as exc:
        if "job_id" in locals():
            shutil.rmtree(job_dir(job_id), ignore_errors=True)
        if request.headers.get("x-requested-with") == "XMLHttpRequest":
            return JSONResponse({"error": str(exc)}, status_code=400)
        return render_home(request, error=str(exc), form_values=form_values)

    return templates.TemplateResponse(
        request,
        "submitted.html",
        {
            "job_id": job_id,
            "project_name": clean_project_name,
            "retention_days": settings.retention_days,
        },
    )


@app.get("/analyzer/view", response_class=HTMLResponse)
def view_form(request: Request, error: str | None = None) -> HTMLResponse:
    if error == "missing":
        error = "Please enter a job ID."
    return templates.TemplateResponse(request, "view.html", {"error": error})


@app.post("/analyzer/view")
async def view_submit(projectid: str = Form("")) -> RedirectResponse:
    job_id = projectid.strip()
    if not job_id:
        return RedirectResponse(url="/analyzer/view?error=missing", status_code=303)
    return RedirectResponse(url=f"/analyzer/results/{job_id}", status_code=303)


@app.get("/analyzer/results/{job_id}", response_class=HTMLResponse)
def result_page(request: Request, job_id: str) -> HTMLResponse:
    cleanup_expired_jobs()
    row = get_job(job_id)
    if row is None:
        return templates.TemplateResponse(
            request,
            "view.html",
            {"error": f"No job was found for ID {job_id}."},
            status_code=404,
        )
    job = row_to_job(row)
    manifest_path = outputs_dir(job_id) / "job.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8")) if manifest_path.exists() else None
    available_files = [
        name
        for name in DOWNLOADABLE_FILES
        if (outputs_dir(job_id) / name).exists()
    ]
    feature_rows = read_marker_feature_rows(job_id, job["params"]) if job["status"] == "succeeded" else []
    show_replicate_hint = any(row.get("has_missing_replicate_metrics") for row in feature_rows)
    return templates.TemplateResponse(
        request,
        "result.html",
        {
            "job": job,
            "display_params": display_job_params(job["params"]),
            "runtime_estimate": estimate_job_runtime(job),
            "manifest": manifest,
            "available_files": sorted(available_files),
            "download_files": downloadable_result_files(sorted(available_files)),
            "feature_rows": feature_rows,
            "show_replicate_hint": show_replicate_hint,
            "na_diagnostic_explanation": NA_DIAGNOSTIC_EXPLANATION,
            "candidate_marker_explanation": CANDIDATE_MARKER_EXPLANATION,
            "primer_available_explanation": PRIMER_AVAILABLE_EXPLANATION,
            "refresh": job["status"] in {"queued", "running"},
        },
    )


def available_result_files(job_id: str) -> list[str]:
    return sorted(
        name
        for name in DOWNLOADABLE_FILES
        if (outputs_dir(job_id) / name).exists()
    )


def downloadable_result_files(available: list[str]) -> list[str]:
    return [name for name in available if not name.endswith(DOWNLOAD_LIST_HIDDEN_SUFFIXES)]


def demo_job_from_manifest(job_id: str, manifest: dict[str, Any] | None) -> dict[str, Any] | None:
    row = get_job(job_id) if settings.db_path.exists() else None
    if row is not None:
        return row_to_job(row)
    if manifest is None:
        return None
    return {
        "job_id": job_id,
        "status": manifest.get("status", "succeeded"),
        "project_name": manifest.get("project_name", "MarkerSeek demo"),
        "created_at": manifest.get("created_at", "permanent demo"),
        "updated_at": manifest.get("completed_at", manifest.get("created_at", "permanent demo")),
        "started_at": None,
        "completed_at": manifest.get("completed_at"),
        "expires_at": "permanent",
        "params": manifest.get("params", {}),
        "input_files": manifest.get("input_files", []),
        "reference_file": manifest.get("reference_file"),
        "error": None,
        "permanent": True,
    }


EXAMPLE_FEATURE_LABEL = "rpl32-trnL-UAG"


def _select_example_feature_row(feature_rows: list[dict[str, Any]]) -> dict[str, Any] | None:
    for row in feature_rows:
        if (row.get("label_name") or "").strip() == EXAMPLE_FEATURE_LABEL:
            return row
    return feature_rows[0] if feature_rows else None


def example_result_context(job_id: str) -> dict[str, Any]:
    manifest_path = outputs_dir(job_id) / "job.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8")) if manifest_path.exists() else None
    job = demo_job_from_manifest(job_id, manifest)
    available_files = available_result_files(job_id)
    feature_rows = read_marker_feature_rows(job_id, job.get("params", {}) if job else {})
    example_feature = None
    example_metadata_rows: list[tuple[str, str]] = []
    example_alignment_html = ""
    selected_row = _select_example_feature_row(feature_rows)
    if selected_row is not None:
        payload_path = feature_payload_dir(job_id) / f"{feature_id_safe(str(selected_row.get('feature_id', '')))}.json"
        if payload_path.exists():
            payload = json.loads(payload_path.read_text(encoding="utf-8"))
            example_feature = payload.get("feature")
            if example_feature:
                example_metadata_rows = feature_metadata_rows(example_feature)
            example_alignment_html = render_alignment_svg_block(payload.get("aligned_block", {}))
    return {
        "job": job,
        "manifest": manifest,
        "available_files": available_files,
        "download_files": downloadable_result_files(available_files),
        "feature_rows": feature_rows,
        "show_replicate_hint": any(row.get("has_missing_replicate_metrics") for row in feature_rows),
        "na_diagnostic_explanation": NA_DIAGNOSTIC_EXPLANATION,
        "candidate_marker_explanation": CANDIDATE_MARKER_EXPLANATION,
        "primer_available_explanation": PRIMER_AVAILABLE_EXPLANATION,
        "example_feature": example_feature,
        "example_metadata_rows": example_metadata_rows,
        "example_alignment_html": example_alignment_html,
    }


@app.get("/analyzer/example", response_class=HTMLResponse)
def example_page(request: Request) -> HTMLResponse:
    context = example_result_context(DEMO_JOB_ID)
    context["demo_job_id"] = DEMO_JOB_ID
    return templates.TemplateResponse(request, "example.html", context)


@app.get("/analyzer/results/{job_id}/feature/{feature_id:path}/data.json")
def feature_data_json(job_id: str, feature_id: str) -> JSONResponse:
    _job, payload = load_feature_payload(job_id, feature_id)
    return JSONResponse({key: payload.get(key) for key in FEATURE_DATA_KEYS})


@app.get("/analyzer/results/{job_id}/feature/{feature_id:path}/download/{kind}")
def feature_download(job_id: str, feature_id: str, kind: str) -> Response:
    if kind not in {"fasta", "csv"}:
        raise HTTPException(status_code=404, detail="Feature download is not available.")
    _job, payload = load_feature_payload(job_id, feature_id)
    safe_id = feature_id_safe(payload["feature"]["feature_id"])
    if kind == "fasta":
        content = feature_fasta(payload)
        media_type = "text/plain"
        filename = f"{safe_id}.fasta"
    else:
        content = feature_csv(payload)
        media_type = "text/csv"
        filename = f"{safe_id}.csv"
    return Response(
        content,
        media_type=media_type,
        headers={"Content-Disposition": f'attachment; filename="{filename}"'},
    )


@app.get("/analyzer/results/{job_id}/feature/{feature_id:path}", response_class=HTMLResponse)
def feature_detail_page(request: Request, job_id: str, feature_id: str) -> HTMLResponse:
    job, payload = load_feature_payload(job_id, feature_id)
    alignment_html = render_alignment_svg_block(payload.get("aligned_block", {}))
    base_path = f"/analyzer/results/{job_id}/feature/{quote(feature_id, safe='')}"
    return templates.TemplateResponse(
        request,
        "feature_detail.html",
        {
            "job": job,
            "feature": payload["feature"],
            "feature_metadata_rows": feature_metadata_rows(payload["feature"]),
            "primers": payload.get("primers", []),
            "alignment_html": alignment_html,
            "feature_base_url": f"{base_path}/",
            "back_url": f"/analyzer/results/{job_id}",
            "back_label": "Back to results",
            "download_fasta_url": f"{base_path}/download/fasta",
            "download_csv_url": f"{base_path}/download/csv",
            "show_downloads": True,
        },
    )


@app.get("/analyzer/results/{job_id}/download/all")
def download_all(job_id: str) -> FileResponse:
    row = get_job(job_id)
    if row is None or row["status"] != "succeeded":
        raise HTTPException(status_code=404, detail="Results are not available.")
    path = outputs_dir(job_id) / "results.zip"
    if not path.exists():
        raise HTTPException(status_code=404, detail="Results archive is missing.")
    return FileResponse(path, filename=f"{job_id}_results.zip", media_type="application/zip")


@app.get("/analyzer/results/{job_id}/download/{filename}")
def download_file(job_id: str, filename: str) -> FileResponse:
    row = get_job(job_id)
    if row is None or row["status"] != "succeeded":
        raise HTTPException(status_code=404, detail="Results are not available.")
    safe_name = safe_filename(filename)
    if safe_name not in DOWNLOADABLE_FILES:
        raise HTTPException(status_code=404, detail="File is not downloadable.")
    path = (outputs_dir(job_id) / safe_name).resolve()
    if outputs_dir(job_id).resolve() not in path.parents or not path.exists():
        raise HTTPException(status_code=404, detail="File was not found.")
    media_type = "application/octet-stream"
    if safe_name.endswith(".png"):
        media_type = "image/png"
    elif safe_name.endswith(".pdf"):
        media_type = "application/pdf"
    elif safe_name.endswith(".tsv"):
        media_type = "text/plain"
    elif safe_name.endswith(".zip"):
        media_type = "application/zip"
    return FileResponse(path, filename=safe_name, media_type=media_type)


@app.get("/analyzer/help", response_class=HTMLResponse)
def help_page(request: Request) -> HTMLResponse:
    return templates.TemplateResponse(
        request,
        "help.html",
        {
            "retention_days": settings.retention_days,
            "max_files": MAX_FILES,
            "max_upload_mb": settings.max_upload_bytes // (1024 * 1024),
            "demo_job_id": DEMO_JOB_ID,
            "na_diagnostic_explanation": NA_DIAGNOSTIC_EXPLANATION,
        },
    )


DATABASE_RESULT_FILES = set(RESULT_FILENAMES) | {"results.zip"}
DATABASE_FIGURE_FILES = ("pi_plot.png", "similarity_plot.png")
DATABASE_DOWNLOAD_FILES = ("candidate_marker_features.tsv", "primers.tsv")
DATABASE_PRIMER_DISPLAY_LIMIT = 30
DATABASE_EXTRA_MARKER_LIMIT = 10
DATABASE_MARKER_DOWNLOAD_DROP = INTRA_DEPENDENT_METRICS


def _genus_row_to_dict(row: sqlite3.Row) -> dict[str, Any]:
    summary = json.loads(row["summary_json"]) if row["summary_json"] else {}
    return {
        "genus": row["genus"],
        "family": row["family"],
        "taxonomic_order": row["taxonomic_order"],
        "species_count": row["species_count"],
        "ref_dir": row["ref_dir"],
        "results_dir": row["results_dir"],
        "has_primers": bool(row["has_primers"]),
        "genome_length_min": row["genome_length_min"],
        "genome_length_max": row["genome_length_max"],
        "top_marker_label": row["top_marker_label"],
        "top_marker_score": row["top_marker_score"],
        "marker_count": int(row["marker_count"] or 0),
        "primer_count": int(row["primer_count"] or 0),
        "indexed_at": row["indexed_at"],
        "summary": summary,
        "has_results": bool(row["results_dir"]) and int(row["marker_count"] or 0) > 0,
    }


def _genus_or_404(genus: str) -> dict[str, Any]:
    row = db_catalog.get_genus(settings.database_db_path, genus)
    if row is None:
        raise HTTPException(status_code=404, detail=f"Genus '{genus}' is not indexed.")
    return _genus_row_to_dict(row)


def _read_database_marker_rows(results_dir: Path) -> list[dict[str, Any]]:
    path = results_dir / "candidate_marker_features.tsv"
    if not path.exists():
        return []
    payload_dir = results_dir / "feature_payload"
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    enriched: list[dict[str, Any]] = []
    for row in rows:
        feature_id = row.get("feature_id", "")
        row["has_detail"] = (payload_dir / f"{feature_id_safe(feature_id)}.json").exists()
        row["has_missing_replicate_metrics"] = any(
            is_na_value(row.get(metric)) for metric in INTRA_DEPENDENT_METRICS
        )
        enriched.append(row)
    return enriched


def _read_hotspot_labels(results_dir: Path) -> set[str]:
    """Return the distinct hotspot ``label_name`` values from ``pi_windows.tsv``."""

    path = results_dir / "pi_windows.tsv"
    if not path.exists():
        return set()
    labels: set[str] = set()
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if (row.get("is_hotspot") or "").strip().lower() == "yes":
                label = (row.get("label_name") or "").strip()
                if label:
                    labels.add(label)
    return labels


def _select_database_marker_rows(
    rows: list[dict[str, Any]], hotspot_labels: set[str]
) -> tuple[list[dict[str, Any]], int, int]:
    """Mirror the analyzer/example selection: keep rows whose feature_payload was
    persisted by the analysis (i.e. those marked is_hotspot during the run),
    show hotspot-label features in genomic order, then the score-promoted rest
    by MarkerSeek score. This keeps the database genus page identical to
    /analyzer/example and /analyzer/results/.
    """

    def _start(row: dict[str, Any]) -> int:
        try:
            return int(row.get("start") or 0)
        except (TypeError, ValueError):
            return 0

    def _score(row: dict[str, Any]) -> float:
        try:
            return float(row.get("markerseek_score") or 0)
        except (TypeError, ValueError):
            return 0.0

    visible = [r for r in rows if r.get("has_detail")]
    hotspot_rows = [r for r in visible if (r.get("label_name") or "").strip() in hotspot_labels]
    hotspot_rows.sort(key=_start)
    selected_ids = {r.get("feature_id") for r in hotspot_rows}
    remaining = [r for r in visible if r.get("feature_id") not in selected_ids]
    remaining.sort(key=_score, reverse=True)
    return hotspot_rows + remaining, len(hotspot_rows), len(remaining)


def _read_database_primer_rows(results_dir: Path) -> list[dict[str, Any]]:
    path = results_dir / "primers.tsv"
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _load_database_feature_payload(genus_record: dict[str, Any], feature_id: str) -> dict[str, Any]:
    results_dir = genus_record.get("results_dir")
    if not results_dir:
        raise HTTPException(status_code=404, detail="No analysis results for this genus.")
    payload_path = Path(results_dir) / "feature_payload" / f"{feature_id_safe(feature_id)}.json"
    if not payload_path.exists():
        raise HTTPException(status_code=404, detail="Feature detail is not available.")
    return json.loads(payload_path.read_text(encoding="utf-8"))


@app.get("/database", response_class=HTMLResponse)
def database_home(request: Request) -> HTMLResponse:
    db_path = settings.database_db_path
    stats = db_catalog.database_stats(db_path)
    facets = db_catalog.list_taxonomy_facets(db_path) if db_path.exists() else {"families": [], "orders": []}
    alphabet = db_catalog.alphabet_index(db_path)
    return templates.TemplateResponse(
        request,
        "database_home.html",
        {
            "stats": stats,
            "facets": facets,
            "alphabet": alphabet,
            "database_root": str(settings.database_root),
            "is_empty": stats["genus_count"] == 0,
        },
    )


@app.get("/database/browse", response_class=HTMLResponse)
def database_browse(
    request: Request,
    q: str = "",
    family: str = "",
    order: str = "",
    page: int = 1,
) -> HTMLResponse:
    db_path = settings.database_db_path
    rows, total = db_catalog.list_genera(
        db_path,
        query=q.strip() or None,
        family=family.strip() or None,
        taxonomic_order=order.strip() or None,
        page=page,
    )
    page_size = 50
    total_pages = max(1, (total + page_size - 1) // page_size)
    facets = db_catalog.list_taxonomy_facets(db_path) if db_path.exists() else {"families": [], "orders": []}
    return templates.TemplateResponse(
        request,
        "database_browse.html",
        {
            "rows": [_genus_row_to_dict(row) for row in rows],
            "total": total,
            "page": page,
            "total_pages": total_pages,
            "page_size": page_size,
            "q": q,
            "family": family,
            "order": order,
            "facets": facets,
        },
    )


@app.get("/database/api/search")
def database_search_api(q: str = "", limit: int = 10) -> JSONResponse:
    if not q.strip():
        return JSONResponse({"results": []})
    rows, _total = db_catalog.list_genera(
        settings.database_db_path,
        query=q.strip(),
        page=1,
        page_size=max(1, min(limit, 50)),
    )
    return JSONResponse(
        {
            "results": [
                {
                    "genus": row["genus"],
                    "species_count": row["species_count"],
                    "has_results": bool(row["results_dir"]) and int(row["marker_count"] or 0) > 0,
                }
                for row in rows
            ]
        }
    )


@app.get("/database/genus/{genus}", response_class=HTMLResponse)
def database_genus_page(request: Request, genus: str) -> HTMLResponse:
    record = _genus_or_404(genus)
    species: list[db_catalog.SpeciesEntry] = []
    if record["ref_dir"]:
        species = db_catalog.list_species_entries(Path(record["ref_dir"]))
    selected_markers: list[dict[str, Any]] = []
    primer_rows: list[dict[str, Any]] = []
    figure_files: list[str] = []
    available_files: list[str] = []
    marker_total = 0
    hotspot_count = 0
    extra_count = 0
    if record["results_dir"]:
        results_dir = Path(record["results_dir"])
        marker_rows = _read_database_marker_rows(results_dir)
        marker_total = len(marker_rows)
        hotspot_labels = _read_hotspot_labels(results_dir)
        selected_markers, hotspot_count, extra_count = _select_database_marker_rows(
            marker_rows, hotspot_labels
        )
        primer_rows = _read_database_primer_rows(results_dir)
        figure_files = [name for name in DATABASE_FIGURE_FILES if (results_dir / name).exists()]
        available_files = [name for name in DATABASE_DOWNLOAD_FILES if (results_dir / name).exists()]
    return templates.TemplateResponse(
        request,
        "database_genus.html",
        {
            "record": record,
            "species": species,
            "species_total": len(species),
            "marker_rows": selected_markers,
            "marker_total": marker_total,
            "hotspot_count": hotspot_count,
            "extra_count": extra_count,
            "primer_rows": primer_rows[:DATABASE_PRIMER_DISPLAY_LIMIT],
            "primer_total": len(primer_rows),
            "figure_files": figure_files,
            "available_files": available_files,
            "candidate_marker_explanation": CANDIDATE_MARKER_EXPLANATION,
            "primer_available_explanation": PRIMER_AVAILABLE_EXPLANATION,
            "na_diagnostic_explanation": NA_DIAGNOSTIC_EXPLANATION,
        },
    )


@app.get("/database/genus/{genus}/feature/{feature_id:path}/data.json")
def database_feature_data(genus: str, feature_id: str) -> JSONResponse:
    record = _genus_or_404(genus)
    payload = _load_database_feature_payload(record, feature_id)
    return JSONResponse({key: payload.get(key) for key in FEATURE_DATA_KEYS})


@app.get("/database/genus/{genus}/feature/{feature_id:path}/download/{kind}")
def database_feature_download(genus: str, feature_id: str, kind: str) -> Response:
    record = _genus_or_404(genus)
    if kind not in {"fasta", "csv"}:
        raise HTTPException(status_code=404, detail="Feature download is not available.")
    payload = _load_database_feature_payload(record, feature_id)
    safe_id = feature_id_safe(payload["feature"]["feature_id"])
    if kind == "fasta":
        content = feature_fasta(payload)
        media_type = "text/plain"
        filename = f"{genus}_{safe_id}.fasta"
    else:
        content = feature_csv(payload)
        media_type = "text/csv"
        filename = f"{genus}_{safe_id}.csv"
    return Response(
        content,
        media_type=media_type,
        headers={"Content-Disposition": f'attachment; filename="{filename}"'},
    )


@app.get("/database/genus/{genus}/feature/{feature_id:path}", response_class=HTMLResponse)
def database_feature_page(request: Request, genus: str, feature_id: str) -> HTMLResponse:
    record = _genus_or_404(genus)
    payload = _load_database_feature_payload(record, feature_id)
    alignment_html = render_alignment_svg_block(payload.get("aligned_block", {}))
    base_path = f"/database/genus/{genus}/feature/{quote(feature_id, safe='')}"
    fake_job = {"job_id": genus, "project_name": f"{genus} reference catalogue"}
    return templates.TemplateResponse(
        request,
        "feature_detail.html",
        {
            "job": fake_job,
            "feature": payload["feature"],
            "feature_metadata_rows": feature_metadata_rows(payload["feature"]),
            "primers": payload.get("primers", []),
            "alignment_html": alignment_html,
            "feature_base_url": f"{base_path}/",
            "back_url": f"/database/genus/{genus}",
            "back_label": f"Back to {genus}",
            "download_fasta_url": f"{base_path}/download/fasta",
            "download_csv_url": f"{base_path}/download/csv",
            "show_downloads": False,
        },
    )


@app.get("/database/genus/{genus}/download/{filename}")
def database_genus_download(genus: str, filename: str) -> Response:
    record = _genus_or_404(genus)
    results_dir = record.get("results_dir")
    if not results_dir:
        raise HTTPException(status_code=404, detail="No analysis results for this genus.")
    safe_name = safe_filename(filename)
    if safe_name not in DATABASE_RESULT_FILES:
        raise HTTPException(status_code=404, detail="File is not downloadable.")
    base = Path(results_dir).resolve()
    path = (base / safe_name).resolve()
    if base not in path.parents or not path.exists():
        raise HTTPException(status_code=404, detail="File was not found.")
    if safe_name == "candidate_marker_features.tsv":
        content = _filtered_candidate_marker_tsv(path)
        return Response(
            content,
            media_type="text/plain",
            headers={"Content-Disposition": f'attachment; filename="{genus}_{safe_name}"'},
        )
    media_type = "application/octet-stream"
    if safe_name.endswith(".png"):
        media_type = "image/png"
    elif safe_name.endswith(".pdf"):
        media_type = "application/pdf"
    elif safe_name.endswith(".tsv"):
        media_type = "text/plain"
    elif safe_name.endswith(".zip"):
        media_type = "application/zip"
    elif safe_name.endswith(".fasta"):
        media_type = "text/plain"
    return FileResponse(path, filename=f"{genus}_{safe_name}", media_type=media_type)


def _filtered_candidate_marker_tsv(path: Path) -> str:
    """Return the contents of candidate_marker_features.tsv with replicate-only metric columns dropped."""

    buffer = StringIO()
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            return ""
        keep_idx = [i for i, name in enumerate(header) if name not in DATABASE_MARKER_DOWNLOAD_DROP]
        writer = csv.writer(buffer, delimiter="\t", lineterminator="\n")
        writer.writerow([header[i] for i in keep_idx])
        for row in reader:
            writer.writerow([row[i] if i < len(row) else "" for i in keep_idx])
    return buffer.getvalue()


@app.get("/database/genus/{genus}/species/{filename}")
def database_species_download(genus: str, filename: str) -> FileResponse:
    record = _genus_or_404(genus)
    ref_dir = record.get("ref_dir")
    if not ref_dir:
        raise HTTPException(status_code=404, detail="No reference set for this genus.")
    safe_name = safe_filename(filename)
    if Path(safe_name).suffix.lower() not in db_catalog.GENBANK_SUFFIXES:
        raise HTTPException(status_code=404, detail="Only GenBank files are downloadable.")
    base = Path(ref_dir).resolve()
    path = (base / safe_name).resolve()
    if base not in path.parents or not path.exists():
        raise HTTPException(status_code=404, detail="File was not found.")
    return FileResponse(path, filename=safe_name, media_type="text/plain")


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="markerseek-web", description="Run the MarkerSeek web server.")
    parser.add_argument("--host", default=os.environ.get("MARKERSEEK_WEB_HOST", "127.0.0.1"))
    parser.add_argument("--port", type=int, default=int(os.environ.get("MARKERSEEK_WEB_PORT", "8000")))
    parser.add_argument("--reload", action="store_true", help="Enable uvicorn reload for development.")
    return parser


def main(argv: list[str] | None = None) -> None:
    args = build_arg_parser().parse_args(argv)
    uvicorn.run("markerseek.web.app:app", host=args.host, port=args.port, reload=args.reload, workers=1)


if __name__ == "__main__":
    main()
