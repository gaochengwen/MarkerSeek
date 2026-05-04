"""FastAPI web application for MarkerSeek."""

from __future__ import annotations

import argparse
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
from pathlib import Path
from secrets import randbelow
from statistics import median
from typing import Any

import uvicorn
from fastapi import FastAPI, File, Form, HTTPException, Request, UploadFile
from fastapi.responses import FileResponse, HTMLResponse, JSONResponse, RedirectResponse, Response
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates

from markerseek.analysis import MarkerSeekError, run_analysis
from markerseek.output import RESULT_FILENAMES, create_results_zip, write_analysis_outputs, write_json

PACKAGE_DIR = Path(__file__).resolve().parent
TEMPLATE_DIR = PACKAGE_DIR / "templates"
STATIC_DIR = PACKAGE_DIR / "static"

SUPPORTED_UPLOAD_SUFFIXES = {".gb", ".gbk", ".genbank"}
JOB_STATUS_ORDER = {"queued", "running", "succeeded", "failed", "expired"}
DOWNLOADABLE_FILES = set(RESULT_FILENAMES) | {"results.zip"}

MAX_FILES = 30
MIN_FILES = 2
DEFAULT_RETENTION_DAYS = 7
DEFAULT_MAX_UPLOAD_BYTES = 20 * 1024 * 1024
DEFAULT_MAFFT_THREADS = min(4, os.cpu_count() or 1)
SAFE_FILENAME_RE = re.compile(r"[^A-Za-z0-9_.-]+")


class WebSettings:
    def __init__(self) -> None:
        self.data_dir = Path(os.environ.get("MARKERSEEK_WEB_DATA", "markerseek_jobs")).resolve()
        self.retention_days = int(os.environ.get("MARKERSEEK_RETENTION_DAYS", str(DEFAULT_RETENTION_DAYS)))
        self.max_upload_bytes = int(os.environ.get("MARKERSEEK_MAX_UPLOAD_BYTES", str(DEFAULT_MAX_UPLOAD_BYTES)))
        self.mafft_bin = os.environ.get("MARKERSEEK_MAFFT_BIN", "mafft")

    @property
    def jobs_dir(self) -> Path:
        return self.data_dir / "jobs"

    @property
    def db_path(self) -> Path:
        return self.data_dir / "jobs.sqlite3"


settings = WebSettings()
templates = Jinja2Templates(directory=str(TEMPLATE_DIR))

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
                error TEXT
            )
            """
        )
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
) -> None:
    now = utcnow()
    expires_at = now + timedelta(days=settings.retention_days)
    db_execute(
        """
        INSERT INTO jobs (
            job_id, status, project_name, created_at, updated_at, expires_at,
            params_json, input_files_json, reference_file
        )
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            job_id,
            "queued",
            project_name,
            isoformat(now),
            isoformat(now),
            isoformat(expires_at),
            json.dumps(params, sort_keys=True),
            json.dumps(input_files),
            reference_file,
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
    hidden = {"label_max", "label_min_distance", "label_mode", "mafft_threads", "similarity_floor"}
    labels = {"window": "hotspot_window", "step": "hotspot_step"}
    priority = [
        "hotspot_mode",
        "hotspot_value",
        "window",
        "step",
        "similarity_window",
        "similarity_step",
        "no_similarity_plot",
    ]
    ordered_keys = [key for key in priority if key in params]
    ordered_keys.extend(key for key in params if key not in hidden and key not in ordered_keys)
    return [(labels.get(key, key), params[key]) for key in ordered_keys if key not in hidden]


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
    rows = db_fetchall("SELECT job_id, status, expires_at FROM jobs WHERE status != 'expired'")
    for row in rows:
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
    cleanup_expired_jobs()
    recover_unfinished_jobs()
    if not WORKER_STARTED:
        thread = threading.Thread(target=worker_loop, name="markerseek-worker", daemon=True)
        thread.start()
        WORKER_STARTED = True


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
    return templates.TemplateResponse(
        request,
        "result.html",
        {
            "job": job,
            "display_params": display_job_params(job["params"]),
            "runtime_estimate": estimate_job_runtime(job),
            "manifest": manifest,
            "available_files": sorted(available_files),
            "refresh": job["status"] in {"queued", "running"},
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
        },
    )


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
