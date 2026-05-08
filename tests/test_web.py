from __future__ import annotations

import sqlite3
import re
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from markerseek.web import app as web_app


def test_parse_job_params_validates_mafft_threads(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(web_app.os, "cpu_count", lambda: 4)
    form = web_app.default_form_values()
    form["mafft_threads"] = "5"

    with pytest.raises(ValueError, match="cannot exceed"):
        web_app.parse_job_params(form)


def test_default_similarity_values_match_web_defaults() -> None:
    values = web_app.default_form_values()

    assert values["similarity_window"] == "200"
    assert values["similarity_step"] == "60"


def test_primer_design_form_field_round_trip() -> None:
    form = web_app.default_form_values()

    enabled = web_app.parse_job_params(form)
    form["skip_primer_design"] = "on"
    disabled = web_app.parse_job_params(form)

    assert enabled["primer_design"] is True
    assert disabled["primer_design"] is False


def test_save_filename_strips_path_and_unsafe_characters() -> None:
    assert web_app.safe_filename("../../bad name.gb") == "bad_name.gb"
    assert web_app.safe_filename("") == "input.gb"


def test_generate_job_id_uses_compact_public_format(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    settings = web_app.WebSettings()
    settings.data_dir = tmp_path
    monkeypatch.setattr(web_app, "settings", settings)
    web_app.init_db()

    job_id = web_app.generate_job_id()

    assert re.fullmatch(r"MSK-\d{15}", job_id)


def test_download_file_rejects_non_whitelisted_output(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    settings = web_app.WebSettings()
    settings.data_dir = tmp_path
    monkeypatch.setattr(web_app, "settings", settings)
    web_app.init_db()

    job_id = "MSK-20260501-ABCDEF12"
    outdir = web_app.outputs_dir(job_id)
    outdir.mkdir(parents=True)
    (outdir / "secret.txt").write_text("private", encoding="utf-8")
    (outdir / "job.json").write_text("{}", encoding="utf-8")
    (outdir / "run.log").write_text("log", encoding="utf-8")
    now = web_app.isoformat(web_app.utcnow())
    with sqlite3.connect(settings.db_path) as conn:
        conn.execute(
            """
            INSERT INTO jobs (
                job_id, status, project_name, created_at, updated_at, expires_at,
                params_json, input_files_json
            )
            VALUES (?, 'succeeded', 'test', ?, ?, ?, '{}', '[]')
            """,
            (job_id, now, now, now),
        )
        conn.commit()

    client = TestClient(web_app.app)
    response = client.get(f"/analyzer/results/{job_id}/download/secret.txt")
    job_json_response = client.get(f"/analyzer/results/{job_id}/download/job.json")
    run_log_response = client.get(f"/analyzer/results/{job_id}/download/run.log")

    assert response.status_code == 404
    assert job_json_response.status_code == 404
    assert run_log_response.status_code == 404


def test_download_candidate_marker_features_tsv_whitelisted(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    settings = web_app.WebSettings()
    settings.data_dir = tmp_path
    monkeypatch.setattr(web_app, "settings", settings)
    web_app.init_db()

    job_id = "MSK-202605010000001"
    outdir = web_app.outputs_dir(job_id)
    outdir.mkdir(parents=True)
    (outdir / "candidate_marker_features.tsv").write_text("feature_id\n", encoding="utf-8")
    (outdir / "Marker_features.tsv").write_text("feature_id\n", encoding="utf-8")
    now = web_app.isoformat(web_app.utcnow())
    with sqlite3.connect(settings.db_path) as conn:
        conn.execute(
            """
            INSERT INTO jobs (
                job_id, status, project_name, created_at, updated_at, expires_at,
                params_json, input_files_json
            )
            VALUES (?, 'succeeded', 'test', ?, ?, ?, '{}', '[]')
            """,
            (job_id, now, now, now),
        )
        conn.commit()

    client = TestClient(web_app.app)
    marker_response = client.get(f"/analyzer/results/{job_id}/download/candidate_marker_features.tsv")
    old_response = client.get(f"/analyzer/results/{job_id}/download/Marker_features.tsv")

    assert marker_response.status_code == 200
    assert old_response.status_code == 404


def test_candidate_marker_rows_keep_hotspots_then_score_candidates(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    settings = web_app.WebSettings()
    settings.data_dir = tmp_path
    monkeypatch.setattr(web_app, "settings", settings)

    job_id = "MSK-202605010000002"
    outdir = web_app.outputs_dir(job_id)
    payload_dir = web_app.feature_payload_dir(job_id)
    payload_dir.mkdir(parents=True)
    (outdir / "candidate_marker_features.tsv").write_text(
        "\n".join(
            [
                "feature_id\tlabel_name\tstart\tend\tpi\tmarkerseek_score",
                "f1\tpi_rank_1\t10\t19\t0.900000\t10.000000",
                "f2\tscore_rank_2\t20\t29\t0.100000\t90.000000",
                "f3\tpi_rank_2\t5\t9\t0.800000\t50.000000",
                "f4\tscore_rank_3\t40\t49\t0.200000\t80.000000",
                "f5\thidden\t50\t59\t0.300000\t70.000000",
                "f6\tscore_rank_1\t60\t69\t0.400000\t95.000000",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    for feature_id in ("f1", "f2", "f3", "f4", "f6"):
        (payload_dir / f"{feature_id}.json").write_text("{}", encoding="utf-8")

    rows = web_app.read_marker_feature_rows(job_id, {"hotspot_mode": "top-n", "hotspot_value": 2})

    assert [row["feature_id"] for row in rows] == ["f3", "f1", "f6", "f2", "f4"]


def test_candidate_marker_rows_use_pi_plot_labels_before_score_candidates(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    settings = web_app.WebSettings()
    settings.data_dir = tmp_path
    monkeypatch.setattr(web_app, "settings", settings)

    job_id = "MSK-202605010000003"
    outdir = web_app.outputs_dir(job_id)
    payload_dir = web_app.feature_payload_dir(job_id)
    payload_dir.mkdir(parents=True)
    (outdir / "candidate_marker_features.tsv").write_text(
        "\n".join(
            [
                "feature_id\tlabel_name\tstart\tend\tpi\tmarkerseek_score",
                "f1\tycf1\t100\t159\t0.100000\t10.000000",
                "f2\tscore_rank_2\t200\t259\t0.900000\t90.000000",
                "f3\tndhF\t300\t359\t0.200000\t50.000000",
                "f4\tscore_rank_3\t400\t459\t0.800000\t80.000000",
                "f5\thidden\t500\t559\t0.700000\t70.000000",
                "f6\tscore_rank_1\t600\t659\t0.600000\t95.000000",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (outdir / "pi_windows.tsv").write_text(
        "\n".join(
            [
                "window_id\tstart\tend\tmidpoint\tpi\tvalid_sites\tregion\tlabel_name\tis_hotspot",
                "W0001\t105\t154\t130\t0.040000\t50\tSSC\tycf1\tyes",
                "W0002\t305\t354\t330\t0.050000\t50\tSSC\tndhF\tyes",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    for feature_id in ("f1", "f2", "f3", "f4", "f6"):
        (payload_dir / f"{feature_id}.json").write_text("{}", encoding="utf-8")

    rows = web_app.read_marker_feature_rows(
        job_id,
        {"hotspot_mode": "top-n", "hotspot_value": 2, "label_mode": "all"},
    )

    assert [row["feature_id"] for row in rows] == ["f1", "f3", "f6", "f2", "f4"]


def test_example_page_renders() -> None:
    client = TestClient(web_app.app)

    response = client.get("/analyzer/example")

    assert response.status_code == 200
    assert "Demo dataset" in response.text
    assert "MSK-EXAMPLE-DEMO" in response.text


def test_example_demo_job_route_returns_200(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    settings = web_app.WebSettings()
    settings.data_dir = tmp_path
    monkeypatch.setattr(web_app, "settings", settings)
    web_app.init_db()

    job_id = "MSK-EXAMPLE-DEMO"
    outdir = web_app.outputs_dir(job_id)
    outdir.mkdir(parents=True)
    (outdir / "candidate_marker_features.tsv").write_text("feature_id\n", encoding="utf-8")
    now = web_app.isoformat(web_app.utcnow())
    (outdir / "job.json").write_text(
        """{
  "job_id": "MSK-EXAMPLE-DEMO",
  "status": "succeeded",
  "permanent": true,
  "expires_at": null,
  "params": {},
  "input_files": [],
  "reference_file": "Salvia_chinensis.gb",
  "sample_count": 18,
  "genome_length": 151000,
  "window_count": 1,
  "hotspot_window_count": 1,
  "feature_count": 1,
  "primer_pair_count": 0
}""",
        encoding="utf-8",
    )
    with sqlite3.connect(settings.db_path) as conn:
        conn.execute(
            """
            INSERT INTO jobs (
                job_id, status, project_name, created_at, updated_at, completed_at,
                expires_at, params_json, input_files_json, reference_file, permanent
            )
            VALUES (?, 'succeeded', 'MarkerSeek demo', ?, ?, ?, 'permanent', '{}', '[]', 'Salvia_chinensis.gb', 1)
            """,
            (job_id, now, now, now),
        )
        conn.commit()

    client = TestClient(web_app.app)
    response = client.get(f"/analyzer/results/{job_id}")

    assert response.status_code == 200
    assert "MSK-EXAMPLE-DEMO" in response.text


def test_runtime_estimate_uses_uploaded_size_and_queue(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    settings = web_app.WebSettings()
    settings.data_dir = tmp_path
    monkeypatch.setattr(web_app, "settings", settings)
    monkeypatch.setattr(web_app.os, "cpu_count", lambda: 4)
    monkeypatch.setattr(web_app.os, "getloadavg", lambda: (0.0, 0.0, 0.0))
    web_app.init_db()

    params = web_app.parse_job_params(web_app.default_form_values())
    web_app.create_job_record(
        job_id="MSK-202605020000002",
        project_name="test",
        params=params,
        input_files=["a.gb", "b.gb"],
        reference_file=None,
    )
    input_dir = web_app.inputs_dir("MSK-202605020000002")
    input_dir.mkdir(parents=True)
    (input_dir / "a.gb").write_bytes(b"A" * 1024)
    (input_dir / "b.gb").write_bytes(b"C" * 1024)

    job = web_app.row_to_job(web_app.get_job("MSK-202605020000002"))
    estimate = web_app.estimate_job_runtime(job)

    assert estimate["primary_label"] == "Estimated time to completion"
    assert estimate["data_text"] == "2 files, 2.0 KB uploaded"
    assert estimate["queue_text"] == "0 jobs ahead"
    assert estimate["server_text"] == "0 running, 1 queued; 4 CPU cores, normal load"


def test_result_page_hides_internal_parameters(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    settings = web_app.WebSettings()
    settings.data_dir = tmp_path
    monkeypatch.setattr(web_app, "settings", settings)
    web_app.init_db()

    params = {
        "window": 600,
        "step": 200,
        "hotspot_mode": "top-percent",
        "hotspot_value": 3.0,
        "label_mode": "peak-only",
        "label_max": None,
        "label_min_distance": 0,
        "similarity_window": 200,
        "similarity_step": 60,
        "similarity_floor": 0.5,
        "no_similarity_plot": False,
        "mafft_threads": 4,
    }
    web_app.create_job_record(
        job_id="MSK-202605020000001",
        project_name="test",
        params=params,
        input_files=["a.gb", "b.gb"],
        reference_file=None,
    )

    client = TestClient(web_app.app)
    response = client.get("/analyzer/results/MSK-202605020000001")

    assert response.status_code == 200
    assert "Estimated Runtime" in response.text
    assert response.text.index("hotspot_value") < response.text.index("hotspot_window")
    assert response.text.index("hotspot_window") < response.text.index("hotspot_step")
    assert "<dt>window</dt>" not in response.text
    assert "<dt>step</dt>" not in response.text
    assert "similarity_window" in response.text
    assert "similarity_step" in response.text
    assert "label_max" not in response.text
    assert "label_min_distance" not in response.text
    assert "label_mode" not in response.text
    assert "mafft_threads" not in response.text
    assert "similarity_floor" not in response.text
