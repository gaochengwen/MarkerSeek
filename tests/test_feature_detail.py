from __future__ import annotations

import json
import sqlite3
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from markerseek.web import app as web_app


def test_feature_detail_route_renders_plotly_skeleton(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    job_id = _write_succeeded_feature_job(tmp_path, monkeypatch)
    client = TestClient(web_app.app)

    response = client.get(f"/analyzer/results/{job_id}/feature/IGS_001")

    assert response.status_code == 200
    assert "cdn.plot.ly" in response.text
    assert 'id="pi-curve"' in response.text
    assert 'id="haplotype-net"' in response.text
    assert "amplicon_size" in response.text
    assert "100-102 bp" in response.text
    assert "primer3_penalty" not in response.text


def test_feature_data_json_endpoint_keys(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    job_id = _write_succeeded_feature_job(tmp_path, monkeypatch)
    client = TestClient(web_app.app)

    response = client.get(f"/analyzer/results/{job_id}/feature/IGS_001/data.json")

    assert response.status_code == 200
    assert set(response.json()) == {
        "pi_curve",
        "snp_positions",
        "indel_positions",
        "haplotype_network",
        "species_pca",
        "primers",
    }
    assert "aligned_block" not in response.json()


def test_feature_detail_404_for_non_hotspot(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    job_id = _write_succeeded_feature_job(tmp_path, monkeypatch)
    client = TestClient(web_app.app)

    response = client.get(f"/analyzer/results/{job_id}/feature/IGS_999")

    assert response.status_code == 404


def test_feature_download_fasta_returns_aligned_block(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    job_id = _write_succeeded_feature_job(tmp_path, monkeypatch)
    client = TestClient(web_app.app)

    response = client.get(f"/analyzer/results/{job_id}/feature/IGS_001/download/fasta")

    assert response.status_code == 200
    assert "attachment" in response.headers["Content-Disposition"]
    assert response.text.count(">") == 2


def _write_succeeded_feature_job(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> str:
    settings = web_app.WebSettings()
    settings.data_dir = tmp_path
    monkeypatch.setattr(web_app, "settings", settings)
    web_app.init_db()

    job_id = "MSK-202605040000001"
    outdir = web_app.outputs_dir(job_id)
    payload_dir = outdir / "feature_payload"
    payload_dir.mkdir(parents=True)
    payload = _feature_payload()
    (payload_dir / "IGS_001.json").write_text(json.dumps(payload), encoding="utf-8")
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
    return job_id


def _feature_payload() -> dict:
    return {
        "feature": {
            "feature_id": "IGS_001",
            "feature_type": "IGS",
            "parent_gene": "geneA|geneB",
            "label_name": "geneA-geneB",
            "start": 1,
            "end": 4,
            "strand": 0,
            "length_bp": 4,
            "region": "LSC",
            "pi": 0.125,
            "variable_sites": 1,
            "indel_sites": 1,
            "conserved_left_bp": 50,
            "conserved_right_bp": 50,
            "primer_available": "yes",
            "species_resolution": 1.0,
            "unique_haplotype_count": 2,
            "species_specific_haplotype_ratio": 1.0,
            "interspecific_divergence": 0.1,
            "intraspecific_divergence": 0.0,
            "nearest_neighbor_discrimination": 1.0,
            "barcoding_gap": 0.1,
            "misclassification_risk": 0.0,
            "alignment_reliability": 1.0,
            "markerseek_score": 88.5,
        },
        "pi_curve": {"positions": [1, 2, 3, 4], "pi": [0.0, 0.0, 0.25, 0.25]},
        "snp_positions": [3],
        "indel_positions": [4],
        "aligned_block": {"sample_a": "ACGT", "sample_b": "AC-T"},
        "haplotypes": {"sample_a": "H001", "sample_b": "H002"},
        "species_map": {"sample_a": "Species a", "sample_b": "Species b"},
        "haplotype_network": {
            "nodes": [
                {"id": "H001", "frequency": 1, "species": ["Species a"], "x": 0.0, "y": 0.0},
                {"id": "H002", "frequency": 1, "species": ["Species b"], "x": 1.0, "y": 0.0},
            ],
            "edges": [{"source": "H001", "target": "H002", "weight": 1}],
        },
        "species_pca": {
            "samples": [
                {"sample_name": "sample_a", "species": "Species a", "pc1": -0.1, "pc2": 0.0},
                {"sample_name": "sample_b", "species": "Species b", "pc1": 0.1, "pc2": 0.0},
            ],
            "explained_variance": [1.0, 0.0],
        },
        "primers": [
            {
                "primer_id": "geneA-geneB_p1",
                "label_name": "geneA-geneB",
                "rank": 1,
                "fwd_seq": "ACGTACGTACGTACGTAC",
                "rev_seq": "TGCATGCATGCATGCATG",
                "fwd_len": 18,
                "rev_len": 18,
                "fwd_gc": 50.0,
                "rev_gc": 50.0,
                "fwd_tm": 58.0,
                "rev_tm": 58.0,
                "fwd_self_any_th": 0.0,
                "rev_self_any_th": 0.0,
                "primer3_penalty": 0.1,
                "target_start": 1,
                "target_end": 4,
                "amplicon_min_len": 100,
                "amplicon_max_len": 102,
                "amplicon_mean_len": 101.0,
                "cross_species_success_rate": 1.0,
                "amplicon_variable_sites": 1,
                "amplicon_indel_sites": 0,
                "primer_score": 95.0,
            }
        ],
    }
