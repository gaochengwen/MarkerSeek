from __future__ import annotations

import json
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from markerseek import database as db_catalog
from markerseek.web import app as web_app


GENBANK_FIXTURE = """\
LOCUS       NC_TEST01              {length} bp    DNA     circular PLN 01-JAN-2000
DEFINITION  {organism} chloroplast, complete genome.
ACCESSION   NC_TEST01
ORGANISM    {organism}
            Cellular organisms; Eukaryota; Viridiplantae.
FEATURES             Location/Qualifiers
     source          1..{length}
ORIGIN
        1 {sequence}
//
"""


def _write_fixture_genbank(path: Path, organism: str, length: int = 60, gc_ratio: float = 0.5) -> None:
    gc_bases = int(round(length * gc_ratio))
    at_bases = length - gc_bases
    sequence = ("gc" * (gc_bases // 2 + 1))[:gc_bases] + ("at" * (at_bases // 2 + 1))[:at_bases]
    sequence = sequence[:length]
    path.write_text(
        GENBANK_FIXTURE.format(organism=organism, length=length, sequence=sequence),
        encoding="utf-8",
    )


def _build_root(root: Path, *, with_results: bool = False) -> None:
    astragalus_ref = root / "Astragalus_ref"
    astragalus_ref.mkdir(parents=True)
    _write_fixture_genbank(astragalus_ref / "Astragalus_alpinus.gb", "Astragalus alpinus", length=120, gc_ratio=0.34)
    _write_fixture_genbank(astragalus_ref / "Astragalus_membranaceus.gb", "Astragalus membranaceus", length=130, gc_ratio=0.34)

    bupleurum_ref = root / "Bupleurum_ref"
    bupleurum_ref.mkdir(parents=True)
    _write_fixture_genbank(bupleurum_ref / "Bupleurum_chinense.gb", "Bupleurum chinense", length=160, gc_ratio=0.38)

    if with_results:
        results = root / "Astragalus_results"
        results.mkdir()
        marker_tsv = results / "candidate_marker_features.tsv"
        marker_tsv.write_text(
            "feature_id\tfeature_type\tparent_gene\tlabel_name\tstart\tend\tstrand\tlength_bp\t"
            "region\tpi\tvariable_sites\tindel_sites\tconserved_left_bp\tconserved_right_bp\t"
            "primer_available\tspecies_resolution\tunique_haplotype_count\tspecies_specific_haplotype_ratio\t"
            "interspecific_divergence\tintraspecific_divergence\tnearest_neighbor_discrimination\t"
            "barcoding_gap\tmisclassification_risk\talignment_reliability\tmarkerseek_score\n"
            "matK\tgene\tmatK\tmatK\t100\t1500\t1\t1400\tLSC\t0.012\t30\t0\t200\t200\t"
            "yes\t0.85\t12\t0.4\t0.04\tNA\tNA\tNA\tNA\t0.95\t72.5\n"
            "rbcL\tgene\trbcL\trbcL\t1600\t3000\t1\t1400\tLSC\t0.005\t10\t0\t200\t200\t"
            "no\t0.50\t6\t0.2\t0.02\tNA\tNA\tNA\tNA\t0.80\t40.0\n",
            encoding="utf-8",
        )
        primers_tsv = results / "primers.tsv"
        primers_tsv.write_text(
            "primer_id\tlabel_name\trank\tfwd_seq\trev_seq\tfwd_len\trev_len\tfwd_gc\trev_gc\t"
            "fwd_tm\trev_tm\tfwd_self_any_th\trev_self_any_th\tprimer3_penalty\ttarget_start\t"
            "target_end\tamplicon_min_len\tamplicon_max_len\tamplicon_mean_len\t"
            "cross_species_success_rate\tamplicon_variable_sites\tamplicon_indel_sites\tprimer_score\n"
            "matK_p1\tmatK\t1\tACGTACGTACGTACGT\tTGCATGCATGCATGCA\t16\t16\t50.0\t50.0\t"
            "58.0\t58.0\t0.0\t0.0\t1.5\t100\t1500\t1400\t1400\t1400.0\t0.92\t30\t0\t78.0\n",
            encoding="utf-8",
        )


def test_reindex_indexes_two_genera_with_correct_summaries(tmp_path: Path) -> None:
    _build_root(tmp_path)
    db_path = tmp_path / "catalogue.sqlite3"

    result = db_catalog.reindex(root=tmp_path, db_path=db_path)

    assert result == {"indexed": 2, "skipped": 0, "total_seen": 2}
    rows, total = db_catalog.list_genera(db_path)
    assert total == 2
    by_name = {row["genus"]: row for row in rows}
    assert by_name["Astragalus"]["species_count"] == 2
    assert by_name["Bupleurum"]["species_count"] == 1
    assert by_name["Astragalus"]["genome_length_min"] == 120
    assert by_name["Astragalus"]["genome_length_max"] == 130
    assert by_name["Astragalus"]["marker_count"] == 0
    assert by_name["Astragalus"]["primer_count"] == 0
    assert by_name["Astragalus"]["family"] is None
    assert by_name["Astragalus"]["taxonomic_order"] is None


def test_reindex_picks_up_results_drop(tmp_path: Path) -> None:
    _build_root(tmp_path, with_results=True)
    db_path = tmp_path / "catalogue.sqlite3"

    db_catalog.reindex(root=tmp_path, db_path=db_path)

    row = db_catalog.get_genus(db_path, "Astragalus")
    assert row is not None
    assert row["marker_count"] == 2
    assert row["primer_count"] == 1
    assert row["has_primers"] == 1
    assert row["top_marker_label"] == "matK"
    assert row["top_marker_score"] == pytest.approx(72.5)


def test_taxonomy_csv_populates_family_and_order(tmp_path: Path) -> None:
    _build_root(tmp_path)
    (tmp_path / "taxonomy.csv").write_text(
        "genus,family,order\nAstragalus,Fabaceae,Fabales\nBupleurum,Apiaceae,Apiales\n",
        encoding="utf-8",
    )
    db_path = tmp_path / "catalogue.sqlite3"

    db_catalog.reindex(root=tmp_path, db_path=db_path)

    facets = db_catalog.list_taxonomy_facets(db_path)
    assert facets == {"families": ["Apiaceae", "Fabaceae"], "orders": ["Apiales", "Fabales"]}
    row = db_catalog.get_genus(db_path, "Astragalus")
    assert row["family"] == "Fabaceae"
    assert row["taxonomic_order"] == "Fabales"


def test_database_routes_validate_genus_and_block_path_traversal(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _build_root(tmp_path)
    monkeypatch.setenv(db_catalog.DEFAULT_DB_ROOT_ENV, str(tmp_path))
    monkeypatch.setenv(db_catalog.DEFAULT_DB_PATH_ENV, str(tmp_path / "catalogue.sqlite3"))
    fresh_settings = web_app.WebSettings()
    monkeypatch.setattr(web_app, "settings", fresh_settings)
    db_catalog.reindex(root=tmp_path, db_path=fresh_settings.database_db_path)

    client = TestClient(web_app.app)

    home = client.get("/database")
    assert home.status_code == 200
    assert "Genera indexed" in home.text

    browse = client.get("/database/browse?q=Astr")
    assert browse.status_code == 200
    assert "Astragalus" in browse.text

    genus = client.get("/database/genus/Astragalus")
    assert genus.status_code == 200
    assert "Reference only" in genus.text or "Markers + reference" in genus.text

    missing = client.get("/database/genus/DoesNotExist")
    assert missing.status_code == 404

    traversal = client.get("/database/genus/Astragalus/species/..%2F..%2Fetc%2Fpasswd")
    assert traversal.status_code == 404

    species_dl = client.get("/database/genus/Astragalus/species/Astragalus_alpinus.gb")
    assert species_dl.status_code == 200
    assert "LOCUS" in species_dl.text


def test_search_api_returns_prefix_matches(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    _build_root(tmp_path)
    monkeypatch.setenv(db_catalog.DEFAULT_DB_ROOT_ENV, str(tmp_path))
    monkeypatch.setenv(db_catalog.DEFAULT_DB_PATH_ENV, str(tmp_path / "catalogue.sqlite3"))
    fresh_settings = web_app.WebSettings()
    monkeypatch.setattr(web_app, "settings", fresh_settings)
    db_catalog.reindex(root=tmp_path, db_path=fresh_settings.database_db_path)

    client = TestClient(web_app.app)
    payload = client.get("/database/api/search?q=Bup").json()

    assert payload["results"] == [
        {"genus": "Bupleurum", "species_count": 1, "has_results": False}
    ]
