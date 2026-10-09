"""
Tests for run_factorbook_logo.py: choosing which TF name a sample is looked up
under, matching it against the Factorbook catalog, and picking the motif
with the most sites.
"""

import pytest

import run_factorbook_logo as m


TSV = (
    "dataset_accession\taccession\ttarget\tbiosample\tassembly\tname\n"
    "ENCSR1\tF1\tCTCF\tH1\tGRCh38\tCCGCG\n"
    "ENCSR2\tF2\tCTCF\tK562\tGRCh38\tCCACT\n"
    "ENCSR3\tF3\tMYC\tK562\tGRCh38\tCACGTG\n"
)

MEME = """\
MEME version 5

MOTIF ENCSR1_CCGCG
letter-probability matrix: alength= 4 w= 2 nsites= 40 E= 0
 0.1 0.7 0.1 0.1
 0.1 0.1 0.7 0.1

MOTIF ENCSR2_CCACT
letter-probability matrix: alength= 4 w= 2 nsites= 90 E= 0
 0.7 0.1 0.1 0.1
 0.1 0.1 0.1 0.7

MOTIF ENCSR3_CACGTG
letter-probability matrix: alength= 4 w= 1 nsites= 10 E= 0
 0.25 0.25 0.25 0.25
"""


@pytest.fixture
def catalog(tmp_path):
    tsv, meme = tmp_path / "fb.tsv", tmp_path / "fb.meme"
    tsv.write_text(TSV)
    meme.write_text(MEME)
    return str(tsv), str(meme)


# ── which names a sample is looked up under ──────────────────────────────────

def test_explicit_tf_is_the_only_candidate():
    assert m.tf_candidates("TF_A", "ctcf") == ["CTCF"]


def test_without_tf_try_full_name_then_prefix():
    assert m.tf_candidates("CTCF_rep1") == ["CTCF_REP1", "CTCF"]
    assert m.tf_candidates("myc-2") == ["MYC-2", "MYC"]
    assert m.tf_candidates("Nfkb.r1") == ["NFKB.R1", "NFKB"]


def test_plain_name_has_one_candidate():
    assert m.tf_candidates("CTCF") == ["CTCF"]


# ── catalog lookup ────────────────────────────────────────────────────────────

def test_first_matching_candidate_wins(catalog):
    tsv, _ = catalog
    name, ids = m.find_motif_ids(tsv, ["CTCF_REP1", "CTCF"])
    assert name == "CTCF"
    assert ids == ["ENCSR1_CCGCG", "ENCSR2_CCACT"]


def test_no_match_returns_none(catalog):
    tsv, _ = catalog
    assert m.find_motif_ids(tsv, ["AT1G01010"]) == (None, [])


def test_missing_catalog_column_names_the_column(tmp_path):
    bad = tmp_path / "bad.tsv"
    bad.write_text("dataset_accession\tname\nENCSR1\tCCGCG\n")
    with pytest.raises(ValueError, match="target"):
        m.find_motif_ids(str(bad), ["CTCF"])


def test_best_motif_has_most_sites(catalog):
    _, meme = catalog
    ppm = m.best_motif_ppm(meme, ["ENCSR1_CCGCG", "ENCSR2_CCACT"])
    assert ppm == [[0.7, 0.1, 0.1, 0.1], [0.1, 0.1, 0.1, 0.7]]   # ENCSR2: 90 sites


def test_best_motif_none_when_ids_absent_from_meme(catalog):
    _, meme = catalog
    assert m.best_motif_ppm(meme, ["ENCSR9_NOPE"]) is None
