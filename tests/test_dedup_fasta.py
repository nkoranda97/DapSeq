"""
Tests for scripts/dedup_fasta.py

Logic: first occurrence of each FASTA name is kept; later duplicates are
dropped.  Output is written with one unbroken sequence line per record.
"""

from pathlib import Path

import dedup_fasta as m

RESULTS = Path(__file__).parent.parent / "results" / "fasta"

FASTA_WITH_DUPS = """\
>seq1
ACGT
>seq2
TTTT
>seq1
GGGG
>seq3
CCCC
"""


# ── helpers ───────────────────────────────────────────────────────────────────

def _headers(path: Path) -> list[str]:
    return [l for l in path.read_text().splitlines() if l.startswith(">")]


# ── unit tests ────────────────────────────────────────────────────────────────

def test_duplicate_names_removed(tmp_path):
    inp = tmp_path / "in.fasta"
    out = tmp_path / "out.fasta"
    inp.write_text(FASTA_WITH_DUPS)

    m.dedup_fasta(str(inp), str(out))

    assert _headers(out) == [">seq1", ">seq2", ">seq3"]


def test_first_occurrence_kept(tmp_path):
    """When a name appears twice, the first sequence is written (not the second)."""
    inp = tmp_path / "in.fasta"
    out = tmp_path / "out.fasta"
    inp.write_text(FASTA_WITH_DUPS)

    m.dedup_fasta(str(inp), str(out))

    lines = out.read_text().splitlines()
    seq1_seq = lines[lines.index(">seq1") + 1]
    assert seq1_seq == "ACGT"  # first occurrence, not the duplicate "GGGG"


def test_no_duplicates_unchanged(tmp_path):
    fasta = ">a\nACGT\n>b\nTTTT\n"
    inp = tmp_path / "in.fasta"
    out = tmp_path / "out.fasta"
    inp.write_text(fasta)

    m.dedup_fasta(str(inp), str(out))

    assert _headers(out) == [">a", ">b"]


def test_empty_input(tmp_path):
    inp = tmp_path / "in.fasta"
    out = tmp_path / "out.fasta"
    inp.write_text("")

    m.dedup_fasta(str(inp), str(out))

    assert out.read_text() == ""


# ── real-data regression ──────────────────────────────────────────────────────

def test_real_tf1_no_duplicates_present(tmp_path):
    """TF1.fasta has no duplicate headers — output should be identical."""
    out = tmp_path / "out.nodup"
    m.dedup_fasta(str(RESULTS / "TF1.fasta"), str(out))

    assert out.read_text() == (RESULTS / "TF1.fasta.nodup").read_text()
