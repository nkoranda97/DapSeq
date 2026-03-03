"""
Tests for scripts/tandem_filter.py

Logic: remove FASTA records that have >= k_max occurrences of their most
frequent k-mer, or that contain more than one N.  The rest are written to
the output file.
"""

from pathlib import Path

from Bio import SeqIO

import tandem_filter as m

RESULTS = Path(__file__).parent.parent / "results" / "fasta"


# ── helpers ───────────────────────────────────────────────────────────────────

def _write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n{seq}\n")


def _read_names(path: Path) -> list[str]:
    with open(path) as fh:
        return [rec.name for rec in SeqIO.parse(fh, "fasta")]


# ── unit tests ────────────────────────────────────────────────────────────────

def test_tandem_repeat_removed(tmp_path):
    """A sequence dominated by one k-mer (count >= k_max) is filtered out."""
    # "ACGTAC..." repeated → ACGT appears 20 times; k_max=5 is well below that.
    # "clean" sequence's top 4-mer ("ATCG") appears only 3 times → kept.
    tandem = "ACGT" * 20
    clean  = "AAACCCGGGTTTATCGATCGATCG"
    inp = tmp_path / "in.fasta"
    out = tmp_path / "out.fasta"
    _write_fasta(inp, [("tandem", tandem), ("clean", clean)])

    m.tandem_filter(str(inp), str(out), k=4, k_max=5)

    kept = _read_names(out)
    assert "clean" in kept
    assert "tandem" not in kept


def test_clean_sequences_kept(tmp_path):
    """Sequences with no dominant k-mer and no Ns are all kept."""
    inp = tmp_path / "in.fasta"
    out = tmp_path / "out.fasta"
    _write_fasta(inp, [
        ("s1", "ATCGATCGATCGATCG"),
        ("s2", "GCTAGCTAGCTAGCTA"),
    ])

    m.tandem_filter(str(inp), str(out), k=4, k_max=10)  # very permissive

    assert _read_names(out) == ["s1", "s2"]


def test_n_run_removed(tmp_path):
    """A sequence with > 1 N is removed regardless of k-mer counts."""
    n_seq = "ACGTNNNNACGT"
    ok_seq = "ACGTACGTACGT"
    inp = tmp_path / "in.fasta"
    out = tmp_path / "out.fasta"
    _write_fasta(inp, [("n_rich", n_seq), ("ok", ok_seq)])

    m.tandem_filter(str(inp), str(out), k=4, k_max=100)  # k_max won't trigger

    kept = _read_names(out)
    assert "ok" in kept
    assert "n_rich" not in kept


def test_single_n_allowed(tmp_path):
    """Exactly one N is below the threshold (len(ns) > 1 is the condition)."""
    one_n = "ACGTNACGT"
    inp = tmp_path / "in.fasta"
    out = tmp_path / "out.fasta"
    _write_fasta(inp, [("one_n", one_n)])

    m.tandem_filter(str(inp), str(out), k=4, k_max=100)

    assert _read_names(out) == ["one_n"]


def test_empty_input(tmp_path):
    inp = tmp_path / "in.fasta"
    out = tmp_path / "out.fasta"
    inp.write_text("")

    m.tandem_filter(str(inp), str(out), k=4, k_max=3)

    assert _read_names(out) == []


# ── real-data regression ──────────────────────────────────────────────────────

def test_real_tf1_matches_known_good(tmp_path):
    """Filtering TF1.fasta.nodup with pipeline defaults reproduces the known output."""
    out = tmp_path / "filtered.fasta"
    m.tandem_filter(
        str(RESULTS / "TF1.fasta.nodup"),
        str(out),
        k=6,
        k_max=4,  # matches Snakefile params
    )
    assert _read_names(out) == _read_names(RESULTS / "TF1.fasta.filtered.fasta")
