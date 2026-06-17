"""
Tests for narrow_peak_to_fasta._triplet_entropy and the integrated
low-complexity filter inside narrow_peak_to_fasta().

conftest.py adds workflow/scripts/ to sys.path so the module imports directly.
"""

from pathlib import Path

import pytest
from Bio import SeqIO

import narrow_peak_to_fasta as m


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _write_genome(path: Path, records: list[tuple[str, str]]) -> None:
    """Write a minimal genome FASTA."""
    with open(path, "w") as fh:
        for chrom, seq in records:
            fh.write(f">{chrom}\n{seq}\n")


def _write_narrowpeak(path: Path, rows: list[dict]) -> None:
    """Write narrowPeak rows (10-column BED)."""
    with open(path, "w") as fh:
        for r in rows:
            fh.write(
                f"{r['chr']}\t{r['start']}\t{r['stop']}\t{r.get('name','peak')}\t"
                f"{r.get('score',100)}\t.\t{r.get('fc',5.0)}\t"
                f"{r.get('qval',10.0)}\t{r.get('pval',12.0)}\t"
                f"{r.get('summit', (r['stop']-r['start'])//2)}\n"
            )


def _read_fasta_names(path: Path) -> list[str]:
    with open(path) as fh:
        return [rec.id for rec in SeqIO.parse(fh, "fasta")]


def _read_fasta_seqs(path: Path) -> list[str]:
    with open(path) as fh:
        return [str(rec.seq) for rec in SeqIO.parse(fh, "fasta")]


# ---------------------------------------------------------------------------
# _triplet_entropy unit tests
# ---------------------------------------------------------------------------


def test_triplet_entropy_flags_perfect_repeat_as_low_complexity():
    """A perfect repeat uses few 3-mers and has low entropy."""
    assert m._triplet_entropy("ACGTAC" * 10) < 3.0


def test_triplet_entropy_clean_sequence_is_higher_complexity():
    """A diverse sequence has clearly higher 3-mer entropy."""
    clean = "ATCGATCGGCTAGCTAGCTAGCTAGCTAGCATGCATGCATGCATGCATGC"
    assert m._triplet_entropy(clean) >= 3.0


def test_triplet_entropy_degenerate_repeats_are_low_complexity():
    """Fuzzy GA/GGA repeats still score low by triplet entropy."""
    assert m._triplet_entropy("GA" * 60) < 2.0
    assert m._triplet_entropy("GGA" * 40) < 2.0


def test_triplet_entropy_all_n_window_scores_zero():
    """Triplets containing Ns are ignored; an all-N window scores 0."""
    assert m._triplet_entropy("N" * 100) == 0.0


def test_triplet_entropy_ignores_n_without_hard_n_count_filter():
    """Multiple Ns do not automatically fail when remaining triplets are diverse."""
    clean = "AGTCTCGATACTTTTCCATAGAGGCACAAAGTGGGGTCACCTAAGCATGGGCCGAGGTAAATCAACGTCCGGGAAACACGGACGTGCCTGTATACCGTCTTGCGGACTCGCTACTGCTAC"
    clean_with_ns = clean[:40] + "NN" + clean[42:]
    assert m._triplet_entropy(clean_with_ns) >= 3.0


def test_triplet_entropy_case_insensitive():
    """Lowercase input is normalised to uppercase before checking."""
    assert m._triplet_entropy(("acgtac" * 10).lower()) < 3.0


# ---------------------------------------------------------------------------
# narrow_peak_to_fasta integration tests with complexity filter
# ---------------------------------------------------------------------------

# Shared genome: chr1 has ~444 bp.
# First 102 bp  → low-complexity ("ACGTAC" * 17)
# Next  102 bp  → also low-complexity
# Next  120 bp  → truly non-repetitive (SHA256-derived; verified max 6-mer count = 1)
# Next  120 bp  → truly non-repetitive (SHA256-derived; verified max 6-mer count = 1)

_LOW_COMPLEXITY_SEQ = "ACGTAC" * 17
_CLEAN_SEQ_A = "AGTCTCGATACTTTTCCATAGAGGCACAAAGTGGGGTCACCTAAGCATGGGCCGAGGTAAATCAACGTCCGGGAAACACGGACGTGCCTGTATACCGTCTTGCGGACTCGCTACTGCTAC"
_CLEAN_SEQ_B = "ATCGACCAATACAAAGGTTGATAAATGCCAGTTGTAGGTTAGCTGTTCCGCCACCCAGGATTCTCGTGACTATCAAACTCTAGTGCACACAATCAATAAACATTATCGTGTGTGTAAGTT"

_GENOME = [("chr1", _LOW_COMPLEXITY_SEQ + _LOW_COMPLEXITY_SEQ + _CLEAN_SEQ_A + _CLEAN_SEQ_B)]


def _make_peaks_and_genome(tmp_path):
    """Return (genome_path, narrowpeak_path) for a 4-peak test set.

    Peak layout on chr1 (extend_bp=50 around summit):
      peak1: summit=51   → coords [1, 101]   → low-complexity segment (fc=10, highest)
      peak2: summit=153  → coords [103, 203] → low-complexity segment (fc=9)
      peak3: summit=255  → coords [205, 305] → clean segment A (fc=8)
      peak4: summit=357  → coords [307, 407] → clean segment B (fc=7)
    """
    genome = tmp_path / "genome.fa"
    _write_genome(genome, _GENOME)

    peaks = tmp_path / "peaks.narrowPeak"
    _write_narrowpeak(peaks, [
        {"chr": "chr1", "start": 0,   "stop": 102,  "name": "peak1", "fc": 10.0, "summit": 51},
        {"chr": "chr1", "start": 102, "stop": 204,  "name": "peak2", "fc": 9.0,  "summit": 51},
        {"chr": "chr1", "start": 204, "stop": 324,  "name": "peak3", "fc": 8.0,  "summit": 51},
        {"chr": "chr1", "start": 324, "stop": 444,  "name": "peak4", "fc": 7.0,  "summit": 51},
    ])
    return str(genome), str(peaks)


def test_complexity_filter_fills_to_maxpeaks(tmp_path):
    """With complexity filter enabled, skipped repeats are replaced by subsequent clean peaks.

    The top 2 peaks are low-complexity repeats; peaks 3 and 4 are clean.
    With maxpeaks=2 and filter enabled, output should contain exactly 2 clean sequences.
    """
    genome, peaks = _make_peaks_and_genome(tmp_path)
    out = str(tmp_path / "out.fasta")

    m.narrow_peak_to_fasta(
        peaks, genome, out,
        maxpeaks=2, extend_bp=50, fimocoords=False,
        complexity_filter_enabled=True, min_entropy=3.0,
    )

    names = _read_fasta_names(Path(out))
    assert len(names) == 2
    assert all("peak3" in n or "peak4" in n for n in names), (
        f"Expected only clean peaks in output, got: {names}"
    )
    assert not any("peak1" in n or "peak2" in n for n in names)


def test_complexity_filter_disabled_keeps_top_peaks(tmp_path):
    """With complexity filter disabled, the top maxpeaks by fold-change are always taken."""
    genome, peaks = _make_peaks_and_genome(tmp_path)
    out = str(tmp_path / "out.fasta")

    m.narrow_peak_to_fasta(
        peaks, genome, out,
        maxpeaks=2, extend_bp=50, fimocoords=False,
        complexity_filter_enabled=False,
    )

    names = _read_fasta_names(Path(out))
    assert len(names) == 2
    assert all("peak1" in n or "peak2" in n for n in names), (
        f"Expected top-2 peaks (low-complexity), got: {names}"
    )


def test_complexity_filter_all_low_complexity_produces_empty_fasta(tmp_path):
    """If all peaks are low-complexity repeats, output is an empty FASTA with no crash."""
    genome = tmp_path / "genome.fa"
    _write_genome(genome, [("chr1", _LOW_COMPLEXITY_SEQ * 4)])

    peaks = tmp_path / "peaks.narrowPeak"
    _write_narrowpeak(peaks, [
        {"chr": "chr1", "start": 0,   "stop": 102, "name": "p1", "fc": 5.0, "summit": 51},
        {"chr": "chr1", "start": 102, "stop": 204, "name": "p2", "fc": 4.0, "summit": 51},
    ])
    out = str(tmp_path / "out.fasta")

    m.narrow_peak_to_fasta(
        str(peaks), str(genome), out,
        maxpeaks=2, extend_bp=50, fimocoords=False,
        complexity_filter_enabled=True, min_entropy=3.0,
    )

    assert _read_fasta_names(Path(out)) == []


def test_complexity_filter_no_maxpeaks_cap(tmp_path):
    """With maxpeaks=None and complexity filter enabled, all clean peaks are kept."""
    genome, peaks = _make_peaks_and_genome(tmp_path)
    out = str(tmp_path / "out.fasta")

    m.narrow_peak_to_fasta(
        peaks, genome, out,
        maxpeaks=None, extend_bp=50, fimocoords=False,
        complexity_filter_enabled=True, min_entropy=3.0,
    )

    names = _read_fasta_names(Path(out))
    # Only the 2 clean peaks should be in output
    assert len(names) == 2
    assert all("peak3" in n or "peak4" in n for n in names)
