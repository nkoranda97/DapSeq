"""
Tests for narrow_peak_to_fasta.narrow_peak_to_fasta(): summit extension,
full-peak mode and its sub-peak dedup, maxpeaks ranking, chromosome-edge
clipping, chrom_filter, FIMO headers, and error handling.

conftest.py adds workflow/scripts/ to sys.path so the module imports directly.
"""

import random
from pathlib import Path

import pytest
from Bio import SeqIO

import narrow_peak_to_fasta as m


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

random.seed(7)
CHR1 = "".join(random.choice("ACGT") for _ in range(500))
CHR2 = "".join(random.choice("ACGT") for _ in range(300))


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


def _read_fasta(path: Path) -> list[tuple[str, str]]:
    with open(path) as fh:
        return [(rec.description, str(rec.seq)) for rec in SeqIO.parse(fh, "fasta")]


def _run(tmp_path, rows, maxpeaks=None, extend_bp=10, fimocoords=False,
         filter_chroms=None, genome=(("chr1", CHR1), ("chr2", CHR2))):
    peaks, fasta, out = tmp_path / "p.narrowPeak", tmp_path / "g.fa", tmp_path / "o.fa"
    _write_narrowpeak(peaks, rows)
    _write_genome(fasta, list(genome))
    m.narrow_peak_to_fasta(str(peaks), str(fasta), str(out), maxpeaks, extend_bp,
                           fimocoords, filter_chroms)
    return _read_fasta(out)


# ---------------------------------------------------------------------------
# Summit extension
# ---------------------------------------------------------------------------

def test_summit_mode_takes_extend_bp_either_side_of_summit(tmp_path):
    rows = [{"chr": "chr1", "start": 100, "stop": 200, "summit": 50, "name": "p1", "fc": 4.5, "qval": 7.25}]
    [(header, seq)] = _run(tmp_path, rows, extend_bp=10)
    assert seq == CHR1[140:160]
    assert header == "p1_foldch=4.5_qscore=7.25_loc=chr1:141-160"


def test_summit_window_is_clipped_at_chromosome_start(tmp_path):
    rows = [{"chr": "chr1", "start": 0, "stop": 50, "summit": 5}]
    [(header, seq)] = _run(tmp_path, rows, extend_bp=20)
    assert seq == CHR1[0:25]
    assert header.endswith("loc=chr1:1-25")


def test_summit_window_is_clipped_at_chromosome_end(tmp_path):
    rows = [{"chr": "chr2", "start": 250, "stop": 300, "summit": 45}]
    [(_, seq)] = _run(tmp_path, rows, extend_bp=20)
    assert seq == CHR2[275:300]


def test_interval_entirely_past_chromosome_end_is_skipped(tmp_path, capsys):
    rows = [{"chr": "chr2", "start": 400, "stop": 450, "summit": 10, "name": "off"}]
    assert _run(tmp_path, rows, extend_bp=5) == []
    assert "Skipping invalid interval for off" in capsys.readouterr().err


# ---------------------------------------------------------------------------
# Full-peak mode
# ---------------------------------------------------------------------------

def test_all_mode_uses_full_peak(tmp_path):
    rows = [{"chr": "chr1", "start": 100, "stop": 180, "summit": 10}]
    [(header, seq)] = _run(tmp_path, rows, extend_bp="all")
    assert seq == CHR1[100:180]
    assert header.endswith("loc=chr1:101-180")


def test_all_mode_keeps_only_highest_fold_subpeak_per_interval(tmp_path):
    # --call-summits writes one row per sub-peak, all sharing the parent interval.
    rows = [
        {"chr": "chr1", "start": 100, "stop": 180, "summit": 10, "name": "sub_a", "fc": 3.0},
        {"chr": "chr1", "start": 100, "stop": 180, "summit": 60, "name": "sub_b", "fc": 8.0},
    ]
    out = _run(tmp_path, rows, extend_bp="all")
    assert [h.split("_foldch")[0] for h, _ in out] == ["sub_b"]


def test_summit_mode_keeps_every_subpeak(tmp_path):
    rows = [
        {"chr": "chr1", "start": 100, "stop": 180, "summit": 10, "name": "sub_a", "fc": 3.0},
        {"chr": "chr1", "start": 100, "stop": 180, "summit": 60, "name": "sub_b", "fc": 8.0},
    ]
    assert len(_run(tmp_path, rows, extend_bp=5)) == 2


# ---------------------------------------------------------------------------
# Ranking and filtering
# ---------------------------------------------------------------------------

def test_peaks_are_written_by_descending_fold_change_and_capped_by_maxpeaks(tmp_path):
    rows = [
        {"chr": "chr1", "start": 50, "stop": 90, "name": "low", "fc": 2.0},
        {"chr": "chr1", "start": 200, "stop": 240, "name": "high", "fc": 9.0},
        {"chr": "chr2", "start": 50, "stop": 90, "name": "mid", "fc": 5.0},
    ]
    names = [h.split("_foldch")[0] for h, _ in _run(tmp_path, rows, maxpeaks=2)]
    assert names == ["high", "mid"]


def test_no_maxpeaks_keeps_all_peaks(tmp_path):
    rows = [{"chr": "chr1", "start": s, "stop": s + 40, "name": f"p{s}"} for s in (50, 150, 250)]
    assert len(_run(tmp_path, rows, maxpeaks=None)) == 3


def test_filter_chroms_drops_peaks_and_logs_counts(tmp_path, capsys):
    rows = [
        {"chr": "chr1", "start": 50, "stop": 90, "name": "keep"},
        {"chr": "chr2", "start": 50, "stop": 90, "name": "drop1"},
        {"chr": "chr2", "start": 150, "stop": 190, "name": "drop2"},
    ]
    out = _run(tmp_path, rows, filter_chroms=["chr2", "chrEBV"])
    assert [h.split("_foldch")[0] for h, _ in out] == ["keep"]
    err = capsys.readouterr().err
    assert "Filtered 2 peaks on chromosome chr2" in err
    assert "Filtered 0 peaks on chromosome chrEBV" in err


def test_filtering_every_peak_writes_empty_fasta(tmp_path):
    rows = [{"chr": "chr2", "start": 50, "stop": 90}]
    out_path = tmp_path / "o.fa"
    assert _run(tmp_path, rows, filter_chroms=["chr2"]) == []
    assert out_path.read_text() == ""


# ---------------------------------------------------------------------------
# Headers and errors
# ---------------------------------------------------------------------------

def test_fimocoords_header_is_one_based_closed_interval(tmp_path):
    rows = [{"chr": "chr1", "start": 100, "stop": 200, "summit": 50}]
    [(header, _)] = _run(tmp_path, rows, extend_bp=10, fimocoords=True)
    assert header == "chr1:141-160"


def test_chromosome_missing_from_genome_raises_naming_hint(tmp_path):
    rows = [{"chr": "1", "start": 50, "stop": 90}]
    with pytest.raises(KeyError, match="chromosome naming"):
        _run(tmp_path, rows)
