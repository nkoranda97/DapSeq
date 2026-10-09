"""
Tests for the collect_stats.py helpers that read tool output: bbduk trim logs,
narrowPeak scores, FIMO hits and bamPEFragmentSize output.
"""

import subprocess
from types import SimpleNamespace

import collect_stats as cs
import fimo_fixtures as ff


def _write(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(text)
    return str(p)


def test_bbduk_log_result_line_gives_trimmed_reads(tmp_path):
    log = _write(tmp_path, "trim.log",
                 "Input:                  \t1000 reads \t\t100000 bases.\n"
                 "Result:                 \t987 reads (98.70%) \t98123 bases (98.12%)\n")
    assert cs._parse_bbduk_log(log) == "987"


def test_bbduk_log_without_result_is_na(tmp_path):
    assert cs._parse_bbduk_log(_write(tmp_path, "trim.log", "crashed\n")) == "NA"


def _peak(signal):
    return f"chr1\t10\t50\tp\t100\t.\t{signal}\t5.0\t6.0\t20\n"


def test_narrowpeak_max_score_is_largest_signal_value(tmp_path):
    np = _write(tmp_path, "p.narrowPeak", _peak(3.5) + _peak(12.25) + _peak(7))
    assert cs._narrowpeak_max_score(np) == "12.25"


def test_narrowpeak_max_score_skips_short_rows_and_empty_files(tmp_path):
    assert cs._narrowpeak_max_score(_write(tmp_path, "e.narrowPeak", "")) == "NA"
    np = _write(tmp_path, "s.narrowPeak", "chr1\t1\t2\n" + _peak(4))
    assert cs._narrowpeak_max_score(np) == "4.0"


def test_fimo_counts_peaks_not_chromosomes(tmp_path):
    # Three peaks have hits, two of them on chr1; rows are interleaved by p-value.
    tsv = _write(tmp_path, "fimo.tsv", ff.FIMO_PEAKS_THREE_HITS)
    assert cs._fimo_motif_peaks(tsv) == "3"


def test_fimo_peak_hit_by_two_motifs_counts_once(tmp_path):
    tsv = _write(tmp_path, "fimo.tsv", ff.FIMO_PEAK_TWO_MOTIFS)
    assert cs._fimo_motif_peaks(tsv) == "1"


def test_fimo_ran_with_no_hits_is_zero_not_na(tmp_path):
    tsv = _write(tmp_path, "fimo.tsv", ff.FIMO_NO_HITS)
    assert cs._fimo_motif_peaks(tsv) == "0"


def test_fimo_old_coordinate_mode_output_is_na(tmp_path):
    # Hits named by chromosome alone cannot be counted as peaks.
    assert cs._fimo_motif_peaks(_write(tmp_path, "a.tsv", ff.FIMO_OLD_MODE_HITS)) == "NA"
    assert cs._fimo_motif_peaks(_write(tmp_path, "b.tsv", ff.FIMO_OLD_MODE_NO_HITS)) == "NA"


def test_fimo_missing_or_empty_output_is_na(tmp_path):
    assert cs._fimo_motif_peaks(str(tmp_path / "missing.tsv")) == "NA"
    assert cs._fimo_motif_peaks(_write(tmp_path, "empty.tsv", "")) == "NA"


def test_bampe_median_frag_parses_median(monkeypatch):
    out = "Fragment lengths:\nMin.: 50.0\nMedian: 182.0\nMax.: 600.0\n"
    monkeypatch.setattr(subprocess, "run", lambda *a, **k: SimpleNamespace(stdout=out, returncode=0))
    assert cs._bampe_median_frag("x.bam") == "182"


def test_bampe_median_frag_without_median_is_na(monkeypatch):
    monkeypatch.setattr(subprocess, "run", lambda *a, **k: SimpleNamespace(stdout="error\n", returncode=1))
    assert cs._bampe_median_frag("x.bam") == "NA"
