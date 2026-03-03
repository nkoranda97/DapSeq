"""
Tests for scripts/combine_peaks.py

Tests cover the pure functions (Peak.to_bed_row, parse_macs_summits,
parse_gem_events, write_bed) with synthetic data and real-data smoke tests.
There is no golden combined output to diff against (the old script was silently
corrupting it), so integration is validated structurally.
"""

from pathlib import Path

import pytest

import combine_peaks as m

RESULTS = Path(__file__).parent.parent / "results"


# ── helpers ───────────────────────────────────────────────────────────────────

def _write_macs_summits(path: Path, rows: list[tuple]) -> None:
    """rows: (chrom, start, end, name, score)"""
    with open(path, "w") as fh:
        for chrom, start, end, name, score in rows:
            fh.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\n")


def _write_gem_events(path: Path, rows: list[tuple], has_control: bool = True) -> None:
    """rows: (position_str, ip, fold_or_none)"""
    with open(path, "w") as fh:
        fh.write(
            "Position\tIP\tControl\tFold\tExpectd\tQ_-lg10\tP_-lg10"
            "\tP_poiss\tIPvsEMP\tNoise\tKmerGroup\tMotifId\tKG_score\tStrand\n"
        )
        for pos, ip, fold in rows:
            fold_val = str(fold) if (has_control and fold is not None) else "NaN"
            ctrl_val = "10.0" if has_control else "NaN"
            fh.write(
                f"{pos}\t{ip}\t{ctrl_val}\t{fold_val}"
                "\t1.0\t999\tNaN\t999\t0.0\t0.0\t----\t-1\t0.0\t*\n"
            )


# ── Peak.to_bed_row ───────────────────────────────────────────────────────────

def test_to_bed_row_basic():
    peak = m.Peak(chrom="chr1", pos=1000, name="p1", score=5.0)
    row = peak.to_bed_row(half_win=80)
    assert row == "chr1\t920\t1080\tp1\t5.0\t.\n"


def test_to_bed_row_clamps_negative_start():
    """pos - half_win < 0 must be clamped to 0."""
    peak = m.Peak(chrom="chrX", pos=10, name="p1", score=1.0)
    fields = peak.to_bed_row(half_win=80).strip().split("\t")
    assert int(fields[1]) == 0


def test_to_bed_row_fields():
    peak = m.Peak(chrom="2", pos=500, name="GEM_7", score=3.14)
    fields = peak.to_bed_row(half_win=50).strip().split("\t")
    assert fields[0] == "2"
    assert int(fields[1]) == 450
    assert int(fields[2]) == 550
    assert fields[3] == "GEM_7"
    assert fields[5] == "."


# ── parse_macs_summits ────────────────────────────────────────────────────────

def test_parse_macs_summits_basic(tmp_path):
    f = tmp_path / "summits.bed"
    _write_macs_summits(f, [
        ("chr1", 500, 501, "peak1", 5.0),
        ("chr2", 200, 201, "peak2", 2.0),
    ])
    peaks = m.parse_macs_summits(str(f), min_score=1.0)
    assert len(peaks) == 2
    assert peaks[0].chrom == "chr1"
    assert peaks[0].pos == 500


def test_parse_macs_summits_score_filter(tmp_path):
    f = tmp_path / "summits.bed"
    _write_macs_summits(f, [
        ("chr1", 500, 501, "peak1", 5.0),
        ("chr2", 200, 201, "peak2", 0.5),  # below threshold
    ])
    peaks = m.parse_macs_summits(str(f), min_score=1.0)
    assert len(peaks) == 1
    assert peaks[0].chrom == "chr1"


def test_parse_macs_summits_names_prefixed(tmp_path):
    """Output peaks should be named MACS_<index>."""
    f = tmp_path / "summits.bed"
    _write_macs_summits(f, [("chr1", 100, 101, "p1", 5.0)])
    peaks = m.parse_macs_summits(str(f), min_score=0.0)
    assert peaks[0].name.startswith("MACS_")


# ── parse_gem_events ──────────────────────────────────────────────────────────

def test_parse_gem_events_with_control(tmp_path):
    f = tmp_path / "events.txt"
    _write_gem_events(f, [
        ("chr1:1000", 100.0, 5.0),
        ("chr2:2000",  50.0, 0.5),  # below threshold
    ], has_control=True)

    peaks = m.parse_gem_events(str(f), min_score=1.0)
    assert len(peaks) == 1
    assert peaks[0].chrom == "chr1"
    assert peaks[0].pos == 1000
    assert peaks[0].score == pytest.approx(5.0)


def test_parse_gem_events_no_control_keeps_all(tmp_path):
    """When Fold is all NaN, all events are kept regardless of min_score."""
    f = tmp_path / "events.txt"
    _write_gem_events(f, [
        ("chr1:1000", 100.0, None),
        ("chr2:2000",  50.0, None),
    ], has_control=False)

    peaks = m.parse_gem_events(str(f), min_score=999.0)
    assert len(peaks) == 2


def test_parse_gem_events_no_control_uses_ip_score(tmp_path):
    """Without a control, peak score is the IP read count."""
    f = tmp_path / "events.txt"
    _write_gem_events(f, [("chr1:500", 123.4, None)], has_control=False)

    peaks = m.parse_gem_events(str(f), min_score=0.0)
    assert peaks[0].score == pytest.approx(123.4, rel=1e-3)


def test_parse_gem_events_position_parsed(tmp_path):
    """chrom and pos are split correctly from 'chr:pos' format."""
    f = tmp_path / "events.txt"
    _write_gem_events(f, [("9:3024502", 500.0, None)], has_control=False)

    peaks = m.parse_gem_events(str(f), min_score=0.0)
    assert peaks[0].chrom == "9"
    assert peaks[0].pos == 3024502


# ── write_bed ─────────────────────────────────────────────────────────────────

def test_write_bed_row_count(tmp_path):
    peaks = [m.Peak("chr1", 1000, "p1", 5.0), m.Peak("chr2", 2000, "p2", 3.0)]
    out = tmp_path / "out.bed"
    m.write_bed(peaks, str(out), half_win=80)
    lines = out.read_text().splitlines()
    assert len(lines) == 2


def test_write_bed_coordinates(tmp_path):
    peaks = [m.Peak("chr1", 1000, "p1", 5.0)]
    out = tmp_path / "out.bed"
    m.write_bed(peaks, str(out), half_win=80)
    fields = out.read_text().strip().split("\t")
    assert fields[:3] == ["chr1", "920", "1080"]


def test_write_bed_empty(tmp_path):
    out = tmp_path / "out.bed"
    m.write_bed([], str(out), half_win=80)
    assert out.read_text() == ""


# ── real-data smoke tests ─────────────────────────────────────────────────────

def test_real_macs_summits_parsed():
    """TF1 MACS summits parse without error and produce well-formed Peaks."""
    summits = RESULTS / "MACS" / "TF1_summits.bed"
    peaks = m.parse_macs_summits(str(summits), min_score=1.0)
    assert len(peaks) > 0
    for p in peaks:
        assert p.chrom != ""
        assert p.pos >= 0
        assert p.score > 1.0


def test_real_gem_no_control_all_kept():
    """TF1 GEM events have no control — all are kept and scored by IP."""
    events = RESULTS / "GEM" / "TF1" / "TF1.GEM_events.txt"
    peaks = m.parse_gem_events(str(events), min_score=1.0)
    assert len(peaks) > 0
    assert all(p.score > 0 for p in peaks)


def test_real_macs_bed_matches_known_good(tmp_path):
    """Per-tool MACS BED written from real summits must match the known output."""
    summits = RESULTS / "MACS" / "TF1_summits.bed"
    expected = RESULTS / "compare_bed" / "TF1.MACS.bed"
    out = tmp_path / "MACS.bed"

    peaks = m.parse_macs_summits(str(summits), min_score=1.0)
    m.write_bed(peaks, str(out), half_win=40)  # window_size=80 → half=40

    assert out.read_text() == expected.read_text()


def test_real_gem_bed_structural(tmp_path):
    """GEM BED rows must be 6-column tab-separated with valid coordinates.

    No golden output exists for GEM.bed — the old script was silently writing
    empty files.  We validate structure instead.
    """
    events = RESULTS / "GEM" / "TF1" / "TF1.GEM_events.txt"
    out = tmp_path / "GEM.bed"

    peaks = m.parse_gem_events(str(events), min_score=1.0)
    m.write_bed(peaks, str(out), half_win=40)

    lines = out.read_text().splitlines()
    assert len(lines) > 0
    for line in lines:
        fields = line.split("\t")
        assert len(fields) == 6, f"Expected 6 BED columns, got {len(fields)}"
        start, end = int(fields[1]), int(fields[2])
        assert start >= 0
        assert end > start
