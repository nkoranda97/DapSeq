"""
Tests for mapping_pct computation added to collect_stats._safe_pct
and the mapping_pct column surfaced through report.make_cols / build_row.
"""

import csv
import os
import sys
from pathlib import Path

import pytest

import collect_stats as cs
import report as rp


# ---------------------------------------------------------------------------
# _safe_pct unit tests
# ---------------------------------------------------------------------------

def test_safe_pct_normal():
    assert cs._safe_pct("90000", "100000") == "90.0"


def test_safe_pct_rounds_to_2dp():
    # 1/3 * 100 = 33.33333... → "33.33"
    assert cs._safe_pct("1", "3") == "33.33"


def test_safe_pct_mapped_equals_trimmed():
    assert cs._safe_pct("50000", "50000") == "100.0"


def test_safe_pct_numerator_na():
    assert cs._safe_pct("NA", "100000") == "NA"


def test_safe_pct_denominator_na():
    assert cs._safe_pct("90000", "NA") == "NA"


def test_safe_pct_both_na():
    assert cs._safe_pct("NA", "NA") == "NA"


def test_safe_pct_zero_denominator():
    assert cs._safe_pct("1000", "0") == "NA"


# ---------------------------------------------------------------------------
# report.make_cols includes mapping_pct
# ---------------------------------------------------------------------------

def test_make_cols_contains_mapping_pct():
    assert "mapping_pct" in rp.make_cols()


def test_make_cols_does_not_contain_mapping_rate():
    assert "mapping_rate" not in rp.make_cols()


def test_columns_does_not_contain_mapping_rate():
    assert "mapping_rate" not in cs.COLUMNS


def test_make_cols_order():
    assert rp.make_cols() == [
        "sample",
        "total_reads",
        "subsampled_frags",
        "trimmed_reads",
        "mapped_reads",
        "mapping_pct",
        "median_frag_size",
        "num_peaks",
        "num_peaks_filt",
        "num_peaks_bl",
        "num_peaks_rmsk",
        "reads_in_peaks",
        "reads_in_peaks_filt",
        "frip",
        "frip_filt",
        "max_peak_score",
        "motif_peaks",
    ]


def test_mapping_pct_before_reads_in_peaks():
    cols = rp.make_cols()
    assert cols.index("mapping_pct") < cols.index("reads_in_peaks")


def test_renamed_columns_present():
    cols = rp.make_cols()
    assert "reads_in_peaks_filt" in cols
    assert "frip_filt" in cols
    assert "num_peaks_filt" in cols
    # legacy / per-fold-level names must be gone from the reporting surface
    assert "reads_in_peaks_5fold" not in cols
    for legacy in (
        "foldch_level1", "foldch_level2", "foldch_level3",
        "num_peaks_fold1", "num_peaks_fold2", "num_peaks_fold3",
        "reads_in_peaks_fold1", "reads_in_peaks_fold2", "reads_in_peaks_fold3",
        "frip_fold1", "frip_fold2", "frip_fold3",
    ):
        assert legacy not in cols


# ---------------------------------------------------------------------------
# collect_stats._gate_subsampled
# ---------------------------------------------------------------------------

def test_gate_subsampled_happy_path():
    assert cs._gate_subsampled("5000000", "5000000") == "5000000"


def test_gate_subsampled_max_frags_unset():
    assert cs._gate_subsampled("5000000", None) == "NA"


def test_gate_subsampled_subsampled_na():
    assert cs._gate_subsampled("NA", "5000000") == "NA"


def test_gate_subsampled_max_frags_zero():
    assert cs._gate_subsampled("5000000", 0) == "NA"


# ---------------------------------------------------------------------------
# collect_stats._parse_subsample_log
# ---------------------------------------------------------------------------

def _write_subsample_log(tmp_path, input_reads, output_reads, output_pct="100.00"):
    """Write a minimal reformat.sh-style subsample log and return its path."""
    log_path = tmp_path / "subsample.log"
    log_path.write_text(
        "Set INTERLEAVED to true\n"
        "Input is being processed as paired\n"
        f"Input:                  \t{input_reads} reads          \t0 bases\n"
        f"Output:                 \t{output_reads} reads ({output_pct}%) \t0 bases ({output_pct}%)\n"
    )
    return str(log_path)


def test_parse_subsample_log_pe_no_subsampling(tmp_path):
    # Real CXXC01.subsample.log: Input/Output both 125671976 reads (R1+R2 combined).
    log = _write_subsample_log(tmp_path, 125671976, 125671976)
    total, subsampled = cs._parse_subsample_log(log, is_pe=True)
    assert total == "125671976"
    assert subsampled == "62835988"


def test_parse_subsample_log_se_no_subsampling(tmp_path):
    # Real DEL1.subsample.log: Input/Output both 3536279 reads.
    log = _write_subsample_log(tmp_path, 3536279, 3536279)
    total, subsampled = cs._parse_subsample_log(log, is_pe=False)
    assert total == "3536279"
    assert subsampled == "3536279"


def test_parse_subsample_log_pe_with_subsampling(tmp_path):
    log = _write_subsample_log(tmp_path, 125671976, 10000000, output_pct="7.96")
    total, subsampled = cs._parse_subsample_log(log, is_pe=True)
    assert total == "125671976"
    assert subsampled == "5000000"


def test_parse_subsample_log_total_not_halved_for_pe(tmp_path):
    # Encodes the bug: total_reads must match the raw Input: count, not Input: // 2,
    # so it stays in the same units as trimmed_reads/mapped_reads.
    log = _write_subsample_log(tmp_path, 125671976, 125671976)
    total, _ = cs._parse_subsample_log(log, is_pe=True)
    assert total == "125671976"
    assert total != "62835988"


# ---------------------------------------------------------------------------
# report.build_row passes mapping_pct through from stats CSV
# ---------------------------------------------------------------------------

def _write_stats_csv(tmp_path, sample, extras=None):
    """Write a minimal {sample}.stats.csv and return the stats dir."""
    row = {c: "NA" for c in cs.COLUMNS}
    row["sample"] = sample
    row["trimmed_reads"] = "100000"
    row["mapped_reads"]  = "90000"
    row["mapping_pct"]   = "90.0"
    if extras:
        row.update(extras)
    stats_dir = tmp_path / "stats"
    stats_dir.mkdir(exist_ok=True)
    path = stats_dir / f"{sample}.stats.csv"
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=cs.COLUMNS)
        writer.writeheader()
        writer.writerow(row)
    return str(stats_dir)


def test_build_row_mapping_pct_present(tmp_path):
    stats_dir = _write_stats_csv(tmp_path, "s1")
    row = rp.build_row("s1", stats_dir)
    assert row["mapping_pct"] == "90.0"


def test_build_row_mapping_pct_na_when_missing(tmp_path):
    """If stats CSV has mapping_pct=NA, build_row passes it through as NA."""
    stats_dir = _write_stats_csv(tmp_path, "s2", extras={"mapping_pct": "NA"})
    row = rp.build_row("s2", stats_dir)
    assert row["mapping_pct"] == "NA"


def test_run_csv_includes_mapping_pct(tmp_path):
    stats_dir = _write_stats_csv(tmp_path, "s1")
    csv_out = str(tmp_path / "report.csv")
    rp.run_csv(["s1"], stats_dir, csv_out)
    with open(csv_out) as fh:
        reader = csv.DictReader(fh)
        rows = list(reader)
    assert "mapping_pct" in reader.fieldnames
    assert rows[0]["mapping_pct"] == "90.0"


def test_run_csv_excludes_mapping_rate_and_includes_renamed_frip(tmp_path):
    stats_dir = _write_stats_csv(
        tmp_path, "s1", extras={"reads_in_peaks_filt": "9000"}
    )
    csv_out = str(tmp_path / "report.csv")
    rp.run_csv(["s1"], stats_dir, csv_out)
    with open(csv_out) as fh:
        reader = csv.DictReader(fh)
        rows = list(reader)
    assert "mapping_rate" not in reader.fieldnames
    assert "frip_filt" in reader.fieldnames
    assert rows[0]["frip_filt"] == "10.0"


# ---------------------------------------------------------------------------
# num_peaks_bl — blacklist filter stats
# ---------------------------------------------------------------------------

def test_make_cols_contains_num_peaks_bl():
    assert "num_peaks_bl" in rp.make_cols()


def test_num_peaks_bl_in_int_cols():
    assert "num_peaks_bl" in rp.INT_COLS


def test_collect_stats_columns_contains_num_peaks_bl():
    assert "num_peaks_bl" in cs.COLUMNS


def test_num_peaks_bl_after_num_peaks_filt_in_make_cols():
    cols = rp.make_cols()
    assert cols.index("num_peaks_bl") == cols.index("num_peaks_filt") + 1


def test_num_peaks_bl_na_by_default(tmp_path):
    """When blacklist filtering is disabled, num_peaks_bl should be NA in the report."""
    stats_dir = _write_stats_csv(tmp_path, "s1")
    row = rp.build_row("s1", stats_dir)
    assert row.get("num_peaks_bl", "NA") == "NA"


def test_num_peaks_bl_counted_when_enabled(tmp_path):
    """When blacklist filtering is enabled, num_peaks_bl reflects the post-filter count."""
    stats_dir = _write_stats_csv(tmp_path, "s1", extras={"num_peaks_bl": "42"})
    row = rp.build_row("s1", stats_dir)
    assert row["num_peaks_bl"] == "42"


def test_count_lines_empty_file(tmp_path):
    f = tmp_path / "empty.narrowPeak"
    f.write_text("")
    assert cs._count_lines(str(f)) == "0"


def test_count_lines_ten_peaks(tmp_path):
    f = tmp_path / "peaks.narrowPeak"
    f.write_text("\n".join([f"chr1\t{i}\t{i+100}" for i in range(10)]) + "\n")
    assert cs._count_lines(str(f)) == "10"


# ---------------------------------------------------------------------------
# collect_stats.COLUMNS is simplified to a single filtered set
# ---------------------------------------------------------------------------

def test_collect_stats_columns_simplified():
    assert cs.COLUMNS == [
        "sample",
        "total_reads",
        "subsampled_frags",
        "trimmed_reads",
        "mapped_reads",
        "mapping_pct",
        "median_frag_size",
        "num_peaks",
        "num_peaks_filt",
        "num_peaks_bl",
        "num_peaks_rmsk",
        "reads_in_peaks",
        "reads_in_peaks_filt",
        "max_peak_score",
        "motif_peaks",
    ]


def test_collect_stats_columns_drop_per_fold_level():
    for legacy in (
        "foldch_level1", "num_peaks_fold1", "num_peaks_fold3",
        "reads_in_peaks_fold1", "reads_in_peaks_fold2", "reads_in_peaks_fold3",
    ):
        assert legacy not in cs.COLUMNS


# ---------------------------------------------------------------------------
# report.build_row derives frip_filt (and no per-fold frip)
# ---------------------------------------------------------------------------

def test_build_row_frip_filt_derived(tmp_path):
    stats_dir = _write_stats_csv(
        tmp_path, "s1", extras={"reads_in_peaks_filt": "9000"}
    )
    row = rp.build_row("s1", stats_dir)
    assert row["frip_filt"] == 10.0  # 9000 / 90000 * 100
    for legacy in ("frip_fold1", "frip_fold2", "frip_fold3"):
        assert legacy not in row


# ---------------------------------------------------------------------------
# HTML header renders experiment metadata + filter note only when provided
# ---------------------------------------------------------------------------

def _write_report_html(tmp_path, **kwargs):
    out = tmp_path / "report.html"
    rp.write_html(
        rows=[{"sample": "s1"}],
        cols=rp.make_cols(),
        logo_b64_map={},
        logo_rc_b64_map={},
        factorbook_logo_map={},
        out_path=str(out),
        **kwargs,
    )
    return out.read_text()


def test_html_header_shows_experiment_metadata(tmp_path):
    html = _write_report_html(
        tmp_path, experiment_date="2026-07-14", gdna_batch="batch-07"
    )
    assert "Experiment date" in html and "2026-07-14" in html
    assert "gDNA batch" in html and "batch-07" in html


def test_html_header_omits_metadata_when_unset(tmp_path):
    html = _write_report_html(tmp_path)
    assert "Experiment date" not in html
    assert "gDNA batch" not in html


def test_html_header_shows_filter_foldch(tmp_path):
    html = _write_report_html(tmp_path, filter_foldch=5)
    assert "Filtered peaks" in html
    assert "5" in html


# ---------------------------------------------------------------------------
# report.read_run_metadata parses the config snapshot
# ---------------------------------------------------------------------------

def test_read_run_metadata_from_snapshot(tmp_path):
    (tmp_path / "config_used.yaml").write_text(
        "experiment_date: '2026-07-14'\n"
        "gdna_batch: batch-07\n"
        "macs3:\n"
        "  foldch_levels: [2, 5, 15]\n"
        "  meme_foldch_level: 2\n"
    )
    filt, exp, batch = rp.read_run_metadata(str(tmp_path))
    assert filt == 5           # foldch_levels[meme_foldch_level - 1]
    assert exp == "2026-07-14"
    assert batch == "batch-07"


def test_read_run_metadata_missing_snapshot(tmp_path):
    assert rp.read_run_metadata(str(tmp_path)) == (None, None, None)
