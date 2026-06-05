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


def test_mapping_pct_after_mapping_rate():
    cols = rp.make_cols()
    assert cols.index("mapping_pct") == cols.index("mapping_rate") + 1


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
