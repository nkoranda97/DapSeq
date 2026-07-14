"""
Tests for workflow/scripts/update_db.py

Covers idempotency (re-run replaces rows) and concurrent-write safety
(two processes writing different output_dirs at the same time both land
in the database without data loss). Tests cover both pipeline_runs and
run_metadata tables.
"""

import sqlite3
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pytest

import update_db as m


# ── helpers — pipeline_runs ───────────────────────────────────────────────────

def _make_rows(output_dir, n=2):
    """Return n minimal row tuples ordered by COLS."""
    base = {c: "" for c in m.COLS}
    rows = []
    for i in range(n):
        r = dict(base)
        r["output_dir"] = output_dir
        r["sample"] = f"s{i}"
        r["run_date"] = "2026-01-01T00:00:00"
        rows.append(tuple(r[c] for c in m.COLS))
    return rows


def _row_count(db_path):
    con = sqlite3.connect(str(db_path))
    n = con.execute("SELECT COUNT(*) FROM pipeline_runs").fetchone()[0]
    con.close()
    return n


def _rows_for(db_path, output_dir):
    con = sqlite3.connect(str(db_path))
    rows = con.execute(
        "SELECT * FROM pipeline_runs WHERE output_dir = ?", (output_dir,)
    ).fetchall()
    con.close()
    return rows


# ── helpers — run_metadata ────────────────────────────────────────────────────

def _make_meta_rows(output_dir, n=2, author="testuser"):
    """Return n minimal meta row tuples ordered by META_COLS."""
    base = {c: "" for c in m.META_COLS}
    rows = []
    for i in range(n):
        r = dict(base)
        r["output_dir"] = output_dir
        r["sample"]     = f"s{i}"
        r["author"]     = author
        r["run_date"]   = "2026-01-01T00:00:00"
        rows.append(tuple(r[c] for c in m.META_COLS))
    return rows


def _meta_row_count(db_path):
    con = sqlite3.connect(str(db_path))
    n = con.execute("SELECT COUNT(*) FROM run_metadata").fetchone()[0]
    con.close()
    return n


def _meta_rows_for(db_path, output_dir):
    con = sqlite3.connect(str(db_path))
    rows = con.execute(
        "SELECT * FROM run_metadata WHERE output_dir = ?", (output_dir,)
    ).fetchall()
    con.close()
    return rows


# ── idempotency — pipeline_runs ───────────────────────────────────────────────

def test_first_write_creates_rows(tmp_path):
    db = tmp_path / "pipeline_db.db"
    m.write_rows(db, "/out/A", _make_rows("/out/A", n=3))
    assert _row_count(db) == 3


def test_rerun_replaces_rows(tmp_path):
    db = tmp_path / "pipeline_db.db"
    m.write_rows(db, "/out/A", _make_rows("/out/A", n=3))
    m.write_rows(db, "/out/A", _make_rows("/out/A", n=2))
    assert _row_count(db) == 2


def test_different_output_dirs_accumulate(tmp_path):
    db = tmp_path / "pipeline_db.db"
    m.write_rows(db, "/out/A", _make_rows("/out/A", n=2))
    m.write_rows(db, "/out/B", _make_rows("/out/B", n=3))
    assert _row_count(db) == 5
    assert len(_rows_for(db, "/out/A")) == 2
    assert len(_rows_for(db, "/out/B")) == 3


def test_rerun_does_not_affect_other_dirs(tmp_path):
    db = tmp_path / "pipeline_db.db"
    m.write_rows(db, "/out/A", _make_rows("/out/A", n=2))
    m.write_rows(db, "/out/B", _make_rows("/out/B", n=3))
    m.write_rows(db, "/out/A", _make_rows("/out/A", n=1))  # re-run A
    assert len(_rows_for(db, "/out/A")) == 1
    assert len(_rows_for(db, "/out/B")) == 3              # B untouched


# ── idempotency — run_metadata ────────────────────────────────────────────────

def test_meta_first_write_creates_rows(tmp_path):
    db = tmp_path / "pipeline_db.db"
    m.write_meta_rows(db, "/out/A", _make_meta_rows("/out/A", n=3))
    assert _meta_row_count(db) == 3


def test_meta_rerun_replaces_rows(tmp_path):
    db = tmp_path / "pipeline_db.db"
    m.write_meta_rows(db, "/out/A", _make_meta_rows("/out/A", n=3))
    m.write_meta_rows(db, "/out/A", _make_meta_rows("/out/A", n=2))
    assert _meta_row_count(db) == 2


def test_meta_different_output_dirs_accumulate(tmp_path):
    db = tmp_path / "pipeline_db.db"
    m.write_meta_rows(db, "/out/A", _make_meta_rows("/out/A", n=2))
    m.write_meta_rows(db, "/out/B", _make_meta_rows("/out/B", n=3))
    assert _meta_row_count(db) == 5
    assert len(_meta_rows_for(db, "/out/A")) == 2
    assert len(_meta_rows_for(db, "/out/B")) == 3


def test_meta_rerun_does_not_affect_other_dirs(tmp_path):
    db = tmp_path / "pipeline_db.db"
    m.write_meta_rows(db, "/out/A", _make_meta_rows("/out/A", n=2))
    m.write_meta_rows(db, "/out/B", _make_meta_rows("/out/B", n=3))
    m.write_meta_rows(db, "/out/A", _make_meta_rows("/out/A", n=1))
    assert len(_meta_rows_for(db, "/out/A")) == 1
    assert len(_meta_rows_for(db, "/out/B")) == 3


# ── path builder unit tests ───────────────────────────────────────────────────

# peaks_filt_narrowpeak is intentionally excluded: it is empty from
# _build_meta_paths and filled dynamically in main() from meme_foldch_level.
TREATMENT_ONLY_COLS = [
    "peaks_narrowpeak", "summits_bed",
    "summits_fasta", "peaks_fasta",
    "meme_summits_dir", "meme_peaks_dir",
    "fimo_summits_tsv", "fimo_peaks_tsv",
    "homer_annotations",
]


def test_peaks_filt_narrowpeak_empty_from_builder():
    """_build_meta_paths leaves peaks_filt_narrowpeak empty for both sample
    kinds; main() fills the treatment value once meme_foldch_level is known."""
    for is_treatment in (True, False):
        paths = m._build_meta_paths("/out/test", "s1", is_treatment=is_treatment)
        assert paths["peaks_filt_narrowpeak"] == ""


def test_meta_control_paths_are_empty():
    paths = m._build_meta_paths("/out/test", "ctrl", is_treatment=False)
    for col in TREATMENT_ONLY_COLS:
        assert paths[col] == "", f"{col} should be empty for control"


def test_meta_treatment_paths_are_nonempty():
    paths = m._build_meta_paths("/out/test", "treat1", is_treatment=True)
    for col in TREATMENT_ONLY_COLS:
        assert paths[col] != "", f"{col} should be non-empty for treatment"


def test_meta_shared_paths_always_present():
    for is_treatment in (True, False):
        paths = m._build_meta_paths("/out/test", "s1", is_treatment=is_treatment)
        for col in ["bam", "bigwig", "trimmed_r1", "trimmed_r2"]:
            assert paths[col] != "", f"{col} should always be set (is_treatment={is_treatment})"


def test_meta_paths_contain_sample_and_output_dir():
    paths = m._build_meta_paths("/my/output", "sample_A", is_treatment=True)
    assert "/my/output" in paths["bam"]
    assert "sample_A" in paths["bam"]
    assert "/my/output" in paths["peaks_narrowpeak"]
    assert "sample_A" in paths["peaks_narrowpeak"]


# ── concurrent writes ─────────────────────────────────────────────────────────

def _worker(db_path_str, output_dir, n):
    """Top-level function so ProcessPoolExecutor can pickle it."""
    import update_db as m  # re-import in subprocess
    m.write_rows(db_path_str, output_dir, m._make_rows_for_test(output_dir, n))
    m.write_meta_rows(db_path_str, output_dir, m._make_meta_rows_for_test(output_dir, n))
    return output_dir


# Expose picklable helpers for the worker (ProcessPoolExecutor needs top-level fns)
def _make_rows_for_test(output_dir, n):
    return _make_rows(output_dir, n)


def _make_meta_rows_for_test(output_dir, n):
    return _make_meta_rows(output_dir, n)


m._make_rows_for_test = _make_rows_for_test          # monkey-patch onto module for worker
m._make_meta_rows_for_test = _make_meta_rows_for_test


def test_concurrent_writes_no_data_loss(tmp_path):
    """Two processes writing different output_dirs simultaneously must both land."""
    db = tmp_path / "pipeline_db.db"

    # Pre-create both tables so workers don't race on CREATE TABLE
    m.write_rows(db, "/init", [])
    m.write_meta_rows(db, "/init", [])

    tasks = {"/out/P": 2, "/out/Q": 3}
    with ProcessPoolExecutor(max_workers=2) as ex:
        futures = [
            ex.submit(_worker, str(db), od, n)
            for od, n in tasks.items()
        ]
        results = [f.result() for f in as_completed(futures)]

    assert set(results) == {"/out/P", "/out/Q"}

    # pipeline_runs
    assert len(_rows_for(db, "/out/P")) == 2
    assert len(_rows_for(db, "/out/Q")) == 3
    assert _row_count(db) == 5  # init rows were 0

    # run_metadata
    assert len(_meta_rows_for(db, "/out/P")) == 2
    assert len(_meta_rows_for(db, "/out/Q")) == 3
    assert _meta_row_count(db) == 5


# ── simplified columns + experiment metadata ──────────────────────────────────

def test_cols_use_single_filtered_set():
    assert "num_peaks_filt" in m.COLS
    assert "reads_in_peaks_filt" in m.COLS
    for legacy in (
        "num_peaks_fold1", "num_peaks_fold2", "num_peaks_fold3",
        "reads_in_peaks_fold1", "reads_in_peaks_fold2", "reads_in_peaks_fold3",
    ):
        assert legacy not in m.COLS


def test_meta_cols_include_experiment_metadata():
    assert "experiment_date" in m.META_COLS
    assert "gdna_batch" in m.META_COLS


def test_experiment_metadata_lands_in_run_metadata(tmp_path):
    db = tmp_path / "pipeline_db.db"
    base = {c: "" for c in m.META_COLS}
    base.update({
        "output_dir": "/out/A", "sample": "s0",
        "experiment_date": "2026-07-14", "gdna_batch": "batch-07",
    })
    m.write_meta_rows(db, "/out/A", [tuple(base[c] for c in m.META_COLS)])

    con = sqlite3.connect(str(db))
    row = dict(zip(
        [d[0] for d in con.execute("SELECT * FROM run_metadata").description],
        con.execute("SELECT * FROM run_metadata").fetchone(),
    ))
    con.close()
    assert row["experiment_date"] == "2026-07-14"
    assert row["gdna_batch"] == "batch-07"


def test_write_rows_migrates_legacy_schema(tmp_path):
    """A DB whose pipeline_runs predates the new columns must accept new rows.

    _ensure_columns adds the missing columns; legacy fold columns are left in
    place (unpopulated) — no destructive migration.
    """
    db = tmp_path / "pipeline_db.db"
    con = sqlite3.connect(str(db))
    con.execute(
        'CREATE TABLE pipeline_runs '
        '(id INTEGER PRIMARY KEY AUTOINCREMENT, '
        '"output_dir" TEXT, "sample" TEXT, "reads_in_peaks_fold1" TEXT)'
    )
    con.commit()
    con.close()

    m.write_rows(db, "/out/A", _make_rows("/out/A", n=2))

    con = sqlite3.connect(str(db))
    cols = {r[1] for r in con.execute("PRAGMA table_info(pipeline_runs)")}
    con.close()
    assert "num_peaks_filt" in cols            # newly added
    assert "reads_in_peaks_filt" in cols       # newly added
    assert "reads_in_peaks_fold1" in cols      # legacy column preserved
    assert len(_rows_for(db, "/out/A")) == 2
