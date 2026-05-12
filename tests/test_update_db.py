"""
Tests for workflow/scripts/update_db.py

Covers idempotency (re-run replaces rows) and concurrent-write safety
(two processes writing different output_dirs at the same time both land
in the database without data loss).
"""

import sqlite3
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pytest

import update_db as m


# ── helpers ───────────────────────────────────────────────────────────────────

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


# ── idempotency ───────────────────────────────────────────────────────────────

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


# ── concurrent writes ─────────────────────────────────────────────────────────

def _worker(db_path_str, output_dir, n):
    """Top-level function so ProcessPoolExecutor can pickle it."""
    import update_db as m  # re-import in subprocess
    m.write_rows(db_path_str, output_dir, m._make_rows_for_test(output_dir, n))
    return output_dir


# Expose a picklable helper for the worker (ProcessPoolExecutor needs top-level fns)
def _make_rows_for_test(output_dir, n):
    return _make_rows(output_dir, n)


m._make_rows_for_test = _make_rows_for_test  # monkey-patch onto module for worker


def test_concurrent_writes_no_data_loss(tmp_path):
    """Two processes writing different output_dirs simultaneously must both land."""
    db = tmp_path / "pipeline_db.db"

    # Pre-create the DB so workers don't race on CREATE TABLE
    m.write_rows(db, "/init", [])

    tasks = {"/out/P": 2, "/out/Q": 3}
    with ProcessPoolExecutor(max_workers=2) as ex:
        futures = [
            ex.submit(_worker, str(db), od, n)
            for od, n in tasks.items()
        ]
        results = [f.result() for f in as_completed(futures)]

    assert set(results) == {"/out/P", "/out/Q"}
    assert len(_rows_for(db, "/out/P")) == 2
    assert len(_rows_for(db, "/out/Q")) == 3
    assert _row_count(db) == 5  # init rows were 0
