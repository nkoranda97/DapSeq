"""
Append per-sample pipeline run metadata to a shared SQLite database.

Called via Snakemake's script: directive; uses the snakemake object for I/O.
Each run writes one row per sample. Re-running the same output_dir replaces
its rows (idempotent).

Concurrency note: SQLite uses fcntl POSIX record locks internally, which work
on Lustre/GPFS via NLM. The 60-second timeout serialises the rare case of two
runs finishing simultaneously without failing. WAL mode is intentionally NOT
used — WAL's shared-memory coordination requires both writers to be on the
same physical node, which SLURM does not guarantee.
"""

import csv
import sqlite3
from datetime import datetime
from pathlib import Path


COLS = [
    "run_date",
    "output_dir",
    "genome_ref",
    "genome_size",
    "input_control",
    "threads",
    "mapq",
    "max_frags",
    "macs3_format",
    "macs3_min_foldch",
    "meme_nmotifs",
    "meme_minw",
    "meme_maxw",
    "meme_maxpeaks",
    "fimo_thresh",
    "sample",
    "r1",
    "r2",
    "is_treatment",
    "total_frags",
    "clean_reads",
    "align_rate%",
    "filtered_reads",
    "peak#",
    "min5fold_peak#",
    "max_peak_score",
    "peak_reads#",
    "FRiP_score",
]

_col_defs  = ", ".join(f'"{c}" TEXT' for c in COLS)
_col_names = ", ".join(f'"{c}"' for c in COLS)
_col_ph    = ", ".join("?" * len(COLS))

_CREATE = (
    f"CREATE TABLE IF NOT EXISTS pipeline_runs "
    f"(id INTEGER PRIMARY KEY AUTOINCREMENT, {_col_defs})"
)
_INSERT = f'INSERT INTO pipeline_runs ({_col_names}) VALUES ({_col_ph})'


def read_report(path):
    """Return a dict of {sample: {col: value}} from the report TSV."""
    result = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            result[row["sample"]] = dict(row)
    return result


def get_r1(sample_cfg):
    r1 = sample_cfg.get("r1")
    if r1 is None:
        return None
    return r1[0] if isinstance(r1, list) else r1


def get_r2(sample_cfg):
    r2 = sample_cfg.get("r2")
    if r2 is None:
        return None
    return r2[0] if isinstance(r2, list) else r2


def write_rows(db_path, output_dir, rows):
    """Write rows (list of tuples ordered by COLS) to the SQLite database.

    Replaces all existing rows for output_dir atomically. Safe to call
    concurrently from multiple processes on shared filesystems.
    """
    con = sqlite3.connect(str(db_path), timeout=60)
    con.execute("PRAGMA journal_mode=DELETE")
    con.execute(_CREATE)
    with con:
        con.execute("DELETE FROM pipeline_runs WHERE output_dir = ?", (output_dir,))
        con.executemany(_INSERT, rows)
    con.close()


def main():
    sm = snakemake  # noqa: F821 — injected by Snakemake

    db_path          = sm.params.db_path
    output_dir       = sm.params.output_dir
    samples_cfg      = sm.params.samples_cfg
    treatment_set    = set(sm.params.treatment_samples)
    run_date         = datetime.now().isoformat(timespec="seconds")

    qc_stats = read_report(sm.input.report)

    shared = {
        "run_date":          run_date,
        "output_dir":        output_dir,
        "genome_ref":        sm.params.genome_ref,
        "genome_size":       sm.params.genome_size,
        "input_control":     sm.params.input_control or "",
        "threads":           sm.params.threads,
        "mapq":              sm.params.mapq,
        "max_frags":         sm.params.max_frags or "",
        "macs3_format":      sm.params.macs3_format,
        "macs3_min_foldch":  sm.params.macs3_min_foldch,
        "meme_nmotifs":      sm.params.meme_nmotifs,
        "meme_minw":         sm.params.meme_minw,
        "meme_maxw":         sm.params.meme_maxw,
        "meme_maxpeaks":     sm.params.meme_maxpeaks,
        "fimo_thresh":       sm.params.fimo_thresh,
    }

    new_rows = []
    for sample, scfg in samples_cfg.items():
        r1 = get_r1(scfg)
        if r1 is None:
            continue

        stats = qc_stats.get(sample, {})
        row = dict(shared)
        row["sample"]       = sample
        row["r1"]           = r1
        row["r2"]           = get_r2(scfg) or ""
        row["is_treatment"] = sample in treatment_set

        row["total_frags"]    = stats.get("total_frags", "NA")
        row["clean_reads"]    = stats.get("clean_reads", "NA")
        row["align_rate%"]    = stats.get("align_rate%", "NA")
        row["filtered_reads"] = stats.get("filtered_reads", "NA")
        row["peak#"]          = stats.get("peak#", "NA")
        row["min5fold_peak#"] = stats.get("min5fold peak#", "NA")
        row["max_peak_score"] = stats.get("max peak score", "NA")
        row["peak_reads#"]    = stats.get("peak reads#", "NA")
        row["FRiP_score"]     = stats.get("FRiP_score", "NA")

        new_rows.append(tuple(row.get(c, "") for c in COLS))

    write_rows(db_path, output_dir, new_rows)

    Path(sm.output.flag).touch()


if "snakemake" in dir():
    main()
