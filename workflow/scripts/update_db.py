"""
Append per-sample pipeline run metadata to a shared CSV database.

Called via Snakemake's script: directive; uses the snakemake object for I/O.
Each run writes one row per sample. Re-running the same output_dir replaces
its rows (idempotent).
"""

import csv
import os
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


def read_qc_summary(path):
    """Return a dict of {sample: {col: value}} from the qc_summary TSV."""
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


def main():
    sm = snakemake  # noqa: F821 — injected by Snakemake

    db_path          = sm.params.db_path
    output_dir       = sm.params.output_dir
    samples_cfg      = sm.params.samples_cfg
    treatment_set    = set(sm.params.treatment_samples)
    run_date         = datetime.now().isoformat(timespec="seconds")

    qc_stats = read_qc_summary(sm.input.qc_summary)

    # Build shared param fields
    shared = {
        "run_date":        run_date,
        "output_dir":      output_dir,
        "genome_ref":      sm.params.genome_ref,
        "genome_size":     sm.params.genome_size,
        "input_control":   sm.params.input_control or "",
        "threads":         sm.params.threads,
        "mapq":            sm.params.mapq,
        "max_frags":       sm.params.max_frags or "",
        "macs3_format":    sm.params.macs3_format,
        "macs3_min_foldch": sm.params.macs3_min_foldch,
        "meme_nmotifs":    sm.params.meme_nmotifs,
        "meme_minw":       sm.params.meme_minw,
        "meme_maxw":       sm.params.meme_maxw,
        "meme_maxpeaks":   sm.params.meme_maxpeaks,
        "fimo_thresh":     sm.params.fimo_thresh,
    }

    new_rows = []
    for sample, scfg in samples_cfg.items():
        r1 = get_r1(scfg)
        if r1 is None:
            continue  # sample has no data (template placeholder)

        stats = qc_stats.get(sample, {})
        row = dict(shared)
        row["sample"]       = sample
        row["r1"]           = r1
        row["r2"]           = get_r2(scfg) or ""
        row["is_treatment"] = sample in treatment_set

        # QC summary columns (rename to match our COLS)
        row["total_frags"]    = stats.get("total_frags", "NA")
        row["clean_reads"]    = stats.get("clean_reads", "NA")
        row["align_rate%"]    = stats.get("align_rate%", "NA")
        row["filtered_reads"] = stats.get("filtered_reads", "NA")
        row["peak#"]          = stats.get("peak#", "NA")
        row["min5fold_peak#"] = stats.get("min5fold peak#", "NA")
        row["max_peak_score"] = stats.get("max peak score", "NA")
        row["peak_reads#"]    = stats.get("peak reads#", "NA")
        row["FRiP_score"]     = stats.get("FRiP_score", "NA")

        new_rows.append(row)

    # Load existing CSV and drop rows from this output_dir (idempotent re-run)
    existing_rows = []
    db_file = Path(db_path)
    if db_file.exists():
        with open(db_file) as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                if row.get("output_dir") != output_dir:
                    existing_rows.append(row)

    all_rows = existing_rows + new_rows

    with open(db_file, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(all_rows)

    # Touch the sentinel flag
    Path(sm.output.flag).touch()


if "snakemake" in dir():
    main()
