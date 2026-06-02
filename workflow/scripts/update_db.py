"""
Append per-sample pipeline run metadata to a shared SQLite database.

Called via Snakemake's script: directive; uses the snakemake object for I/O.
Each run writes one row per sample to two tables:
  - pipeline_runs  : QC metrics and parameters
  - run_metadata   : file paths and provenance (author, dates, all I/O paths)

Both tables are joined via (output_dir, sample). Re-running the same
output_dir replaces its rows in both tables (idempotent).

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
    "total_reads",
    "trimmed_reads",
    "mapped_reads",
    "alignment_rate",
    "reads_in_peaks",
    "reads_5fold",
    "reads_nfold",
    "max_peak_score",
    "mapping_pct",
    "motif_peaks",
    "subsampled_frags",
    "median_frag_size",
    "num_peaks",
    "num_peaks_filt",
]

META_COLS = [
    # join keys
    "output_dir",
    "sample",
    # provenance
    "author",
    "run_date",
    # reference paths (from config)
    "genome_ref",
    "gene_annotation",
    # raw input paths
    "r1",
    "r2",
    # generated paths — trimming
    "trimmed_r1",
    "trimmed_r2",
    # generated paths — alignment
    "bam",
    "bigwig",
    # generated paths — peaks (treatment only)
    "peaks_narrowpeak",
    "peaks_filt_narrowpeak",
    "summits_bed",
    # generated paths — fasta (treatment only)
    "summits_fasta",
    "peaks_fasta",
    # generated paths — meme (treatment only)
    "meme_summits_dir",
    "meme_peaks_dir",
    # generated paths — fimo (treatment only)
    "fimo_summits_tsv",
    "fimo_peaks_tsv",
    # generated paths — homer (treatment only, if gene_annotation set)
    "homer_annotations",
    # generated paths — report
    "report_csv",
]

_col_defs  = ", ".join(f'"{c}" TEXT' for c in COLS)
_col_names = ", ".join(f'"{c}"' for c in COLS)
_col_ph    = ", ".join("?" * len(COLS))

_CREATE = (
    f"CREATE TABLE IF NOT EXISTS pipeline_runs "
    f"(id INTEGER PRIMARY KEY AUTOINCREMENT, {_col_defs})"
)
_INSERT = f'INSERT INTO pipeline_runs ({_col_names}) VALUES ({_col_ph})'

_meta_col_defs  = ", ".join(f'"{c}" TEXT' for c in META_COLS)
_meta_col_names = ", ".join(f'"{c}"' for c in META_COLS)
_meta_col_ph    = ", ".join("?" * len(META_COLS))

_CREATE_META = (
    f"CREATE TABLE IF NOT EXISTS run_metadata "
    f"(id INTEGER PRIMARY KEY AUTOINCREMENT, {_meta_col_defs})"
)
_INSERT_META = f'INSERT INTO run_metadata ({_meta_col_names}) VALUES ({_meta_col_ph})'


def _ensure_columns(con, table, cols):
    """Add any columns in *cols* that are missing from *table*.

    ALTER TABLE must run outside a transaction (SQLite DDL implicitly commits),
    so callers must invoke this before opening a 'with con:' block.
    """
    existing = {row[1] for row in con.execute(f"PRAGMA table_info({table})")}
    for col in cols:
        if col not in existing:
            con.execute(f'ALTER TABLE "{table}" ADD COLUMN "{col}" TEXT')


def read_report(path):
    """Return a dict of {sample: {col: value}} from the report CSV."""
    result = {}
    with open(path) as fh:
        reader = csv.DictReader(fh)
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
    _ensure_columns(con, "pipeline_runs", COLS)  # DDL outside transaction
    with con:
        con.execute("DELETE FROM pipeline_runs WHERE output_dir = ?", (output_dir,))
        con.executemany(_INSERT, rows)
    con.close()


def write_meta_rows(db_path, output_dir, rows):
    """Write rows (list of tuples ordered by META_COLS) to run_metadata.

    Same idempotency and concurrency guarantees as write_rows().
    """
    con = sqlite3.connect(str(db_path), timeout=60)
    con.execute("PRAGMA journal_mode=DELETE")
    con.execute(_CREATE_META)
    with con:
        con.execute("DELETE FROM run_metadata WHERE output_dir = ?", (output_dir,))
        con.executemany(_INSERT_META, rows)
    con.close()


def _build_meta_paths(output_dir, sample, is_treatment):
    """Return dict of generated file paths for one sample.

    Treatment-only steps (peaks, MEME, FIMO, HOMER) are empty strings
    for control samples since those rules are not run on controls.
    trimmed_r2 is recorded as a path string even for SE samples — the file
    may not exist, but the expected path is stored for reference.
    """
    o = output_dir.rstrip("/")
    paths = {
        "trimmed_r1": f"{o}/trimmed/{sample}.R1.fastq.gz",
        "trimmed_r2": f"{o}/trimmed/{sample}.R2.fastq.gz",
        "bam":    f"{o}/bam/{sample}.bam",
        "bigwig": f"{o}/bigWig/{sample}.bw",
    }
    if is_treatment:
        paths.update({
            "peaks_narrowpeak":      f"{o}/MACS/{sample}_peaks.narrowPeak",
            "peaks_filt_narrowpeak": f"{o}/MACS/{sample}_peaks_filt.narrowPeak",
            "summits_bed":           f"{o}/MACS/{sample}_summits.bed",
            "summits_fasta":         f"{o}/fasta/{sample}.summits.fasta",
            "peaks_fasta":           f"{o}/fasta/{sample}.peaks.fasta",
            "meme_summits_dir":      f"{o}/meme/{sample}/summits",
            "meme_peaks_dir":        f"{o}/meme/{sample}/peaks",
            "fimo_summits_tsv":      f"{o}/fimo/{sample}/summits/fimo.tsv",
            "fimo_peaks_tsv":        f"{o}/fimo/{sample}/peaks/fimo.tsv",
            "homer_annotations":     f"{o}/annotations/{sample}.peak_annotations.txt",
        })
    else:
        paths.update({k: "" for k in [
            "peaks_narrowpeak", "peaks_filt_narrowpeak", "summits_bed",
            "summits_fasta", "peaks_fasta",
            "meme_summits_dir", "meme_peaks_dir",
            "fimo_summits_tsv", "fimo_peaks_tsv",
            "homer_annotations",
        ]})
    paths["report_csv"] = f"{o}/stats/report.csv"
    return paths


def main():
    sm = snakemake  # noqa: F821 — injected by Snakemake

    db_path          = sm.params.db_path
    output_dir       = sm.params.output_dir
    samples_cfg      = sm.params.samples_cfg
    treatment_set    = set(sm.params.treatment_samples)
    run_date         = datetime.now().isoformat(timespec="seconds")
    author           = sm.params.author
    gene_annotation  = sm.params.gene_annotation or ""

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
    meta_rows = []
    for sample, scfg in samples_cfg.items():
        r1 = get_r1(scfg)
        if r1 is None:
            continue

        is_treatment = sample in treatment_set

        stats = qc_stats.get(sample, {})
        row = dict(shared)
        row["sample"]       = sample
        row["r1"]           = r1
        row["r2"]           = get_r2(scfg) or ""
        row["is_treatment"] = is_treatment
        row["total_reads"]   = stats.get("total_reads", "NA")
        row["trimmed_reads"] = stats.get("trimmed_reads", "NA")
        row["mapped_reads"]  = stats.get("mapped_reads", "NA")
        row["alignment_rate"]       = stats.get("alignment_rate", "NA")
        row["reads_in_peaks"]  = stats.get("reads_in_peaks", "NA")
        row["reads_5fold"]     = stats.get("reads_5fold", "NA")
        row["reads_nfold"]     = stats.get("reads_nfold", "NA")
        row["max_peak_score"]    = stats.get("max_peak_score", "NA")
        row["mapping_pct"]       = stats.get("mapping_pct", "NA")
        row["motif_peaks"]       = stats.get("motif_peaks", "NA")
        row["subsampled_frags"]  = stats.get("subsampled_frags", "NA")
        row["median_frag_size"]  = stats.get("median_frag_size", "NA")
        row["num_peaks"]         = stats.get("num_peaks", "NA")
        row["num_peaks_filt"]    = stats.get("num_peaks_filt", "NA")

        new_rows.append(tuple(row.get(c, "") for c in COLS))

        generated = _build_meta_paths(output_dir, sample, is_treatment)
        meta = {
            "output_dir":      output_dir,
            "sample":          sample,
            "author":          author,
            "run_date":        run_date,
            "genome_ref":      sm.params.genome_ref,
            "gene_annotation": gene_annotation,
            "r1":              r1,
            "r2":              get_r2(scfg) or "",
        }
        meta.update(generated)
        meta_rows.append(tuple(meta.get(c, "") for c in META_COLS))

    write_rows(db_path, output_dir, new_rows)
    write_meta_rows(db_path, output_dir, meta_rows)

    Path(sm.output.flag).touch()


if "snakemake" in dir():
    main()
