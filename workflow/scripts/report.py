"""
Aggregate per-sample QC statistics into a CSV summary table.

One row per treatment sample. Peak stats are reported for both bam types
(full = all contigs, chr = canonical chromosomes). FRiP always uses
reads_on_chromosomes as denominator. Called via Snakemake's script: directive
or directly from the CLI.
"""

import argparse
import csv
import os
import re
import sys


def parse_kv_file(path):
    """Read a key<TAB>value file, skipping blank lines and '#' comments."""
    result = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t", 1)
            if len(parts) == 2:
                result[parts[0]] = parts[1]
    return result


def parse_bbduk_log(path):
    """Return output read count from a bbduk trim log, or 'NA'."""
    with open(path) as fh:
        content = fh.read()
    m = re.search(r"^Result:\s+(\d+)\s+reads", content, re.MULTILINE)
    return int(m.group(1)) if m else "NA"


def parse_narrowpeak_max(path):
    """Return max signalValue (column 6, 0-indexed) from a narrowPeak file, or 'NA'."""
    max_signal = None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            fields = line.split("\t")
            signal = float(fields[6])
            if max_signal is None or signal > max_signal:
                max_signal = signal
    return round(max_signal, 4) if max_signal is not None else "NA"


def parse_fimo(path):
    """
    Count unique peak names (sequence_name column) in a FIMO TSV.
    FIMO output has a header line and optional '#' comment lines.
    Returns integer count or 'NA' if file is missing/empty.
    """
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return "NA"
    seen = set()
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if parts[0] == "motif_id":
                continue
            if len(parts) >= 3:
                seen.add(parts[2])
    return len(seen) if seen else "NA"


def make_cols():
    return [
        # ── pre-processing ──────────────────────────────────────────────
        "sample",
        "reads_original",
        "reads_trimmed",
        # ── alignment ───────────────────────────────────────────────────
        "reads_mapped",           # full.bam (all contigs, MAPQ-filtered)
        "reads_on_chromosomes",   # chr.bam (canonical chromosomes only)
        # ── full.bam peak stats (pipeline order) ────────────────────────
        "reads_in_peaks_full",
        "max_peak_score_full",
        "reads_5fold_full",
        "reads_nfold_full",
        "mapping_pct_full",       # reads_in_peaks_full / reads_on_chromosomes × 100
        "motif_peaks_full",
        # ── chr.bam peak stats (pipeline order) ─────────────────────────
        "reads_in_peaks_chr",
        "max_peak_score_chr",
        "reads_5fold_chr",
        "reads_nfold_chr",
        "mapping_pct_chr",        # reads_in_peaks_chr / reads_on_chromosomes × 100
        "motif_peaks_chr",
    ]


def _safe_frip(reads_in_peaks, reads_chr):
    if reads_in_peaks == "NA" or reads_chr == "NA":
        return "NA"
    denom = int(reads_chr)
    return round(int(reads_in_peaks) / denom * 100, 2) if denom > 0 else "NA"


def build_row(sample, trim_log_dir, stats_dir, macs_dir, fimo_dir):
    row = {"sample": sample}

    total_frags_path = os.path.join(stats_dir, f"{sample}.total_frags.txt")
    row["reads_original"] = open(total_frags_path).read().strip() if os.path.exists(total_frags_path) else "NA"

    trim_log_path = os.path.join(trim_log_dir, f"{sample}.trim.log")
    row["reads_trimmed"] = parse_bbduk_log(trim_log_path) if os.path.exists(trim_log_path) else "NA"

    full_stats = parse_kv_file(os.path.join(stats_dir, f"{sample}.full.filtered_stats.txt"))
    row["reads_mapped"] = full_stats.get("filtered_reads", "NA")

    chr_stats = parse_kv_file(os.path.join(stats_dir, f"{sample}.chr.filtered_stats.txt"))
    reads_chr = chr_stats.get("filtered_reads", "NA")
    row["reads_on_chromosomes"] = reads_chr

    for bt in ("full", "chr"):
        frip_data = parse_kv_file(os.path.join(stats_dir, f"{sample}.{bt}.frip_macs.txt"))
        reads_in_peaks = frip_data.get("reads_in_peaks_macs", "NA")
        row[f"reads_in_peaks_{bt}"] = reads_in_peaks

        frip_5fold_data = parse_kv_file(os.path.join(stats_dir, f"{sample}.{bt}.frip_macs_5fold.txt"))
        row[f"reads_5fold_{bt}"] = frip_5fold_data.get("reads_in_peaks_macs_5fold", "NA")

        frip_filt_data = parse_kv_file(os.path.join(stats_dir, f"{sample}.{bt}.frip_macs_filt.txt"))
        row[f"reads_nfold_{bt}"] = frip_filt_data.get("reads_in_peaks_macs_filt", "NA")

        np_path = os.path.join(macs_dir, f"{sample}.{bt}_peaks.narrowPeak")
        row[f"max_peak_score_{bt}"] = parse_narrowpeak_max(np_path)

        row[f"mapping_pct_{bt}"] = _safe_frip(reads_in_peaks, reads_chr)

        fimo_path = os.path.join(fimo_dir, f"{sample}.{bt}", "peaks", "fimo.tsv")
        row[f"motif_peaks_{bt}"] = parse_fimo(fimo_path)

    return row


def run(treatment_samples, trim_log_dir, stats_dir, macs_dir, fimo_dir, csv_out):
    cols = make_cols()
    rows = [build_row(s, trim_log_dir, stats_dir, macs_dir, fimo_dir) for s in treatment_samples]

    with open(csv_out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=cols)
        writer.writeheader()
        for row in rows:
            writer.writerow({c: row.get(c, "NA") for c in cols})


def main():
    sm = snakemake  # noqa: F821 — injected by Snakemake
    run(
        treatment_samples = list(sm.params.treatment_samples),
        trim_log_dir      = sm.params.trim_log_dir,
        macs_dir          = sm.params.macs_dir,
        stats_dir         = sm.params.stats_dir,
        fimo_dir          = sm.params.fimo_dir,
        csv_out           = sm.output.csv,
    )


def cli_main():
    parser = argparse.ArgumentParser(
        description="Regenerate the DAP-seq stats CSV from completed pipeline output."
    )
    parser.add_argument("out_dir", help="Pipeline output directory (contains stats/, MACS/, fimo/, logs/bbduk/)")
    args = parser.parse_args()

    out_dir      = os.path.abspath(args.out_dir)
    stats_dir    = os.path.join(out_dir, "stats")
    macs_dir     = os.path.join(out_dir, "MACS")
    fimo_dir     = os.path.join(out_dir, "fimo")
    trim_log_dir = os.path.join(out_dir, "logs", "bbduk")
    csv_out      = os.path.join(stats_dir, "report.csv")

    samples = set()
    for fname in os.listdir(stats_dir):
        if fname.endswith(".chr.frip_macs.txt"):
            samples.add(fname[: -len(".chr.frip_macs.txt")])
    if not samples:
        sys.exit(f"No *.chr.frip_macs.txt files found in {stats_dir} — is the pipeline complete?")

    run(
        treatment_samples = sorted(samples),
        trim_log_dir      = trim_log_dir,
        stats_dir         = stats_dir,
        macs_dir          = macs_dir,
        fimo_dir          = fimo_dir,
        csv_out           = csv_out,
    )
    print(f"Wrote {csv_out}")


if "snakemake" in dir():
    main()
elif __name__ == "__main__":
    cli_main()
