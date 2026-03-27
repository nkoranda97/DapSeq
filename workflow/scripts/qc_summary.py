"""
Aggregate per-sample QC statistics into a single TSV summary table.

Called via Snakemake's script: directive; uses the snakemake object for I/O.
"""

import os
import re


def parse_kv_file(path):
    """
    Read a key<TAB>value file, skipping blank lines and lines starting with '#'.
    Returns a dict of {key: value}.
    """
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
    """
    Parse a bbduk trim log and return the output read count (clean reads after trimming).
    """
    with open(path) as fh:
        content = fh.read()
    m = re.search(r"^Result:\s+(\d+)\s+reads", content, re.MULTILINE)
    return int(m.group(1)) if m else "NA"


def parse_narrowpeak(path):
    """
    Return (total_peaks, min5fold_peaks, max_signal) from a narrowPeak file.
    Column 7 (0-indexed: 6) is signalValue (fold enrichment).
    """
    total = 0
    min5 = 0
    max_signal = None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            fields = line.split("\t")
            signal = float(fields[6])
            total += 1
            if signal >= 5:
                min5 += 1
            if max_signal is None or signal > max_signal:
                max_signal = signal
    if max_signal is None:
        return "NA", "NA", "NA"
    return total, min5, round(max_signal, 4)


def parse_top_motif(meme_txt_path):
    """Return the consensus sequence of the top MEME motif, or 'NA'."""
    with open(meme_txt_path) as fh:
        for line in fh:
            if line.startswith("MOTIF "):
                return line.split()[1]
    return "NA"


COLS = [
    "sample",
    "total_frags",
    "clean_reads",
    "align_rate%",
    "filtered_reads",
    "peak#",
    "min5fold peak#",
    "max peak score",
    "peak reads#",
    "FRiP_score",
    "motif",
]


def build_row(sample, is_treatment, trim_log_dir, stats_dir, macs_dir, meme_dir):
    total_frags_path = os.path.join(stats_dir, f"{sample}.total_frags.txt")
    total_frags = open(total_frags_path).read().strip() if os.path.exists(total_frags_path) else "NA"

    trim_log_path = os.path.join(trim_log_dir, f"{sample}.trim.log")
    clean_reads = parse_bbduk_log(trim_log_path) if os.path.exists(trim_log_path) else "NA"

    align_rate_path = os.path.join(stats_dir, f"{sample}.align_rate.txt")
    align_rate = open(align_rate_path).read().strip() if os.path.exists(align_rate_path) else "NA"

    filt_stats = parse_kv_file(os.path.join(stats_dir, f"{sample}.filtered_stats.txt"))
    filtered_reads = filt_stats.get("filtered_reads", "NA")

    if is_treatment:
        np_path = os.path.join(macs_dir, f"{sample}_peaks_filt.narrowPeak")
        peak_total, min5fold, max_score = parse_narrowpeak(np_path)

        frip_macs = parse_kv_file(os.path.join(stats_dir, f"{sample}.frip_macs.txt"))
        peak_reads = frip_macs.get("reads_in_peaks_macs", "NA")
        if peak_reads != "NA" and filtered_reads != "NA":
            frip = round(int(peak_reads) / int(filtered_reads) * 100, 2)
        else:
            frip = "NA"

        meme_path = os.path.join(meme_dir, sample, "summits", "meme.txt")
        motif = parse_top_motif(meme_path) if os.path.exists(meme_path) else "NA"
    else:
        peak_total = min5fold = max_score = peak_reads = frip = motif = "NA"

    return {
        "sample":          sample,
        "total_frags":     total_frags,
        "clean_reads":     clean_reads,
        "align_rate%":     align_rate,
        "filtered_reads":  filtered_reads,
        "peak#":           peak_total,
        "min5fold peak#":  min5fold,
        "max peak score":  max_score,
        "peak reads#":     peak_reads,
        "FRiP_score":      frip,
        "motif":           motif,
    }


def main():
    sm = snakemake  # noqa: F821 — injected by Snakemake

    samples           = list(sm.params.samples)
    treatment_samples = set(sm.params.treatment_samples)
    trim_log_dir      = sm.params.trim_log_dir
    macs_dir          = sm.params.macs_dir
    meme_dir          = sm.params.meme_dir
    stats_dir         = sm.params.stats_dir

    rows = [
        build_row(s, s in treatment_samples, trim_log_dir, stats_dir, macs_dir, meme_dir)
        for s in samples
    ]

    with open(sm.output[0], "w") as out:
        out.write("\t".join(COLS) + "\n")
        for row in rows:
            out.write("\t".join(str(row.get(c, "NA")) for c in COLS) + "\n")


if "snakemake" in dir():
    main()
