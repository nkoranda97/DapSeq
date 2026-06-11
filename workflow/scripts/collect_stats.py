"""
Compute all per-sample QC statistics and write {sample}.stats.csv.

Called via Snakemake's script: directive. Uses snakemake.input, .output,
and .params (is_pe: bool, is_treatment: bool).
"""

import csv
import os
import re
import subprocess


NA = "NA"


# ---------------------------------------------------------------------------
# Log parsing helpers
# ---------------------------------------------------------------------------

def _parse_subsample_log(path, is_pe):
    """Return (total_reads, subsampled_frags) from a reformat.sh subsample log.

    reformat.sh logs 'Input: N reads' and 'Output: N reads', counting
    individual reads (R1+R2 combined for PE) -- the same units as
    trimmed_reads and mapped_reads. total_reads is reported as-is so all
    three columns share units. subsampled_frags represents fragment/pair
    counts (compared against the fragment-based max_frags config), so it
    is halved for PE.
    """
    total = NA
    subsampled = NA
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            m = re.match(r"^Input:\s+(\d[\d,]*)\s+reads", line)
            if m:
                n = int(m.group(1).replace(",", ""))
                total = str(n)
            m = re.match(r"^Output:\s+(\d[\d,]*)\s+reads", line)
            if m:
                n = int(m.group(1).replace(",", ""))
                subsampled = str(n // 2 if is_pe else n)
    return total, subsampled


def _parse_bbduk_log(path):
    """Return trimmed read count from a bbduk trim log ('Result: N reads')."""
    with open(path) as fh:
        content = fh.read()
    m = re.search(r"^Result:\s+(\d+)\s+reads", content, re.MULTILINE)
    return str(m.group(1)) if m else NA


# ---------------------------------------------------------------------------
# Subprocess helpers
# ---------------------------------------------------------------------------

def _safe_pct(numerator, denominator):
    """Return numerator/denominator*100 rounded to 2dp, or NA on bad inputs."""
    if numerator == NA or denominator == NA:
        return NA
    try:
        denom = int(denominator)
        return str(round(int(numerator) / denom * 100, 2)) if denom > 0 else NA
    except (ValueError, ZeroDivisionError):
        return NA


def _gate_subsampled(subsampled, max_frags):
    """Return subsampled count if max_frags is set, else NA."""
    return subsampled if max_frags else NA


def _samtools_count(bam):
    """Return filtered read count via samtools view -c -F 2308."""
    result = subprocess.run(
        ["samtools", "view", "-c", "-F", "2308", bam],
        capture_output=True, text=True, check=True,
    )
    return result.stdout.strip()


def _bedtools_intersect_count(bam, peaks):
    """Count reads overlapping peaks: bedtools intersect | samtools view -c."""
    intersect = subprocess.Popen(
        ["bedtools", "intersect", "-abam", bam, "-b", peaks, "-u"],
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
    )
    count = subprocess.run(
        ["samtools", "view", "-c"],
        stdin=intersect.stdout, capture_output=True, text=True, check=True,
    )
    intersect.stdout.close()
    rc = intersect.wait()
    if rc != 0:
        raise subprocess.CalledProcessError(rc, "bedtools intersect")
    return count.stdout.strip()


def _bampe_median_frag(bam):
    """Return median fragment size from bamPEFragmentSize."""
    result = subprocess.run(
        ["bamPEFragmentSize", "-b", bam],
        capture_output=True, text=True,
    )
    for line in result.stdout.splitlines():
        m = re.search(r"Median:\s*([\d.]+)", line)
        if m:
            return str(int(float(m.group(1))))
    return NA


# ---------------------------------------------------------------------------
# File-read helpers
# ---------------------------------------------------------------------------

def _count_lines(path):
    """Return line count of a file (for peak files)."""
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return "0"
    with open(path) as fh:
        return str(sum(1 for _ in fh))


def _narrowpeak_max_score(path):
    """Return max signalValue (col 6, 0-indexed) from a narrowPeak file."""
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return NA
    max_signal = None
    with open(path) as fh:
        for line in fh:
            fields = line.strip().split("\t")
            if len(fields) < 10:
                continue
            try:
                v = float(fields[6])
                if max_signal is None or v > max_signal:
                    max_signal = v
            except ValueError:
                continue
    return str(round(max_signal, 4)) if max_signal is not None else NA


def _fimo_motif_peaks(path):
    """Count unique peak names (sequence_name column) in a FIMO TSV."""
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return NA
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
    return str(len(seen)) if seen else NA


# ---------------------------------------------------------------------------
# Column definition
# ---------------------------------------------------------------------------

COLUMNS = [
    "sample",
    "total_reads",
    "subsampled_frags",
    "trimmed_reads",
    "mapped_reads",
    "mapping_pct",
    "median_frag_size",
    "num_peaks",
    "num_peaks_filt",
    "reads_in_peaks",
    "reads_in_peaks_5fold",
    "reads_in_peaks_filt",
    "max_peak_score",
    "motif_peaks",
]


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def main():
    sm = snakemake  # noqa: F821 — injected by Snakemake

    is_pe         = bool(sm.params.is_pe)
    is_treatment  = bool(sm.params.is_treatment)
    sample        = sm.wildcards.sample

    row = {c: NA for c in COLUMNS}
    row["sample"] = sample

    # --- Log-parsed stats (all samples) ---
    total_reads, subsampled = _parse_subsample_log(sm.input.subsample_log, is_pe)
    row["total_reads"] = total_reads
    row["subsampled_frags"] = _gate_subsampled(subsampled, sm.params.max_frags)
    row["trimmed_reads"] = _parse_bbduk_log(sm.input.trim_log)

    # --- Samtools count (all samples) ---
    row["mapped_reads"] = _samtools_count(sm.input.bam)
    row["mapping_pct"] = _safe_pct(row["mapped_reads"], row["trimmed_reads"])

    # --- PE-only stats ---
    if is_pe:
        row["median_frag_size"] = _bampe_median_frag(sm.input.bam)

    # --- Treatment-only stats ---
    if is_treatment:
        row["reads_in_peaks"] = _bedtools_intersect_count(
            sm.input.bam, sm.input.peaks
        )
        row["reads_in_peaks_5fold"] = _bedtools_intersect_count(
            sm.input.bam, sm.input.peaks_5fold
        )
        row["reads_in_peaks_filt"] = _bedtools_intersect_count(
            sm.input.bam, sm.input.peaks_filt
        )
        row["num_peaks"]      = _count_lines(sm.input.peaks)
        row["num_peaks_filt"] = _count_lines(sm.input.peaks_filt)
        row["max_peak_score"] = _narrowpeak_max_score(sm.input.peaks_filt)
        row["motif_peaks"]    = _fimo_motif_peaks(sm.input.fimo)

    with open(sm.output.stats_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerow(row)


if "snakemake" in dir():
    main()
