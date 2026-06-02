"""
Aggregate per-sample QC statistics into a CSV summary table.

One row per treatment sample. Reads from {sample}.stats.csv files produced
by collect_stats.py. Called via Snakemake's script: directive or from CLI.
"""

import argparse
import csv
import os
import sys


def make_cols():
    return [
        "sample",
        # ── reads ──────────────────────────────────────────────────────────
        "total_reads",
        "trimmed_reads",
        "mapped_reads",
        "reads_in_peaks",
        "reads_5fold",
        "reads_nfold",
        # ── alignment / mapping ────────────────────────────────────────────
        "alignment_rate",
        "mapping_pct",
        # ── peak quality ───────────────────────────────────────────────────
        "max_peak_score",
        "motif_peaks",
        # ── additional QC ──────────────────────────────────────────────────
        "subsampled_frags",
        "median_frag_size",
        "num_peaks",
        "num_peaks_filt",
    ]


def _safe_mapping_pct(reads_in_peaks, mapped_reads):
    if reads_in_peaks == "NA" or mapped_reads == "NA":
        return "NA"
    try:
        denom = int(mapped_reads)
        return round(int(reads_in_peaks) / denom * 100, 2) if denom > 0 else "NA"
    except (ValueError, ZeroDivisionError):
        return "NA"


def build_row(sample, stats_dir):
    stats_path = os.path.join(stats_dir, f"{sample}.stats.csv")
    if not os.path.exists(stats_path):
        return {"sample": sample}

    with open(stats_path, newline="") as fh:
        reader = csv.DictReader(fh)
        rows = list(reader)

    if not rows:
        return {"sample": sample}

    data = rows[0]
    row = dict(data)
    row["mapping_pct"] = _safe_mapping_pct(
        data.get("reads_in_peaks", "NA"),
        data.get("mapped_reads", "NA"),
    )
    return row


def run(treatment_samples, stats_dir, csv_out):
    cols = make_cols()
    rows = [build_row(s, stats_dir) for s in treatment_samples]

    with open(csv_out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({c: row.get(c, "NA") for c in cols})


def main():
    sm = snakemake  # noqa: F821 — injected by Snakemake
    run(
        treatment_samples = list(sm.params.treatment_samples),
        stats_dir         = sm.params.stats_dir,
        csv_out           = sm.output.csv,
    )


def cli_main():
    parser = argparse.ArgumentParser(
        description="Regenerate the DAP-seq stats CSV from completed pipeline output."
    )
    parser.add_argument(
        "out_dir",
        help="Pipeline output directory (contains stats/ with *.stats.csv files)",
    )
    args = parser.parse_args()

    out_dir   = os.path.abspath(args.out_dir)
    stats_dir = os.path.join(out_dir, "stats")
    csv_out   = os.path.join(stats_dir, "report.csv")

    samples = sorted(
        f[: -len(".stats.csv")]
        for f in os.listdir(stats_dir)
        if f.endswith(".stats.csv")
    )
    if not samples:
        sys.exit(f"No *.stats.csv files found in {stats_dir} — is the pipeline complete?")

    run(
        treatment_samples=samples,
        stats_dir=stats_dir,
        csv_out=csv_out,
    )
    print(f"Wrote {csv_out}")


if "snakemake" in dir():
    main()
elif __name__ == "__main__":
    cli_main()
