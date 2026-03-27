"""
Standalone script to regenerate qc_summary.tsv from an existing results directory.

Usage:
    uv run python workflow/scripts/make_report.py --results-dir results/niu_a_thaliana_bowtie2
    uv run python workflow/scripts/make_report.py --results-dir results/niu_a_thaliana_bowtie2 --output my_report.tsv

Sample discovery:
    - All samples:       stats/{sample}.filtered_stats.txt
    - Treatment samples: MACS/{sample}_peaks_filt.narrowPeak (those that have a narrowPeak file)
"""

import argparse
import glob
import os
import sys

# Reuse helpers from qc_summary.py
sys.path.insert(0, os.path.dirname(__file__))
from qc_summary import COLS, build_row


def discover_samples(results_dir):
    pattern = os.path.join(results_dir, "stats", "*.filtered_stats.txt")
    paths = glob.glob(pattern)
    if not paths:
        sys.exit(f"No filtered_stats.txt files found under {results_dir}/stats/")
    samples = sorted(
        os.path.basename(p).replace(".filtered_stats.txt", "") for p in paths
    )
    treatment = [
        s for s in samples
        if os.path.exists(os.path.join(results_dir, "MACS", f"{s}_peaks_filt.narrowPeak"))
    ]
    return samples, treatment


def main():
    parser = argparse.ArgumentParser(description="Regenerate qc_summary.tsv from existing results.")
    parser.add_argument("--results-dir", required=True, help="Path to pipeline results directory")
    parser.add_argument("--output", help="Output TSV path (default: <results-dir>/stats/qc_summary.tsv)")
    args = parser.parse_args()

    results_dir = args.results_dir
    out_path = args.output or os.path.join(results_dir, "stats", "qc_summary.tsv")

    samples, treatment_samples = discover_samples(results_dir)
    treatment_set = set(treatment_samples)

    stats_dir    = os.path.join(results_dir, "stats")
    trim_log_dir = os.path.join(results_dir, "logs", "bbduk")
    macs_dir     = os.path.join(results_dir, "MACS")
    meme_dir     = os.path.join(results_dir, "meme")

    rows = [
        build_row(s, s in treatment_set, trim_log_dir, stats_dir, macs_dir, meme_dir)
        for s in samples
    ]

    with open(out_path, "w") as fh:
        fh.write("\t".join(COLS) + "\n")
        for row in rows:
            fh.write("\t".join(str(row.get(c, "NA")) for c in COLS) + "\n")

    print(f"Written: {out_path}")
    for row in rows:
        print("\t".join(str(row.get(c, "NA")) for c in COLS))


if __name__ == "__main__":
    main()
