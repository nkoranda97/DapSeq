"""
Standalone script to regenerate qc_summary.tsv and qc_summary.html from an
existing results directory.

Usage:
    uv run python workflow/scripts/make_report.py --results-dir results/my_run
    uv run python workflow/scripts/make_report.py --results-dir results/my_run --output my_report.tsv

Sample discovery:
    - Treatment samples: MACS/{sample}_peaks_filt.narrowPeak (those that have a narrowPeak file)
    - Control samples are excluded (consistent with pipeline behaviour).
"""

import argparse
import glob
import os
import sys

# Reuse helpers from qc_summary.py
sys.path.insert(0, os.path.dirname(__file__))
from qc_summary import COLS, build_row, logo_to_base64, write_html


def discover_treatment_samples(results_dir):
    pattern = os.path.join(results_dir, "MACS", "*_peaks_filt.narrowPeak")
    paths = glob.glob(pattern)
    if not paths:
        sys.exit(f"No *_peaks_filt.narrowPeak files found under {results_dir}/MACS/")
    return sorted(
        os.path.basename(p).replace("_peaks_filt.narrowPeak", "") for p in paths
    )


def main():
    parser = argparse.ArgumentParser(description="Regenerate qc_summary from existing results.")
    parser.add_argument("--results-dir", required=True, help="Path to pipeline results directory")
    parser.add_argument("--output", help="Output TSV path (default: <results-dir>/stats/qc_summary.tsv)")
    args = parser.parse_args()

    results_dir = args.results_dir
    tsv_path  = args.output or os.path.join(results_dir, "stats", "qc_summary.tsv")
    html_path = tsv_path.replace(".tsv", ".html")

    treatment_samples = discover_treatment_samples(results_dir)

    stats_dir    = os.path.join(results_dir, "stats")
    trim_log_dir = os.path.join(results_dir, "logs", "bbduk")
    macs_dir     = os.path.join(results_dir, "MACS")
    meme_dir     = os.path.join(results_dir, "meme")

    rows = [
        build_row(s, trim_log_dir, stats_dir, macs_dir)
        for s in treatment_samples
    ]

    with open(tsv_path, "w") as fh:
        fh.write("\t".join(COLS) + "\n")
        for row in rows:
            fh.write("\t".join(str(row.get(c, "NA")) for c in COLS) + "\n")

    logo_b64_map = {
        s: logo_to_base64(os.path.join(meme_dir, s, "summits", "logo1.png"))
        for s in treatment_samples
    }
    write_html(rows, logo_b64_map, html_path)

    print(f"Written: {tsv_path}")
    print(f"Written: {html_path}")
    for row in rows:
        print("\t".join(str(row.get(c, "NA")) for c in COLS))


if __name__ == "__main__":
    main()
