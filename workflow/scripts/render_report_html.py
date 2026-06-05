"""
Render report.html from an existing report.csv and logo PNG directories.

Called via Snakemake's script: directive (report_html rule) or directly from
the CLI to regenerate the HTML without recomputing per-sample statistics.
"""

import argparse
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
from report import logo_to_base64, make_cols, write_html  # noqa: E402


def render(treatment_samples, report_csv, meme_dir, factorbook_dir, html_out):
    cols = make_cols()

    with open(report_csv, newline="") as fh:
        rows = list(csv.DictReader(fh))

    logo_b64_map = {
        s: logo_to_base64(os.path.join(meme_dir, s, "summits", "logo1.png"))
        for s in treatment_samples
    }
    logo_rc_b64_map = {
        s: logo_to_base64(os.path.join(meme_dir, s, "summits", "logo_rc1.png"))
        for s in treatment_samples
    }
    factorbook_logo_map = {
        s: logo_to_base64(os.path.join(factorbook_dir, f"{s}.logo.png"))
        for s in treatment_samples
    }
    write_html(rows, cols, logo_b64_map, logo_rc_b64_map, factorbook_logo_map, html_out)


def main():
    """Snakemake entry point — writes report.html from report.csv + logos."""
    sm = snakemake  # noqa: F821 — injected by Snakemake
    render(
        treatment_samples=list(sm.params.treatment_samples),
        report_csv=str(sm.input.report_csv),
        meme_dir=sm.params.meme_dir,
        factorbook_dir=sm.params.factorbook_dir,
        html_out=str(sm.output.html),
    )


def cli_main():
    """CLI entry point — regenerates report.html from an existing out directory."""
    parser = argparse.ArgumentParser(
        description="Re-render report.html from an existing report.csv and logo PNGs."
    )
    parser.add_argument("out_dir", help="Pipeline output directory (contains stats/report.csv, meme/, factorbook/)")
    args = parser.parse_args()

    out_dir        = os.path.abspath(args.out_dir)
    stats_dir      = os.path.join(out_dir, "stats")
    report_csv     = os.path.join(stats_dir, "report.csv")
    meme_dir       = os.path.join(out_dir, "meme")
    factorbook_dir = os.path.join(out_dir, "factorbook")
    html_out       = os.path.join(stats_dir, "report.html")

    if not os.path.exists(report_csv):
        sys.exit(f"report.csv not found at {report_csv} — run the pipeline first.")

    with open(report_csv, newline="") as fh:
        reader = csv.DictReader(fh)
        samples = [row["sample"] for row in reader if row.get("sample")]

    if not samples:
        sys.exit(f"No samples found in {report_csv}.")

    render(
        treatment_samples=samples,
        report_csv=report_csv,
        meme_dir=meme_dir,
        factorbook_dir=factorbook_dir,
        html_out=html_out,
    )
    print(f"Wrote {html_out}")


if "snakemake" in dir():
    main()
elif __name__ == "__main__":
    cli_main()
