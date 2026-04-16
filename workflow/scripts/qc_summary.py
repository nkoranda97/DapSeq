"""
Aggregate per-sample QC statistics into a TSV summary table and a companion
HTML report with embedded motif logo images.

Only treatment samples are included — control samples are excluded.
Called via Snakemake's script: directive; uses the snakemake object for I/O.
"""

import base64
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


def logo_to_base64(png_path):
    """Return base64-encoded PNG data URI, or None if file is empty/missing."""
    if not os.path.exists(png_path) or os.path.getsize(png_path) == 0:
        return None
    with open(png_path, "rb") as fh:
        data = fh.read()
    return base64.b64encode(data).decode("ascii")


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
]


def build_row(sample, trim_log_dir, stats_dir, macs_dir):
    total_frags_path = os.path.join(stats_dir, f"{sample}.total_frags.txt")
    total_frags = open(total_frags_path).read().strip() if os.path.exists(total_frags_path) else "NA"

    trim_log_path = os.path.join(trim_log_dir, f"{sample}.trim.log")
    clean_reads = parse_bbduk_log(trim_log_path) if os.path.exists(trim_log_path) else "NA"

    align_rate_path = os.path.join(stats_dir, f"{sample}.align_rate.txt")
    align_rate = open(align_rate_path).read().strip() if os.path.exists(align_rate_path) else "NA"

    filt_stats = parse_kv_file(os.path.join(stats_dir, f"{sample}.filtered_stats.txt"))
    filtered_reads = filt_stats.get("filtered_reads", "NA")

    np_path = os.path.join(macs_dir, f"{sample}_peaks_filt.narrowPeak")
    peak_total, min5fold, max_score = parse_narrowpeak(np_path)

    frip_macs = parse_kv_file(os.path.join(stats_dir, f"{sample}.frip_macs.txt"))
    peak_reads = frip_macs.get("reads_in_peaks_macs", "NA")
    if peak_reads != "NA" and filtered_reads != "NA":
        frip = round(int(peak_reads) / int(filtered_reads) * 100, 2)
    else:
        frip = "NA"

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
    }


_HTML_STYLE = """
<style>
  body { font-family: sans-serif; font-size: 13px; margin: 20px; }
  table { border-collapse: collapse; width: 100%; }
  th { background: #2c3e50; color: white; padding: 8px 10px; text-align: left; }
  td { padding: 6px 10px; border-bottom: 1px solid #ddd; vertical-align: middle; }
  tr:nth-child(even) { background: #f7f7f7; }
  tr:hover { background: #eaf4fb; }
  img.motif { max-width: 280px; max-height: 100px; height: auto; display: block; }
  .na { color: #aaa; font-style: italic; }
</style>
"""

_INT_COLS = {"total_frags", "clean_reads", "filtered_reads",
             "peak#", "min5fold peak#", "peak reads#"}
_PCT_COLS = {"align_rate%", "FRiP_score"}


def _fmt_html(col, val):
    if val == "NA":
        return "NA"
    if col in _INT_COLS:
        return f"{int(val):,}"
    if col in _PCT_COLS:
        return f"{val}%"
    return str(val)

_HTML_COLS = COLS + ["top motif"]


def write_html(rows, logo_b64_map, out_path):
    """Write a self-contained HTML summary table with embedded motif logos."""
    lines = [
        "<!DOCTYPE html>",
        "<html><head><meta charset='utf-8'>",
        "<title>DAP-seq QC Summary</title>",
        _HTML_STYLE,
        "</head><body>",
        "<h2>DAP-seq QC Summary</h2>",
        "<table>",
        "<thead><tr>",
    ]
    for col in _HTML_COLS:
        lines.append(f"  <th>{col}</th>")
    lines.append("</tr></thead><tbody>")

    for row in rows:
        lines.append("<tr>")
        for col in COLS:
            val = row.get(col, "NA")
            css = ' class="na"' if val == "NA" else ""
            lines.append(f"  <td{css}>{_fmt_html(col, val)}</td>")
        # Motif logo column
        b64 = logo_b64_map.get(row["sample"])
        if b64:
            lines.append(f'  <td><img class="motif" src="data:image/png;base64,{b64}" alt="motif logo"/></td>')
        else:
            lines.append('  <td class="na">NA</td>')
        lines.append("</tr>")

    lines += ["</tbody></table>", "</body></html>"]

    with open(out_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def main():
    sm = snakemake  # noqa: F821 — injected by Snakemake

    treatment_samples = list(sm.params.treatment_samples)
    trim_log_dir      = sm.params.trim_log_dir
    macs_dir          = sm.params.macs_dir
    meme_dir          = sm.params.meme_dir
    stats_dir         = sm.params.stats_dir

    rows = [
        build_row(s, trim_log_dir, stats_dir, macs_dir)
        for s in treatment_samples
    ]

    # TSV output
    with open(sm.output.tsv, "w") as out:
        out.write("\t".join(COLS) + "\n")
        for row in rows:
            out.write("\t".join(str(row.get(c, "NA")) for c in COLS) + "\n")

    # HTML output — embed motif logos as base64
    logo_b64_map = {
        s: logo_to_base64(os.path.join(meme_dir, s, "summits", "logo1.png"))
        for s in treatment_samples
    }
    write_html(rows, logo_b64_map, sm.output.html)


if "snakemake" in dir():
    main()
