"""
Aggregate per-sample QC statistics into a CSV summary table.

One row per treatment sample. Reads from {sample}.stats.csv files produced
by collect_stats.py. Called via Snakemake's script: directive (CSV only) or
from the CLI (CSV + HTML together).
"""

import argparse
import base64
import csv
import io
import os
import sys


INT_COLS = {
    "total_reads", "trimmed_reads", "mapped_reads",
    "reads_in_peaks", "reads_5fold", "reads_nfold",
    "motif_peaks", "subsampled_frags", "num_peaks", "num_peaks_filt",
}
PCT_COLS = {"frip", "frip_top_n_fold", "mapping_pct"}


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
        # ── alignment / FRiP ───────────────────────────────────────────────
        "mapping_rate",
        "mapping_pct",
        "frip",
        "frip_top_n_fold",
        # ── peak quality ───────────────────────────────────────────────────
        "max_peak_score",
        "motif_peaks",
        # ── additional QC ──────────────────────────────────────────────────
        "subsampled_frags",
        "median_frag_size",
        "num_peaks",
        "num_peaks_filt",
    ]


def _safe_frip(reads_in_peaks, mapped_reads):
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
    row["frip"] = _safe_frip(
        data.get("reads_in_peaks", "NA"),
        data.get("mapped_reads", "NA"),
    )
    row["frip_top_n_fold"] = _safe_frip(
        data.get("reads_nfold", "NA"),
        data.get("mapped_reads", "NA"),
    )
    return row


# ---------------------------------------------------------------------------
# HTML report helpers
# ---------------------------------------------------------------------------

_HTML_STYLE = """
<style>
  body { font-family: sans-serif; font-size: 13px; margin: 20px; }
  table { border-collapse: collapse; width: 100%; }
  th { background: #2c3e50; color: white; padding: 8px 10px; text-align: left; cursor: pointer; user-select: none; }
  th.no-sort { cursor: default; }
  th.sort-asc::after  { content: " ▲"; font-size: 10px; }
  th.sort-desc::after { content: " ▼"; font-size: 10px; }
  td { padding: 6px 10px; border-bottom: 1px solid #ddd; vertical-align: middle; }
  tr:nth-child(even) { background: #f7f7f7; }
  tr:hover { background: #eaf4fb; }
  img.motif { max-width: 280px; max-height: 100px; height: auto; display: block; }
  .na { color: #aaa; font-style: italic; }
</style>
"""

_HTML_SCRIPT = """
<script>
(function () {
  var table = document.querySelector('table');
  var headers = table.querySelectorAll('thead th');
  var tbody = table.querySelector('tbody');
  var sortCol = -1;
  var sortAsc = true;

  headers.forEach(function (th, idx) {
    if (th.classList.contains('no-sort')) return;
    th.addEventListener('click', function () {
      var isNumeric = th.dataset.colType === 'numeric';
      if (sortCol === idx) {
        sortAsc = !sortAsc;
      } else {
        sortCol = idx;
        sortAsc = true;
      }
      headers.forEach(function (h) { h.classList.remove('sort-asc', 'sort-desc'); });
      th.classList.add(sortAsc ? 'sort-asc' : 'sort-desc');

      var rows = Array.from(tbody.querySelectorAll('tr'));
      rows.sort(function (a, b) {
        var av = a.cells[idx] ? a.cells[idx].textContent.trim() : '';
        var bv = b.cells[idx] ? b.cells[idx].textContent.trim() : '';
        var aNa = av === 'NA';
        var bNa = bv === 'NA';
        if (aNa && bNa) return 0;
        if (aNa) return 1;
        if (bNa) return -1;
        if (isNumeric) {
          var an = parseFloat(av.replace(/,/g, '').replace('%', ''));
          var bn = parseFloat(bv.replace(/,/g, '').replace('%', ''));
          return sortAsc ? an - bn : bn - an;
        }
        return sortAsc ? av.localeCompare(bv) : bv.localeCompare(av);
      });
      rows.forEach(function (r) { tbody.appendChild(r); });
    });
  });
})();
</script>
"""


def logo_to_base64(png_path):
    """Return base64-encoded PNG, or None if file is empty/missing."""
    if not os.path.exists(png_path) or os.path.getsize(png_path) == 0:
        return None
    from PIL import Image, ImageChops
    img = Image.open(png_path).convert("RGB")
    bg = Image.new("RGB", img.size, (255, 255, 255))
    bbox = ImageChops.difference(img, bg).getbbox()
    if bbox:
        pad = 6
        w, h = img.size
        bbox = (max(0, bbox[0] - pad), max(0, bbox[1] - pad),
                min(w, bbox[2] + pad), min(h, bbox[3] + pad))
        img = img.crop(bbox)
    buf = io.BytesIO()
    img.save(buf, format="PNG")
    return base64.b64encode(buf.getvalue()).decode("ascii")


def _fmt_html(col, val):
    if val == "NA":
        return "NA"
    if col in INT_COLS:
        return f"{int(val):,}"
    if col in PCT_COLS:
        return f"{val}%"
    return str(val)


def write_html(rows, cols, logo_b64_map, logo_rc_b64_map, factorbook_logo_map, out_path):
    """Write a self-contained HTML report table with embedded motif logos."""
    html_cols = cols + ["top motif", "top motif RC", "factorbook motif"]
    lines = [
        "<!DOCTYPE html>",
        "<html><head><meta charset='utf-8'>",
        "<title>DAP-seq Report</title>",
        _HTML_STYLE,
        "</head><body>",
        "<h2>DAP-seq Report</h2>",
        "<table>",
        "<thead><tr>",
    ]
    for col in html_cols:
        if col in ("top motif", "top motif RC", "factorbook motif"):
            lines.append(f'  <th class="no-sort">{col}</th>')
        elif col == "sample":
            lines.append(f'  <th data-col-type="text">{col}</th>')
        else:
            lines.append(f'  <th data-col-type="numeric">{col}</th>')
    lines.append("</tr></thead><tbody>")

    for row in rows:
        sample = row.get("sample", "")
        lines.append("<tr>")
        for col in cols:
            val = row.get(col, "NA")
            css = ' class="na"' if val == "NA" else ""
            lines.append(f"  <td{css}>{_fmt_html(col, val)}</td>")
        for logo_map, alt in (
            (logo_b64_map,        "motif logo"),
            (logo_rc_b64_map,     "motif RC logo"),
            (factorbook_logo_map, "factorbook motif"),
        ):
            b64 = logo_map.get(sample)
            if b64:
                lines.append(f'  <td><img class="motif" src="data:image/png;base64,{b64}" alt="{alt}"/></td>')
            else:
                lines.append('  <td class="na">NA</td>')
        lines.append("</tr>")

    lines += ["</tbody></table>", _HTML_SCRIPT, "</body></html>"]
    with open(out_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


# ---------------------------------------------------------------------------
# Run functions
# ---------------------------------------------------------------------------

def run_csv(treatment_samples, stats_dir, csv_out):
    """Build per-sample rows and write report.csv."""
    cols = make_cols()
    rows = [build_row(s, stats_dir) for s in treatment_samples]
    with open(csv_out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({c: row.get(c, "NA") for c in cols})
    return rows


def run_html(treatment_samples, stats_dir, meme_dir, factorbook_dir, csv_out, html_out):
    """Write report.csv then render report.html from the same rows + logos."""
    rows = run_csv(treatment_samples, stats_dir, csv_out)
    cols = make_cols()
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


# ---------------------------------------------------------------------------
# Entry points
# ---------------------------------------------------------------------------

def main():
    """Snakemake entry point — writes report.csv only."""
    sm = snakemake  # noqa: F821 — injected by Snakemake
    run_csv(
        treatment_samples=list(sm.params.treatment_samples),
        stats_dir=sm.params.stats_dir,
        csv_out=sm.output.csv,
    )


def cli_main():
    """CLI entry point — writes both report.csv and report.html."""
    parser = argparse.ArgumentParser(
        description="Regenerate the DAP-seq stats report from completed pipeline output."
    )
    parser.add_argument("out_dir", help="Pipeline output directory (contains stats/, meme/, factorbook/)")
    args = parser.parse_args()

    out_dir        = os.path.abspath(args.out_dir)
    stats_dir      = os.path.join(out_dir, "stats")
    meme_dir       = os.path.join(out_dir, "meme")
    factorbook_dir = os.path.join(out_dir, "factorbook")
    csv_out        = os.path.join(stats_dir, "report.csv")
    html_out       = os.path.join(stats_dir, "report.html")

    samples = sorted(
        f[: -len(".stats.csv")]
        for f in os.listdir(stats_dir)
        if f.endswith(".stats.csv")
    )
    if not samples:
        sys.exit(f"No *.stats.csv files found in {stats_dir} — is the pipeline complete?")

    run_html(
        treatment_samples=samples,
        stats_dir=stats_dir,
        meme_dir=meme_dir,
        factorbook_dir=factorbook_dir,
        csv_out=csv_out,
        html_out=html_out,
    )
    print(f"Wrote {csv_out}")
    print(f"Wrote {html_out}")


if "snakemake" in dir():
    main()
elif __name__ == "__main__":
    cli_main()
