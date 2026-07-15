"""
Aggregate per-sample QC statistics into a CSV summary table.

One row per treatment sample. Reads from {sample}.stats.csv files produced
by collect_stats.py. Called via Snakemake's script: directive (CSV only) or
from the CLI (CSV + HTML together).
"""

import argparse
import base64
import csv
import html
import io
import os
import sys


INT_COLS = {
    "total_reads", "trimmed_reads", "mapped_reads",
    "reads_in_peaks", "reads_in_peaks_filt",
    "motif_peaks", "subsampled_frags",
    "num_peaks", "num_peaks_filt", "num_peaks_bl", "num_peaks_rmsk",
}
PCT_COLS = {"frip", "frip_filt", "mapping_pct"}
# Columns rendered (and sorted) as free-form text rather than numeric.
TEXT_COLS = {"sample", "experiment_date", "gdna_batch"}


def make_cols():
    return [
        "sample",
        "experiment_date",
        "gdna_batch",
        "total_reads",
        "subsampled_frags",
        "trimmed_reads",
        "mapped_reads",
        "mapping_pct",
        "median_frag_size",
        "num_peaks",
        "num_peaks_filt",
        "num_peaks_bl",
        "num_peaks_rmsk",
        "reads_in_peaks",
        "reads_in_peaks_filt",
        "frip",
        "frip_filt",
        "max_peak_score",
        "motif_peaks",
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
    row["frip_filt"] = _safe_frip(
        data.get("reads_in_peaks_filt", "NA"),
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
  p.legend { color: #555; font-size: 12px; }
  /* Control-row highlight. Placed last so these win specificity ties
     against tr:nth-child(even) and tr:hover. */
  tr.control-row { background: #fff4e5; }
  tr.control-row:hover { background: #ffe9cc; }
  .control-badge { display: inline-block; background: #e67e22; color: white;
    font-size: 10px; padding: 1px 6px; border-radius: 3px; margin-left: 6px;
    vertical-align: middle; }
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
    return html.escape(str(val))


def _report_header_html(filter_foldch=None):
    """Return the header <p> block describing the _filt fold threshold.

    Experiment metadata (experiment_date, gdna_batch) is per-sample and lives
    in the table columns, not this header. The filter note is emitted only when
    a fold value is supplied.
    """
    if filter_foldch is None:
        return ""
    return (
        f"<p><strong>Filtered peaks</strong> "
        f"(<code>num_peaks_filt</code>, <code>reads_in_peaks_filt</code>, "
        f"<code>frip_filt</code>) use fold-change&nbsp;&ge;&nbsp;"
        f"{filter_foldch}&times; — the peak set fed to MEME/FIMO.</p>"
    )


def write_html(rows, cols, logo_b64_map, logo_rc_b64_map, factorbook_logo_map, out_path,
               filter_foldch=None, control_samples=None):
    """Write a self-contained HTML report table with embedded motif logos.

    Rows whose sample is in *control_samples* are flagged (a highlighted row and
    a "control" badge) — those controls were peak-called against themselves.
    """
    control_set = set(control_samples or [])
    html_cols = cols + ["top motif", "top motif RC", "factorbook motif"]

    header = _report_header_html(filter_foldch=filter_foldch)
    if control_set:
        header += (
            '\n<p class="legend">Rows marked '
            '<span class="control-badge">control</span> were peak-called against '
            'themselves (expect few or no peaks).</p>'
        )

    lines = [
        "<!DOCTYPE html>",
        "<html><head><meta charset='utf-8'>",
        "<title>DAP-seq Report</title>",
        _HTML_STYLE,
        "</head><body>",
        "<h2>DAP-seq Report</h2>",
        header,
        "<table>",
        "<thead><tr>",
    ]
    for col in html_cols:
        if col in ("top motif", "top motif RC", "factorbook motif"):
            lines.append(f'  <th class="no-sort">{col}</th>')
        elif col in TEXT_COLS:
            lines.append(f'  <th data-col-type="text">{col}</th>')
        else:
            lines.append(f'  <th data-col-type="numeric">{col}</th>')
    lines.append("</tr></thead><tbody>")

    for row in rows:
        sample = row.get("sample", "")
        is_ctrl = sample in control_set
        lines.append('<tr class="control-row">' if is_ctrl else "<tr>")
        for col in cols:
            val = row.get(col, "NA")
            css = ' class="na"' if val == "NA" else ""
            if col == "sample" and is_ctrl:
                lines.append(
                    f'  <td{css}>{_fmt_html(col, val)} '
                    f'<span class="control-badge">control</span></td>'
                )
            else:
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

def run_csv(samples, stats_dir, csv_out):
    """Build per-sample rows and write report.csv."""
    cols = make_cols()
    rows = [build_row(s, stats_dir) for s in samples]
    with open(csv_out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({c: row.get(c, "NA") for c in cols})
    return rows


def _load_config_snapshot(out_dir):
    """Best-effort load of ``{out_dir}/config_used.yaml`` (written by the
    pipeline's onstart hook). Returns the parsed dict, or {} when the snapshot
    is missing or unparseable."""
    path = os.path.join(out_dir, "config_used.yaml")
    if not os.path.exists(path):
        return {}
    try:
        import yaml
        with open(path) as fh:
            return yaml.safe_load(fh) or {}
    except Exception:
        return {}


def read_run_metadata(out_dir):
    """Best-effort read of the _filt fold value from the config snapshot.

    Returns the fold-change threshold that the reported ``_filt`` columns
    represent, or None when the snapshot is missing, unparseable, or lacks the
    value. Experiment metadata is per-sample and comes from the report CSV, not
    this snapshot.
    """
    cfg    = _load_config_snapshot(out_dir)
    macs3  = cfg.get("macs3") or {}
    levels = macs3.get("foldch_levels")
    idx    = macs3.get("meme_foldch_level", 2)
    if isinstance(levels, (list, tuple)) and isinstance(idx, int) and 1 <= idx <= len(levels):
        return levels[idx - 1]
    return None


def read_control(out_dir):
    """Return the list of control sample names from the config snapshot.

    Normalizes the ``control`` field (scalar / list / null) to a list; returns
    [] when the snapshot is missing, unparseable, or has no control.
    """
    ctrl = _load_config_snapshot(out_dir).get("control")
    if not ctrl:
        return []
    return list(ctrl) if isinstance(ctrl, (list, tuple)) else [ctrl]


def run_html(samples, stats_dir, meme_dir, factorbook_dir, csv_out, html_out,
             filter_foldch=None, control_samples=None):
    """Write report.csv then render report.html from the same rows + logos."""
    rows = run_csv(samples, stats_dir, csv_out)
    cols = make_cols()
    logo_b64_map = {
        s: logo_to_base64(os.path.join(meme_dir, s, "summits", "logo1.png"))
        for s in samples
    }
    logo_rc_b64_map = {
        s: logo_to_base64(os.path.join(meme_dir, s, "summits", "logo_rc1.png"))
        for s in samples
    }
    factorbook_logo_map = {
        s: logo_to_base64(os.path.join(factorbook_dir, f"{s}.logo.png"))
        for s in samples
    }
    write_html(rows, cols, logo_b64_map, logo_rc_b64_map, factorbook_logo_map, html_out,
               filter_foldch=filter_foldch, control_samples=control_samples)


# ---------------------------------------------------------------------------
# Entry points
# ---------------------------------------------------------------------------

def main():
    """Snakemake entry point — writes report.csv only."""
    sm = snakemake  # noqa: F821 — injected by Snakemake
    run_csv(
        samples=list(sm.params.report_samples),
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

    filter_foldch = read_run_metadata(out_dir)

    run_html(
        samples=samples,
        stats_dir=stats_dir,
        meme_dir=meme_dir,
        factorbook_dir=factorbook_dir,
        csv_out=csv_out,
        html_out=html_out,
        filter_foldch=filter_foldch,
        control_samples=read_control(out_dir),
    )
    print(f"Wrote {csv_out}")
    print(f"Wrote {html_out}")


if "snakemake" in dir():
    main()
elif __name__ == "__main__":
    cli_main()
