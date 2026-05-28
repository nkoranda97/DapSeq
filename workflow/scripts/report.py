"""
Aggregate per-sample QC statistics into a TSV summary table and a companion
HTML report with embedded motif logo images.

Only treatment samples are included — control samples are excluded.
Called via Snakemake's script: directive; uses the snakemake object for I/O.
"""

import argparse
import base64
import os
import re
import sys

import subprocess

_SCRIPT_DIR      = os.path.dirname(os.path.abspath(__file__))
_DEFAULT_FB_TSV  = os.path.normpath(os.path.join(_SCRIPT_DIR, '..', '..', 'factorbook', 'factorbook_chip_seq_meme_motifs.tsv'))
_DEFAULT_FB_MEME = os.path.normpath(os.path.join(_SCRIPT_DIR, '..', '..', 'factorbook', 'factorbook_chip_seq_meme_motif_catalog.meme'))


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
    Return (total_peaks, max_signal) from a narrowPeak file.
    Column 7 (0-indexed: 6) is signalValue (fold enrichment).
    """
    total = 0
    max_signal = None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            fields = line.split("\t")
            signal = float(fields[6])
            total += 1
            if max_signal is None or signal > max_signal:
                max_signal = signal
    if max_signal is None:
        return "NA", "NA"
    return total, round(max_signal, 4)


def logo_to_base64(png_path):
    """Return base64-encoded PNG data URI, or None if file is empty/missing."""
    if not os.path.exists(png_path) or os.path.getsize(png_path) == 0:
        return None
    from PIL import Image, ImageChops
    import io
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



def _render_logo_with_logomaker(ppm, out_path):
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import logomaker

    df    = pd.DataFrame(ppm, columns=['A', 'C', 'G', 'T'])
    ic_df = logomaker.transform_matrix(df, from_type='probability', to_type='information')

    fig, ax = plt.subplots(figsize=(max(2, len(ppm) * 0.35), 2.5))
    logomaker.Logo(ic_df, ax=ax, color_scheme='classic')
    ax.set_xticks(range(len(ppm)))
    ax.set_xticklabels(range(1, len(ppm) + 1), fontsize=7)
    ax.set_ylabel('bits', fontsize=8)
    ax.spines[['top', 'right']].set_visible(False)
    plt.tight_layout()
    fig.savefig(out_path, format='png', dpi=150, bbox_inches='tight')
    plt.close(fig)


def _make_factorbook_logos(samples, tsv_path, meme_path, factorbook_dir):
    """Generate factorbook logo PNGs for any samples that don't already have one."""
    os.makedirs(factorbook_dir, exist_ok=True)
    for sample in samples:
        out_png = os.path.join(factorbook_dir, f"{sample}.logo.png")
        if os.path.exists(out_png) and os.path.getsize(out_png) > 0:
            continue

        tf = sample.upper()
        motif_id = None
        with open(tsv_path) as fh:
            header = fh.readline().rstrip('\n').split('\t')
            ti = header.index('target')
            ai = header.index('dataset_accession')
            ni = header.index('name')
            for line in fh:
                parts = line.rstrip('\n').split('\t')
                if parts[ti].upper() == tf:
                    motif_id = f"{parts[ai]}_{parts[ni]}"
                    break

        if motif_id is None:
            open(out_png, 'w').close()
            continue

        header_lines, motif_lines = [], []
        ppm = []
        in_header  = True
        in_motif   = False
        in_matrix  = False
        with open(meme_path) as fh:
            for line in fh:
                if in_motif:
                    if line.startswith('MOTIF '):
                        break
                    motif_lines.append(line)
                    if 'letter-probability matrix' in line:
                        in_matrix = True
                    elif in_matrix:
                        try:
                            vals = [float(x) for x in line.split()]
                            if len(vals) == 4:
                                ppm.append(vals)
                            else:
                                in_matrix = False
                        except ValueError:
                            in_matrix = False
                elif line.startswith('MOTIF '):
                    in_header = False
                    if line.split()[1] == motif_id:
                        in_motif = True
                        motif_lines.append(line)
                elif in_header:
                    header_lines.append(line)

        tmp_meme = out_png + ".meme"
        tmp_eps  = out_png + ".eps"
        try:
            with open(tmp_meme, 'w') as fh:
                fh.writelines(header_lines)
                fh.writelines(motif_lines)
            r = subprocess.run(['ceqlogo', '-i1', tmp_meme, '-o', tmp_eps, '-f', 'EPS'],
                               capture_output=True)
            if r.returncode == 0 and os.path.exists(tmp_eps):
                subprocess.run(['gs', '-dNOPAUSE', '-dBATCH', '-sDEVICE=png16m', '-r150',
                                f'-sOutputFile={out_png}', tmp_eps],
                               capture_output=True)
        except FileNotFoundError:
            pass  # ceqlogo not on PATH
        finally:
            for f in (tmp_meme, tmp_eps):
                if os.path.exists(f):
                    os.remove(f)
        if not os.path.exists(out_png) or os.path.getsize(out_png) == 0:
            if ppm:
                _render_logo_with_logomaker(ppm, out_png)
            else:
                open(out_png, 'w').close()


def _detect_samples_and_bam_types(stats_dir):
    """Scan stats_dir for *.frip_macs.txt to infer treatment samples and bam_types."""
    samples, bam_types = set(), set()
    for fname in os.listdir(stats_dir):
        if fname.endswith('.frip_macs.txt'):
            stem = fname[:-len('.frip_macs.txt')]
            parts = stem.rsplit('.', 1)
            if len(parts) == 2:
                samples.add(parts[0])
                bam_types.add(parts[1])
    if not samples:
        sys.exit(f"No *.frip_macs.txt files found in {stats_dir} — is the pipeline complete?")
    return sorted(samples), sorted(bam_types)


def _fmt_n(min_foldch):
    """Format fold-change threshold as int when whole, else as-is (e.g. 2.0 -> '2', 2.5 -> '2.5')."""
    return str(int(min_foldch)) if min_foldch == int(min_foldch) else str(min_foldch)


def make_cols(n):
    return [
        "sample",
        "total_frags",
        "clean_reads",
        "align_rate%",
        "filtered_reads",
        "peak#",
        f"peaks_{n}fold#",
        "max peak score",
        "peak reads#",
        "peak reads filt#",
        "FRiP_score",
        "FRiP_top_n_fold",
    ]


def int_cols(n):
    return {"total_frags", "clean_reads", "filtered_reads", "peak#", f"peaks_{n}fold#", "peak reads#", "peak reads filt#"}


PCT_COLS = {"align_rate%", "FRiP_score", "FRiP_top_n_fold"}


def build_row(sample, bam_type, trim_log_dir, stats_dir, macs_dir, min_foldch):
    n = _fmt_n(min_foldch)

    # Pre-filter stats — keyed by sample only (same for both bam_types)
    total_frags_path = os.path.join(stats_dir, f"{sample}.total_frags.txt")
    total_frags = open(total_frags_path).read().strip() if os.path.exists(total_frags_path) else "NA"

    trim_log_path = os.path.join(trim_log_dir, f"{sample}.trim.log")
    clean_reads = parse_bbduk_log(trim_log_path) if os.path.exists(trim_log_path) else "NA"

    align_rate_path = os.path.join(stats_dir, f"{sample}.align_rate.txt")
    align_rate = open(align_rate_path).read().strip() if os.path.exists(align_rate_path) else "NA"

    # Post-filter stats — keyed by sample.bam_type
    filt_stats = parse_kv_file(os.path.join(stats_dir, f"{sample}.{bam_type}.filtered_stats.txt"))
    filtered_reads = filt_stats.get("filtered_reads", "NA")

    # peak# from unfiltered file; peaks_nfold# from filtered file
    unfilt_np = os.path.join(macs_dir, f"{sample}.{bam_type}_peaks.narrowPeak")
    peak_total, max_score = parse_narrowpeak(unfilt_np)

    filt_np = os.path.join(macs_dir, f"{sample}.{bam_type}_peaks_filt.narrowPeak")
    peaks_nfold, _ = parse_narrowpeak(filt_np)

    frip_macs_data = parse_kv_file(os.path.join(stats_dir, f"{sample}.{bam_type}.frip_macs.txt"))
    peak_reads = frip_macs_data.get("reads_in_peaks_macs", "NA")
    frip_val   = frip_macs_data.get("frip_macs", "NA")
    frip       = round(float(frip_val) * 100, 2) if frip_val != "NA" else "NA"

    frip_macs_filt_data = parse_kv_file(os.path.join(stats_dir, f"{sample}.{bam_type}.frip_macs_filt.txt"))
    peak_reads_filt = frip_macs_filt_data.get("reads_in_peaks_macs_filt", "NA")
    frip_filt_val   = frip_macs_filt_data.get("frip_macs_filt", "NA")
    frip_top_n      = round(float(frip_filt_val) * 100, 2) if frip_filt_val != "NA" else "NA"

    return {
        "sample":          f"{sample}.{bam_type}",
        "total_frags":     total_frags,
        "clean_reads":     clean_reads,
        "align_rate%":     align_rate,
        "filtered_reads":  filtered_reads,
        "peak#":           peak_total,
        f"peaks_{n}fold#": peaks_nfold,
        "max peak score":  max_score,
        "peak reads#":      peak_reads,
        "peak reads filt#": peak_reads_filt,
        "FRiP_score":       frip,
        "FRiP_top_n_fold":  frip_top_n,
    }


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


def _fmt_html(col, val, i_cols, p_cols):
    if val == "NA":
        return "NA"
    if col in i_cols:
        return f"{int(val):,}"
    if col in p_cols:
        return f"{val}%"
    return str(val)


def write_html(rows, cols, i_cols, p_cols, logo_b64_map, logo_rc_b64_map, factorbook_logo_map, out_path):
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
        lines.append("<tr>")
        for col in cols:
            val = row.get(col, "NA")
            css = ' class="na"' if val == "NA" else ""
            lines.append(f"  <td{css}>{_fmt_html(col, val, i_cols, p_cols)}</td>")
        b64 = logo_b64_map.get(row["sample"])
        if b64:
            lines.append(f'  <td><img class="motif" src="data:image/png;base64,{b64}" alt="motif logo"/></td>')
        else:
            lines.append('  <td class="na">NA</td>')
        b64_rc = logo_rc_b64_map.get(row["sample"])
        if b64_rc:
            lines.append(f'  <td><img class="motif" src="data:image/png;base64,{b64_rc}" alt="motif RC logo"/></td>')
        else:
            lines.append('  <td class="na">NA</td>')
        tf = row["sample"].rsplit(".", 1)[0]
        fb_b64 = factorbook_logo_map.get(tf)
        if fb_b64:
            lines.append(f'  <td><img class="motif" src="data:image/png;base64,{fb_b64}" alt="factorbook motif"/></td>')
        else:
            lines.append('  <td class="na">NA</td>')
        lines.append("</tr>")

    lines += ["</tbody></table>", _HTML_SCRIPT, "</body></html>"]

    with open(out_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def run(treatment_samples, bam_types, trim_log_dir, stats_dir, macs_dir, meme_dir,
        min_foldch, factorbook_dir, tsv_out, html_out):
    n      = _fmt_n(min_foldch)
    cols   = make_cols(n)
    i_cols = int_cols(n)
    p_cols = PCT_COLS

    rows = [
        build_row(s, bt, trim_log_dir, stats_dir, macs_dir, min_foldch)
        for s in treatment_samples
        for bt in bam_types
    ]

    with open(tsv_out, "w") as out:
        out.write("\t".join(cols) + "\n")
        for row in rows:
            out.write("\t".join(str(row.get(c, "NA")) for c in cols) + "\n")

    logo_b64_map = {
        f"{s}.{bt}": logo_to_base64(os.path.join(meme_dir, f"{s}.{bt}", "summits", "logo1.png"))
        for s in treatment_samples
        for bt in bam_types
    }
    logo_rc_b64_map = {
        f"{s}.{bt}": logo_to_base64(os.path.join(meme_dir, f"{s}.{bt}", "summits", "logo_rc1.png"))
        for s in treatment_samples
        for bt in bam_types
    }
    factorbook_logo_map = {
        s: logo_to_base64(os.path.join(factorbook_dir, f"{s}.logo.png"))
        for s in treatment_samples
    }
    write_html(rows, cols, i_cols, p_cols, logo_b64_map, logo_rc_b64_map, factorbook_logo_map, html_out)


def main():
    sm = snakemake  # noqa: F821 — injected by Snakemake
    run(
        treatment_samples = list(sm.params.treatment_samples),
        bam_types         = list(sm.params.bam_types),
        trim_log_dir      = sm.params.trim_log_dir,
        macs_dir          = sm.params.macs_dir,
        meme_dir          = sm.params.meme_dir,
        stats_dir         = sm.params.stats_dir,
        min_foldch        = float(sm.params.min_foldch),
        factorbook_dir    = sm.params.factorbook_dir,
        tsv_out           = sm.output.tsv,
        html_out          = sm.output.html,
    )


def cli_main():
    parser = argparse.ArgumentParser(
        description="Regenerate the DAP-seq QC report from completed pipeline output."
    )
    parser.add_argument("out_dir", help="Pipeline output directory (contains stats/, MACS/, meme/, logs/bbduk/)")
    args = parser.parse_args()

    out_dir      = os.path.abspath(args.out_dir)
    stats_dir    = os.path.join(out_dir, "stats")
    macs_dir     = os.path.join(out_dir, "MACS")
    meme_dir     = os.path.join(out_dir, "meme")
    trim_log_dir = os.path.join(out_dir, "logs", "bbduk")
    tsv_out      = os.path.join(stats_dir, "report.tsv")
    html_out     = os.path.join(stats_dir, "report.html")

    samples, bam_types = _detect_samples_and_bam_types(stats_dir)

    if os.path.exists(_DEFAULT_FB_TSV) and os.path.exists(_DEFAULT_FB_MEME):
        _make_factorbook_logos(samples, _DEFAULT_FB_TSV, _DEFAULT_FB_MEME,
                               os.path.join(out_dir, "factorbook"))

    # Try to read min_foldch from config; fall back to 2.0
    min_foldch = 2.0
    config_path = os.path.join(_SCRIPT_DIR, '..', '..', 'config', 'config.yaml')
    if os.path.exists(config_path):
        with open(config_path) as fh:
            for line in fh:
                if 'min_foldch' in line:
                    try:
                        min_foldch = float(line.split(':')[1].strip())
                    except (IndexError, ValueError):
                        pass
                    break

    run(
        treatment_samples = samples,
        bam_types         = bam_types,
        trim_log_dir      = trim_log_dir,
        stats_dir         = stats_dir,
        macs_dir          = macs_dir,
        meme_dir          = meme_dir,
        min_foldch        = min_foldch,
        factorbook_dir    = os.path.join(out_dir, "factorbook"),
        tsv_out           = tsv_out,
        html_out          = html_out,
    )
    print(f"Wrote {tsv_out}")
    print(f"Wrote {html_out}")


if "snakemake" in dir():
    main()
elif __name__ == "__main__":
    cli_main()
