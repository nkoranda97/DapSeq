"""
generate_report.py — Standalone DAP-seq HTML report generator.

Usage:
    uv run python scripts/generate_report.py \
        --results /path/to/results \
        --output  dapseq_report.html \
        --control CONTROL_SAMPLE_NAME   # optional
"""

import argparse
import base64
import io
import math
import re
import subprocess
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import pysam
from matplotlib_venn import venn3
from PIL import Image, ImageOps


# ── Helpers ───────────────────────────────────────────────────────────────────

MIME = {
    '.png': 'image/png',
    '.jpg': 'image/jpeg',
    '.gif': 'image/gif',
    '.svg': 'image/svg+xml',
}


def fig_to_base64(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight", dpi=150)
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


def img_tag(encoded: str, style: str = "") -> str:
    return f'<img src="data:image/png;base64,{encoded}" style="{style}">'


def iframe_from_file(path: Path, height: str = "700px") -> str:
    if not path.exists():
        return f'<p class="missing">File not found: {path}</p>'
    encoded = base64.b64encode(path.read_bytes()).decode()
    return (
        f'<iframe src="data:text/html;base64,{encoded}" '
        f'width="100%" height="{height}" style="border:none;"></iframe>'
    )


def load_logo(path: Path):
    """Load a motif logo PNG with auto-cropped whitespace. Returns np.ndarray."""
    import numpy as np
    img = Image.open(str(path)).convert("RGB")
    inv = ImageOps.invert(img)
    bbox = inv.getbbox()
    if bbox:
        pad = 8
        bbox = (
            max(0, bbox[0] - pad), max(0, bbox[1] - pad),
            min(img.width, bbox[2] + pad), min(img.height, bbox[3] + pad),
        )
        img = img.crop(bbox)
    return img


def logo_to_base64(path: Path) -> str | None:
    if not path.exists():
        return None
    img = load_logo(path)
    buf = io.BytesIO()
    img.save(buf, format="PNG")
    return base64.b64encode(buf.getvalue()).decode()


def gem_self_contained(results: Path, sample: str) -> str | None:
    """Read GEM results HTML, follow the meta-refresh redirect, and inline all
    relative <img src="..."> as base64 data URIs so the page renders correctly
    inside a data: or blob: iframe."""
    top = results / "GEM" / sample / f"{sample}.results.htm"
    if not top.exists():
        return None
    m = re.search(
        r'content=["\']0;\s*url=([^"\']+)["\']',
        top.read_text(errors='replace'),
        re.I,
    )
    if m:
        content_path = top.parent / m.group(1)
        base_dir = content_path.parent
    else:
        content_path, base_dir = top, top.parent
    html = content_path.read_text(errors='replace')

    def inline_src(match):
        src = match.group(1)
        if src.startswith(('http', 'data:', '//')):
            return match.group(0)
        img = (base_dir / src).resolve()
        if img.exists():
            mime = MIME.get(img.suffix.lower(), 'image/png')
            enc = base64.b64encode(img.read_bytes()).decode()
            return f'src="data:{mime};base64,{enc}"'
        return match.group(0)

    return re.sub(r"""src=["']([^"']+)["']""", inline_src, html)


def blob_iframe(path: Path, iframe_id: str, height: str = "700px") -> str:
    """Embed an HTML file as a Blob URL iframe so JavaScript executes freely.
    Chrome 60+ blocks JS in data: URI iframes; Blob URLs are same-origin."""
    if not path.exists():
        return f'<p class="missing">File not found: {path}</p>'
    encoded = base64.b64encode(path.read_bytes()).decode()
    return (
        f'<iframe id="{iframe_id}" width="100%" height="{height}" style="border:none;"></iframe>\n'
        f'<script>(function(){{\n'
        f'  var b64="{encoded}";\n'
        f'  var bin=atob(b64);\n'
        f'  var arr=new Uint8Array(bin.length);\n'
        f'  for(var i=0;i<bin.length;i++) arr[i]=bin.charCodeAt(i);\n'
        f'  var blob=new Blob([arr],{{type:"text/html"}});\n'
        f'  document.getElementById("{iframe_id}").src=URL.createObjectURL(blob);\n'
        f'}})();</script>'
    )


def ensure_logo_png(logo_eps: Path) -> Path:
    logo_png = logo_eps.with_suffix(".png")
    if logo_eps.exists() and not logo_png.exists():
        subprocess.run(
            ["gs", "-dNOPAUSE", "-dBATCH", "-sDEVICE=png16m", "-r600",
             f"-sOutputFile={logo_png}", str(logo_eps)],
            capture_output=True,
        )
    return logo_png


# ── Data collection ───────────────────────────────────────────────────────────

def parse_trimmomatic(results: Path, sample: str):
    log = results / "logs" / "trimmomatic" / f"{sample}.log"
    if not log.exists():
        return None, None
    text = log.read_text()
    # PE format
    m_raw = re.search(r"Input Read Pairs: (\d+)", text)
    m_clean = re.search(r"Both Surviving: (\d+)", text)
    if m_raw and m_clean:
        return int(m_raw.group(1)), int(m_clean.group(1))
    # SE format
    m = re.search(r"Input Reads: (\d+) Surviving: (\d+)", text)
    if m:
        return int(m.group(1)), int(m.group(2))
    return None, None


def parse_bowtie2_log(results: Path, sample: str):
    """Returns (total_reads, unique_aligned, overall_rate_pct) from bowtie2 log."""
    log = results / "logs" / "bowtie2" / f"{sample}.log"
    if not log.exists():
        return None, None, None
    text = log.read_text()
    # Total reads (first number on "reads; of these:" line)
    m_total = re.search(r"^(\d+) reads; of these:", text, re.MULTILINE)
    # SE: "aligned exactly 1 time" — PE: "aligned concordantly exactly 1 time"
    m_unique = re.search(r"(\d+) \(\d+\.\d+%\) aligned (?:concordantly )?exactly 1 time", text)
    # Overall rate
    m_rate = re.search(r"([\d.]+)% overall alignment rate", text)
    total = int(m_total.group(1)) if m_total else None
    unique = int(m_unique.group(1)) if m_unique else None
    rate = float(m_rate.group(1)) if m_rate else None
    return total, unique, rate


def parse_bwa_log(results: Path, sample: str):
    """Returns (total_reads, unique_aligned, overall_rate_pct) from bwa mem log.
    BWA mem doesn't report alignment rate; returns None for those fields."""
    log = results / "logs" / "bwa" / f"{sample}.log"
    if not log.exists():
        return None, None, None
    text = log.read_text()
    m_total = re.search(r"\[M::mem_process_seqs\] Processed (\d+) reads", text)
    total = int(m_total.group(1)) if m_total else None
    return total, None, None


def detect_aligner(results: Path, sample: str) -> str | None:
    """Detect which aligner was used by checking which log directory has a file."""
    for aligner in ("bowtie2", "bwa"):
        if (results / "logs" / aligner / f"{sample}.log").exists():
            return aligner
    return None


def parse_alignment_log(results: Path, sample: str):
    """Dispatch to the correct aligner log parser based on auto-detection."""
    aligner = detect_aligner(results, sample)
    if aligner == "bowtie2":
        return parse_bowtie2_log(results, sample)
    if aligner == "bwa":
        return parse_bwa_log(results, sample)
    return None, None, None


def get_alignment_stats(results: Path, sample: str):
    bam = results / "bam" / f"{sample}.bam"
    if not bam.exists():
        return None, None
    stats = pysam.flagstat(str(bam))
    m_mapped = re.search(r"(\d+) \+ \d+ mapped", stats)
    m_paired = re.search(r"(\d+) \+ \d+ with itself and mate mapped", stats)
    mapped = int(m_mapped.group(1)) if m_mapped else None
    unique_pairs = int(m_paired.group(1)) // 2 if m_paired else None
    return mapped, unique_pairs


def get_peak_stats(results: Path, sample: str):
    peak_file = results / "MACS" / f"{sample}_peaks.narrowPeak"
    if not peak_file.exists():
        return None, None, None
    folds = [float(line.split("\t")[6]) for line in open(peak_file)]
    return (
        len(folds),
        sum(1 for f in folds if f > 5),
        round(max(folds), 4) if folds else None,
    )


def get_frip(results: Path, sample: str, mapped_reads):
    bam_path = results / "bam" / f"{sample}.bam"
    peak_file = results / "MACS" / f"{sample}_peaks.narrowPeak"
    if not bam_path.exists() or not peak_file.exists() or not mapped_reads:
        return None, None
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    reads_in_peaks = sum(
        bam.count(chrom, int(start), int(end))
        for line in open(peak_file)
        for chrom, start, end, *_ in [line.strip().split("\t")]
    )
    bam.close()
    return reads_in_peaks, round(100 * reads_in_peaks / mapped_reads, 2)


def build_stats_df(results: Path, samples: list[str]) -> pd.DataFrame:
    rows = []
    for sample in samples:
        gem_file = results / "GEM" / sample / f"{sample}.GEM_events.txt"
        raw_reads, _ = parse_trimmomatic(results, sample)
        clean_reads, unique_aligned, mapping_pct = parse_alignment_log(results, sample)
        mapped, _ = get_alignment_stats(results, sample)
        macs_peaks, min5fold, max_score = get_peak_stats(results, sample)
        gem_count = (
            max(0, sum(1 for _ in open(gem_file)) - 1) if gem_file.exists() else None
        )
        peak_reads, frip = get_frip(results, sample, mapped)
        rows.append({
            "sample":          sample,
            "raw_reads":       raw_reads,
            "clean_reads":     clean_reads,
            "mapping_ratio%":  mapping_pct,
            "unique_aligned":  unique_aligned,
            "MACS_peaks":      macs_peaks,
            "min5fold_peaks":  min5fold,
            "max_peak_score":  max_score,
            "GEM_peaks":       gem_count,
            "peak_reads":      peak_reads,
            "FRiP%":           frip,
        })
    return pd.DataFrame(rows)


def compute_venn_sets(results: Path, samples: list[str]) -> dict:
    venn_set = {}
    for sample in samples:
        try:
            def _read(p):
                return pd.read_table(p, header=None)
            df_MACS      = _read(results / "compare_bed" / f"{sample}.MACS.bed")
            df_GEM       = _read(results / "compare_bed" / f"{sample}.GEM.bed")
            df_Motif     = _read(results / "fimo" / sample / "fimo.bed")
            df_MACS_peak = _read(results / "compare_bed" / f"{sample}.MACS_peak.bed")
            df_GEM_peak  = _read(results / "compare_bed" / f"{sample}.GEM_peak.bed")
            df_compare   = _read(results / "compare_bed" / f"{sample}.compare_peak.bed")

            All_overlay     = len(set(df_MACS_peak[3]).intersection(df_GEM_peak[3]))
            MACS_only       = len(df_MACS) - len(df_compare) - len(df_MACS_peak) + All_overlay
            GEM_only        = len(df_GEM) - len(df_compare) - len(df_GEM_peak) + All_overlay
            MOTIF_only      = len(df_Motif) - len(df_MACS_peak) - len(df_GEM_peak) + All_overlay
            GEM_MACS_only   = len(df_compare) - All_overlay
            MACS_MOTIF_only = len(df_MACS_peak) - All_overlay
            GEM_MOTIF_only  = len(df_GEM_peak) - All_overlay

            venn_set[sample] = [
                MACS_only, GEM_only, GEM_MACS_only,
                MOTIF_only, MACS_MOTIF_only, GEM_MOTIF_only, All_overlay,
            ]
        except Exception as e:
            print(f"  [warn] Venn for {sample}: {e}")
    return venn_set


# ── Section builders ──────────────────────────────────────────────────────────

def section_stats(results: Path, samples: list[str]) -> str:
    print("  Building sample statistics table…")
    df = build_stats_df(results, samples)

    int_cols = [
        "raw_reads", "clean_reads", "unique_aligned",
        "MACS_peaks", "min5fold_peaks", "GEM_peaks", "peak_reads",
    ]

    # Convert EPS logos if needed
    for sample in samples:
        ensure_logo_png(results / "meme" / sample / "logo1.eps")

    def _logo_cell(sample):
        enc = logo_to_base64(results / "meme" / sample / "logo1.png")
        if enc:
            return img_tag(enc, "max-width:180px;")
        return "—"

    df_disp = df.copy()
    for col in int_cols:
        if col in df_disp:
            df_disp[col] = df_disp[col].apply(
                lambda x: f"{int(x):,}" if pd.notna(x) and x is not None else "—"
            )
    df_disp["motif"] = [_logo_cell(s) for s in df_disp["sample"]]

    table_html = df_disp.to_html(escape=False, index=False, classes="stats-table")
    return f'<section id="stats"><h2>1. Sample Statistics</h2>{table_html}</section>'


def section_multiqc(results: Path) -> str:
    print("  Embedding MultiQC report…")
    path = results / "multiqc_report.html"
    if not path.exists():
        iframe = f'<p class="missing">File not found: {path}</p>'
    else:
        html = path.read_bytes()
        # Firefox throws NS_ERROR_FAILURE when blob:null-origin pages access
        # localStorage/sessionStorage (null-origin restriction). MultiQC's dark
        # mode toggler calls localStorage.getItem() without a try-catch, which
        # crashes the script and leaves window.bootstrap undefined. Inject an
        # in-memory polyfill before any other script runs.
        polyfill = (
            b'<script>(function(){'
            b'function _ms(){var s={};return{'
            b'getItem:function(k){return s.hasOwnProperty(k)?s[k]:null;},'
            b'setItem:function(k,v){s[k]=String(v);},'
            b'removeItem:function(k){delete s[k];},'
            b'clear:function(){s={};},'
            b'get length(){return Object.keys(s).length;},'
            b'key:function(i){return Object.keys(s)[i]||null;}}};'
            b'try{localStorage.getItem("__t__");}catch(e){'
            b'Object.defineProperty(window,"localStorage",'
            b'{value:_ms(),configurable:true,writable:true});}'
            b'try{sessionStorage.getItem("__t__");}catch(e){'
            b'Object.defineProperty(window,"sessionStorage",'
            b'{value:_ms(),configurable:true,writable:true});}'
            b'})();</script>'
        )
        html = html.replace(b'<head>', b'<head>' + polyfill, 1)
        encoded = base64.b64encode(html).decode()
        iframe = (
            '<iframe id="multiqc-iframe" width="100%" height="800px" style="border:none;"></iframe>\n'
            f'<script>(function(){{'
            f'var b64="{encoded}";'
            f'var bin=atob(b64);'
            f'var arr=new Uint8Array(bin.length);'
            f'for(var i=0;i<bin.length;i++)arr[i]=bin.charCodeAt(i);'
            f'var blob=new Blob([arr],{{type:"text/html"}});'
            f'document.getElementById("multiqc-iframe").src=URL.createObjectURL(blob);'
            f'}})();</script>'
        )
    return (
        '<section id="multiqc"><h2>2. QC — MultiQC</h2>'
        + iframe
        + "</section>"
    )


def _venn_figure(results: Path, samples: list[str], venn_set: dict, logo_key_fn, title: str) -> str:
    """Render a grid of logo + venn pairs, return base64 PNG."""
    valid = [s for s in samples if s in venn_set]
    if not valid:
        return "<p class='missing'>No data available.</p>"
    nrows = math.ceil(len(valid) / 2)
    fig, axes = plt.subplots(nrows, 4, figsize=(28, nrows * 8), squeeze=False)
    axes_flat = axes.ravel()
    for ax in axes_flat:
        ax.set_axis_off()
    for ix, sample in enumerate(valid):
        logo_path = logo_key_fn(sample)
        ax_logo = axes_flat[2 * ix]
        ax_logo.set_title(sample)
        if logo_path.exists():
            ax_logo.imshow(load_logo(logo_path))
            ax_logo.set_axis_on()
            ax_logo.set_xticks([])
            ax_logo.set_yticks([])
            ax_logo.spines[:].set_visible(False)

        ax_venn = axes_flat[2 * ix + 1]
        ax_venn.set_title(sample)
        ax_venn.set_axis_on()
        venn3(subsets=venn_set[sample], set_labels=("MACS", "GEM", "Motif"), ax=ax_venn)
    fig.tight_layout()
    return img_tag(fig_to_base64(fig), "max-width:100%;")


def section_venn_first(results: Path, samples: list[str], venn_set: dict) -> str:
    print("  Building first-round Venn diagrams…")
    content = _venn_figure(
        results, samples, venn_set,
        logo_key_fn=lambda s: results / "meme" / s / "logo1.png",
        title="Peak Caller & Motif Overlap",
    )
    return (
        '<section id="venn-first">'
        "<h2>3. Peak Caller &amp; Motif Overlap</h2>"
        + content
        + "</section>"
    )


def section_venn_intersection(results: Path, samples: list[str], venn_set: dict) -> str:
    print("  Building intersection Venn diagrams…")
    # Convert EPS logos for intersection meme dirs
    for sample in samples:
        ensure_logo_png(results / "meme" / f"{sample}-intersection" / "logo1.eps")
    content = _venn_figure(
        results, samples, venn_set,
        logo_key_fn=lambda s: results / "meme" / f"{s}-intersection" / "logo1.png",
        title="Refined Motifs (Intersection)",
    )
    return (
        '<section id="venn-intersect">'
        "<h2>4. Refined Motifs (Intersection)</h2>"
        + content
        + "</section>"
    )


def section_per_sample(results: Path, samples: list[str]) -> str:
    print("  Embedding per-sample reports…")
    parts = ['<section id="per-sample"><h2>5. Per-Sample Detail</h2>']
    for sample in samples:
        anchor = f"sample-{sample.replace(' ', '_')}"
        parts.append(f'<div id="{anchor}" class="sample-block">')
        parts.append(f"<h3>{sample}</h3>")

        parts.append("<h4>MEME Motif Report</h4>")
        parts.append(blob_iframe(results / "meme" / sample / "meme.html", f"meme-{anchor}"))

        parts.append("<h4>FIMO Motif Scan</h4>")
        parts.append(blob_iframe(results / "fimo" / sample / "fimo.html", f"fimo-{anchor}"))

        parts.append("<h4>GEM Results</h4>")
        gem_html = gem_self_contained(results, sample)
        if gem_html:
            enc = base64.b64encode(gem_html.encode()).decode()
            parts.append(
                f'<iframe src="data:text/html;base64,{enc}" '
                f'width="100%" height="700px" style="border:none;"></iframe>'
            )
        else:
            parts.append('<p class="missing">GEM results not found.</p>')

        parts.append("</div>")
    parts.append("</section>")
    return "\n".join(parts)


# ── TOC ───────────────────────────────────────────────────────────────────────

def build_toc(samples: list[str]) -> str:
    items = [
        ('<a href="#stats">1. Sample Statistics</a>', ""),
        ('<a href="#multiqc">2. QC — MultiQC</a>', ""),
        ('<a href="#venn-first">3. Peak Caller &amp; Motif Overlap</a>', ""),
        ('<a href="#venn-intersect">4. Refined Motifs (Intersection)</a>', ""),
        ('<a href="#per-sample">5. Per-Sample Detail</a>', ""),
    ]
    sub_items = "".join(
        f'<li><a href="#sample-{s.replace(" ", "_")}">{s}</a></li>'
        for s in samples
    )
    li = "".join(f"<li>{link}</li>" for link, _ in items)
    li = li.replace(
        "<li><a href=\"#per-sample\">5. Per-Sample Detail</a></li>",
        f'<li><a href="#per-sample">5. Per-Sample Detail</a><ul>{sub_items}</ul></li>',
    )
    return f'<nav id="toc"><h2>Contents</h2><ul>{li}</ul></nav>'


# ── HTML assembly ─────────────────────────────────────────────────────────────

CSS = """
* { box-sizing: border-box; }
body { margin: 0; font-family: sans-serif; display: flex; }
#toc {
    position: sticky; top: 0; height: 100vh; overflow-y: auto;
    width: 240px; min-width: 220px; padding: 16px;
    background: #f4f4f4; border-right: 1px solid #ccc;
    font-size: 13px;
}
#toc h2 { font-size: 14px; margin-top: 0; }
#toc ul { padding-left: 16px; margin: 4px 0; }
#toc li { margin: 3px 0; }
#toc a { text-decoration: none; color: #2255aa; }
#toc a:hover { text-decoration: underline; }
#content { flex: 1; padding: 24px 32px; max-width: 1200px; overflow-x: auto; }
h2 { border-bottom: 2px solid #2255aa; padding-bottom: 4px; }
h3 { color: #444; }
h4 { color: #666; margin-top: 20px; }
.sample-block { margin-bottom: 48px; border-top: 1px solid #ddd; padding-top: 16px; }
.stats-table { border-collapse: collapse; width: 100%; font-size: 13px; }
.stats-table th, .stats-table td {
    border: 1px solid #ccc; padding: 6px 10px; text-align: right;
}
.stats-table th { background: #e8eef8; text-align: center; }
.stats-table td:first-child { text-align: left; }
.missing { color: #999; font-style: italic; }
section { margin-bottom: 48px; }
"""


def assemble_html(toc: str, sections: list[str], generated_at: str) -> str:
    body = "\n".join(sections)
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>DAP-seq Report</title>
<style>{CSS}</style>
</head>
<body>
{toc}
<div id="content">
<h1>DAP-seq Analysis Report</h1>
<p style="color:#888;font-size:13px;">Generated: {generated_at}</p>
{body}
</div>
</body>
</html>"""


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Generate a standalone HTML DAP-seq report."
    )
    parser.add_argument("--results", required=True, help="Path to Snakemake results directory")
    parser.add_argument("--output", default="dapseq_report.html", help="Output HTML file")
    parser.add_argument("--control", default=None, help="Control sample name to exclude")
    args = parser.parse_args()

    results = Path(args.results)
    if not results.is_dir():
        raise SystemExit(f"ERROR: results directory not found: {results}")

    # Discover samples
    all_samples = sorted(
        p.stem.replace("_summits", "")
        for p in (results / "MACS").glob("*_summits.bed")
    )
    exclude = [args.control, "TF1", "TF2"]
    samples = [s for s in all_samples if s not in exclude]
    print(f"Found {len(samples)} treatment samples: {samples}")

    if not samples:
        raise SystemExit("ERROR: No treatment samples found — check --results and --control.")

    print("Building report sections…")
    venn_set = compute_venn_sets(results, samples)

    sections = [
        section_stats(results, samples),
        section_multiqc(results),
        section_venn_first(results, samples, venn_set),
        section_venn_intersection(results, samples, venn_set),
        section_per_sample(results, samples),
    ]

    toc = build_toc(samples)
    html = assemble_html(toc, sections, datetime.now().strftime("%Y-%m-%d %H:%M"))

    out = Path(args.output)
    out.write_text(html, encoding="utf-8")
    size_mb = out.stat().st_size / 1_048_576
    print(f"Report written to: {out}  ({size_mb:.1f} MB)")


if __name__ == "__main__":
    main()
