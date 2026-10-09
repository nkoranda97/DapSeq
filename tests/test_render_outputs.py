"""
Tests for rendering outputs from finished runs: render_report_html.render()
(embeds logos into report.html, also used to re-render from the CLI) and
render_meme_logo.render() (forward and reverse-complement MEME logos).
"""

import csv
import os

import report as rp
import render_meme_logo as rml
import render_report_html as rrh


MEME_TXT = """\
MEME version 5.5.5

MOTIF AAGT MEME-1	width =   2
letter-probability matrix: alength= 4 w= 2 nsites= 10 E= 1.0e-010
 0.900000  0.033333  0.033333  0.033333
 0.025000  0.025000  0.025000  0.925000
"""


def _report_csv(tmp_path, rows):
    path = tmp_path / "report.csv"
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=rp.make_cols())
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "NA") for c in rp.make_cols()})
    return str(path)


# ── render_meme_logo ─────────────────────────────────────────────────────────

def test_meme_logo_renders_forward_and_reverse_complement(tmp_path):
    txt = tmp_path / "meme.txt"
    txt.write_text(MEME_TXT)
    fwd, rc = tmp_path / "logo1.png", tmp_path / "logo_rc1.png"
    rml.render(str(txt), str(fwd), str(rc))
    assert fwd.stat().st_size > 0 and rc.stat().st_size > 0
    assert fwd.read_bytes() != rc.read_bytes()


def test_meme_logo_writes_empty_files_without_motifs(tmp_path):
    fwd, rc = tmp_path / "logo1.png", tmp_path / "logo_rc1.png"
    rml.render(str(tmp_path / "missing.txt"), str(fwd), str(rc))
    assert fwd.stat().st_size == 0 and rc.stat().st_size == 0


# ── render_report_html ───────────────────────────────────────────────────────

def test_report_embeds_logo_only_for_samples_that_have_one(tmp_path):
    meme_dir = tmp_path / "meme"
    (meme_dir / "TF1" / "summits").mkdir(parents=True)
    txt = tmp_path / "meme.txt"
    txt.write_text(MEME_TXT)
    rml.render(str(txt), str(meme_dir / "TF1" / "summits" / "logo1.png"),
               str(meme_dir / "TF1" / "summits" / "logo_rc1.png"))

    csv_path = _report_csv(tmp_path, [{"sample": "TF1"}, {"sample": "TF2"}])
    html_out = tmp_path / "report.html"
    rrh.render(["TF1", "TF2"], csv_path, str(meme_dir), str(tmp_path / "factorbook"), str(html_out))
    html = html_out.read_text()
    assert html.count("data:image/png;base64,") == 2      # TF1 forward + RC logo
    assert "TF2" in html


def test_report_rerender_tolerates_blank_count_cells(tmp_path):
    # An older or hand-edited report.csv can have blank numeric cells.
    csv_path = _report_csv(tmp_path, [{"sample": "TF1", "num_peaks": "", "mapped_reads": ""}])
    html_out = tmp_path / "report.html"
    rrh.render(["TF1"], csv_path, str(tmp_path / "meme"), str(tmp_path / "fb"), str(html_out))
    assert "TF1" in html_out.read_text()
