"""
Tests for scripts/fimo_to_bed.py

Logic: parse a FIMO TSV (skipping header, blank lines, and # comments) and
emit a 6-column BED: chrom, start, stop, name(=motif_id), score, strand.
"""

from pathlib import Path

import fimo_to_bed as m

RESULTS_FIMO = Path(__file__).parent.parent / "results" / "fimo"

# Minimal TSV with one header, two data rows, a blank line, and a footer comment
TSV_CONTENT = (
    "motif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tscore\tp-value\tq-value\tmatched_sequence\n"
    "MOTIF1\tALT1\tchr1\t100\t110\t+\t15.5\t1e-6\t0.1\tACGTACGTAC\n"
    "MOTIF1\tALT1\tchr2\t200\t210\t-\t12.3\t5e-5\t0.5\tGGCCTTAAGG\n"
    "\n"
    "# trailing comment — should be ignored\n"
)

EXPECTED_BED = (
    "chr1\t100\t110\tMOTIF1\t15.5\t+\n"
    "chr2\t200\t210\tMOTIF1\t12.3\t-\n"
)


# ── unit tests ────────────────────────────────────────────────────────────────

def test_basic_conversion(tmp_path):
    tsv = tmp_path / "fimo.tsv"
    bed = tmp_path / "fimo.bed"
    tsv.write_text(TSV_CONTENT)

    m.fimo_to_bed(str(tsv), str(bed))

    assert bed.read_text() == EXPECTED_BED


def test_skips_header_and_comments(tmp_path):
    tsv = tmp_path / "fimo.tsv"
    bed = tmp_path / "fimo.bed"
    tsv.write_text(TSV_CONTENT)

    m.fimo_to_bed(str(tsv), str(bed))

    lines = bed.read_text().splitlines()
    assert all(not l.startswith("#") for l in lines)
    assert all(not l.startswith("motif_id") for l in lines)


def test_empty_tsv(tmp_path):
    """Only a header + blank lines → empty BED."""
    tsv = tmp_path / "fimo.tsv"
    bed = tmp_path / "fimo.bed"
    tsv.write_text("motif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tscore\n\n")

    m.fimo_to_bed(str(tsv), str(bed))

    assert bed.read_text() == ""


def test_column_mapping(tmp_path):
    """Verify each BED column comes from the correct TSV column."""
    tsv = tmp_path / "fimo.tsv"
    bed = tmp_path / "fimo.bed"
    tsv.write_text(
        "motif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tscore\tp-value\tq-value\tmatched\n"
        "MY_MOTIF\tALT\tscaffold_7\t42\t55\t-\t9.99\t1e-4\t0.9\tAAAA\n"
    )

    m.fimo_to_bed(str(tsv), str(bed))

    fields = bed.read_text().strip().split("\t")
    assert fields == ["scaffold_7", "42", "55", "MY_MOTIF", "9.99", "-"]


# ── real-data regression ──────────────────────────────────────────────────────

def test_real_tf1_matches_known_good(tmp_path):
    """Running the script on TF1 fimo.tsv must reproduce TF1 fimo.bed exactly."""
    out = tmp_path / "fimo.bed"
    m.fimo_to_bed(
        str(RESULTS_FIMO / "TF1" / "fimo.tsv"),
        str(out),
    )
    assert out.read_text() == (RESULTS_FIMO / "TF1" / "fimo.bed").read_text()
