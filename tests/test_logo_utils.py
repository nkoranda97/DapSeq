"""
Tests for workflow/scripts/logo_utils.py: reading a motif's probability
matrix out of MEME text output, and reverse-complementing it.
"""

import logo_utils as m


MEME_TXT = """\
MEME version 5.5.5

MOTIF AAGT MEME-1	width =   2  sites =  10  llr = 50  E-value = 1.0e-010
--------------------------------------------------------------------------------
	Motif AAGT MEME-1 position-specific probability matrix
--------------------------------------------------------------------------------
letter-probability matrix: alength= 4 w= 2 nsites= 10 E= 1.0e-010
 0.700000  0.100000  0.100000  0.100000
 0.000000  0.000000  0.900000  0.100000
--------------------------------------------------------------------------------

MOTIF CCTG MEME-2	width =   3  sites =   8  llr = 40  E-value = 2.0e-005
letter-probability matrix: alength= 4 w= 3 nsites= 8 E= 2.0e-005
 0.100000  0.800000  0.050000  0.050000
 0.000000  0.000000  0.000000  1.000000
 0.250000  0.250000  0.250000  0.250000

Time  1.23 secs.
"""


def _meme(tmp_path):
    p = tmp_path / "meme.txt"
    p.write_text(MEME_TXT)
    return str(p)


def test_parse_meme_ppm_reads_first_motif_and_stops_at_separator(tmp_path):
    assert m.parse_meme_ppm(_meme(tmp_path)) == [
        [0.7, 0.1, 0.1, 0.1],
        [0.0, 0.0, 0.9, 0.1],
    ]


def test_parse_meme_ppm_selects_motif_by_index(tmp_path):
    ppm = m.parse_meme_ppm(_meme(tmp_path), motif_index=1)
    assert ppm == [
        [0.1, 0.8, 0.05, 0.05],
        [0.0, 0.0, 0.0, 1.0],
        [0.25, 0.25, 0.25, 0.25],
    ]


def test_parse_meme_ppm_returns_none_for_missing_motif(tmp_path):
    assert m.parse_meme_ppm(_meme(tmp_path), motif_index=5) is None


def test_parse_meme_ppm_returns_none_when_meme_found_no_motifs(tmp_path):
    p = tmp_path / "empty.txt"
    p.write_text("MEME version 5.5.5\n\nTime 0.01 secs.\n")
    assert m.parse_meme_ppm(str(p)) is None


def test_rc_ppm_reverses_positions_and_swaps_complementary_bases():
    ppm = [[0.7, 0.1, 0.15, 0.05], [0.0, 0.2, 0.3, 0.5]]
    assert m.rc_ppm(ppm) == [[0.5, 0.3, 0.2, 0.0], [0.05, 0.15, 0.1, 0.7]]
