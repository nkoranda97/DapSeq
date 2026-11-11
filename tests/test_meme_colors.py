"""
Tests for workflow/scripts/meme_colors.py: meme.base_colors validation and the
MEME custom alphabet built from it.
"""

import pytest
import yaml

import meme_colors as m


def _colors(text):
    """base_colors() of a YAML snippet, as a user's config would parse."""
    return m.base_colors(yaml.safe_load(text))


def test_unset_means_no_custom_colours():
    assert m.base_colors(None) == {}
    assert m.base_colors({}) == {}


def test_quoted_colours_are_normalised():
    assert _colors('{A: "#00aa00", C: "0000FF"}') == {"A": "00AA00", "C": "0000FF"}


@pytest.mark.parametrize("text, shown", [
    ("{A: 000000}", "0"),           # YAML reads unquoted digits as a number
    ("{A: 123456}", "123456"),
    ("A: #008000\n", "None"),       # an unquoted # starts a comment
    ('{A: "#0080"}', "'#0080'"),    # too short
    ('{A: "green"}', "'green'"),
])
def test_unusable_colour_values_are_rejected(text, shown):
    with pytest.raises(ValueError, match=r"meme\.base_colors\.A must be a quoted") as e:
        _colors(text)
    assert f"got {shown}" in str(e.value)


def test_unknown_base_is_rejected():
    with pytest.raises(ValueError, match=r"unknown base\(s\) \['U'\]"):
        _colors('{U: "FF0000"}')


def test_non_mapping_is_rejected():
    with pytest.raises(ValueError, match="must map bases"):
        _colors('"FF0000"')


def test_alphabet_uses_given_colours_and_meme_defaults_for_the_rest():
    text = m.meme_alphabet({"A": "111111"})
    lines = text.splitlines()
    assert lines[0] == 'ALPHABET "DNA" DNA-LIKE'
    assert 'A "Adenine" 111111 ~ T "Thymine" FF0000' in lines
    assert 'C "Cytosine" 0000FF ~ G "Guanine" FFB300' in lines
    # MEME 5.5.9 rejects an END ALPHABET line in this format
    assert not any(line.startswith("END") for line in lines)


def test_alphabet_keeps_the_standard_dna_wildcards():
    lines = m.meme_alphabet({}).splitlines()
    for symbol in ("N = ACGT", "R = AG", "Y = CT", "U = T"):
        assert symbol in lines
