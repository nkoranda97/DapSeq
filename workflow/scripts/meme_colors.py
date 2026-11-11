"""meme.base_colors: validation and the MEME custom alphabet built from it.

Imported by workflow/rules/common.smk and motifs.smk. Kept free of Snakemake
objects so the logic is unit-testable.
"""

import re

BASES = ("A", "C", "G", "T")
_HEX = re.compile(r"#?[0-9A-Fa-f]{6}")

# MEME's standard DNA alphabet (doc/alphabet-format.html, MEME 5.5.9) with the
# core colours left as placeholders. Only upper-case symbols are defined, so
# MEME reads the alphabet case-insensitively, the same as -dna.
_ALPHABET = (
    'ALPHABET "DNA" DNA-LIKE\n'
    'A "Adenine" {A} ~ T "Thymine" {T}\n'
    'C "Cytosine" {C} ~ G "Guanine" {G}\n'
    "U = T\nR = AG\nY = CT\nK = GT\nM = AC\nS = CG\nW = AT\n"
    "B = CGT\nD = GAT\nH = ACT\nV = ACG\nN = ACGT\nX = ACGT\n"
)
MEME_DEFAULT_COLORS = {"A": "008000", "C": "0000FF", "G": "FFB300", "T": "FF0000"}


def base_colors(value):
    """meme.base_colors as {base: "RRGGBB"} for the bases it sets ({} if unset).

    Raises ValueError for anything that is not a mapping of A/C/G/T to a
    quoted six-digit hex colour. YAML turns unquoted colours into other
    values: 000000 becomes the number 0, and #008000 starts a comment, so the
    key is left empty.
    """
    if value is None:
        return {}
    if not isinstance(value, dict):
        raise ValueError(
            f"meme.base_colors must map bases (A, C, G, T) to colours, got {value!r}"
        )
    unknown = sorted(str(k) for k in value if k not in BASES)
    if unknown:
        raise ValueError(
            f"meme.base_colors has unknown base(s) {unknown}; use A, C, G and T"
        )
    colors = {}
    for base, color in value.items():
        if not isinstance(color, str) or not _HEX.fullmatch(color):
            raise ValueError(
                f"meme.base_colors.{base} must be a quoted six-digit hex colour such as "
                f'"008000" or "#008000", got {color!r}. Unquoted, YAML reads 000000 as '
                "a number and #008000 as a comment."
            )
        colors[base] = color.lstrip("#").upper()
    return colors


def meme_alphabet(colors):
    """MEME custom alphabet text: the standard DNA alphabet with *colors*
    (from base_colors) replacing MEME's defaults for the bases they set."""
    return _ALPHABET.format(**{**MEME_DEFAULT_COLORS, **colors})
