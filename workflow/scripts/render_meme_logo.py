"""
Render the top MEME motif logo (forward and RC) from a completed meme.txt output file.

Called via Snakemake's script: directive after the meme shell rule has run.
Produces empty files when meme.txt is absent or contains no motifs.
"""

import os
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(__file__))
from logo_utils import _render_logo_with_logomaker, parse_meme_ppm, rc_ppm  # noqa: E402


def render(txt_path, logo_path, logo_rc_path, base_colors=None):
    """Write forward and reverse-complement logos of the first motif, or two
    empty files when there is no motif."""
    if os.path.exists(txt_path) and os.path.getsize(txt_path) > 0:
        ppm = parse_meme_ppm(txt_path)
        if ppm:
            _render_logo_with_logomaker(ppm, logo_path, base_colors)
            _render_logo_with_logomaker(rc_ppm(ppm), logo_rc_path, base_colors)
            return
    Path(logo_path).touch()
    Path(logo_rc_path).touch()


if "snakemake" in dir():
    render(
        str(snakemake.input[0]),           # noqa: F821
        str(snakemake.output.logo),        # noqa: F821
        str(snakemake.output.logo_rc),     # noqa: F821
        snakemake.params.base_colors,      # noqa: F821
    )
