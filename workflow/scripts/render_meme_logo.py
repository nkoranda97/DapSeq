"""
Render the top MEME motif logo from a completed meme.txt output file.

Called via Snakemake's script: directive after the meme shell rule has run.
Produces an empty logo file when meme.txt is absent or contains no motifs.
"""

import os
import sys
from pathlib import Path

sm = snakemake  # noqa: F821

sys.path.insert(0, os.path.dirname(__file__))
from logo_utils import _render_logo_with_logomaker, parse_meme_ppm  # noqa: E402

txt_path = str(sm.input[0])
logo_path = str(sm.output[0])

if os.path.exists(txt_path) and os.path.getsize(txt_path) > 0:
    ppm = parse_meme_ppm(txt_path)
    if ppm:
        _render_logo_with_logomaker(ppm, logo_path, sm.params.base_colors)
        sys.exit(0)

Path(logo_path).touch()
