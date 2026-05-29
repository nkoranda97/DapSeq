"""
Run MEME motif discovery on a FASTA file and render the top motif logo.

Called via Snakemake's script: directive; uses the snakemake object for I/O.
Produces empty output files when the input FASTA is absent or empty so that
downstream rules can depend on the outputs unconditionally.
"""

import os
import sys
from pathlib import Path

sm = snakemake  # noqa: F821


def _touch(*paths):
    for p in paths:
        Path(str(p)).touch()


if not (os.path.exists(str(sm.input[0])) and os.path.getsize(str(sm.input[0])) > 0):
    _touch(sm.output.txt, sm.output.xml, sm.output.logo)
else:
    base_colors = sm.params.base_colors
    active = {k: v for k, v in base_colors.items() if v}

    if active:
        defaults = {"A": "008000", "C": "0000FF", "G": "FFB300", "T": "FF0000"}

        def _hex(k):
            return active.get(k, defaults[k]).lstrip("#")

        alph_file = os.path.join(sm.params.outdir, "custom_alph.txt")
        os.makedirs(sm.params.outdir, exist_ok=True)
        with open(alph_file, "w") as fh:
            fh.write('ALPHABET "DNA"\n')
            fh.write(f'A ADENINE {_hex("A")} ~ T THYMINE {_hex("T")}\n')
            fh.write(f'C CYTOSINE {_hex("C")} ~ G GUANINE {_hex("G")}\n')
            fh.write("END ALPHABET\n")
        alph_arg = f"-alph {alph_file}"
    else:
        alph_arg = "-dna"

    shell_cmd = (
        f"meme {sm.input[0]} -oc {sm.params.outdir} "
        f"{alph_arg} -revcomp "
        f"-mod {sm.params.mod} -nmotifs {sm.params.nmotifs} "
        f"-minw {sm.params.minw} -maxw {sm.params.maxw} "
        f"-maxsize {sm.params.maxsize} -p {sm.threads} -nostatus {sm.params.extra} "
        f"2>{sm.log[0]}"
    )
    os.system(shell_cmd)

    sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
    from logo_utils import _render_logo_with_logomaker, parse_meme_ppm

    ppm = parse_meme_ppm(sm.output.txt) if os.path.getsize(sm.output.txt) > 0 else None
    if ppm:
        _render_logo_with_logomaker(ppm, str(sm.output.logo), base_colors)
    else:
        Path(str(sm.output.logo)).touch()
