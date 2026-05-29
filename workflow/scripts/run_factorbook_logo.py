"""
Look up a TF in the Factorbook ChIP-seq MEME catalog and render its best motif logo.

Called via Snakemake's script: directive; uses the snakemake object for I/O.
Produces an empty PNG when no matching motif is found so downstream rules can
depend on the output unconditionally.
"""

import os
import re
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(__file__))
from logo_utils import _render_logo_with_logomaker  # noqa: E402

sm = snakemake  # noqa: F821

tf = sm.wildcards.sample.upper()

tf_motif_ids = []
seen = set()
with open(sm.input.tsv) as fh:
    header = fh.readline().rstrip("\n").split("\t")
    ti = header.index("target")
    ai = header.index("dataset_accession")
    ni = header.index("name")
    for line in fh:
        parts = line.rstrip("\n").split("\t")
        if parts[ti].upper() == tf:
            mid = f"{parts[ai]}_{parts[ni]}"
            if mid not in seen:
                seen.add(mid)
                tf_motif_ids.append(mid)

if not tf_motif_ids:
    Path(str(sm.output[0])).touch()
else:
    tf_id_set = set(tf_motif_ids)
    motif_ppms = {}
    motif_nsites = {}
    current = None
    in_matrix = False
    with open(sm.input.meme) as fh:
        for line in fh:
            if line.startswith("MOTIF "):
                mid = line.split()[1]
                current = mid if mid in tf_id_set else None
                in_matrix = False
                if current:
                    motif_ppms[current] = []
                    motif_nsites[current] = 0
            elif current is not None:
                if "letter-probability matrix" in line:
                    m = re.search(r"nsites= *(\d+)", line)
                    if m:
                        motif_nsites[current] = int(m.group(1))
                    in_matrix = True
                elif in_matrix:
                    try:
                        vals = [float(x) for x in line.split()]
                        if len(vals) == 4:
                            motif_ppms[current].append(vals)
                        else:
                            in_matrix = False
                    except ValueError:
                        in_matrix = False

    motif_id = max(
        (mid for mid in tf_motif_ids if motif_ppms.get(mid)),
        key=lambda mid: motif_nsites.get(mid, 0),
        default=None,
    )

    os.makedirs(os.path.dirname(str(sm.output[0])), exist_ok=True)
    ppm = motif_ppms.get(motif_id) if motif_id else None
    if ppm:
        _render_logo_with_logomaker(ppm, str(sm.output[0]), sm.params.base_colors)
    else:
        Path(str(sm.output[0])).touch()
