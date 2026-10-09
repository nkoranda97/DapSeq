"""
Look up a TF in the Factorbook ChIP-seq MEME catalog and render its best motif logo.

Called via Snakemake's script: directive; uses the snakemake object for I/O.
Produces an empty PNG when no matching motif is found so downstream rules can
depend on the output unconditionally.

The TF is the sample's optional `tf:` config key. Without it, the sample name
is tried as-is, then the part before the first "_" or "-" (CTCF_rep1 ->
CTCF). Matching against the catalog's `target` column is case-insensitive.
"""

import os
import re
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(__file__))
from logo_utils import _render_logo_with_logomaker  # noqa: E402

_REQUIRED_COLUMNS = ("target", "dataset_accession", "name")


def tf_candidates(sample, tf=None):
    """Upper-cased TF names to look up, in order of preference."""
    if tf:
        return [str(tf).upper()]
    full = sample.upper()
    prefix = re.split(r"[_\-]", full, maxsplit=1)[0]
    return [full] if prefix in ("", full) else [full, prefix]


def find_motif_ids(tsv_path, candidates):
    """Return (matched TF name, motif IDs) for the first candidate present in
    the catalog's `target` column, or (None, []) when none is."""
    with open(tsv_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        missing = [c for c in _REQUIRED_COLUMNS if c not in header]
        if missing:
            raise ValueError(
                f"Factorbook TSV {tsv_path} is missing column(s) {missing}; "
                f"expected {list(_REQUIRED_COLUMNS)} in its header"
            )
        ti, ai, ni = (header.index(c) for c in _REQUIRED_COLUMNS)
        ids_by_target = {}
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            ids = ids_by_target.setdefault(parts[ti].upper(), [])
            mid = f"{parts[ai]}_{parts[ni]}"
            if mid not in ids:
                ids.append(mid)

    for name in candidates:
        if ids_by_target.get(name):
            return name, ids_by_target[name]
    return None, []


def best_motif_ppm(meme_path, motif_ids):
    """Probability matrix of the motif in *motif_ids* with the most sites, or None."""
    wanted = set(motif_ids)
    ppms, nsites = {}, {}
    current = None
    in_matrix = False
    with open(meme_path) as fh:
        for line in fh:
            if line.startswith("MOTIF "):
                mid = line.split()[1]
                current = mid if mid in wanted else None
                in_matrix = False
                if current:
                    ppms[current] = []
                    nsites[current] = 0
            elif current is not None:
                if "letter-probability matrix" in line:
                    match = re.search(r"nsites= *(\d+)", line)
                    if match:
                        nsites[current] = int(match.group(1))
                    in_matrix = True
                elif in_matrix:
                    try:
                        vals = [float(x) for x in line.split()]
                    except ValueError:
                        in_matrix = False
                        continue
                    if len(vals) == 4:
                        ppms[current].append(vals)
                    else:
                        in_matrix = False

    best = max(
        (mid for mid in motif_ids if ppms.get(mid)),
        key=lambda mid: nsites.get(mid, 0),
        default=None,
    )
    return ppms[best] if best else None


def main(sm):
    out = str(sm.output[0])
    os.makedirs(os.path.dirname(out), exist_ok=True)
    candidates = tf_candidates(sm.wildcards.sample, sm.params.tf)
    with open(sm.log[0], "w") as log:
        name, ids = find_motif_ids(sm.input.tsv, candidates)
        ppm = best_motif_ppm(sm.input.meme, ids) if ids else None
        if ppm:
            print(f"Matched Factorbook target {name} (tried {candidates})", file=log)
            _render_logo_with_logomaker(ppm, out, sm.params.base_colors)
        else:
            print(f"No Factorbook motif for {candidates}; writing an empty logo", file=log)
            Path(out).touch()


if "snakemake" in dir():
    main(snakemake)  # noqa: F821
