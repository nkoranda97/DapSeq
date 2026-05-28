"""Shared utilities for rendering motif logo PNGs with logomaker."""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import logomaker

_MEME_DEFAULT_COLORS = {"A": "008000", "C": "0000FF", "G": "FFB300", "T": "FF0000"}


def _render_logo_with_logomaker(ppm, out_path, base_colors=None):
    raw = {k: v for k, v in (base_colors or {}).items() if v}
    color_scheme = {b: '#' + raw.get(b, _MEME_DEFAULT_COLORS[b]).lstrip('#') for b in 'ACGT'}

    df    = pd.DataFrame(ppm, columns=['A', 'C', 'G', 'T'])
    ic_df = logomaker.transform_matrix(df, from_type='probability', to_type='information')

    fig, ax = plt.subplots(figsize=(max(2, len(ppm) * 0.35), 2.5))
    logomaker.Logo(ic_df, ax=ax, color_scheme=color_scheme)
    ax.set_xticks(range(len(ppm)))
    ax.set_xticklabels(range(1, len(ppm) + 1), fontsize=7)
    ax.set_ylabel('bits', fontsize=8)
    ax.spines[['top', 'right']].set_visible(False)
    plt.tight_layout()
    fig.savefig(out_path, format='png', dpi=150, bbox_inches='tight')
    plt.close(fig)


def parse_meme_ppm(path, motif_index=0):
    """Return PPM (list of [A,C,G,T] rows) for the nth motif (0-based) in a MEME text file."""
    motifs_seen = -1
    in_target = False
    in_matrix = False
    ppm = []

    with open(path) as fh:
        for line in fh:
            if line.startswith('MOTIF '):
                motifs_seen += 1
                if motifs_seen > motif_index:
                    break
                in_target = (motifs_seen == motif_index)
                in_matrix = False
            elif in_target and 'letter-probability matrix' in line:
                in_matrix = True
            elif in_target and in_matrix:
                try:
                    vals = [float(x) for x in line.split()]
                    if len(vals) == 4:
                        ppm.append(vals)
                    else:
                        in_matrix = False
                except ValueError:
                    in_matrix = False

    return ppm if ppm else None
