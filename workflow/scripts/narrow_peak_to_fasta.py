"""
narrow_peak_to_fasta.py

Make a fasta file for motif search within peak sequences.

Snakemake inputs:
    snakemake.input.narrowpeak  — filtered narrowPeak file
    snakemake.input.genome      — genome FASTA

Snakemake output:
    snakemake.output[0]         — output FASTA file

Snakemake params:
    snakemake.params.maxpeaks   — max peaks to use (prioritizes highest fold-change)
    snakemake.params.extend_bp  — bp around summit (int), or "all" for full peak
    snakemake.params.fimocoords — bool; write chr:start-end headers for FIMO
"""

from os import path

import pandas as pd
from Bio import SeqIO


def get_peak_seq(start, stop, record):
    return str(record[start:stop].seq)


def narrow_peak_to_fasta(narrowpeak, genome, outfile, maxpeaks, extend_bp, fimocoords):
    chr_records = SeqIO.to_dict(SeqIO.parse(genome, "fasta"))

    peak_col_dtypes = {
        "chromosome": str,
        "start": int,
        "stop": int,
        "peak_name": str,
        "display_score": float,
        "unused": object,
        "fold-change": float,
        "-log10qvalue": float,
        "-log10pvalue": float,
        "relative_summit": int,
    }

    peaks = pd.read_csv(
        narrowpeak, sep="\t", header=None,
        names=list(peak_col_dtypes.keys()), dtype=peak_col_dtypes
    )
    peaks = peaks.sort_values(by=["fold-change"], ascending=False)

    if peaks.empty:
        open(outfile, "w").close()
        return

    if extend_bp == "all":
        foldch_maxes = peaks.groupby(["chromosome", "start", "stop"])["fold-change"].transform(max)
        peaks = peaks.loc[peaks["fold-change"] == foldch_maxes]

    if maxpeaks is not None:
        peaks = peaks.head(maxpeaks)

    with open(outfile, "w", encoding="utf-8") as f:
        for i in peaks.index:
            peak_chr  = peaks.loc[i, "chromosome"]
            peak_name = peaks.loc[i, "peak_name"]
            start_bp  = peaks.loc[i, "start"]
            stop_bp   = peaks.loc[i, "stop"]
            fold_ch   = round(peaks.loc[i, "fold-change"], 2)
            q_score   = round(peaks.loc[i, "-log10qvalue"], 2)

            if extend_bp != "all":
                abs_summit = start_bp + peaks.loc[i, "relative_summit"]
                start_bp   = abs_summit - int(extend_bp)
                stop_bp    = abs_summit + int(extend_bp)

            try:
                chr_record = chr_records[peak_chr]
            except KeyError:
                raise KeyError(
                    f"Chromosome '{peak_chr}' not found in genome FASTA. "
                    "Check that chromosome naming matches between peaks and genome "
                    "(e.g. 'chr1' vs '1')."
                )
            SEQ = get_peak_seq(start=start_bp, stop=stop_bp, record=chr_record)

            if fimocoords:
                header = f">{peak_chr}:{start_bp + 1}-{stop_bp}"
            else:
                header = f">{peak_name}_foldch={fold_ch}_qscore={q_score}_loc={peak_chr}:{start_bp+1}-{stop_bp}"

            f.write(header + "\n" + SEQ + "\n")


if "snakemake" in dir():
    narrow_peak_to_fasta(
        snakemake.input.narrowpeak,    # noqa: F821
        snakemake.input.genome,        # noqa: F821
        snakemake.output[0],           # noqa: F821
        snakemake.params.maxpeaks,     # noqa: F821
        snakemake.params.extend_bp,    # noqa: F821
        snakemake.params.fimocoords,   # noqa: F821
    )
