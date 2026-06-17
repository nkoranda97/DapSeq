"""
narrow_peak_to_fasta.py

Make a fasta file for motif search within peak sequences.

Snakemake inputs:
    snakemake.input.narrowpeak  — filtered narrowPeak file
    snakemake.input.genome      — genome FASTA

Snakemake output:
    snakemake.output[0]         — output FASTA file

Snakemake params:
    snakemake.params.maxpeaks                  — max peaks to use (prioritizes highest fold-change)
    snakemake.params.extend_bp                 — bp around summit (int), or "all" for full peak
    snakemake.params.fimocoords                — bool; write chr:start-end headers for FIMO
    snakemake.params.complexity_filter_enabled — bool; skip low-complexity sequences
    snakemake.params.min_entropy               — minimum 3-mer Shannon entropy in bits
"""

import sys
from collections import Counter
from math import log2

import pandas as pd
from Bio import SeqIO


def get_peak_seq(start, stop, record):
    return str(record[start:stop].seq)


def _triplet_entropy(seq: str) -> float:
    """Shannon entropy (bits) of the 3-mer distribution.

    Length-stable and bounded [0, log2(64)=6]. Low values indicate
    low-complexity / repetitive sequence (perfect OR degenerate),
    because few distinct triplets dominate. Triplets containing
    non-ACGT characters (e.g. N) are ignored; an all-N window
    therefore scores 0 and is filtered.
    """
    seq = seq.upper()
    counts = Counter(
        seq[i:i + 3]
        for i in range(len(seq) - 2)
        if set(seq[i:i + 3]) <= set("ACGT")
    )
    total = sum(counts.values())
    if total == 0:
        return 0.0
    return -sum((c / total) * log2(c / total) for c in counts.values())


def narrow_peak_to_fasta(
    narrowpeak,
    genome,
    outfile,
    maxpeaks,
    extend_bp,
    fimocoords,
    filter_chroms=None,
    complexity_filter_enabled=False,
    min_entropy=3.0,
):
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
    if filter_chroms:
        for chrom in filter_chroms:
            count = (peaks["chromosome"] == chrom).sum()
            print(f"Filtered {count} peaks on chromosome {chrom}", file=sys.stderr)
        peaks = peaks[~peaks["chromosome"].isin(filter_chroms)]

    peaks = peaks.sort_values(by=["fold-change"], ascending=False)

    if peaks.empty:
        open(outfile, "w").close()
        return

    if extend_bp == "all":
        foldch_maxes = peaks.groupby(["chromosome", "start", "stop"])["fold-change"].transform(max)
        peaks = peaks.loc[peaks["fold-change"] == foldch_maxes]

    # When complexity filtering is disabled, pre-limit to maxpeaks as before.
    # When enabled, defer limiting so we can fill up to maxpeaks after skipping repeats.
    if not complexity_filter_enabled and maxpeaks is not None:
        peaks = peaks.head(maxpeaks)

    accepted = 0
    low_complexity_removed = 0

    with open(outfile, "w", encoding="utf-8") as f:
        for i in peaks.index:
            if maxpeaks is not None and accepted >= maxpeaks:
                break

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

            start_bp = max(0, start_bp)
            stop_bp = min(len(chr_record), stop_bp)
            if stop_bp <= start_bp:
                print(
                    f"Skipping invalid interval for {peak_name}: "
                    f"{peak_chr}:{start_bp + 1}-{stop_bp}",
                    file=sys.stderr,
                )
                continue

            SEQ = get_peak_seq(start=start_bp, stop=stop_bp, record=chr_record)

            if complexity_filter_enabled and _triplet_entropy(SEQ) < min_entropy:
                low_complexity_removed += 1
                print(f"Low-complexity filtered: {peak_name}", file=sys.stderr)
                continue

            if fimocoords:
                header = f">{peak_chr}:{start_bp + 1}-{stop_bp}"
            else:
                header = f">{peak_name}_foldch={fold_ch}_qscore={q_score}_loc={peak_chr}:{start_bp+1}-{stop_bp}"

            f.write(header + "\n" + SEQ + "\n")
            accepted += 1

    if complexity_filter_enabled:
        print(
            f"Complexity filter: {low_complexity_removed} removed, "
            f"{accepted} kept",
            file=sys.stderr,
        )


if "snakemake" in dir():
    sys.stderr = open(snakemake.log[0], "w")  # noqa: F821
    narrow_peak_to_fasta(
        snakemake.input.narrowpeak,                   # noqa: F821
        snakemake.input.genome,                       # noqa: F821
        snakemake.output[0],                          # noqa: F821
        snakemake.params.maxpeaks,                    # noqa: F821
        snakemake.params.extend_bp,                   # noqa: F821
        snakemake.params.fimocoords,                  # noqa: F821
        snakemake.params.filter_chroms,               # noqa: F821
        snakemake.params.complexity_filter_enabled,   # noqa: F821
        snakemake.params.min_entropy,                 # noqa: F821
    )
