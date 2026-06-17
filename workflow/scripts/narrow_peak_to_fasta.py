"""
narrow_peak_to_fasta.py

Make a fasta file for motif search within peak sequences.

Snakemake inputs:
    snakemake.input.narrowpeak  — filtered narrowPeak file
    snakemake.input.genome      — genome FASTA

Snakemake output:
    snakemake.output[0]         — output FASTA file

Snakemake params:
    snakemake.params.maxpeaks            — max peaks to use (prioritizes highest fold-change)
    snakemake.params.extend_bp           — bp around summit (int), or "all" for full peak
    snakemake.params.fimocoords          — bool; write chr:start-end headers for FIMO
    snakemake.params.tandem_filter_enabled — bool; skip sequences with tandem repeats
    snakemake.params.tandem_k            — k-mer length for tandem repeat detection (default 6)
    snakemake.params.tandem_k_max        — max occurrences of top k-mer before flagging (default 3)
"""

import itertools
import re
import sys

import pandas as pd
from Bio import SeqIO


def get_peak_seq(start, stop, record):
    return str(record[start:stop].seq)


def _is_tandem_repeat(seq: str, k: int = 6, k_max: int = 3) -> bool:
    """Return True if seq is dominated by a single k-mer or contains multiple Ns."""
    seq = seq.upper()
    kmers = ["".join(c) for c in itertools.product("ACGT", repeat=k)]
    counts = [seq.count(kmer) for kmer in kmers]
    top_kmer = kmers[counts.index(max(counts))]
    n_occurrences = len(re.findall(top_kmer, seq))
    return n_occurrences >= k_max or seq.count("N") > 1


def narrow_peak_to_fasta(
    narrowpeak,
    genome,
    outfile,
    maxpeaks,
    extend_bp,
    fimocoords,
    filter_chroms=None,
    tandem_filter_enabled=False,
    tandem_k=6,
    tandem_k_max=3,
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

    # When tandem filtering is disabled, pre-limit to maxpeaks as before.
    # When enabled, defer limiting so we can fill up to maxpeaks after skipping repeats.
    if not tandem_filter_enabled and maxpeaks is not None:
        peaks = peaks.head(maxpeaks)

    accepted = 0
    tandem_removed = 0

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

            if tandem_filter_enabled and _is_tandem_repeat(SEQ, tandem_k, tandem_k_max):
                tandem_removed += 1
                print(f"Tandem repeat filtered: {peak_name}", file=sys.stderr)
                continue

            if fimocoords:
                header = f">{peak_chr}:{start_bp + 1}-{stop_bp}"
            else:
                header = f">{peak_name}_foldch={fold_ch}_qscore={q_score}_loc={peak_chr}:{start_bp+1}-{stop_bp}"

            f.write(header + "\n" + SEQ + "\n")
            accepted += 1

    if tandem_filter_enabled:
        print(f"Tandem filter: {tandem_removed} removed, {accepted} kept", file=sys.stderr)


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
        snakemake.params.tandem_filter_enabled,       # noqa: F821
        snakemake.params.tandem_k,                    # noqa: F821
        snakemake.params.tandem_k_max,                # noqa: F821
    )
