import re
import itertools
from Bio import SeqIO


def tandem_filter(fasta_file: str, out_file: str, k: int = 6, k_max: int = 3) -> None:
    """
    Remove sequences with tandem repeats or runs of Ns.

    A sequence is removed if:
      - The most frequent k-mer appears >= k_max times, or
      - There is more than 1 N in the sequence
    """
    kmers = ["".join(c) for c in itertools.product("ACGT", repeat=k)]

    kept = 0
    removed = 0

    with open(fasta_file) as f_in, open(out_file, "w") as f_out:
        for record in SeqIO.parse(f_in, "fasta"):
            record = record.upper()
            seq = str(record.seq)

            counts  = [seq.count(kmer) for kmer in kmers]
            n_peaks = max(counts)
            n_ns    = seq.count("N")

            if n_peaks >= k_max or n_ns > 1:
                removed += 1
            else:
                SeqIO.write(record, f_out, "fasta")
                kept += 1

    print(f"tandem_filter: kept {kept}, removed {removed}", flush=True)


if "snakemake" in dir():  # noqa: F821
    tandem_filter(
        snakemake.input[0],   # noqa: F821
        snakemake.output[0],  # noqa: F821
        k=snakemake.params.k,
        k_max=snakemake.params.k_max,
    )