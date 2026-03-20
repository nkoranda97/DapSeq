"""
dedup_fasta.py — Step 7.2

Removes duplicate FASTA records (same sequence name = same genomic region).
First occurrence of each name is kept; subsequent duplicates are discarded.

Snakemake input:  fasta/{sample}.fasta
Snakemake output: fasta/{sample}.fasta.nodup
"""

from Bio import SeqIO


def dedup_fasta(fasta_in: str, fasta_out: str) -> None:
    """
    Copy records from fasta_in to fasta_out, skipping any whose name has
    already been seen.  Each sequence is written as a single unbroken line
    to preserve the original format.
    """
    seen: set[str] = set()
    with open(fasta_in) as handle, open(fasta_out, "w") as out_handle:
        for rec in SeqIO.parse(handle, "fasta"):
            if rec.name not in seen:
                seen.add(rec.name)
                out_handle.write(">" + rec.name + "\n" + str(rec.seq) + "\n")


if "snakemake" in dir():  # noqa: F821
    dedup_fasta(snakemake.input[0], snakemake.output[0])  # noqa: F821
