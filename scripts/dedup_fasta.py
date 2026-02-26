"""
dedup_fasta.py — Step 7.2
Removes duplicate FASTA records (same sequence name = same genomic region).

Snakemake input:  fasta/{sample}.fasta
Snakemake output: fasta/{sample}.fasta.nodup
"""

from Bio import SeqIO

fasta_in  = snakemake.input[0]
fasta_out = snakemake.output[0]

seen = set()
with open(fasta_in, "r") as handle, open(fasta_out, "w") as out_handle:
    for rec in SeqIO.parse(handle, "fasta"):
        if rec.name not in seen:
            seen.add(rec.name)
            out_handle.write(">" + rec.name + "\n" + str(rec.seq) + "\n")
