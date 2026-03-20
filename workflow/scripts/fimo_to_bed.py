"""
fimo_to_bed.py — Step 8.1

Converts a FIMO TSV output file to a standard 6-column BED file.

Snakemake input:  fimo/{sample}/fimo.tsv
Snakemake output: fimo/{sample}/fimo.bed

FIMO TSV column layout (0-indexed):
    0  motif_id       → BED name
    1  motif_alt_id   (discarded)
    2  sequence_name  → BED chrom
    3  start          → BED start
    4  stop           → BED end
    5  strand         → BED strand
    6  score          → BED score
    7+ p-value, q-value, matched_sequence  (discarded)
"""


def fimo_to_bed(tsv_path: str, bed_path: str) -> None:
    """
    Parse a FIMO TSV file and write a 6-column BED file.
    Skips the header row, blank lines, and trailing comment lines (# …).
    """
    with open(tsv_path) as fin, open(bed_path, "w") as fout:
        for line in fin:
            if line.startswith("#") or line.startswith("motif_id") or line.strip() == "":
                continue
            fields = line.strip().split("\t")
            chrom  = fields[2]
            start  = fields[3]
            stop   = fields[4]
            name   = fields[0]
            score  = fields[6]
            strand = fields[5]
            fout.write(f"{chrom}\t{start}\t{stop}\t{name}\t{score}\t{strand}\n")


if "snakemake" in dir():  # noqa: F821
    fimo_to_bed(snakemake.input[0], snakemake.output[0])  # noqa: F821
