with open(snakemake.input[0]) as fin, open(snakemake.output[0], "w") as fout:
    for line in fin:
        if line.startswith("#") or line.startswith("motif_id") or line.strip() == "":
            continue
        fields = line.strip().split("\t")
        chrom = fields[2]
        start = fields[3]
        stop = fields[4]
        name = fields[0]
        score = fields[6]
        strand = fields[5]
        fout.write(f"{chrom}\t{start}\t{stop}\t{name}\t{score}\t{strand}\n")