import re
import itertools
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import SeqIO


def tandem_filter(fasta_file, out_file, k=6, k_max=3):
    """Remove sequences with tandem repeats or runs of Ns."""
    kmers = ["".join(c) for c in itertools.product("ACGT", repeat=k)]

    with open(fasta_file) as f:
        fasta_dict = SeqIO.to_dict(rec.upper() for rec in SeqIO.parse(f, "fasta"))

    kept = []
    removed = []
    peak_locs = []
    top_kmers = []
    n_locs = []

    for record in fasta_dict.values():
        seq = str(record.seq)
        counts = [seq.count(kmer) for kmer in kmers]
        top = kmers[counts.index(max(counts))]
        peaks = [m.start() for m in re.finditer(top, seq)]
        ns = [m.start() for m in re.finditer("N", seq)]
        if len(peaks) >= k_max or len(ns) > 1:
            removed.append(record)
            peak_locs.append(peaks)
            top_kmers.append(top)
            n_locs.append(ns)
        else:
            kept.append(record)

    # Plot removed sequences for QC
    n_rows = len(removed) // 3 + 1
    fig, axes = plt.subplots(n_rows, 3, figsize=(12, max(3, len(removed) // 10)))
    for ax in axes.ravel():
        ax.set_axis_off()
    axes = axes.ravel()
    plt.suptitle(fasta_file.split("/")[-1].split(".")[0])
    for i, record in enumerate(removed):
        axes[i].plot([0, len(record)], [0, 0], "lightblue")
        if len(n_locs[i]) > 1:
            axes[i].plot([n_locs[i][0], n_locs[i][-1]], [0, 0], "grey")
        for x in peak_locs[i]:
            axes[i].plot([x, x + k], [0, 0], "r")
        if k_max > 2:
            axes[i].set_title(f"{top_kmers[i]}={len(peak_locs[i])}", fontsize=10)

    with open(out_file, "w") as f:
        SeqIO.write(kept, f, "fasta")


tandem_filter(
    snakemake.input[0],
    snakemake.output[0],
    k=snakemake.params.k,
    k_max=snakemake.params.k_max,
)