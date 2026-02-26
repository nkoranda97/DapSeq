"""
combine_peaks.py — Step 7
Combines MACS2 summits and GEM events into a single BED file,
and writes per-tool BED files for downstream comparison.

Snakemake inputs:
    macs_summits  — MACS/{sample}_summits.bed
    gem_events    — GEM_BED/{sample}/{sample}.GEM_events.txt

Snakemake outputs:
    combined_bed  — Combined_BED/{sample}.combined.bed
    macs_bed      — compare_bed/{sample}.MACS.bed
    gem_bed       — compare_bed/{sample}.GEM.bed

Snakemake params:
    window_size   — bp window centered on summit (default 200)
    min_score     — minimum fold-enrichment score (default 1)
"""

import pandas as pd

macs_summits = snakemake.input.macs_summits
gem_events   = snakemake.input.gem_events
combined_out = snakemake.output.combined_bed
macs_out     = snakemake.output.macs_bed
gem_out      = snakemake.output.gem_bed
window_size  = snakemake.params.window_size
min_score    = snakemake.params.min_score
half_win     = window_size // 2

# ── MACS2 summits ──────────────────────────────────────────────────────────
# Columns: chr, start, end, name, score (fold-enrichment)
df_bed  = pd.read_table(macs_summits, header=None)
df_macs = df_bed[df_bed[4] > min_score]

macs_chrs   = list(df_macs[0])
macs_peaks  = list(df_macs[1])
macs_names  = ["MACS_" + str(i) for i in range(len(df_macs))]
macs_scores = list(df_macs[4])

with open(macs_out, "w") as f:
    for chrom, pos, name, score in zip(macs_chrs, macs_peaks, macs_names, macs_scores):
        f.write(f"{chrom}\t{pos - half_win}\t{pos + half_win}\t{name}\t{score}\t.\n")

# ── GEM events ─────────────────────────────────────────────────────────────
# Columns include: Position (chrN:pos), Fold, IP, ...
df_GEM      = pd.read_table(gem_events)
df_GEM_filt = df_GEM[df_GEM["Fold"] > min_score]

gem_chrs   = ["chr" + pos.split(":")[0] for pos in df_GEM_filt["Position"]]
gem_peaks  = [int(pos.split(":")[1])    for pos in df_GEM_filt["Position"]]
gem_names  = ["GEM_" + str(i) for i in range(len(df_GEM_filt))]
gem_scores = list(df_GEM_filt["Fold"])

with open(gem_out, "w") as f:
    for chrom, pos, name, score in zip(gem_chrs, gem_peaks, gem_names, gem_scores):
        f.write(f"{chrom}\t{pos - half_win}\t{pos + half_win}\t{name}\t{score}\t.\n")

# ── Combined BED ───────────────────────────────────────────────────────────
all_chrs   = macs_chrs   + gem_chrs
all_peaks  = macs_peaks  + gem_peaks
all_names  = macs_names  + gem_names
all_scores = macs_scores + gem_scores

with open(combined_out, "w") as f:
    for chrom, pos, name, score in zip(all_chrs, all_peaks, all_names, all_scores):
        f.write(f"{chrom}\t{pos - half_win}\t{pos + half_win}\t{name}\t{score}\t.\n")
