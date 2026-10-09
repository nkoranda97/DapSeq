"""
Real FIMO 5.5.9 output (the version pinned in apptainer_build/pixi.lock),
copied verbatim for the motif_peaks tests.

Generated with apptainer_build/.pixi/envs/default/bin/fimo, scanning two
8-bp motifs (M1 ACGTACGT, M2 GGATCCAA) against FASTA records whose headers
use the pipeline's chr:start-end form. The "new" outputs pass --no-pgc, as
the fimo rule does; the "old mode" outputs pass --parse-genomic-coord, as the
rule did before, which makes FIMO name each hit by chromosome only.
"""

# Hits in three peaks: two on chr1, one on chr2 (a fourth peak, chr3, has none).
FIMO_PEAKS_THREE_HITS = (
    'motif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tscore\tp-value\tq-value\tmatched_sequence\n'
    'M1\tACGTACGT\tchr2:1001-1040\t6\t13\t+\t16\t1.47e-05\t0.000617\tACGTACGT\n'
    'M1\tACGTACGT\tchr1:101-140\t11\t18\t+\t16\t1.47e-05\t0.000617\tACGTACGT\n'
    'M1\tACGTACGT\tchr2:1001-1040\t6\t13\t-\t16\t1.47e-05\t0.000617\tACGTACGT\n'
    'M1\tACGTACGT\tchr1:101-140\t11\t18\t-\t16\t1.47e-05\t0.000617\tACGTACGT\n'
    'M1\tACGTACGT\tchr1:501-540\t21\t28\t+\t16\t1.47e-05\t0.000617\tACGTACGT\n'
    'M1\tACGTACGT\tchr1:501-540\t21\t28\t-\t16\t1.47e-05\t0.000617\tACGTACGT\n'
    '\n'
    '# FIMO (Find Individual Motif Occurrences): Version 5.5.9 compiled on Nov 24 2025 at 04:06:59\n'
    '# The format of this file is described at https://meme-suite.org/meme/doc/fimo-output-format.html#tsv_results.\n'
    '# fimo --thresh 1e-4 --no-pgc --oc new_three motifs.meme three.fa\n'
)

# One peak hit eight times by two different motifs.
FIMO_PEAK_TWO_MOTIFS = (
    'motif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tscore\tp-value\tq-value\tmatched_sequence\n'
    'M2\tGGATCCAA\tchr1:301-360\t18\t25\t+\t16\t1.47e-05\t0.000376\tGGATCCAA\n'
    'M2\tGGATCCAA\tchr1:301-360\t16\t23\t-\t16\t1.47e-05\t0.000376\tGGATCCAA\n'
    'M2\tGGATCCAA\tchr1:301-360\t43\t50\t+\t16\t1.47e-05\t0.000376\tGGATCCAA\n'
    'M2\tGGATCCAA\tchr1:301-360\t41\t48\t-\t16\t1.47e-05\t0.000376\tGGATCCAA\n'
    'M1\tACGTACGT\tchr1:301-360\t6\t13\t+\t16\t1.47e-05\t0.000354\tACGTACGT\n'
    'M1\tACGTACGT\tchr1:301-360\t6\t13\t-\t16\t1.47e-05\t0.000354\tACGTACGT\n'
    'M1\tACGTACGT\tchr1:301-360\t31\t38\t+\t16\t1.47e-05\t0.000354\tACGTACGT\n'
    'M1\tACGTACGT\tchr1:301-360\t31\t38\t-\t16\t1.47e-05\t0.000354\tACGTACGT\n'
    '\n'
    '# FIMO (Find Individual Motif Occurrences): Version 5.5.9 compiled on Nov 24 2025 at 04:06:59\n'
    '# The format of this file is described at https://meme-suite.org/meme/doc/fimo-output-format.html#tsv_results.\n'
    '# fimo --thresh 1e-4 --no-pgc --oc new_multi motifs.meme multi.fa\n'
)

# FIMO ran and matched nothing: no header row, only the trailer.
FIMO_NO_HITS = (
    '\n'
    '# FIMO (Find Individual Motif Occurrences): Version 5.5.9 compiled on Nov 24 2025 at 04:06:59\n'
    '# The format of this file is described at https://meme-suite.org/meme/doc/fimo-output-format.html#tsv_results.\n'
    '# fimo --thresh 1e-4 --no-pgc --oc new_nohit motifs.meme nohit.fa\n'
)

# The same three-peak scan in the old coordinate mode: hits named by chromosome.
FIMO_OLD_MODE_HITS = (
    'motif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tscore\tp-value\tq-value\tmatched_sequence\n'
    'M1\tACGTACGT\tchr1\t111\t118\t+\t16\t1.47e-05\t0.000617\tACGTACGT\n'
    'M1\tACGTACGT\tchr1\t111\t118\t-\t16\t1.47e-05\t0.000617\tACGTACGT\n'
    'M1\tACGTACGT\tchr1\t521\t528\t+\t16\t1.47e-05\t0.000617\tACGTACGT\n'
    'M1\tACGTACGT\tchr1\t521\t528\t-\t16\t1.47e-05\t0.000617\tACGTACGT\n'
    'M1\tACGTACGT\tchr2\t1006\t1013\t+\t16\t1.47e-05\t0.000617\tACGTACGT\n'
    'M1\tACGTACGT\tchr2\t1006\t1013\t-\t16\t1.47e-05\t0.000617\tACGTACGT\n'
    '\n'
    '# FIMO (Find Individual Motif Occurrences): Version 5.5.9 compiled on Nov 24 2025 at 04:06:59\n'
    '# The format of this file is described at https://meme-suite.org/meme/doc/fimo-output-format.html#tsv_results.\n'
    '# fimo --parse-genomic-coord --thresh 1e-4 --oc old_three motifs.meme three.fa\n'
)

# Old coordinate mode, no hits.
FIMO_OLD_MODE_NO_HITS = (
    '\n'
    '# FIMO (Find Individual Motif Occurrences): Version 5.5.9 compiled on Nov 24 2025 at 04:06:59\n'
    '# The format of this file is described at https://meme-suite.org/meme/doc/fimo-output-format.html#tsv_results.\n'
    '# fimo --parse-genomic-coord --thresh 1e-4 --oc old_nohit motifs.meme nohit.fa\n'
)
