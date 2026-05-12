# DAP-seq Pipeline

This is an implementation of the JGI DAP-seq analysis pipeline configured to run in Snakemake with SLURM.

## Basic Setup

We are using `Pixi` for environment management. Install with:

```sh
curl -fsSL https://pixi.sh/install.sh | sh
```

### Running the pipeline

Because multiple users share this installation, always pass `--directory` pointing to your own output directory. Snakemake stores its working state (`.snakemake/`) in the directory you specify — without this flag it defaults to the pipeline root, which is shared, and a second user trying to run simultaneously will be locked out.

```sh
module load apptainer
pixi run snakemake \
    --snakefile /path/to/pipeline/workflow/Snakefile \
    --configfile /path/to/your/configfile.yaml \
    --directory /your/output/dir
```

`--directory` must match `output_dir` in your config file.

### One-time reference index setup (admin)

The Bowtie2 index must be built once per genome reference before any user runs the pipeline. Once built, Snakemake skips the indexing step automatically for all subsequent runs.

```sh
pixi run snakemake bowtie2_index \
    --snakefile /path/to/pipeline/workflow/Snakefile \
    --configfile /path/to/admin/configfile.yaml \
    --directory /tmp/index_build
```

Run this for each new genome reference. Do not run the pipeline without a pre-built index — two simultaneous first-time runs will race to write the same index files.

## Default Configuration

Additional arguments can be added to any process with `extra: "--parameter example"` syntax. Example:

```yaml
macs3:
 format: BAMPE # default
 min_foldch: 2.0 # default
 extra: "--qvalue 0.01" # extra argument
```

Only modified, required, or most relevant parameters are shown. For program default parameters, refer to its documentation (linked).

Default configs are kept in `./config/config.yaml`. You can override any variable in this config file by declaring it in your experiment's config file.

### bbduk

[BBDuk](https://bbmap.org/tools/bbduk)

- `adaptors: null` — Uses built-in BBTools adapters by default.
- `max_frags: null` — Does not do any sub-sampling of fragments by default.
- `k: 21`: int — Kmer lengths shorter than `k` will not be found.
- `mink: 11`: int — Look for shorter kmers at read tips down to this length.
- `ktrim: r`: f|r|l — Trim reads to remove bases matching reference kmers, plus all bases to the right.
   - `f`: don't trim
   - `r`: trim right
   - `l`: trim left
- `qtrim: r`: rl|f|r|l|w — Trim read ends to remove bases with quality below `trimq`.
   - `rl`: both ends
   - `f`: neither end
   - `r`: right only
   - `l`: left only
   - `w`: sliding window
- `trimq: 6.0`: float — Regions with average quality below this threshold will be trimmed.
- `maq: 10`: int — Remove reads with mapping quality below this value.
- `ow: t`: bool — Grant permission to overwrite files.

### Bowtie2

[Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)

### Samtools view

[samtools view](https://www.htslib.org/doc/samtools-view.html)

- `mapq: 10.0`: float — Quality cutoff value.

### Bamcoverage

[bamCoverage](https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html)

- `normalize_using: RPGC`: RPKM|CPM|BPM|RPGC|None — Normalization method for reads per bin.
   - `RPKM`: Reads Per Kilobase per Million
   - `CPM`: Counts Per Million
   - `BPM`: Bins Per Million
   - `RPGC`: reads per genomic content
   - `None`
- `bin_size: 1`: int — Size of bins in bases for the output bigwig/bedgraph file.
- `ignore_duplicates: true`: bool — If true, reads with the same orientation and start position are counted only once.
- `max_fragment_length: 600`: int — Maximum fragment length for read/pair inclusion.
- `extendReads: true`: bool — If true, each read is extended without exception.

### Bamcompare

[bamCompare](https://deeptools.readthedocs.io/en/develop/content/tools/bamCompare.html)

- `bin_size: 80`: int — Size of bins in bases for the output bigwig/bedgraph file.
- `operation: ratio`: log2|reciprocal_ratio|ratio|subtract|add|mean|first|second — Comparison operation between samples.
   - `log2`: log2 ratio of the two samples
   - `reciprocal_ratio`: negative inverse of ratio when ratio < 0
   - `ratio`: ratio of the two samples
   - `subtract`: difference between samples
   - `add`: sum of two samples
   - `mean`: average of two samples
   - `first`: first sample only
   - `second`: second sample only
- `scaleFactorsMethod: SES`: readCount|SES|None — Method used to scale the samples.
- `n: 1000.0`: float — Number of samplings taken from the genome to compute the scaling factors.

### MACS3

[MACS3](https://macs3-project.github.io/MACS/)

- `format: BAMPE`: ELAND|BED|ELANDMULTI|ELANDEXPORT|SAM|BAM|BOWTIE|BAMPE|BEDPR|FRAG|AUTO — Format of tag file.
- `min_foldch: 2.0`
- `keep_dup: 1`: int|string — Keep this number of duplicate tags at the exact same position. Also accepts `auto` or `all`.

## Scripts

- `maxpeaks: 100`
- `summit_extend: 30`

### MEME

[MEME](https://meme-suite.org/meme/doc/meme.html)

- `nmotifs: 2`: int — MEME will stop searching for motifs after `nmotifs` have been found.
- `minw: 8`: int — Search for motifs with a minimum width of `minw`.
- `maxw: 32`: int — Search for motifs with a maximum width of `maxw`.
- `mod: anr`: oops|zoops|anr — Controls the distribution of motif occurrences.
   - `oops`: exactly one occurrence per sequence
   - `zoops`: zero or one occurrence per sequence
   - `anr`: any number of non-overlapping occurrences per sequence
- `maxsize: 10000000`: int — Maximum allowed dataset size.

### FIMO

[FIMO](https://meme-suite.org/meme/doc/fimo.html?man_type=web)

- `thresh: 1.0e-5` — Only report results with a p-value less than `thresh`.

### HOMER

[HOMER](http://homer.ucsd.edu/homer/ngs/annotation.html)
