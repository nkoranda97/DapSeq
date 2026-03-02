# DAP-seq Pipeline

Snakemake pipeline for DAP-seq analysis. Based on the original notebook pipeline by [ndu-bioinfo](https://github.com/ndu-bioinfo/Dap-Seq-pipeline).

I also bundled the software into a container and call that in the pipeline.

Claude Opus 4.6 helped a lot with the conversion, I have done my best to keep the actual logic as true to the original pipeline as possible, but there could still be bugs.

## Changes (up for scrutiny)

1. **MEME runs on tandem-filtered file (new) vs unfiltered (old)**

   Old step 8 runs MEME on `.fasta.nodup` (before tandem filtering). Now it runs MEME on `.fasta.filtered.fasta` (after tandem filtering). I think this was a bug in the old pipeline. It runs Tandem_filter but never uses its output.

2. **chr prefix in GEM output — intentional organism difference**

   The old pipeline hardcoded "chr". Now it uses whatever chromosome names GEM outputs. It should work without manual modification.

3. **bedtools intersect flag: `-wo` vs `-wa`**

   Old pipeline uses `-wo` which gives output data. Because we aren't using a notebook, there is no reason to do this (I think) and `-wa` is used instead.

4. **BED column count: 5 vs 6**

   Old pipeline writes 5-column BED (no strand). New `combine_peaks.py` writes 6 columns (with `.` for strand). 6-column is apparently the more standard BED format?

## Setup

Edit `config/config.yaml` before running, or create a project-specific config anywhere on the filesystem:

| Key | Description |
|---|---|
| `samples` | Sample names with paths to R1 and R2 fastq files |
| `input_control` | Name of the input/control sample. Set to `null` to run without a control |
| `output_dir` | Where all outputs will be written |
| `genome_ref` | Path to reference genome FASTA |
| `genome_size` | Effective genome size string for MACS2 (e.g. `"1.0e8"`) |
| `threads` | Threads per job |
| `slurm_partition` | SLURM partition |
| `slurm_account` | SLURM account |

## Requirements

You need snakemake and a slurm plugin for snakemake to run the pipeline. A minimal python env will do the trick, but conda can be used as well.

```bash
pip install snakemake snakemake-executor-plugin-slurm

# if using the same environment for analysis
pip install matplotlib pandas
```

## Running

```bash
# uses config/config.yaml
module load singularity
snakemake --profile profile --configfile /path/to/config/file.yaml

# dry run
snakemake --profile profile --configfile /path/to/config/file.yaml -n
```

If a run fails during execution, you can restart it with the `--rerun-incomplete` flag.

To force a run to restart from the beginning use the `--forceall` flag.

```bash
snakemake --profile profile --configfile /path/to/config/file.yaml --rerun-incomplete

snakemake --profile profile --configfile /path/to/config/file.yaml --forceall
```
