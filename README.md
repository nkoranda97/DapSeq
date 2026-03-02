# DAP-seq Pipeline

Snakemake pipeline for DAP-seq analysis. Based on the original notebook pipeline by [ndu-bioinfo](https://github.com/ndu-bioinfo/Dap-Seq-pipeline).

## Setup

Edit `config/config.yaml` before running, or create a project-specific config anywhere on the filesystem:

| Key | Description |
|---|---|
| `samples` | Sample names with paths to R1 and R2 fastq files |
| `input_control` | Name of the input/control sample. Set to `null` to run without a control |
| `output_dir` | Where all outputs will be written |
| `genome_ref` | Path to reference genome FASTA |
| `genome_dir` | Directory of per-chromosome FASTAs (for GEM) |
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

If a run fails during execution, you can restart it with the rerun-incomplete flag

To force a run to restart from the beginning use the --forceall flag

```bash
snakemake --profile profile --configfile /path/to/config/file.yaml --rerun-incomplete

snakemake --profile profile --configfile /path/to/config/file.yaml --forceall
```