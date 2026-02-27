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

```bash
pip install snakemake snakemake-executor-plugin-slurm

# if using sample environment for analysis
pip install matplotlib pandas
```

## Running

```bash
# uses config/config.yaml
bash run

# use a project-specific config outside the pipeline directory
bash run /path/to/project/config.yaml

# dry run
bash run -n
bash run /path/to/project/config.yaml -n
```
