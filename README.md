# DAP-seq Pipeline

Snakemake pipeline for DAP-seq analysis. Based on the original notebook pipeline by [ndu-bioinfo](https://github.com/ndu-bioinfo/Dap-Seq-pipeline).

## Setup

Edit `config/config.yaml` before running:

| Key | Description |
|---|---|
| `samples` | Sample names with paths to R1 and R2 fastq files |
| `input_control` | Name of the input/control sample. Set to `null` to run without a control |
| `genome_ref` | Path to reference genome FASTA |
| `genome_dir` | Directory of per-chromosome FASTAs (for GEM) |
| `genome_size` | Effective genome size string for MACS2 (e.g. `"1.0e8"`) |
| `threads` | Threads per job |
| `slurm_partition` | SLURM partition |
| `slurm_account` | SLURM account |

## Running

```bash
snakemake --profile profile
```
