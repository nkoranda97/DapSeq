# DAP-seq Pipeline

Snakemake implementation of the JGI DAP-seq analysis pipeline, designed to run on SLURM clusters using Apptainer.

## Setup

Install [Pixi](https://pixi.sh):

```sh
curl -fsSL https://pixi.sh/install.sh | sh
```

## Running

```sh
module load apptainer
pixi run snakemake \
    --snakefile /path/to/pipeline/workflow/Snakefile \
    --configfile /path/to/your/config.yaml \
    --directory /your/output/dir
```

`--directory` must match `output_dir` in your config and must be outside the pipeline directory. Snakemake stores its working state (`.snakemake/`) there — without this flag it defaults to the pipeline root, which is shared across users.

## Configuration

Copy `config.yaml` from the repo root and fill in the required fields. Defaults for all tool parameters are in `config/config.yaml`.

### Required fields

```yaml
author: "Your Name"  # optional

samples:
  sample_name:
    r1: /path/to/sample.R1.fastq.gz
    r2: /path/to/sample.R2.fastq.gz  # null for single-end
  control:
    r1: /path/to/control.R1.fastq.gz
    r2: null

input_control: control  # must match a key in samples; null to run without control

output_dir: /path/to/output/
genome_ref: /path/to/genome.fa
genome_size: "3000000000"
gene_annotation: /path/to/annotation.gtf  # null to skip HOMER annotation

slurm_partition: "caslake"
slurm_account: "pi-yourlab"
```

### Notable options

| Option | Default | Notes |
|---|---|---|
| `aligner` | `bowtie2` | `bowtie2` or `bwa_mem2` |
| `samtools.mapq` | `30` | MAPQ filter applied after alignment |
| `chrom_filter` | `[]` | Chromosomes excluded from peaks before MEME (e.g. `[chrEBV]`) |
| `macs3.min_foldch` | `2.0` | Peak fold-change filter |
| `macs3.format` | `BAMPE` | Set to `BAM` for single-end data |
| `meme.base_colors` | unset | Optional hex color overrides for sequence logos |

### Extra arguments

Any tool accepts an `extra` field for additional flags:

```yaml
macs3:
  extra: "--qvalue 0.01"
```

## Index build (admin, one-time per genome)

Build the index once per genome reference before any user runs the pipeline. Snakemake skips this step automatically on subsequent runs.

**Bowtie2:**
```sh
pixi run snakemake bowtie2_index \
    --snakefile /path/to/pipeline/workflow/Snakefile \
    --configfile /path/to/admin/config.yaml \
    --directory /tmp/index_build
```

**bwa-mem2:** *(config must include `aligner: bwa_mem2`)*
```sh
pixi run snakemake bwa_mem2_index \
    --snakefile /path/to/pipeline/workflow/Snakefile \
    --configfile /path/to/admin/config.yaml \
    --directory /tmp/index_build
```

Do not run two simultaneous first-time index builds against the same genome — they will race to write the same index files.
