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

Edit `config/config.yaml` before running, or create a project-specific config anywhere on the filesystem (see the Known bugs note about config merging).

**Required keys**

| Key | Description |
|---|---|
| `samples` | Sample names with paths to R1/R2 fastq files |
| `input_control` | Control sample name (must match a key in `samples`). `null` to run without a control |
| `output_dir` | Root output directory |
| `genome_ref` | Path to reference genome FASTA |
| `genome_size` | Effective genome size for MACS3 — shorthand (`hs`, `mm`, `ce`, `dm`) or a number (`2.7e9`) |
| `threads` | Threads per job (default: `8`) |
| `slurm_partition` | SLURM partition |
| `slurm_account` | SLURM account |

**Tool parameters** (all have defaults matching the original pipeline)

| Key | Default | Description |
|---|---|---|
| `bwa.k` | `60` | BWA MEM minimum seed length (`-k`) |
| `bwa.B` | `7` | BWA MEM mismatch penalty (`-B`) |
| `bwa.O` | `6` | BWA MEM gap open penalty (`-O`) |
| `trimmomatic.minlen` | `36` | Discard reads shorter than this after trimming |
| `macs3.qvalue` | `0.01` | FDR cutoff for peak calling (`-q`) |
| `macs3.min_fold` | `2` | Min fold-enrichment for shifting model (`-m`) |
| `macs3.max_fold` | `50` | Max fold-enrichment for shifting model (`-m`) |
| `gem.k_min` | `6` | GEM min k-mer length (`--k_min`) |
| `gem.k_max` | `20` | GEM max k-mer length (`--k_max`) |
| `gem.k_seqs` | `600` | GEM sequences for motif training (`--k_seqs`) |
| `gem.fold_cutoff` | `1` | GEM IP/control fold-enrichment threshold — keep in sync with `combine_peaks.min_score` |
| `combine_peaks.window_size` | `80` | bp window centered on each summit for sequence extraction |
| `combine_peaks.min_score` | `1` | Minimum MACS3 summit score to include a peak |
| `fimo.thresh` | `1e-5` | FIMO p-value threshold (`--thresh`) |
| `meme.nmotifs` | `1` | Number of motifs to find (`-nmotifs`) |
| `meme.minw` | `4` | Minimum motif width (`-minw`) |
| `meme.maxw` | `12` | Maximum motif width (`-maxw`) |
| `meme.mod` | `oops` | Site distribution model: `oops`, `zoops`, or `anr` (`-mod`) |

**Per-rule resource overrides**

Each rule has a `resources.<rule_name>.mem_mb` and `resources.<rule_name>.runtime` (minutes) key. Override these to tune SLURM job requests without touching the Snakefile. Rule names match exactly — see the `resources:` section of `config.yaml` for the full list.

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

## Docker

The Dockerfile used to build the container is in `docker_build/`, along with the associated `environment.yml`. The image is built for linux/amd64.

The pre-built container (`docker://nkoranda/dapseq:latest`) is pulled automatically by Snakemake. To use a custom build, rebuild from the Dockerfile and update the `container:` directive at the top of the `Snakefile`.

## Known bugs

- GEM insists on putting its log file into the current working directory. Trying to fix this caused a whole new set of errors so for now this little detail is being ignored, but you will find the GEM logs in whatever directory you are in when you start snakemake.
- Snakemake merges external config files with the internal one rather than replacing it. Keys in the external file override matching internal keys, but anything only in the internal config still gets passed into the workflow. This can cause unintended things to bleed in if you have params in both files. Don't delete the internal config file.
