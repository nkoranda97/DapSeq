# DapSeq Pipeline

## Quick Start on midway3
 ```bash
 # install environment
 git clone https://github.com/nkoranda97/DapSeq && cd DapSeq
 module load uv
 uv sync
 
 # run pipeline
 module load singularity
 uv run snakemake --profile profile --configfile /path/to/configfile
 ```