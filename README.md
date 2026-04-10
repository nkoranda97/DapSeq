# DapSeq Pipeline

## Quick Start on midway3
 ```bash
 # install environment
 git clone https://github.com/nkoranda97/DapSeq && cd DapSeq
 module load uv
 uv sync

 # build container (once)
 module load apptainer
 apptainer build apptainer_build/dapseq.sif apptainer_build/dapseq.def

 # run pipeline (profiles/default/ loads automatically)
 uv run snakemake --configfile /path/to/configfile
 ```
