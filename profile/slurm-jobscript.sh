#!/bin/bash
# Snakemake SLURM jobscript — wraps every submitted rule
# Loads Apptainer/Singularity module before the job runs so the
# container is available on compute nodes.

module load singularity

{exec_job}
