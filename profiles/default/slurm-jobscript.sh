#!/bin/bash
# Snakemake SLURM jobscript — wraps every submitted rule
# Loads Apptainer module before the job runs so the container
# is available on compute nodes.

module load apptainer

{exec_job}
