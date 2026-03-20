import os

GEM_JAR       = "/opt/gem/gem.jar"
GEM_READ_DIST = "/opt/gem/Read_Distribution_default.txt"

SAMPLES           = [s for s in config["samples"] if config["samples"][s]["r1"] is not None]
CONTROL           = config.get("input_control") or None
TREATMENT_SAMPLES = [s for s in SAMPLES if s != CONTROL]
OUT               = config["output_dir"]
GENOME_SPLIT_DIR  = OUT + "/genome_split"

PEAK_CALLER = config.get("peak_caller", "both")
USE_MACS    = PEAK_CALLER in ("both", "macs3")
USE_GEM     = PEAK_CALLER in ("both", "gem")
USE_BOTH    = PEAK_CALLER == "both"

ALIGNER = config.get("aligner", "bowtie2")

SE_SAMPLES = {s for s in SAMPLES if config["samples"][s].get("r2") is None}
PE_SAMPLES = {s for s in SAMPLES if config["samples"][s].get("r2") is not None}

SCRIPTS = os.path.join(workflow.basedir, "scripts")


def get_r2(wc):
    if wc.sample in SE_SAMPLES:
        return []
    return [OUT + f"/trimmed/{wc.sample}.R2.fastq.gz"]


def is_pe(wc):
    return "true" if wc.sample in PE_SAMPLES else "false"


wildcard_constraints:
    sample = "[^/.]+",
    read   = "R[12]",


rule mapped:
    input:
        expand(OUT + "/bam/{sample}.bam",     sample=SAMPLES),
        expand(OUT + "/bam/{sample}.bam.bai", sample=SAMPLES),
        expand(OUT + "/stats/{sample}.raw_alignment_stats.txt", sample=SAMPLES),


rule peaked:
    input:
        expand(OUT + "/combined_bed/{sample}.combined.bed", sample=TREATMENT_SAMPLES),


rule motifs_done:
    input:
        expand(OUT + "/meme/{sample}/meme.txt",  sample=TREATMENT_SAMPLES),
        expand(OUT + "/fimo/{sample}/fimo.tsv",  sample=TREATMENT_SAMPLES),


rule qc_done:
    input:
        OUT + "/stats/qc_summary.tsv",
