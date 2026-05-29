import os


def get_r1(sample):
    r1 = config["samples"][sample]["r1"]
    if r1 is None:
        return []
    return r1 if isinstance(r1, list) else [r1]


SAMPLES     = [s for s in config["samples"] if get_r1(s)]

_ctrl_cfg = config.get("input_control") or None
if _ctrl_cfg is None:
    CONTROL_SAMPLES = []
elif isinstance(_ctrl_cfg, list):
    CONTROL_SAMPLES = [s for s in _ctrl_cfg if s in SAMPLES]
else:
    CONTROL_SAMPLES = [_ctrl_cfg] if _ctrl_cfg in SAMPLES else []

CONTROL = (
    "merged_control" if len(CONTROL_SAMPLES) > 1
    else CONTROL_SAMPLES[0] if CONTROL_SAMPLES
    else None
)
TREATMENT_SAMPLES = [s for s in SAMPLES if s not in CONTROL_SAMPLES]
OUT               = config["output_dir"]

SE_SAMPLES = {s for s in SAMPLES if config["samples"][s].get("r2") is None}
PE_SAMPLES = {s for s in SAMPLES if config["samples"][s].get("r2") is not None}
BAM_TYPES  = ["full", "chr"]

SCRIPTS = os.path.join(workflow.basedir, "scripts")


def is_pe(wc):
    return "true" if wc.sample in PE_SAMPLES else "false"


wildcard_constraints:
    sample    = "[^/.]+",
    read      = "R[12]",
    bam_type  = "full|chr",
    peak_type = "summits|peaks",


rule mapped:
    input:
        expand(OUT + "/bam/{sample}.{bam_type}.bam",     sample=SAMPLES, bam_type=BAM_TYPES),
        expand(OUT + "/bam/{sample}.{bam_type}.bam.bai", sample=SAMPLES, bam_type=BAM_TYPES),
        expand(OUT + "/bigWig/{sample}.{bam_type}.bw",   sample=SAMPLES, bam_type=BAM_TYPES),


rule peaked:
    input:
        expand(OUT + "/MACS/{sample}.{bam_type}_peaks_filt.narrowPeak",
               sample=TREATMENT_SAMPLES, bam_type=BAM_TYPES),


rule motifs_done:
    input:
        expand(OUT + "/fimo/{sample}.{bam_type}/summits/fimo.tsv",
               sample=TREATMENT_SAMPLES, bam_type=BAM_TYPES),
        expand(OUT + "/fimo/{sample}.{bam_type}/peaks/fimo.tsv",
               sample=TREATMENT_SAMPLES, bam_type=BAM_TYPES),


rule qc_done:
    input:
        OUT + "/stats/report.csv",
        OUT + "/multiqc_report.html",


if config.get("gene_annotation"):
    rule annotate_done:
        input:
            expand(OUT + "/annotations/{sample}.{bam_type}.peak_annotations.txt",
                   sample=TREATMENT_SAMPLES, bam_type=BAM_TYPES),
