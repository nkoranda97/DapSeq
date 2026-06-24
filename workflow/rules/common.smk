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

SCRIPTS = os.path.join(workflow.basedir, "scripts")

# --- Fold-change filter config validation ---
_fold_levels = config["macs3"]["foldch_levels"]
_meme_idx    = config["macs3"].get("meme_foldch_level", 2)
if len(_fold_levels) != 3:
    raise ValueError(f"macs3.foldch_levels must have exactly 3 values, got {len(_fold_levels)}")
if not (_fold_levels[0] < _fold_levels[1] < _fold_levels[2]):
    raise ValueError(f"macs3.foldch_levels must be strictly increasing, got {_fold_levels}")
if _meme_idx not in (1, 2, 3):
    raise ValueError(f"macs3.meme_foldch_level must be 1, 2, or 3, got {_meme_idx}")
FOLD_LEVELS   = _fold_levels
MEME_FOLD_IDX = _meme_idx


def is_pe(wc):
    return "true" if wc.sample in PE_SAMPLES else "false"


def get_final_filtered_peaks(wc):
    if config.get("blacklist_filter", {}).get("enabled", False):
        return OUT + f"/MACS/{wc.sample}_peaks_fold{MEME_FOLD_IDX}_bl.narrowPeak"
    return OUT + f"/MACS/{wc.sample}_peaks_fold{MEME_FOLD_IDX}.narrowPeak"


wildcard_constraints:
    sample    = "[^/.]+",
    read      = "R[12]",
    peak_type = "summits|peaks",


rule mapped:
    input:
        expand(OUT + "/bam/{sample}.bam",     sample=SAMPLES),
        expand(OUT + "/bam/{sample}.bam.bai", sample=SAMPLES),
        expand(OUT + "/bigWig/{sample}.bw",   sample=SAMPLES),


rule peaked:
    input:
        expand(OUT + "/MACS/{sample}_peaks_fold{level}_bl.narrowPeak",
               sample=TREATMENT_SAMPLES, level=MEME_FOLD_IDX)
        if config.get("blacklist_filter", {}).get("enabled", False)
        else expand(OUT + "/MACS/{sample}_peaks_fold{level}.narrowPeak",
                    sample=TREATMENT_SAMPLES, level=MEME_FOLD_IDX),


rule motifs_done:
    input:
        expand(OUT + "/fimo/{sample}/summits/fimo.tsv",      sample=TREATMENT_SAMPLES),
        expand(OUT + "/fimo/{sample}/peaks/fimo.tsv",        sample=TREATMENT_SAMPLES),
        expand(OUT + "/meme/{sample}/summits/logo1.png",     sample=TREATMENT_SAMPLES),
        expand(OUT + "/meme/{sample}/summits/logo_rc1.png",  sample=TREATMENT_SAMPLES),
        expand(OUT + "/meme/{sample}/peaks/logo1.png",       sample=TREATMENT_SAMPLES),
        expand(OUT + "/meme/{sample}/peaks/logo_rc1.png",    sample=TREATMENT_SAMPLES),
        expand(OUT + "/factorbook/{sample}.logo.png",        sample=TREATMENT_SAMPLES),


rule qc_done:
    input:
        OUT + "/stats/report.csv",
        OUT + "/stats/report.html",
        OUT + "/multiqc_report.html",


if config.get("gene_annotation"):
    rule annotate_done:
        input:
            expand(OUT + "/annotations/{sample}.peak_annotations.txt",
                   sample=TREATMENT_SAMPLES),
