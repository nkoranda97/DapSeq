import os

# Known top-level config keys — union of config.yaml and config/config.yaml,
# plus keys that have no YAML default but are accessed via config.get().
_KNOWN_CONFIG_KEYS = frozenset({
    "aligner", "author", "bamcompare", "bamcoverage", "bbduk",
    "blacklist_filter", "bowtie2", "bwa_mem2", "chrom_filter", "container",
    "control", "db_path", "factorbook", "fastqc", "fimo", "gene_annotation",
    "genome_ref", "genome_size", "homer", "macs3", "meme", "multiqc",
    "output_dir", "resources", "rmsk_filter", "samples", "samtools",
    "threads",
})

_unexpected_keys = set(config.keys()) - _KNOWN_CONFIG_KEYS
if _unexpected_keys:
    _hint = ""
    if "input_control" in _unexpected_keys:
        _hint += " The control field is now named 'control', not 'input_control'."
    if "slurm_partition" in _unexpected_keys or "slurm_account" in _unexpected_keys:
        _hint += (
            " 'slurm_partition' and 'slurm_account' have moved to"
            " profiles/slurm/config.yaml under default-resources:"
            " — remove them from config.yaml."
        )
    raise ValueError(f"Unrecognized config key(s): {sorted(_unexpected_keys)}.{_hint}")


def get_r1(sample):
    r1 = config["samples"][sample]["r1"]
    if r1 is None:
        return []
    return r1 if isinstance(r1, list) else [r1]


SAMPLES     = [s for s in config["samples"] if get_r1(s)]

# Fail early with a clear message if required user config fields are missing.
# This happens when --configfile was not passed and only the pipeline defaults
# (config/config.yaml) are loaded.
_missing_fields = []
if not config.get("output_dir"):
    _missing_fields.append("output_dir")
if not config.get("genome_ref"):
    _missing_fields.append("genome_ref")
if not SAMPLES:
    _missing_fields.append("samples (at least one sample must have r1 set)")
if _missing_fields:
    raise ValueError(
        f"Required config field(s) not set: {', '.join(_missing_fields)}.\n"
        "Pass your experiment config with: --configfile /path/to/your/config.yaml\n"
        "See the README for setup instructions."
    )

_ctrl_cfg = config.get("control") or None
if _ctrl_cfg is None:
    CONTROL_SAMPLES = []
elif isinstance(_ctrl_cfg, list):
    _missing_ctrl = [s for s in _ctrl_cfg if s not in SAMPLES]
    if _missing_ctrl:
        raise ValueError(
            f"control references unknown sample(s) {_missing_ctrl} -- "
            "must match a key under 'samples:' that has a valid r1."
        )
    CONTROL_SAMPLES = list(_ctrl_cfg)
else:
    if _ctrl_cfg not in SAMPLES:
        raise ValueError(
            f"control '{_ctrl_cfg}' does not match any sample under 'samples:' "
            "(or it has no r1 set)."
        )
    CONTROL_SAMPLES = [_ctrl_cfg]

CONTROL = (
    "merged_control" if len(CONTROL_SAMPLES) > 1
    else CONTROL_SAMPLES[0] if CONTROL_SAMPLES
    else None
)
TREATMENT_SAMPLES = [s for s in SAMPLES if s not in CONTROL_SAMPLES]
OUT               = config["output_dir"]

# Regex constraints for wildcard_constraints blocks; "(?!)" never matches (no samples).
CONTROL_SAMPLE_CONSTRAINT   = "|".join(CONTROL_SAMPLES)   if CONTROL_SAMPLES   else "(?!)"
TREATMENT_SAMPLE_CONSTRAINT = "|".join(TREATMENT_SAMPLES) if TREATMENT_SAMPLES else "(?!)"

SE_SAMPLES = {s for s in SAMPLES if config["samples"][s].get("r2") is None}
PE_SAMPLES = {s for s in SAMPLES if config["samples"][s].get("r2") is not None}

SCRIPTS = os.path.join(workflow.basedir, "scripts")

ALIGNER = config.get("aligner", "bowtie2")

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

# Compute combined filter suffix for the MEME-level peaks file
_BL_ENABLED   = config.get("blacklist_filter", {}).get("enabled", False)
_RMSK_ENABLED = config.get("rmsk_filter",     {}).get("enabled", False)
PEAKS_FILTER_SUFFIX = ("_bl" if _BL_ENABLED else "") + ("_rmsk" if _RMSK_ENABLED else "")


def is_pe(wc):
    return "true" if wc.sample in PE_SAMPLES else "false"


def get_final_filtered_peaks(wc):
    return OUT + f"/MACS/{wc.sample}_peaks_fold{MEME_FOLD_IDX}{PEAKS_FILTER_SUFFIX}.narrowPeak"


def get_pre_rmsk_peaks(wc):
    """Peaks file fed into rmsk filtering (post-blacklist if enabled, else post-fold)."""
    pre = "_bl" if _BL_ENABLED else ""
    return OUT + f"/MACS/{wc.sample}_peaks_fold{MEME_FOLD_IDX}{pre}.narrowPeak"


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
        expand(
            OUT + f"/MACS/{{sample}}_peaks_fold{MEME_FOLD_IDX}{PEAKS_FILTER_SUFFIX}.narrowPeak",
            sample=TREATMENT_SAMPLES,
        ),


rule control_peaked:
    input:
        expand(OUT + "/MACS/{sample}_control_peaks.narrowPeak", sample=CONTROL_SAMPLES),


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
        expand(OUT + "/stats/{sample}.idxstats.txt", sample=SAMPLES),


if config.get("gene_annotation"):
    rule annotate_done:
        input:
            expand(OUT + "/annotations/{sample}.peak_annotations.txt",
                   sample=TREATMENT_SAMPLES),
