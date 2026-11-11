"""Startup checks on the shape of the merged config.

Imported by workflow/rules/common.smk. Kept free of Snakemake objects so the
logic is unit-testable.
"""

import copy
import difflib

from layout_utils import parse_positive_whole_number

# Keys a sample may set under samples: <name>:
SAMPLE_KEYS = ("r1", "r2", "control", "tf", "experiment_date", "gdna_batch")

# macs3.format values that fit the pipeline's BAM inputs.
MACS3_FORMATS = ("BAM", "BAMPE")


def fill_empty_sections(config, defaults, _prefix=""):
    """Give every section left empty in the user's config its defaults.

    A section header whose keys are all commented out (``factorbook:``, or
    ``meme:`` under ``resources:``) is null after YAML parsing, and Snakemake's
    merge then replaces the whole default section with None. Treat that as
    "no overrides", at any depth. Modifies *config* in place and returns the
    restored sections as dotted names.
    """
    restored = []
    for key, default in defaults.items():
        if not isinstance(default, dict) or key not in config:
            continue
        if config[key] is None:
            config[key] = copy.deepcopy(default)
            restored.append(_prefix + key)
        elif isinstance(config[key], dict):
            restored += fill_empty_sections(config[key], default, f"{_prefix}{key}.")
    return restored


def max_frags(value):
    """bbduk.max_frags as a positive whole number, or None when subsampling
    is off (null, empty or 0).

    YAML reads 10,000,000 or 10M as a string, which the trim rules' awk would
    turn into 10, so anything that is not a whole number raises ValueError.
    """
    if value is None or value == "" or value == 0:
        return None
    return parse_positive_whole_number(value, "bbduk.max_frags", '10000000 or "1e7"')


def sample_config_errors(samples_cfg):
    """One error line per sample that is not a mapping or sets an unknown key."""
    errors = []
    for name, cfg in samples_cfg.items():
        if not isinstance(cfg, dict):
            errors.append(
                f"  sample '{name}' has no settings; give it at least r1, e.g. "
                f"{name}: {{r1: /path/to/reads.fastq.gz}}"
            )
            continue
        for key in cfg:
            if key in SAMPLE_KEYS:
                continue
            close = difflib.get_close_matches(str(key), SAMPLE_KEYS, n=1)
            hint = f" (did you mean '{close[0]}'?)" if close else ""
            errors.append(
                f"  sample '{name}' sets unknown key '{key}'{hint}; "
                f"allowed keys are {', '.join(SAMPLE_KEYS)}"
            )
    return errors


def macs3_format_override(value, se_samples):
    """macs3.format as "BAM" or "BAMPE", or None when unset.

    Raises ValueError for any other value, and for BAMPE when some samples
    are single-end: MACS3 finds no read pairs in a single-end BAM.
    """
    if value is None or value == "":
        return None
    fmt = str(value).upper()
    if fmt not in MACS3_FORMATS:
        raise ValueError(
            f"macs3.format must be null (choose per peak call), BAM or BAMPE, got {value!r}"
        )
    if fmt == "BAMPE" and se_samples:
        raise ValueError(
            "macs3.format is BAMPE, but these samples are single-end (no r2): "
            f"{', '.join(sorted(se_samples))}. Set it to null to choose per peak call, "
            "or BAM."
        )
    return fmt


def compressed_reference_error(path):
    """Error line when the genome FASTA is gzip- or bgzip-compressed, else None.

    samtools faidx refuses plain gzip, and the MEME FASTA step and HOMER read
    the file as plain text, so a compressed reference fails partway through.
    """
    try:
        with open(path, "rb") as fh:
            magic = fh.read(2)
    except OSError:
        return None  # existence is reported by the path checks
    if magic == b"\x1f\x8b":
        return (
            f"  genome_ref: '{path}' is gzip-compressed; the pipeline needs an "
            "uncompressed FASTA. Decompress it (gunzip -k) and point genome_ref at the .fa"
        )
    return None
