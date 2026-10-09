"""Read-layout decisions derived from the per-sample config.

Imported by workflow/rules/common.smk. Kept free of Snakemake objects so the
logic is unit-testable.
"""

import gzip
import os
from decimal import Decimal, InvalidOperation


def as_list(value):
    """A sample's r1/r2 value (null, a path, or a list of lane paths) as a list."""
    if value is None:
        return []
    return list(value) if isinstance(value, (list, tuple)) else [value]


def lane_count_errors(samples_cfg, pe_samples):
    """One error line per paired-end sample whose r1 and r2 lane lists differ in length."""
    errors = []
    for sample in sorted(pe_samples):
        cfg = samples_cfg[sample]
        n_r1 = len(as_list(cfg.get("r1")))
        n_r2 = len(as_list(cfg.get("r2")))
        if n_r1 != n_r2:
            errors.append(
                f"  sample '{sample}' has {n_r1} r1 file(s) but {n_r2} r2 file(s); "
                "r1 and r2 lanes pair by position"
            )
    return errors


def macs3_format(sample, control, pe_samples, control_samples, override):
    """MACS3 -f value for one peak call.

    An override (macs3.format) wins for every call. Otherwise a control
    sample's self call follows its own layout, and a treatment call is BAMPE
    only when the treatment is paired-end and its control (if any) is too.
    """
    if override:
        return override
    if sample not in pe_samples:
        return "BAM"
    if sample in control_samples or control is None:
        return "BAMPE"
    return "BAMPE" if control in pe_samples else "BAM"


def fallback_pairs(sample_control, pe_samples, override):
    """(treatment, control) pairs where a paired-end treatment is called as BAM
    because its control is single-end. Empty when an override is set."""
    if override:
        return []
    return [
        (t, c) for t, c in sorted(sample_control.items())
        if t in pe_samples and macs3_format(t, c, pe_samples, (), None) == "BAM"
    ]


def parse_genome_size(value):
    """genome_size as a positive whole number. Accepts an int, a digit string,
    or scientific notation such as "2.7e9"; raises ValueError otherwise."""
    if value is None or isinstance(value, bool):
        raise ValueError(f"genome_size must be a positive whole number, got {value!r}")
    try:
        size = Decimal(str(value).strip())
    except InvalidOperation:
        raise ValueError(
            f"genome_size must be a number such as 2700000000 or \"2.7e9\", got {value!r}"
        ) from None
    if not size.is_finite() or size <= 0 or size != size.to_integral_value():
        raise ValueError(f"genome_size must be a positive whole number, got {value!r}")
    return int(size)


def reference_length(fasta_path):
    """Total sequence length of a reference FASTA, or None when it cannot be read.

    Sums <fasta>.fai (samtools faidx) when it exists; otherwise counts the
    FASTA's sequence letters, through gzip for a .gz path. Any read or parse
    error gives None: the caller only warns, and an exception raised in
    onstart would stop the run.
    """
    try:
        fai = fasta_path + ".fai"
        if os.path.exists(fai):
            with open(fai) as fh:
                return sum(int(line.split("\t")[1]) for line in fh if line.strip())
        opener = gzip.open if fasta_path.endswith(".gz") else open
        total = 0
        with opener(fasta_path, "rt") as fh:
            for line in fh:
                if not line.startswith(">"):
                    total += len(line.strip())
        return total
    except (OSError, ValueError, IndexError, EOFError):
        return None


def oversized_genome_size_warning(genome_size, reference_bp):
    """Warning text when genome_size exceeds the reference length, else None.

    The effective genome size can never be larger than the assembly, so such a
    value is certainly wrong (e.g. a human size on a plant genome).
    """
    if reference_bp is None or genome_size <= reference_bp:
        return None
    return (
        f"WARNING: genome_size {genome_size:,} is larger than the reference genome "
        f"({reference_bp:,} bp). genome_size is the effective genome size and cannot "
        "exceed the reference; MACS3 and bamCoverage still use the configured value. "
        'See "Choosing genome_size" in the README.'
    )
