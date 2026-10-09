"""Read-layout decisions derived from the per-sample config.

Imported by workflow/rules/common.smk. Kept free of Snakemake objects so the
logic is unit-testable.
"""

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
