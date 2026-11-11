"""Sample-name rules for the workflow's wildcard matching.

Imported by workflow/rules/common.smk. Kept free of Snakemake objects so the
logic is unit-testable.
"""

import re

# Output paths use "{sample}.{read}" and directories per sample, so a sample
# name cannot contain "." or "/". Names also go unquoted into shell commands,
# so only characters with no shell meaning are allowed.
_FORBIDDEN = re.compile(r"[./\s]")
_ALLOWED = re.compile(r"[A-Za-z0-9_+-]+")


def sample_regex(names):
    """Regex alternation matching exactly these sample names; never matches when empty."""
    names = list(names)
    return "|".join(re.escape(n) for n in names) if names else "(?!)"


def invalid_sample_names(names):
    """One error line per name the workflow cannot use: anything but letters,
    digits, "_", "+" and "-", a leading "-", or a name that is not text."""
    errors = []
    for n in names:
        if not isinstance(n, str):
            errors.append(
                f"  sample {n!r}: names must be text; quote it in the config, e.g. \"{n}\":"
            )
        elif _FORBIDDEN.search(n):
            errors.append(f"  sample {n!r}: names may not contain '.', '/' or whitespace")
        elif not _ALLOWED.fullmatch(n):
            errors.append(
                f"  sample {n!r}: names may only contain letters, digits, '_', '+' and '-'"
            )
        elif n.startswith("-"):
            errors.append(f"  sample {n!r}: names may not start with '-'")
    return errors


def control_name_collisions(names, controls):
    """One error line per sample named "<control>_control", which is the name
    of that control's QC self-call output in MACS/."""
    taken = {f"{c}_control" for c in controls}
    return [
        f"  sample {n!r}: clashes with the QC self-call output of control "
        f"'{n[:-len('_control')]}'; rename the sample"
        for n in names if n in taken
    ]
