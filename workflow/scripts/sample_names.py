"""Sample-name rules for the workflow's wildcard matching.

Imported by workflow/rules/common.smk. Kept free of Snakemake objects so the
logic is unit-testable.
"""

import re

# Output paths use "{sample}.{read}" and directories per sample, so a sample
# name cannot contain "." or "/"; whitespace breaks the shell commands.
_FORBIDDEN = re.compile(r"[./\s]")


def sample_regex(names):
    """Regex alternation matching exactly these sample names; never matches when empty."""
    names = list(names)
    return "|".join(re.escape(n) for n in names) if names else "(?!)"


def invalid_sample_names(names):
    """One error line per name containing ".", "/" or whitespace."""
    return [
        f"  sample {n!r}: names may not contain '.', '/' or whitespace"
        for n in names if _FORBIDDEN.search(n)
    ]


def control_name_collisions(names, controls):
    """One error line per sample named "<control>_control", which is the name
    of that control's QC self-call output in MACS/."""
    taken = {f"{c}_control" for c in controls}
    return [
        f"  sample {n!r}: clashes with the QC self-call output of control "
        f"'{n[:-len('_control')]}'; rename the sample"
        for n in names if n in taken
    ]
