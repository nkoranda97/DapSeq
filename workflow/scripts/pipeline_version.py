"""
Version label for the pipeline checkout, stored with every results-DB row.

Computed on the host when Snakemake parses the workflow (the container has no
git) and passed to update_db as a rule param.
"""

import subprocess

UNKNOWN = "unknown"


def pipeline_version(repo_dir, timeout=10):
    """`git describe --always --dirty` for the checkout at repo_dir, else "unknown".

    safe.directory=* lets lab members run from a shared checkout they do not
    own; git otherwise refuses it as "dubious ownership".
    """
    try:
        result = subprocess.run(
            ["git", "-c", "safe.directory=*", "-C", str(repo_dir),
             "describe", "--always", "--dirty"],
            capture_output=True, text=True, timeout=timeout,
        )
    except (OSError, subprocess.SubprocessError):
        return UNKNOWN
    version = result.stdout.strip()
    return version if result.returncode == 0 and version else UNKNOWN
