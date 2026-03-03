"""
combine_peaks.py — Step 7

Combines MACS3 summits and GEM events into a single BED file
and writes per-tool BED files for downstream comparison.

Snakemake inputs:
    macs_summits  — MACS/{sample}_summits.bed
    gem_events    — GEM/{sample}/{sample}.GEM_events.txt

Snakemake outputs:
    combined_bed  — combined_bed/{sample}.combined.bed
    macs_bed      — compare_bed/{sample}.MACS.bed
    gem_bed       — compare_bed/{sample}.GEM.bed

Snakemake params:
    window_size   — bp half-window centered on summit (default 80)
    min_score     — minimum fold-enrichment score (default 1)
"""


import logging
import sys
from dataclasses import dataclass

import pandas as pd

logging.basicConfig(
    stream=sys.stderr,
    level=logging.INFO,
    format="%(levelname)s: %(message)s",
)
log = logging.getLogger(__name__)


# ── Data model ────────────────────────────────────────────────────────────────

@dataclass
class Peak:
    chrom: str
    pos:   int
    name:  str
    score: float

    def to_bed_row(self, half_win: int) -> str:
        start = max(0, self.pos - half_win)
        end   = self.pos + half_win
        return f"{self.chrom}\t{start}\t{end}\t{self.name}\t{self.score}\t.\n"


# ── Parsers ───────────────────────────────────────────────────────────────────

def parse_macs_summits(path: str, min_score: float) -> list[Peak]:
    """
    Read a MACS3 summits BED file.
    Columns: chr, start, end, name, score (fold-enrichment)
    The 'start' column is the summit position.
    """
    df = pd.read_table(path, header=None, names=["chrom", "start", "end", "name", "score"])
    df = df[df["score"] > min_score].reset_index(drop=True)
    log.info("MACS: %d peaks after score filter (min_score=%.2f)", len(df), min_score)

    return [
        Peak(chrom=row.chrom, pos=row.start, name=f"MACS_{i}", score=row.score)
        for i, row in df.iterrows()
    ]


def parse_gem_events(path: str, min_score: float) -> list[Peak]:
    """
    Read a GEM events file.
    Position column is formatted as 'chr:pos'.
    When no control is provided, Fold is NaN — in that case all events are kept.
    """
    df = pd.read_table(path)

    no_control = df["Fold"].isna().all()
    if no_control:
        log.info("GEM: no control detected (Fold is NaN), keeping all %d events", len(df))
        df_filt = df
    else:
        df_filt = df[df["Fold"] > min_score].reset_index(drop=True)
        log.info("GEM: %d peaks after score filter (min_score=%.2f)", len(df_filt), min_score)

    peaks = []
    for i, row in df_filt.iterrows():
        chrom, pos = row["Position"].split(":")
        peaks.append(Peak(chrom=chrom, pos=int(pos), name=f"GEM_{i}", score=row["Fold"] if not no_control else row["IP"]))

    return peaks


# ── Writer ────────────────────────────────────────────────────────────────────

def write_bed(peaks: list[Peak], path: str, half_win: int) -> None:
    with open(path, "w") as fh:
        for peak in peaks:
            fh.write(peak.to_bed_row(half_win))
    log.info("Wrote %d peaks to %s", len(peaks), path)


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    window_size = snakemake.params.window_size  # noqa: F821
    min_score   = snakemake.params.min_score    # noqa: F821
    half_win    = window_size // 2

    macs_peaks = parse_macs_summits(snakemake.input.macs_summits, min_score)  # noqa: F821
    gem_peaks  = parse_gem_events(snakemake.input.gem_events, min_score)       # noqa: F821

    write_bed(macs_peaks,              snakemake.output.macs_bed,     half_win)  # noqa: F821
    write_bed(gem_peaks,               snakemake.output.gem_bed,      half_win)  # noqa: F821
    write_bed(macs_peaks + gem_peaks,  snakemake.output.combined_bed, half_win)  # noqa: F821


if "snakemake" in dir():  # noqa: F821
    main()