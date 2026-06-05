---
title: "feat: Add chromosome name filtering to narrow_peak_to_fasta"
date: 2026-06-05
sequence: "002"
type: feat
status: active
---

# feat: Add chromosome name filtering to narrow_peak_to_fasta

## Summary

Add a configurable chromosome exclusion list that is applied inside `narrow_peak_to_fasta.py` before MEME runs. Peaks whose chromosome name appears in the list are dropped. Each excluded chromosome is counted and logged (zero when absent, no errors).

## Problem Frame

Certain experiments include peaks on non-target chromosomes such as viral sequences (e.g. `chrEBV`). These peaks break `narrow_peak_to_fasta.py` with a `KeyError` when the chromosome is absent from the genome FASTA, and they pollute MEME motif discovery when present. The fix must be configurable per-run with zero friction for runs that do not need filtering.

## Requirements

- R1: Config accepts a list of chromosome names to exclude (default: empty — no filtering).
- R2: Filtering happens after peak calling and before MEME sequence extraction.
- R3: Every chromosome in the filter list is counted and logged, even when its count is zero.
- R4: Missing chromosomes (zero occurrences) do not raise errors or warnings beyond the log entry.
- R5: An empty filter list is a no-op — output is identical to the current behavior.

## Key Technical Decisions

**Inline filtering in `narrow_peak_to_fasta.py` rather than a separate Snakemake rule.**
Keeps the DAG unchanged, avoids an intermediate file, and the pandas DataFrame is already loaded at that point — filtering is a one-liner. A separate rule would add a file and an extra job for what is logically a preprocessing step of the same operation.

**Top-level `chrom_filter:` config key (list, default empty).**
The filter applies post-MACS3 and is not a MACS3 parameter, so nesting it under `macs3:` would be misleading. A top-level key is consistent with how other cross-step settings (`threads`, `gene_annotation`) are placed.

**Count before filtering, iterate over the config list.**
Iterating the user-supplied list (not the unique chromosomes found in the file) guarantees a log entry for every requested chromosome regardless of presence. This is the intended zero-count-without-error behavior.

## Implementation Units

### U1. Add `chrom_filter` config key

**Goal:** Expose the chromosome exclusion list in the config with a safe default.

**Requirements:** R1, R5

**Dependencies:** none

**Files:**
- `config/config.yaml`

**Approach:** Add `chrom_filter: []` as a top-level key between the `macs3:` and `meme:` blocks. Include a short comment with an example (`# - chrEBV`).

**Test expectation:** none — config change only, no behavioral logic here.

**Verification:** `pixi run snakemake --configfile config/config.yaml --list` completes without error; the key is accessible as `config.get("chrom_filter", [])` in rules.

---

### U2. Pass `filter_chroms` param in `narrow_peak_to_fasta` rule

**Goal:** Wire the config list into the script via `snakemake.params`.

**Requirements:** R1, R2

**Dependencies:** U1

**Files:**
- `workflow/rules/motifs.smk`

**Approach:** Add `filter_chroms = config.get("chrom_filter", [])` to the `params:` block of the `narrow_peak_to_fasta` rule. Use `.get` with `[]` default so runs without `chrom_filter` in the config continue working.

**Test expectation:** none — wiring only.

**Verification:** `snakemake --dryrun` completes; rule params include `filter_chroms`.

---

### U3. Filter and log chromosome counts in `narrow_peak_to_fasta.py`

**Goal:** Apply the exclusion list, log per-chromosome counts, and leave the rest of the script unchanged.

**Requirements:** R1–R5

**Dependencies:** U2

**Files:**
- `workflow/scripts/narrow_peak_to_fasta.py`

**Approach:**
1. Add `filter_chroms` parameter to the `narrow_peak_to_fasta` function signature.
2. After `peaks` is loaded from the narrowPeak file and before the sort, if `filter_chroms` is non-empty:
   - For each chrom in `filter_chroms`, count `(peaks["chromosome"] == chrom).sum()` and print `f"Filtered {count} peaks on chromosome {chrom}"` to `sys.stderr` (which is already redirected to the log file).
   - Drop all rows where `chromosome` is in `filter_chroms`: `peaks = peaks[~peaks["chromosome"].isin(filter_chroms)]`.
3. Update the Snakemake entrypoint at the bottom to pass `snakemake.params.filter_chroms`.

**Patterns to follow:**
- Existing `sys.stderr` log redirect at line 106 — all log output goes to `sys.stderr`, which the rule already redirects.
- Existing `if peaks.empty:` guard — filtering may produce an empty DataFrame; the existing guard handles that case without changes.

**Test scenarios:**
- Empty `filter_chroms` list → `peaks` DataFrame is unchanged; no log output for chrom filtering; FASTA output identical to current behavior.
- `filter_chroms = ["chrEBV"]`, file contains 5 chrEBV peaks → log reads `"Filtered 5 peaks on chromosome chrEBV"`, those 5 rows absent from FASTA.
- `filter_chroms = ["chrEBV"]`, file contains zero chrEBV peaks → log reads `"Filtered 0 peaks on chromosome chrEBV"`, no error, FASTA output unchanged.
- `filter_chroms = ["chrEBV", "chrDecoy"]`, both present → two log lines, both sets removed.
- After filtering, DataFrame is empty (all peaks were on excluded chroms) → `if peaks.empty:` guard produces an empty FASTA file, no crash.

**Verification:** Run `pixi run snakemake` on a sample with `chrom_filter: [chrEBV]` in config; inspect `{OUT}/logs/narrow_peak_to_fasta/{sample}.*.log` for the filter count line; confirm excluded peaks absent from `{OUT}/fasta/{sample}.*.fasta`.

## Scope Boundaries

**In scope:** filtering inside `narrow_peak_to_fasta.py`, config key, rule param wiring.

**Out of scope:**
- Filtering `_peaks_5fold.narrowPeak` or any other MACS output.
- Modifying the `peaked` aggregate target in `common.smk`.
- Pattern/regex matching on chromosome names — exact string match only.

## Sources & Research

- `workflow/scripts/narrow_peak_to_fasta.py` — existing script structure and log redirect pattern.
- `workflow/rules/peaks.smk` — `filter_peaks` rule as a pattern for params and log conventions.
- `workflow/rules/motifs.smk` — `narrow_peak_to_fasta` rule to be modified in U2.
- `config/config.yaml` — existing config structure and top-level key placement convention.
