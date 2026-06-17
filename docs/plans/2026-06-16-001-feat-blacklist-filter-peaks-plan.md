---
title: "feat: Blacklist filtering of narrowPeak files after peak calling"
date: 2026-06-16
sequence: "001"
type: feat
status: active
---

# feat: Blacklist filtering of narrowPeak files after peak calling

## Summary

Add a configurable post-peak-calling step that removes blacklisted genomic regions from `_peaks_filt.narrowPeak` using `bedtools intersect -v`, producing `_peaks_filt_bl.narrowPeak`. When enabled, all downstream consumers (motif analysis, stats) use the blacklist-filtered file. The stats CSV and HTML report gain a `num_peaks_bl` column, giving before/after peak counts flanking the blacklist step. The feature is off by default; users opt in via `blacklist_filter: {enabled: true, bed: <path>}` in their config.

## Problem Frame

Genomic peaks called in known problematic regions (centromeres, telomeres, low-complexity repeats) inflate peak counts and pollute motif discovery. ENCODE blacklists capture these regions for common assemblies (e.g., `hg38-blacklist.v2.bed.gz`). The pipeline has no mechanism to remove them. This adds a clean, optional filtering step that slots between `filter_peaks` and downstream analysis without breaking runs that don't need it.

## Requirements

- R1: Config accepts `blacklist_filter.enabled` (bool, default false) and `blacklist_filter.bed` (path).
- R2: When enabled, `bedtools intersect -v` filters `_peaks_filt.narrowPeak` against the blacklist BED, writing `_peaks_filt_bl.narrowPeak`.
- R3: When disabled, no new file is produced; downstream rules use `_peaks_filt.narrowPeak` unchanged.
- R4: `narrow_peak_to_fasta` (and therefore MEME/FIMO) consumes the bl file when enabled, the filt file when disabled.
- R5: `num_peaks_filt` (pre-blacklist count) and `num_peaks_bl` (post-blacklist count) both appear in the stats CSV and HTML report. `num_peaks_bl` is NA when blacklist filtering is disabled.
- R6: When `enabled=true` and `bed` is null or the file is missing, the rule fails with a clear error.

## Key Technical Decisions

**Separate Snakemake rule rather than inline in `filter_peaks`.**
The blacklist step is logically distinct from the fold-change filter, benefits from its own log, and producing a named output (`_peaks_filt_bl.narrowPeak`) makes before/after counts trivially computable. An inline approach would require passing the log path around to collect the removed-peak count.

**Conditional DAG via a helper function, not a pass-through copy.**
A helper `get_final_filtered_peaks(wc)` in `common.smk` returns the bl file path when enabled, the filt file path otherwise. Downstream rules use this function as their input lambda. This means the `blacklist_filter_peaks` rule is never triggered when disabled — no wasted I/O. A pass-through `cp` approach would produce an extra file every run regardless of the setting.

**`reads_in_peaks_filt` stays anchored to `_peaks_filt.narrowPeak`.**
The user asked for peak counts before/after blacklist. Re-routing the read-depth stat would change the meaning of an existing column silently. The before/after comparison is fully captured by `num_peaks_filt` vs `num_peaks_bl`.

**`blacklist_filter.bed` is a config-supplied path (no pipeline default).**
Blacklist files are assembly-specific. The config templates document the `blacklist/` directory as the conventional location, but the pipeline does not hardcode the path — users must supply it when enabling the filter, matching the pattern for `genome_ref` and `gene_annotation`.

---

## Implementation Units

### U1. Add `blacklist_filter` config keys

**Goal:** Expose the feature flag and bed path in both config files.

**Requirements:** R1, R3, R6

**Dependencies:** none

**Files:**
- `config/config.yaml`
- `config.yaml`

**Approach:**
Add a `blacklist_filter:` block after the `chrom_filter:` line in both files:

```yaml
blacklist_filter:
  enabled: false
  bed: null     # path to blacklist BED/BED.gz; required when enabled=true
```

In `config/config.yaml` the comment can note the conventional location (`blacklist/hg38-blacklist.v2.bed.gz`). In `config.yaml` (the user template) the comment should explain that `bed` must be set when `enabled: true`.

Also add a `blacklist_filter_peaks` resource block to the `resources:` section of both files (mem_mb 32000 / 8000, runtime 10).

**Test expectation:** none — config change only.

**Verification:** `pixi run snakemake --configfile config/config.yaml --list` completes without error. `config.get("blacklist_filter", {}).get("enabled", False)` evaluates to `False`.

---

### U2. `blacklist_filter_peaks` rule and `peaked` target update

**Goal:** Implement the bedtools filtering rule and update the `peaked` aggregate target to reflect the new terminal peak file.

**Requirements:** R2, R3, R6

**Dependencies:** U1

**Files:**
- `workflow/rules/peaks.smk`
- `workflow/rules/common.smk`

**Approach — `peaks.smk`:**
Add the rule after `filter_peaks`:

```
rule blacklist_filter_peaks:
    input:
        peaks    = OUT + "/MACS/{sample}_peaks_filt.narrowPeak",
        blacklist = config["blacklist_filter"]["bed"],
    output:
        OUT + "/MACS/{sample}_peaks_filt_bl.narrowPeak",
    log:
        OUT + "/logs/blacklist_filter/{sample}.log"
    resources: ...
    shell: """
        bedtools intersect -v -a {input.peaks} -b {input.blacklist} \
          > {output} 2>{log}
    """
```

The `config["blacklist_filter"]["bed"]` reference in `input:` means Snakemake will error at parse time if `bed` is null when the rule is referenced. Guard this rule definition inside `if config.get("blacklist_filter", {}).get("enabled", False):` so the rule is never defined (and never parsed for its bed path) when disabled.

**Approach — `common.smk`:**
Add helper function:

```python
def get_final_filtered_peaks(wc):
    if config.get("blacklist_filter", {}).get("enabled", False):
        return OUT + f"/MACS/{wc.sample}_peaks_filt_bl.narrowPeak"
    return OUT + f"/MACS/{wc.sample}_peaks_filt.narrowPeak"
```

Update the `peaked` pseudo-rule:

```python
rule peaked:
    input:
        expand(OUT + "/MACS/{sample}_peaks_filt_bl.narrowPeak", sample=TREATMENT_SAMPLES)
        if config.get("blacklist_filter", {}).get("enabled", False)
        else expand(OUT + "/MACS/{sample}_peaks_filt.narrowPeak",  sample=TREATMENT_SAMPLES)
```

**Patterns to follow:**
- `chrom_filter` conditional-block pattern in `motifs.smk` for the `if config.get(...)` guard.
- Existing `filter_peaks` rule for log/resource/shell structure in `peaks.smk`.

**Test scenarios:**
- `enabled=false` → `blacklist_filter_peaks` rule is not defined; `pixi run snakemake --dryrun` completes; `peaked` target requires `_peaks_filt.narrowPeak`.
- `enabled=true, bed=<valid path>` → dry-run includes `blacklist_filter_peaks` in the DAG; `peaked` target requires `_peaks_filt_bl.narrowPeak`.
- `enabled=true, bed=null` → Snakemake parse error at rule definition, not a silent wrong-file error.

**Verification:** `pixi run snakemake --dryrun` with both config states; confirm job list includes/excludes `blacklist_filter_peaks` appropriately.

---

### U3. Wire downstream consumers

**Goal:** `narrow_peak_to_fasta` and `sample_stats_treatment` use `get_final_filtered_peaks` so they transparently consume whichever file is terminal.

**Requirements:** R4, R5

**Dependencies:** U2

**Files:**
- `workflow/rules/motifs.smk`
- `workflow/rules/stats.smk`

**Approach — `motifs.smk`:**
In `narrow_peak_to_fasta`, change the `narrowpeak` input from a string to the helper:

```python
input:
    narrowpeak = get_final_filtered_peaks,
    genome     = config["genome_ref"],
```

**Approach — `stats.smk`:**
In `sample_stats_treatment`, add an optional `peaks_bl` input via lambda:

```python
peaks_bl = lambda wc: (
    [OUT + f"/MACS/{wc.sample}_peaks_filt_bl.narrowPeak"]
    if config.get("blacklist_filter", {}).get("enabled", False)
    else []
),
```

Pass `blacklist_enabled` as a param so the script can branch without re-reading config:

```python
params:
    ...
    blacklist_enabled = config.get("blacklist_filter", {}).get("enabled", False),
```

**Patterns to follow:**
- Existing `get_r1(sample)` helper in `common.smk` for the function-as-input pattern.
- `control_bam = (OUT + f"/bam/{CONTROL}.bam" if CONTROL else [])` in `peaks.smk` for the conditional-empty-list input pattern.

**Test scenarios:**
- Dry-run with `enabled=false`: `narrow_peak_to_fasta` input resolves to `_peaks_filt.narrowPeak`; `sample_stats_treatment` has no `peaks_bl` input.
- Dry-run with `enabled=true`: `narrow_peak_to_fasta` input resolves to `_peaks_filt_bl.narrowPeak`; `sample_stats_treatment` includes `_peaks_filt_bl.narrowPeak` as a dependency.

**Verification:** `pixi run snakemake --dryrun --dag | grep narrow_peak_to_fasta` shows the correct upstream dependency for each config state.

---

### U4. Stats and report updates

**Goal:** Add `num_peaks_bl` to the per-sample stats CSV and the HTML report table.

**Requirements:** R5

**Dependencies:** U3

**Files:**
- `workflow/scripts/collect_stats.py`
- `workflow/scripts/report.py`
- `tests/test_collect_stats_and_report.py`

**Approach — `collect_stats.py`:**
1. Add `"num_peaks_bl"` to `COLUMNS` after `"num_peaks_filt"`.
2. In `main()`, after the existing peak counts block, add:
   ```python
   if sm.params.blacklist_enabled and sm.input.peaks_bl:
       row["num_peaks_bl"] = _count_lines(sm.input.peaks_bl[0])
   ```
   `sm.input.peaks_bl` is a Snakemake `InputFiles` list — index 0 gives the path.

**Approach — `report.py`:**
1. Add `"num_peaks_bl"` to `INT_COLS`.
2. Add `"num_peaks_bl"` to `make_cols()` after `"num_peaks_filt"`.

**Patterns to follow:**
- Existing `_count_lines` helper for peak file counting.
- Existing `INT_COLS` set and `make_cols()` list in `report.py` for adding new numeric columns.

**Test scenarios:**
- `report.make_cols()` contains `"num_peaks_bl"` after `"num_peaks_filt"`.
- `"num_peaks_bl"` is in `report.INT_COLS`.
- `collect_stats.COLUMNS` contains `"num_peaks_bl"`.
- When `blacklist_enabled=False`, `row["num_peaks_bl"]` remains `"NA"` (default).
- When `blacklist_enabled=True` and `peaks_bl` points to a file with 10 lines, `row["num_peaks_bl"]` is `"10"`.
- When `blacklist_enabled=True` and `peaks_bl` points to an empty file, `row["num_peaks_bl"]` is `"0"` (matches existing `_count_lines` behavior).

**Verification:** `pixi run pytest tests/test_collect_stats_and_report.py` passes.

---

## Scope Boundaries

**In scope:** `_peaks_filt.narrowPeak` → blacklist → `_peaks_filt_bl.narrowPeak`; `narrow_peak_to_fasta` input switching; `num_peaks_bl` in stats CSV and HTML table; both config files updated.

### Deferred to Follow-Up Work
- Filtering `_peaks_5fold.narrowPeak` against the blacklist (not currently consumed downstream; can be added if needed).
- Re-routing `reads_in_peaks_filt` to count against `_peaks_filt_bl.narrowPeak` when enabled.
- Automating the blacklist download/indexing as a pipeline `ref` rule.
- Adding `num_peaks_bl` to the annotation rule's input (`workflow/rules/annotation.smk`).

---

## Open Questions

- **None blocking.** The blacklist file for this project (`blacklist/hg38-blacklist.v2.bed.gz`) is already present. bedtools reads gzipped BED natively with `-b`; no decompression step needed.

---

## Sources & Research

- `workflow/rules/peaks.smk` — `filter_peaks` rule: shell/log/resource pattern to follow for the new rule.
- `workflow/rules/common.smk` — `peaked` pseudo-rule and existing helper functions.
- `workflow/rules/motifs.smk` — `narrow_peak_to_fasta` rule: input to update; `if config.get(...)` guard pattern.
- `workflow/rules/stats.smk` — `sample_stats_treatment`: optional-input lambda pattern.
- `workflow/scripts/collect_stats.py` — `_count_lines`, `COLUMNS`, and `main()` structure.
- `workflow/scripts/report.py` — `INT_COLS`, `make_cols()`: column registration.
- `blacklist/hg38-blacklist.v2.bed.gz` — blacklist file already in repo; gzipped, readable directly by bedtools.
