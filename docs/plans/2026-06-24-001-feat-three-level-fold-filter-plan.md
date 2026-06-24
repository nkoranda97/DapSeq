---
title: "feat: Three-level enrichment fold filter"
type: feat
status: active
date: 2026-06-24
---

# feat: Three-level enrichment fold filter

## Summary

Replace the current two-level post-MACS3 fold-change filter (configurable ~2× and hardcoded 5×) with three configurable levels (defaults 2, 5, 15). All peak-count and reads-in-peaks metrics are computed for every level. MEME continues to run on exactly one level, now configurable (default: the middle level). Blacklist filtering applies only to the MEME-selected level.

---

## Problem Frame

The current pipeline produces two filtered peak sets — a configurable "low" threshold (`_peaks_filt`) and a hardcoded 5× secondary (`_peaks_5fold`). The 5× threshold is buried in a shell command, the two levels are inconsistently named, and MEME always uses the lower threshold regardless of which set is more meaningful for motif discovery. Users also have no way to see the actual threshold values in the QC report.

**Scope boundaries**
- In scope: config schema, filtering rule, blacklist rule, stats collection, report output
- Out of scope: MEME analysis settings (maxpeaks, summit_extend, nmotifs), FIMO, HOMER, the MACS3 calling step itself

---

## Requirements

- R1: Three fold-change thresholds configurable via `macs3.foldch_levels`, defaults `[2, 5, 15]`
- R2: `macs3.meme_foldch_level` (1-based integer, default `2`) selects which level feeds MEME
- R3: Strict-ordering validation at startup: level1 < level2 < level3; valid index 1–3
- R4: `filter_peaks` produces three narrowPeak files, one per level
- R5: Blacklist filtering runs on the MEME-selected level only
- R6: Peak counts and reads-in-peaks computed for all three levels
- R7: Actual threshold values appear in the QC report (CSV columns + HTML summary)

---

## Key Technical Decisions

| Decision | Choice | Rationale |
|---|---|---|
| Output file naming | `_peaks_fold1`, `_peaks_fold2`, `_peaks_fold3` | Generic index names keep filenames stable if thresholds change; embedding the value (e.g. `_peaks_2fold`) would rename outputs every time config changes |
| MEME-level selector | 1-based index `meme_foldch_level: 2` | Intuitive for users; maps directly to "Level 1 / Level 2 / Level 3" as reported in QC |
| Blacklist output filename | `_peaks_fold{N}_bl` where N = meme_foldch_level | Consistent with the fold naming scheme; avoids a static `_filt_bl` suffix that no longer reflects the source |
| Threshold reporting | Add `foldch_level1/2/3` columns to per-sample stats CSV; add a config-summary line above the HTML table | Columns make the CSV self-documenting for archival; the HTML line provides at-a-glance context without cluttering per-sample data |
| `max_peak_score` source | MEME-level peak file | The MEME set is the primary analysis set; reporting its max score is most actionable |
| FRIP columns | `frip_fold1`, `frip_fold2`, `frip_fold3` replace `frip_filt` | Parallel naming with the underlying count columns |

---

## High-Level Technical Design

### Data flow through the reworked filter stage

```
{sample}_peaks.narrowPeak  (raw MACS3 output)
         │
         ├─ awk $7 >= level1 ──► _peaks_fold1.narrowPeak
         ├─ awk $7 >= level2 ──► _peaks_fold2.narrowPeak   ◄─── default MEME input
         └─ awk $7 >= level3 ──► _peaks_fold3.narrowPeak
                                       │  (if blacklist enabled)
                                       └─► _peaks_fold2_bl.narrowPeak  ──► MEME
```

*(Diagram shows the default case where meme_foldch_level = 2; the blacklist and MEME arrow shift to fold1 or fold3 when configured differently.)*

### Config globals resolved in `common.smk`

```
FOLD_LEVELS   = [2, 5, 15]      # from config, after validation
MEME_FOLD_IDX = 2               # 1-based; from config
```

`get_final_filtered_peaks(wc)` returns `_peaks_fold{MEME_FOLD_IDX}[_bl].narrowPeak`.

---

## Implementation Units

### U1. Update config schema

**Goal:** Replace `min_foldch` with `foldch_levels` + `meme_foldch_level` in both config templates.

**Requirements:** R1, R2

**Dependencies:** none

**Files:**
- `config.yaml`
- `config/config.yaml`

**Approach:**
Under the `macs3:` block, remove `min_foldch: 2.0` and add:
```yaml
foldch_levels: [2, 5, 15]   # three enrichment filter thresholds; must be strictly increasing
meme_foldch_level: 2         # 1-based index into foldch_levels for MEME input (1=low, 2=mid, 3=high)
```
Apply the same change to both files. Keep all other `macs3` keys (`format`, `keep_dup`, `extra`) unchanged.

**Test scenarios:**
- Test expectation: none — pure config schema change; validation is tested via U2

---

### U2. Validation and helper update in `common.smk`

**Goal:** Validate the new config at DAG-build time and expose `FOLD_LEVELS` / `MEME_FOLD_IDX` globals; update `get_final_filtered_peaks()` and `rule peaked`.

**Requirements:** R2, R3

**Dependencies:** U1

**Files:**
- `workflow/rules/common.smk`

**Approach:**
After the existing variable definitions (around line 27), add a validation block:

```python
# Validate fold-change filter config
_fold_levels = config["macs3"]["foldch_levels"]
_meme_idx    = config["macs3"].get("meme_foldch_level", 2)
if len(_fold_levels) != 3:
    raise ValueError("macs3.foldch_levels must have exactly 3 values")
if not (_fold_levels[0] < _fold_levels[1] < _fold_levels[2]):
    raise ValueError(f"macs3.foldch_levels must be strictly increasing, got {_fold_levels}")
if _meme_idx not in (1, 2, 3):
    raise ValueError(f"macs3.meme_foldch_level must be 1, 2, or 3, got {_meme_idx}")
FOLD_LEVELS   = _fold_levels
MEME_FOLD_IDX = _meme_idx
```

Update `get_final_filtered_peaks()`:
```python
def get_final_filtered_peaks(wc):
    if config.get("blacklist_filter", {}).get("enabled", False):
        return OUT + f"/MACS/{wc.sample}_peaks_fold{MEME_FOLD_IDX}_bl.narrowPeak"
    return OUT + f"/MACS/{wc.sample}_peaks_fold{MEME_FOLD_IDX}.narrowPeak"
```

Update `rule peaked` to reference `_peaks_fold{MEME_FOLD_IDX}[_bl]` instead of `_peaks_filt[_bl]`.

**Patterns to follow:** existing `CONTROL`/`TREATMENT_SAMPLES` validation style in the same file.

**Test scenarios:**
- Config with `foldch_levels: [2, 5, 15]`, `meme_foldch_level: 2` → no error, `FOLD_LEVELS=[2,5,15]`, `MEME_FOLD_IDX=2`
- Config with `foldch_levels: [5, 2, 15]` → raises `ValueError` mentioning "strictly increasing"
- Config with `foldch_levels: [2, 5]` (only 2 values) → raises `ValueError` mentioning "exactly 3 values"
- Config with `meme_foldch_level: 4` → raises `ValueError` mentioning valid range
- `get_final_filtered_peaks()` with blacklist disabled, `MEME_FOLD_IDX=2` → returns `…_peaks_fold2.narrowPeak`
- `get_final_filtered_peaks()` with blacklist enabled, `MEME_FOLD_IDX=1` → returns `…_peaks_fold1_bl.narrowPeak`

---

### U3. Rework filtering rules in `peaks.smk`

**Goal:** `filter_peaks` produces three fold-filtered output files; `blacklist_filter_peaks` targets the MEME-selected level.

**Requirements:** R4, R5

**Dependencies:** U2

**Files:**
- `workflow/rules/peaks.smk`

**Approach:**

**`filter_peaks` rule** — replace the current two-output rule:
- Outputs: `fold1`, `fold2`, `fold3` (`_peaks_fold1/2/3.narrowPeak`)
- Params: `fch1`, `fch2`, `fch3` from `FOLD_LEVELS[0/1/2]`
- Shell: three separate `awk` invocations, each reading from `{input}` directly (no piping needed); redirect each to its output; `set -euo pipefail`; errors to `{log}`

**`blacklist_filter_peaks` rule** (conditional on `blacklist_filter.enabled`):
- Input: `OUT + f"/MACS/{{sample}}_peaks_fold{MEME_FOLD_IDX}.narrowPeak"`
- Output: `OUT + f"/MACS/{{sample}}_peaks_fold{MEME_FOLD_IDX}_bl.narrowPeak"`
- Shell unchanged (bedtools intersect -v)

**Patterns to follow:** existing awk column-7 filter pattern; `MEME_FOLD_IDX` is available as a module-level global from `common.smk` (included before `peaks.smk`).

**Test scenarios:**
- With default levels `[2, 5, 15]`: `_fold1` contains peaks with score ≥ 2, `_fold2` ≥ 5, `_fold3` ≥ 15; each is a strict subset of the previous
- A peak at score exactly 5.0 appears in `_fold1` and `_fold2` but not `_fold3` (assuming level3=15)
- A peak at score 1.9 appears in none of the three files
- With `meme_foldch_level: 3`, `blacklist_filter_peaks` input is `_peaks_fold3.narrowPeak` and output is `_peaks_fold3_bl.narrowPeak`
- Empty input narrowPeak produces empty (but present) output files for all three levels

---

### U4. Update statistics collection

**Goal:** Compute peak counts and reads-in-peaks for all three fold levels; include actual threshold values as columns; remove the old `_peaks_filt` / `_peaks_5fold` references.

**Requirements:** R6, R7

**Dependencies:** U3

**Files:**
- `workflow/rules/stats.smk`
- `workflow/scripts/collect_stats.py`

**Approach:**

**`stats.smk` — `sample_stats_treatment` rule:**
- Replace `peaks_filt` and `peaks_5fold` inputs with `peaks_fold1`, `peaks_fold2`, `peaks_fold3`
- Pass new params: `foldch_levels = config["macs3"]["foldch_levels"]` and `meme_foldch_level = MEME_FOLD_IDX`
- The `peaks_bl` lambda continues pointing at `_peaks_fold{MEME_FOLD_IDX}_bl.narrowPeak`

**`collect_stats.py`:**

Replace old columns:
```
num_peaks_filt       → num_peaks_fold1, num_peaks_fold2, num_peaks_fold3
reads_in_peaks_5fold → (removed — subsumed by reads_in_peaks_fold2)
reads_in_peaks_filt  → reads_in_peaks_fold1, reads_in_peaks_fold2, reads_in_peaks_fold3
```

New columns added at the front of the "treatment" block:
```
foldch_level1, foldch_level2, foldch_level3   # actual threshold values from params
```

Updated `COLUMNS` list (replace the four old fold columns):
```python
"foldch_level1", "foldch_level2", "foldch_level3",
"num_peaks_fold1", "num_peaks_fold2", "num_peaks_fold3",
"num_peaks_bl",
"reads_in_peaks", "reads_in_peaks_fold1", "reads_in_peaks_fold2", "reads_in_peaks_fold3",
```

In the `if is_treatment:` block:
- Set `foldch_levelN` directly from `sm.params.foldch_levels[N-1]` (not from file content)
- Call `_count_lines` for all three fold files
- Call `_bedtools_intersect_count` for all three fold files
- `max_peak_score` reads from `sm.input[f"peaks_fold{sm.params.meme_foldch_level}"]`
- `num_peaks_bl` logic unchanged (reads from `sm.input.peaks_bl[0]` when blacklist enabled)

**Test scenarios:**
- All three `num_peaks_foldN` values are non-negative integers; fold1 ≥ fold2 ≥ fold3
- `foldch_level1/2/3` columns equal the configured threshold values (not NA), same for all samples
- `reads_in_peaks_fold2` matches what the old `reads_in_peaks_5fold` produced when level2=5
- Empty peak file for any level produces `"0"` for the count and `"0"` for reads, not NA
- Non-treatment samples (control) have NA for all fold columns

---

### U5. Update report

**Goal:** Reflect the three-level naming in report.csv columns and FRIP metrics; display actual threshold values in the HTML report.

**Requirements:** R7

**Dependencies:** U4

**Files:**
- `workflow/scripts/report.py`

**Approach:**

**`make_cols()`:** replace old fold columns with:
```python
"foldch_level1", "foldch_level2", "foldch_level3",
"num_peaks_fold1", "num_peaks_fold2", "num_peaks_fold3",
"num_peaks_bl",
"reads_in_peaks", "reads_in_peaks_fold1", "reads_in_peaks_fold2", "reads_in_peaks_fold3",
"frip", "frip_fold1", "frip_fold2", "frip_fold3",
```

**`INT_COLS`:** replace `reads_in_peaks_5fold`, `reads_in_peaks_filt`, `num_peaks_filt` with the three new count and reads columns.

**`PCT_COLS`:** replace `frip_filt` with `frip_fold1`, `frip_fold2`, `frip_fold3`.

**`build_row()`:** compute three FRIP values:
```python
row["frip_fold1"] = _safe_frip(data.get("reads_in_peaks_fold1"), data.get("mapped_reads"))
row["frip_fold2"] = _safe_frip(data.get("reads_in_peaks_fold2"), data.get("mapped_reads"))
row["frip_fold3"] = _safe_frip(data.get("reads_in_peaks_fold3"), data.get("mapped_reads"))
```

**HTML config summary:** In `write_html()`, add a `<p>` line above the `<table>` tag that reads the threshold values from the first row of data:
```html
<p><strong>Fold-change thresholds:</strong>
  Level 1 = {foldch_level1}×, Level 2 = {foldch_level2}×, Level 3 = {foldch_level3}×
  &nbsp;|&nbsp; MEME input: Level {meme_foldch_level}
</p>
```
The meme level can be inferred as the level whose `frip_foldN` is the "primary" one, or simply read from the first data row if passed as a param. Since `report.py` receives data via `stats_dir`, the simplest path is: read `foldch_level1/2/3` from the first sample row and determine which level has a matching threshold. Alternatively, pass `meme_foldch_level` as a param from `stats.smk`'s `rule report` and thread it through to `run_csv`/`run_html`. The param approach is cleaner — add `meme_foldch_level` to `rule report`'s params and to the function signatures.

**Test scenarios:**
- `make_cols()` returns no references to `num_peaks_filt`, `reads_in_peaks_5fold`, `reads_in_peaks_filt`, or `frip_filt`
- `frip_fold2` is computed correctly from `reads_in_peaks_fold2` / `mapped_reads` * 100
- `frip_fold2` is NA when `reads_in_peaks_fold2` is NA
- HTML output contains the threshold summary line with actual values
- HTML output contains "MEME input: Level 2" (or whichever level is configured)
- CSV output contains `foldch_level1`, `foldch_level2`, `foldch_level3` columns with numeric values

---

## Scope Boundaries

### Deferred to Follow-Up Work
- Renaming the `reads_in_peaks` (raw, unfiltered) column — not needed for this change
- Per-level bigWig tracks or coverage stats at each threshold
- Validation schema via `snakemake --lint` or a schema YAML file

### Non-Goals
- Changing MEME parameters (maxpeaks, nmotifs, summit_extend)
- Adding a fourth fold level or making the number of levels dynamic
- Per-level FIMO or HOMER annotation

---

## Open Questions

| # | Question | Status |
|---|---|---|
| OQ1 | Should `rule report` receive `meme_foldch_level` as a Snakemake param (cleanest) or should `report.py` infer it from the data? | Deferred to implementation — either works; param approach recommended |

---

## Risks & Dependencies

| Risk | Mitigation |
|---|---|
| Old output filenames (`_peaks_filt`, `_peaks_5fold`, `_peaks_filt_bl`) break any cached Snakemake DAG | Users must re-run from the filter step; document in commit message |
| Three bedtools intersect calls per sample increases stats runtime 3× over the old 2× | Intersects are fast for narrowPeak-scale data; the rule's `runtime` resource can be bumped if needed |
| `report.py` CLI mode (`cli_main`) also calls `run_html` — must receive `meme_foldch_level` | Pass it via argparse or default to 2 with a warning when not supplied |

---

## Sources & Research

- Codebase: `workflow/rules/peaks.smk`, `stats.smk`, `common.smk`, `motifs.smk`
- Codebase: `workflow/scripts/collect_stats.py`, `report.py`
- Codebase: `config.yaml`, `config/config.yaml`
- No external research required — all patterns already established in the codebase
