# fix: Clean up report.html column set, ordering, and naming

## Summary

`report.html`/`report.csv` currently have several reporting problems:

1. Two near-identically-named columns — `mapping_rate` (pre-MAPQ-filter
   alignment %) and `mapping_pct` (post-MAPQ-filter
   `mapped_reads / trimmed_reads` %) — read as duplicates but measure
   different things.
2. `subsampled_frags` always shows a populated value equal to `total_reads`
   even when subsampling is disabled (`bbduk.max_frags: null`), because the
   `reformat.sh` pass-through step always runs at `samplerate=1` and its
   "Output: N reads" line gets logged regardless.
3. Column order doesn't follow the order stats are actually computed /
   the logical pipeline flow — e.g. `mapping_pct` appears after
   `reads_in_peaks` even though mapping happens before peak calling.
4. `reads_5fold` and `reads_nfold` are cryptic abbreviations that don't
   describe what they measure, and don't follow the naming convention
   already established by `num_peaks_filt` / `peaks_filt.narrowPeak`.

This plan removes `mapping_rate` (keeping `mapping_pct` as the single
mapped/trimmed retention metric), gates `subsampled_frags` on
`bbduk.max_frags` so it only shows a value when subsampling actually
occurred, reorders `report.py`'s columns (and `collect_stats.py`'s matching
per-sample CSV columns) to follow the pipeline's computation order, and
renames `reads_5fold` / `reads_nfold` / `frip_top_n_fold` to explicit,
self-describing names consistent with the existing `_filt` / `_5fold` naming
used for peak files.

`update_db.py` / the SQLite database sync is explicitly **out of scope** and
deferred — see Scope Boundaries and Risks.

---

## Problem Frame

- `mapping_rate` (~99.7%) and `mapping_pct` (~90%) sit next to each other in
  the report with near-identical names but different denominators/numerators
  (alignment success vs. post-MAPQ-filter retention). `mapping_rate` is a
  rename (commit `40591ca`) of an older `alignment_rate` column that predates
  `mapping_pct` (added in `6839e52`), and is now redundant.
- `subsampled_frags` is sourced from `reformat.sh`'s "Output: N reads" line
  (`collect_stats._parse_subsample_log`). With `bbduk.max_frags: null`,
  `SUB_FRAC=1` and `reformat.sh` passes every read through unchanged, so
  `subsampled_frags == total_reads` for every sample — verified directly
  against the live `report.csv`. The column looks like subsampling occurred
  when it did not.
- The current `report.py` column order is:
  `sample, total_reads, trimmed_reads, mapped_reads, reads_in_peaks,
  reads_5fold, reads_nfold, mapping_rate, mapping_pct, frip,
  frip_top_n_fold, max_peak_score, motif_peaks, subsampled_frags,
  median_frag_size, num_peaks, num_peaks_filt`. This interleaves
  read-counting, peak-counting, and alignment metrics in an order that
  doesn't match either the pipeline's execution order or a reader's mental
  model (e.g. peak-derived columns appear *before* the alignment metrics
  they conceptually depend on).
- `reads_5fold` (reads overlapping peaks with fold-enrichment ≥ 5, a fixed
  threshold from `peaks.smk`) and `reads_nfold` (reads overlapping
  `peaks_filt.narrowPeak`, i.e. fold-enrichment ≥ `macs3.min_foldch`) are
  named after internal MACS thresholds rather than what they represent.
  `frip_top_n_fold` is the FRiP score against that same filtered peak set.
  None of these names make their meaning clear without reading the code.

## Scope Boundaries

**In scope:**
- Remove `mapping_rate` from `collect_stats.py` and `report.py`, and remove
  the now-dead `_parse_align_rate` / `_bwamem2_align_rate` helpers,
  `aligner`-branching logic, and related `align_log`/`flagstat_log` inputs in
  `stats.smk` and `align.smk`.
- Gate `subsampled_frags` on `bbduk.max_frags` in `collect_stats.py`/
  `stats.smk`: `NA` when `max_frags` is not set, the actual subsampled count
  when it is.
- Reorder `report.py`'s `make_cols()` (and `collect_stats.py`'s `COLUMNS`,
  kept in matching order for consistency) to follow the pipeline's
  computation order: read counts → mapping → fragment size → peak calling →
  reads-in-peaks/FRiP → peak quality/motif.
- Rename `reads_5fold` → `reads_in_peaks_5fold`, `reads_nfold` →
  `reads_in_peaks_filt`, and `frip_top_n_fold` → `frip_filt` throughout
  `collect_stats.py` and `report.py`, matching the existing `_filt`/`_5fold`
  naming already used for `num_peaks_filt`/`peaks_filt.narrowPeak` and
  `peaks_5fold.narrowPeak`.
- Update `tests/test_collect_stats_and_report.py` to reflect all of the
  above.

**Out of scope / non-goals:**
- **`update_db.py` / SQLite database sync** — explicitly deferred per user
  direction ("fuck the database right now, we can fix that later"). No
  changes to `update_db.py`'s `COLS`, row population, or schema in this plan.
  See "Deferred to Follow-Up Work" and Risks for the consequences of this.
- Changing `mapping_pct`'s calculation — it stays as the sole mapped/trimmed
  retention metric.
- Changing the bowtie2 alignment pipeline (`align_log` continues to be
  written by bowtie2 as before; it's just no longer read by stats).
- Any change to `bbduk.max_frags` defaults or subsampling behavior itself —
  only how the resulting count is reported.
- Regenerating/backfilling existing `report.csv`/`report.html`/database rows
  for prior pipeline runs.

### Deferred to Follow-Up Work

- **Sync `update_db.py` to the renamed/removed report columns.** Once this
  plan lands, `update_db.py` will reference column names
  (`reads_5fold`, `reads_nfold`, `frip_top_n_fold`) that no longer exist in
  `report.csv`, joining the pre-existing dead `alignment_rate` reference (a
  leftover from before the `alignment_rate` → `mapping_rate` rename, already
  always `NA`). A follow-up should, in one pass:
  - Remove the dead `alignment_rate` column from `update_db.COLS`.
  - Rename `reads_5fold` → `reads_in_peaks_5fold`, `reads_nfold` →
    `reads_in_peaks_filt`, `frip_top_n_fold` → `frip_filt` in
    `update_db.COLS` and `main()`'s row population.
  - Decide whether to migrate/rename existing SQLite columns or simply add
    new ones via `_ensure_columns` and leave old columns stale.

---

## Key Technical Decisions

### KTD1: Drop `mapping_rate` entirely rather than rename

**Decision:** Remove `mapping_rate` from all column lists and stop computing
it, rather than renaming it for clarity.

**Rationale:** `mapping_pct` (mapped_reads / trimmed_reads) already answers
"what fraction of reads ended up in the final BAM" — the metric that matters
for downstream QC. `mapping_rate` (pre-filter alignment %) is consistently
~99-100% across all current samples and adds no differentiating signal.
(User-confirmed.)

### KTD2: `subsampled_frags` reports NA unless `max_frags` is configured

**Decision:** In `collect_stats.py`, only populate `subsampled_frags` from
`_parse_subsample_log` when `bbduk.max_frags` is set (truthy); otherwise
write `NA`.

**Rationale:** When `max_frags` is `null`, `reformat.sh` runs at
`samplerate=1` and passes every read through, so "Output: N reads" is
identical to `total_reads` and conveys no additional information — worse, it
visually implies subsampling happened. Gating on `max_frags` makes the column
mean exactly one thing: "this many fragments remained after subsampling to
`max_frags`." (User-confirmed.)

### KTD3: Full cleanup of now-dead alignment-rate plumbing, including align.smk

**Decision:** Remove `_parse_align_rate`/`_bwamem2_align_rate`, the
`align_log`/`flagstat_log` stats inputs, the `aligner`/`_align_log_dir`/
`_is_bwa_mem2` helpers in `stats.smk`, *and* the
`tee >(samtools flagstat ...)` step in the bwa-mem2 `align_se`/`align_pe`
rules in `align.smk`.

**Rationale:** Once `mapping_rate` is gone, nothing reads
`{sample}.flagstat.log` or bowtie2's `align_log` for stats purposes. Leaving
the `tee` step in place would silently generate an orphaned log file on every
bwa-mem2 run with no consumer. (User-confirmed: full cleanup including
`align.smk`.)

### KTD4: Column order follows pipeline computation order

**Decision:** Reorder `report.py`'s `make_cols()` to:

```
sample
total_reads
subsampled_frags
trimmed_reads
mapped_reads
mapping_pct
median_frag_size
num_peaks
num_peaks_filt
reads_in_peaks
reads_in_peaks_5fold
reads_in_peaks_filt
frip
frip_filt
max_peak_score
motif_peaks
```

and reorder `collect_stats.py`'s `COLUMNS` to match (it currently has a
different order than `make_cols()` already — this brings the per-sample CSV
and the aggregate report into a single consistent order).

**Rationale:** This follows the order stats are actually produced by the
pipeline — read counts (subsample → trim → align/filter → mapping_pct →
fragment size), then peak calling (`num_peaks`/`num_peaks_filt`), then
reads-in-peaks and FRiP derived from those peak sets, then peak-quality and
motif metrics last. `mapping_pct` now appears immediately after
`mapped_reads` and well before any peak-derived columns, addressing the
explicit ordering request that `mapping_pct` precede `reads_in_peaks`.
(User-confirmed.)

### KTD5: Rename `reads_5fold`/`reads_nfold`/`frip_top_n_fold` for clarity

**Decision:** Rename:
- `reads_5fold` → `reads_in_peaks_5fold` (reads overlapping
  `{sample}_peaks_5fold.narrowPeak`, fold-enrichment ≥ 5)
- `reads_nfold` → `reads_in_peaks_filt` (reads overlapping
  `{sample}_peaks_filt.narrowPeak`, fold-enrichment ≥ `macs3.min_foldch`)
- `frip_top_n_fold` → `frip_filt` (FRiP computed against
  `reads_in_peaks_filt`)

**Rationale:** `_5fold` and `_filt` are already established suffixes for the
two MACS-derived peak sets (`peaks_5fold.narrowPeak`, `peaks_filt.narrowPeak`,
`num_peaks_filt`). Extending those suffixes to the corresponding read-count
and FRiP columns makes the relationship between `num_peaks_filt`,
`reads_in_peaks_filt`, and `frip_filt` (and the `_5fold` equivalents)
self-evident from the names alone, without needing to read `peaks.smk`.
(User-confirmed.)

---

## Implementation Units

### U1. Remove `mapping_rate` and dead alignment-rate helpers from `collect_stats.py`

**Goal:** `collect_stats.py` no longer computes or writes `mapping_rate`, and
the alignment-rate parsing helpers and `aligner` branch are removed.

**Requirements:** KTD1, KTD3

**Dependencies:** None

**Files:**
- `workflow/scripts/collect_stats.py` (modify)

**Approach:**
- Remove `_parse_align_rate()` and `_bwamem2_align_rate()` entirely.
- Remove `"mapping_rate"` from `COLUMNS` (this also gets reordered/renamed in
  U4 — remove it here regardless of final ordering).
- In `main()`, remove the `if sm.params.aligner == "bwa_mem2": ... else: ...`
  block that sets `row["mapping_rate"]`, and remove the now-unused
  `sm.params.aligner` read.
- Leave `_safe_pct` and `mapping_pct` computation
  (`mapped_reads / trimmed_reads`) unchanged.

**Patterns to follow:** Existing `COLUMNS` list and `main()` structure in
`workflow/scripts/collect_stats.py`.

**Test scenarios:**
- Test expectation: covered by U6's column-list assertions
  (`main()` itself is exercised via Snakemake, not unit-tested directly).

**Verification:** `grep -n "mapping_rate\|_parse_align_rate\|_bwamem2_align_rate\|sm.params.aligner" workflow/scripts/collect_stats.py` returns nothing.

---

### U2. Gate `subsampled_frags` on `bbduk.max_frags`

**Goal:** `subsampled_frags` is `NA` when subsampling is not configured, and
the real subsampled count when it is.

**Requirements:** KTD2

**Dependencies:** U3 (needs the new `max_frags` param wired through from
`stats.smk`)

**Files:**
- `workflow/scripts/collect_stats.py` (modify)

**Approach:**
- `_parse_subsample_log` keeps parsing both `total_reads` (still always
  reported — it reflects the true input read/pair count) and the raw
  "Output: N reads" value.
- Extract a small pure helper, e.g.:
  ```
  def _gate_subsampled(subsampled, max_frags):
      """Return subsampled count if max_frags is set, else NA."""
      return subsampled if max_frags else NA
  ```
- In `main()`, after computing `total_reads, subsampled =
  _parse_subsample_log(...)`:
  ```
  row["total_reads"] = total_reads
  row["subsampled_frags"] = _gate_subsampled(subsampled, sm.params.max_frags)
  ```
- This keeps `_parse_subsample_log`'s signature/behavior unchanged (still the
  single source for `total_reads`), and isolates the gating in a directly
  unit-testable helper.

**Patterns to follow:** `_safe_pct`'s short-circuit-to-`NA` style.

**Test scenarios:**
- Happy path: `_gate_subsampled("5000000", "5000000")` → `"5000000"`
  (max_frags set to a numeric value).
- Edge case: `_gate_subsampled("5000000", None)` → `"NA"` (max_frags unset,
  matches `config["bbduk"].get("max_frags")` returning `None`).
- Edge case: `_gate_subsampled("NA", "5000000")` → `"NA"` (subsample log had
  no `Output:` line, even though max_frags is set).
- Edge case: `_gate_subsampled("5000000", 0)` → `"NA"` (max_frags falsy
  zero, treated same as unset).

**Verification:** Unit tests in U6 cover all four scenarios via
`_gate_subsampled` directly.

---

### U3. Wire `max_frags` param and remove dead alignment-log inputs in `stats.smk`

**Goal:** `sample_stats_treatment` / `sample_stats_control` pass
`bbduk.max_frags` to `collect_stats.py` and no longer declare `align_log` /
`flagstat_log` inputs or the `aligner` param.

**Requirements:** KTD2, KTD3

**Dependencies:** None

**Files:**
- `workflow/rules/stats.smk` (modify)

**Approach:**
- Remove the module-level `_align_log_dir` and `_is_bwa_mem2` variables (no
  longer referenced anywhere in this file after this change).
- In both `sample_stats_treatment` and `sample_stats_control`:
  - Remove the `align_log = ...` input line.
  - Remove the conditional `**({} if not _is_bwa_mem2 else {"flagstat_log":
    ...})` input block.
  - Remove `aligner = config.get("aligner", "bowtie2")` from `params`.
  - Add `max_frags = config["bbduk"].get("max_frags")` to `params`.
- `idxstats`, `report`, `report_html` rules are unaffected.

**Patterns to follow:** Existing `params:` style of pulling values via
`config[...].get(...)` already used throughout `stats.smk` and `trim.smk`.

**Test scenarios:**
- Test expectation: none -- this is a Snakemake rule wiring change with no
  standalone unit; correctness is verified by U2's `_gate_subsampled` tests
  (which exercise the `max_frags` param contract) and by the grep below.

**Verification:** `grep -n "align_log\|flagstat_log\|_align_log_dir\|_is_bwa_mem2\|aligner" workflow/rules/stats.smk` returns nothing; `grep -n "max_frags" workflow/rules/stats.smk` shows the new param in both stats rules.

---

### U4. Remove bwa-mem2 flagstat tee step from `align.smk`

**Goal:** The bwa-mem2 `align_se` and `align_pe` rules no longer tee their
output through `samtools flagstat`, and no longer declare a `flagstat` log.

**Requirements:** KTD3

**Dependencies:** None (logically grouped with the `mapping_rate` removal
since it shares the same root cause)

**Files:**
- `workflow/rules/align.smk` (modify)

**Approach:**
- In both bwa-mem2 `align_se` and `align_pe` rules:
  - Remove the `flagstat = OUT + "/logs/bwa_mem2/{sample}.flagstat.log"` log
    entry.
  - Simplify the pipeline from:
    ```
    bwa-mem2 mem ... 2>{log.align} \
    | tee >(samtools flagstat - > {log.flagstat}) \
    | samtools view -h -F 4 -q {params.mapq} -b - \
    | samtools sort -@ {threads} -o {output.bam} -
    ```
    to:
    ```
    bwa-mem2 mem ... 2>{log.align} \
    | samtools view -h -F 4 -q {params.mapq} -b - \
    | samtools sort -@ {threads} -o {output.bam} -
    ```
  - `set -euo pipefail` at the top of each shell block is unaffected and
    remains correct without the `tee`/process-substitution.
- Bowtie2 `align_se`/`align_pe` rules are unaffected — their `align_log` was
  already only consumed by the now-removed stats input, and bowtie2 itself
  still writes its own log via `2>{log}` redirection regardless.

**Patterns to follow:** The non-bwa-mem2 (bowtie2) `align_se`/`align_pe`
rules in this same file already use the simpler
`... 2>{log} | samtools view ... | samtools sort ...` shape without `tee` —
mirror that.

**Test scenarios:**
- Test expectation: none -- shell pipeline simplification in a Snakemake
  rule; no existing test executes `align.smk` rules (these require real
  aligners/genome indices). Correctness is verified by reviewing the diff
  against the bowtie2 rules' equivalent (already-working) pipeline shape.

**Verification:** `grep -n "flagstat\|tee >" workflow/rules/align.smk` returns nothing.

---

### U5. Reorder and rename columns in `collect_stats.py` and `report.py`

**Goal:** `collect_stats.COLUMNS` and `report.make_cols()` both follow the
KTD4 order, and `reads_5fold`/`reads_nfold`/`frip_top_n_fold` are renamed per
KTD5.

**Requirements:** KTD4, KTD5

**Dependencies:** U1 (mapping_rate already removed before reordering)

**Files:**
- `workflow/scripts/collect_stats.py` (modify)
- `workflow/scripts/report.py` (modify)

**Approach:**
- `collect_stats.py`:
  - Reorder `COLUMNS` to:
    ```
    sample, total_reads, subsampled_frags, trimmed_reads, mapped_reads,
    mapping_pct, median_frag_size, num_peaks, num_peaks_filt,
    reads_in_peaks, reads_in_peaks_5fold, reads_in_peaks_filt,
    max_peak_score, motif_peaks
    ```
    `frip` and `frip_filt` are computed in `report.py`
    (`build_row`/`_safe_frip`), not in `collect_stats.py`, and are *not*
    part of `collect_stats.COLUMNS` today — `COLUMNS` lists only what
    `collect_stats.py` itself writes, so they're omitted here and added by
    `report.make_cols()` below.
  - Rename `row["reads_5fold"]` → `row["reads_in_peaks_5fold"]` and
    `row["reads_nfold"]` → `row["reads_in_peaks_filt"]` in `main()` (the
    `_bedtools_intersect_count` calls and their `sm.input.peaks_5fold`/
    `sm.input.peaks_filt` arguments are unchanged — only the output dict keys
    change).
- `report.py`:
  - Reorder `make_cols()` to the full KTD4 order (including `frip`/
    `frip_filt`, which `report.py` computes itself).
  - Rename `INT_COLS` entries `"reads_5fold"` → `"reads_in_peaks_5fold"`,
    `"reads_nfold"` → `"reads_in_peaks_filt"`.
  - Rename `PCT_COLS` entry `"frip_top_n_fold"` → `"frip_filt"`.
  - In `_safe_frip`/`build_row`, rename:
    - `data.get("reads_in_peaks", "NA")` stays the same (unchanged column).
    - `data.get("reads_nfold", "NA")` → `data.get("reads_in_peaks_filt",
      "NA")` (this now reads the renamed key written by `collect_stats.py`
      in U1's sibling change above).
    - `row["frip_top_n_fold"]` → `row["frip_filt"]`.

**Patterns to follow:** Existing `COLUMNS`/`make_cols()`/`INT_COLS`/
`PCT_COLS` structure; `_safe_frip` is a pure function and stays unchanged
internally — only its call-site keys change.

**Test scenarios:**
- Covered by U6.

**Verification:** `grep -rn "reads_5fold\|reads_nfold\|frip_top_n_fold\|mapping_rate" workflow/scripts/collect_stats.py workflow/scripts/report.py` returns nothing; `python -c "import sys; sys.path.insert(0,'workflow/scripts'); import report as rp; print(rp.make_cols())"` shows the KTD4 order with renamed columns.

---

### U6. Update `tests/test_collect_stats_and_report.py`

**Goal:** Tests reflect `mapping_rate` removal, `subsampled_frags` gating,
new column order, and renamed columns.

**Requirements:** KTD1, KTD2, KTD4, KTD5

**Dependencies:** U1, U2, U5

**Files:**
- `tests/test_collect_stats_and_report.py` (modify)

**Approach:**
- Remove `test_mapping_pct_after_mapping_rate` (asserts an ordering relative
  to a column that no longer exists).
- Add `test_make_cols_does_not_contain_mapping_rate`: assert
  `"mapping_rate" not in rp.make_cols()`.
- Add `test_columns_does_not_contain_mapping_rate`: assert
  `"mapping_rate" not in cs.COLUMNS`.
- Add `test_make_cols_order`: assert `rp.make_cols()` matches the full KTD4
  ordered list exactly (pin the order so future changes are deliberate).
- Add `test_mapping_pct_before_reads_in_peaks`: assert
  `cols.index("mapping_pct") < cols.index("reads_in_peaks")` (directly pins
  the user's explicit ordering requirement, in addition to the full-list
  pin above).
- Add `test_renamed_columns_present`: assert `"reads_in_peaks_5fold"`,
  `"reads_in_peaks_filt"`, `"frip_filt"` are in `rp.make_cols()` and
  `"reads_5fold"`, `"reads_nfold"`, `"frip_top_n_fold"` are not.
- Add a `_gate_subsampled` test section (4 cases from U2's test scenarios).
- Update `_write_stats_csv`'s hardcoded extras (`row["mapping_pct"] =
  "90.0"`, etc.) and any other references to renamed columns
  (`reads_nfold` → `reads_in_peaks_filt` if used in `_safe_frip`-related
  fixtures) so existing FRiP tests continue to exercise the renamed key.
- `test_run_csv_includes_mapping_pct` and the `_safe_pct`/`_safe_frip` unit
  tests are otherwise unaffected and stay as-is.

**Patterns to follow:** Existing `_write_stats_csv` helper and
`test_make_cols_contains_mapping_pct`-style assertions.

**Test scenarios:**
- Happy path: `rp.make_cols()` does not contain `"mapping_rate"` and does
  contain `"mapping_pct"`.
- Happy path: `cs.COLUMNS` does not contain `"mapping_rate"`.
- Happy path: `rp.make_cols()` equals the exact KTD4-ordered list.
- Happy path: `mapping_pct` index < `reads_in_peaks` index in
  `rp.make_cols()`.
- Happy path: renamed columns (`reads_in_peaks_5fold`, `reads_in_peaks_filt`,
  `frip_filt`) present; old names absent.
- Happy path (U2): `_gate_subsampled("5000000", "5000000")` → `"5000000"`.
- Edge case (U2): `_gate_subsampled("5000000", None)` → `"NA"`.
- Edge case (U2): `_gate_subsampled("NA", "5000000")` → `"NA"`.
- Edge case (U2): `_gate_subsampled("5000000", 0)` → `"NA"`.
- Integration: `_write_stats_csv`-based test confirms `report.run_csv` output
  for a sample never includes a `mapping_rate` column and includes
  `reads_in_peaks_filt`/`frip_filt` with values computed from the renamed
  `reads_in_peaks_filt` input field.

**Verification:** `pytest tests/test_collect_stats_and_report.py -v` passes.

---

## Risks & Dependencies

- **`update_db.py` will reference stale/renamed/removed columns after this
  plan lands** (`alignment_rate` — already dead before this plan;
  `reads_5fold`, `reads_nfold`, `frip_top_n_fold` — newly stale from U5).
  `stats.get(col, "NA")` already defaults missing keys to `"NA"`, so this
  fails *silently* (no exception, just `NA` written to the database for
  those columns) — explicitly accepted per user direction to defer the
  database fix, but flagged here so it's a tracked, known gap rather than a
  rediscovered surprise. See "Deferred to Follow-Up Work".
- **Existing `report.csv`/database rows for prior runs** retain their old
  column names/values (including `subsampled_frags == total_reads` from
  before this fix). Out of scope per Scope Boundaries — only new pipeline
  runs get the corrected/reordered/renamed columns.
- **`align.smk` pipeline edit (U4)** changes a shell pipeline for the
  bwa-mem2 alignment rule. The bowtie2 equivalent already uses the simpler
  no-`tee` shape, so the change is low-risk, but it's the one unit in this
  plan that touches a rule with real runtime side effects (alignment output).
  Recommend running a small smoke sample with `aligner: bwa_mem2` after
  implementation if the environment allows. Removing the `flagstat` tee
  loses a standing diagnostic artifact, but it's cheaply recoverable —
  `samtools flagstat {output.bam}` reproduces the same numbers post-hoc from
  the output BAM if ever needed for debugging.
- **U2 depends on U3** for the `max_frags` param to be available on
  `sm.params`; **U5 depends on U1** for `mapping_rate` to already be gone
  before reordering; **U6 depends on U1, U2, U5**. Suggested implementation
  order: U1 → U3 → U2 → U4 (independent, can interleave) → U5 → U6.
