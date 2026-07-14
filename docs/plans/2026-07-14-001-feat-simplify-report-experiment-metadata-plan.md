---
title: "feat: Simplify report columns and add experiment metadata (date, gDNA batch)"
status: completed
created: 2026-07-14
type: feat
---

## Summary

Two related changes to the QC report and results database:

1. **Add experiment metadata** — two new experiment-wide config fields, `experiment_date` and `gdna_batch`, set once near `author` at the top of the config. They flow into the HTML report as a header line and into the SQLite database (once per run).
2. **Simplify the report/DB columns** — collapse the four per-fold-level column groups (`foldch_level1-3`, `num_peaks_fold1-3`, `reads_in_peaks_fold1-3`, `frip_fold1-3`) down to a single filtered set: base metrics (`num_peaks`, `reads_in_peaks`, `frip`) plus one filtered variant (`num_peaks_filt`, `reads_in_peaks_filt`, `frip_filt`), keeping `num_peaks_bl`/`num_peaks_rmsk` when those filters are enabled.

The three-level fold-filtering **internals stay in place** — `filter.smk` still produces `_fold1/2/3.narrowPeak`, and MEME/FIMO still select their input via `meme_foldch_level`. Only the *reporting* surface (report CSV/HTML + DB columns) is simplified. The single filtered set shown (`_filt`) is the fold level MEME actually uses (`meme_foldch_level`), so the report reflects the peaks carried into motif analysis.

This partially walks back `docs/plans/2026-06-24-001-feat-three-level-fold-filter-plan.md` at the reporting layer only, and realigns the code with the single-filtered-set scheme the README and `tests/test_collect_stats_and_report.py` already describe (4 tests currently fail against the diverged 3-level code).

---

## Problem Frame

The report grew a wide 3-level fold-change breakdown (12 fold-related columns across CSV, HTML, and DB) that is hard to read and rarely used — the biologically meaningful filtered set is the one fed to MEME. At the same time, runs carry no record of *when* the experiment was performed or *which gDNA batch* it came from, so results in the shared database cannot be traced back to bench provenance.

**Success criteria:**
1. A run's `experiment_date` and `gdna_batch` (when set in config) appear in the HTML report header and are stored once per run in the database.
2. The report CSV/HTML and DB expose exactly one filtered peak set (`_filt`) plus the base metrics, with all `*_fold1/2/3` and `foldch_level*` columns removed from the reporting surface.
3. `_filt` columns reflect the `meme_foldch_level`-selected fold peaks (the set carried into motif analysis), and the HTML header states which fold-change threshold that is.
4. The fold-filtering pipeline internals (`filter.smk`, MEME/FIMO input selection, blacklist/rmsk filtering) are unchanged and still functional.
5. `tests/test_collect_stats_and_report.py` passes, updated to the simplified column scheme; new fields and header are covered by tests.
6. Adding the new config keys does not trip the unrecognized-key validation guard in `common.smk`.

---

## Requirements

- **R1** — Add `experiment_date` and `gdna_batch` top-level config fields to both the user template (`config.yaml`) and pipeline defaults (`config/config.yaml`, defaulting to `null`), and register both in `_KNOWN_CONFIG_KEYS` in `workflow/rules/common.smk`.
- **R2** — Both fields are optional and experiment-wide (a single value per run). When unset (`null`), the report omits the corresponding header entry and the DB stores an empty string, matching how `author` behaves.
- **R3** — Define the canonical simplified column set (see KTD 1) and apply it consistently across `collect_stats.py` (`COLUMNS`), `report.py` (`make_cols`, `INT_COLS`, `PCT_COLS`), and `update_db.py` (`COLS`).
- **R4** — Compute `num_peaks_filt`, `reads_in_peaks_filt` (and derive `frip_filt`) from the `meme_foldch_level`-selected fold peak file; remove computation of the `_fold1/2/3` and `foldch_level*` values from the reporting outputs.
- **R5** — Render `experiment_date` and `gdna_batch` in the HTML report header, alongside a one-line summary of which fold-change threshold `_filt` represents.
- **R6** — Store `experiment_date` and `gdna_batch` once per run in the SQLite database (in `run_metadata`, next to `author`/`run_date`), and remove the `*_fold1/2/3` columns from new `pipeline_runs` inserts while adding `num_peaks_filt`/`reads_in_peaks_filt`.
- **R7** — Preserve DB backward compatibility: existing databases must continue to open and accept new rows without manual migration (the `_ensure_columns` helper adds new columns; obsolete columns are left in place, unpopulated).
- **R8** — The standalone HTML regeneration CLI (`render_report_html.py`) must still work when run against a completed output directory; it should source `experiment_date`/`gdna_batch` from the `config_used.yaml` snapshot in that directory when present, and omit them otherwise.
- **R9** — Update `README.md` to document the two new config fields and the simplified report columns.

---

## Key Technical Decisions

**1. Canonical simplified column order.**
One authoritative column list, used everywhere:

```
sample, total_reads, subsampled_frags, trimmed_reads, mapped_reads, mapping_pct,
median_frag_size,
num_peaks, num_peaks_filt, num_peaks_bl, num_peaks_rmsk,
reads_in_peaks, reads_in_peaks_filt,
frip, frip_filt,
max_peak_score, motif_peaks
```

Removed vs. current: `foldch_level1/2/3`, `num_peaks_fold1/2/3`, `reads_in_peaks_fold1/2/3`, `frip_fold1/2/3`. `frip`/`frip_filt` remain derived in `report.py` (not stored by `collect_stats.py`), consistent with today's design. The stale test `test_make_cols_order` lists `reads_in_peaks_5fold` and omits `num_peaks_rmsk`; the tests are **not** authoritative here and will be rewritten to this list (KTD 6).

**2. `_filt` = the `meme_foldch_level`-selected fold peaks.**
`collect_stats.py` already selects `getattr(sm.input, f"peaks_fold{meme_foldch_level}")` for `max_peak_score`. Reuse that same file for `num_peaks_filt` and `reads_in_peaks_filt`. This keeps the report's filtered set identical to the peaks MEME/FIMO consume — the meaningful "filtered" set — without introducing a new config knob.

**3. Keep the 3-level internals untouched (per confirmed scope).**
`filter.smk` still emits `_fold1/2/3.narrowPeak`; `common.smk` fold validation, `foldch_levels`, and `meme_foldch_level` stay. `sample_stats_treatment` keeps `peaks_fold1/2/3` as declared inputs even though `collect_stats.py` now reads only the selected one — leaving the peak-filtering wiring alone is the whole point of the "report/DB only" scope choice, and the extra declared inputs are harmless (the files are produced regardless).

**4. Experiment metadata is experiment-wide, sourced from config params — not per-sample CSV columns.**
Like `author`, `experiment_date`/`gdna_batch` are single values per run. They are **not** added as per-sample columns in the stats/report CSV. They reach the HTML via Snakemake `params` on the `report_html` rule and reach the DB via `params` on the `update_db` rule (both read from `config`). This mirrors the existing `author` path in `update_db.py`.

**5. Store experiment metadata in `run_metadata`, next to `author`/`run_date`.**
`run_metadata` is the provenance table (author, run_date, paths); `experiment_date` and `gdna_batch` are provenance. Add them to `META_COLS`. `pipeline_runs` (`COLS`) is edited only to swap the fold columns for `num_peaks_filt`/`reads_in_peaks_filt`.

**6. DB schema migration is additive and safe.**
`_ensure_columns` only ever *adds* columns, and `_INSERT` names only current `COLS`/`META_COLS`. Dropping `*_fold1/2/3` from `COLS` leaves those columns in existing tables (NULL for new rows) — no destructive migration, no manual step. New columns (`num_peaks_filt`, `reads_in_peaks_filt`, `experiment_date`, `gdna_batch`) are auto-added on first run against an old DB. This is the established pattern (see the file-path-validation plan's deferred note about `_ensure_columns`).

**7. HTML header replaces the 3-level fold summary.**
`report.py:_fold_summary_html` currently prints three thresholds plus the MEME level. Replace it with: (a) an experiment metadata line (`Experiment date`, `gDNA batch`, shown only when provided), and (b) a one-line filter summary stating the fold-change threshold `_filt` represents (the selected `meme_foldch_level` value). The filter summary reads the selected level's value from the config params rather than the removed `foldch_level*` CSV columns.

**8. CLI HTML regeneration reads metadata from `config_used.yaml`.**
`render_report_html.py cli_main` has no `config` object. The pipeline writes `{output_dir}/config_used.yaml` on every run (`onstart` hook). The CLI path parses `experiment_date`/`gdna_batch` from that snapshot when present, and omits the header entries otherwise. This keeps standalone regeneration faithful without new CLI flags.

---

## Scope Boundaries

**In scope:**
- New `experiment_date` / `gdna_batch` config fields (both config files + validation allow-list).
- Report CSV/HTML column simplification to base + single `_filt` set.
- DB column changes (swap fold columns for `_filt`; add experiment metadata to `run_metadata`).
- HTML header rendering of experiment metadata + filter summary.
- Test updates + README documentation for the above.

**Deferred to Follow-Up Work:**
- Collapsing the fold-filter *internals* to a single threshold (`min_foldch`) — the user chose to keep the 3-level machinery; a future cleanup could remove `filter.smk`'s `_fold1/3` outputs and the unused `sample_stats_treatment` inputs.
- Broader `README.md` staleness unrelated to this change (the `slurm_partition`/`slurm_account` required-fields block now rejected by `common.smk`; the documented `complexity_filter`/`min_foldch` options that don't match current config). Only the report-columns and new-fields sections are updated here.
- A DB-backed UI / results browser (tracked in `todo.md`).
- Dropping the now-unpopulated `*_fold1/2/3` columns from existing databases (cosmetic; requires a rebuild).

**Non-goals:**
- No change to peak calling, fold filtering, blacklist/rmsk filtering, MEME, or FIMO behavior.
- No change to `meme_foldch_level` semantics.

---

## Implementation Units

### U1. Add experiment metadata config fields and register them

**Goal:** Introduce `experiment_date` and `gdna_batch` as optional top-level config fields without tripping validation.
**Requirements:** R1, R2
**Dependencies:** none
**Files:**
- `config.yaml` (user template — add near `author`, with an explanatory comment and example format)
- `config/config.yaml` (pipeline defaults — add both keys set to `null`)
- `workflow/rules/common.smk` (add `"experiment_date"`, `"gdna_batch"` to `_KNOWN_CONFIG_KEYS`)

**Approach:** Place the fields immediately after `author` in both configs. In the template, comment them as optional experiment-wide provenance (e.g. `experiment_date: "2026-07-14"  # optional, free-form`, `gdna_batch: null  # optional`). Add both literals to the `_KNOWN_CONFIG_KEYS` frozenset — **critical**, since `common.smk:14-25` raises `ValueError` on any unrecognized top-level key.
**Patterns to follow:** the existing `author` field and the `_KNOWN_CONFIG_KEYS` set in `common.smk`.
**Test scenarios:**
- Config with both fields set loads without a validation error (covered indirectly by U6 DB/HTML tests using a config dict with these keys).
- Test expectation for the YAML files themselves: none (static config) — behavior is exercised through U6.

### U2. Simplify per-sample stats computation (`collect_stats.py`)

**Goal:** Emit the canonical simplified `COLUMNS`; compute `num_peaks_filt`/`reads_in_peaks_filt` from the MEME-selected fold peaks; stop emitting `foldch_level*` and `*_fold1/2/3`.
**Requirements:** R3, R4
**Dependencies:** none (independent of U1)
**Files:**
- `workflow/scripts/collect_stats.py`
- `tests/test_collect_stats_and_report.py` (assertions on `cs.COLUMNS` — see U6 for the full test rewrite)

**Approach:** Replace `COLUMNS` with the KTD 1 list. In `main()`'s treatment block, drop the `foldch_level1/2/3` assignments and the `_fold1/2/3` reads/peaks computations. Set `num_peaks_filt = _count_lines(meme_peak_file)` and `reads_in_peaks_filt = _bedtools_intersect_count(bam, meme_peak_file)`, where `meme_peak_file = getattr(sm.input, f"peaks_fold{meme_foldch_level}")` (the existing selection used for `max_peak_score`). Keep `num_peaks` (raw `_peaks.narrowPeak`), `reads_in_peaks` (raw), `num_peaks_bl`, `num_peaks_rmsk`, `max_peak_score`, `motif_peaks` as-is.
**Patterns to follow:** the existing `meme_peak_file = getattr(sm.input, f"peaks_fold{meme_foldch_level}")` line and the `_count_lines`/`_bedtools_intersect_count` helpers already in the file.
**Test scenarios:**
- `cs.COLUMNS` equals the KTD 1 list exactly (order-sensitive).
- `cs.COLUMNS` contains none of `foldch_level1`, `num_peaks_fold1`, `reads_in_peaks_fold3`, etc.
- `cs.COLUMNS` contains `num_peaks_filt` and `reads_in_peaks_filt`.
- Existing `_count_lines`, `_safe_pct`, `_gate_subsampled`, `_parse_subsample_log` tests continue to pass unchanged (no behavior change to those helpers).

### U3. Simplify report table and add experiment header (`report.py`)

**Goal:** Update the report column set, FRiP derivation, and HTML header to the simplified scheme plus experiment metadata.
**Requirements:** R3, R5
**Dependencies:** U2 (column names must agree)
**Files:**
- `workflow/scripts/report.py`
- `tests/test_collect_stats_and_report.py` (see U6)

**Approach:**
- `make_cols()` → KTD 1 list. `INT_COLS` → `{total_reads, trimmed_reads, mapped_reads, reads_in_peaks, reads_in_peaks_filt, motif_peaks, subsampled_frags, num_peaks, num_peaks_filt, num_peaks_bl, num_peaks_rmsk}`. `PCT_COLS` → `{frip, frip_filt, mapping_pct}`.
- `build_row()` → keep `frip` (from `reads_in_peaks`/`mapped_reads`) and add `frip_filt` (from `reads_in_peaks_filt`/`mapped_reads`); remove `frip_fold1/2/3`.
- Replace `_fold_summary_html()` with a header builder that renders: (a) `Experiment date` / `gDNA batch` entries only when a non-empty value is passed in, and (b) a one-line filter note stating the fold-change threshold that `_filt` represents. Thread `experiment_date`, `gdna_batch`, and the selected fold value through `write_html(...)` (new optional keyword params, defaulting to `None`) and `run_html(...)`.
**Technical design (directional):** `write_html(rows, cols, ..., out_path, filter_foldch=None, experiment_date=None, gdna_batch=None)`; the header `<p>` is assembled from whichever of the three are provided. Keep the change backward-compatible so `render_report_html.render()` can call it with or without the new kwargs.
**Patterns to follow:** the existing `_fold_summary_html` + `<p>` insertion in `write_html`, and the `INT_COLS`/`PCT_COLS`/`_fmt_html` formatting logic.
**Test scenarios:**
- `make_cols()` equals the KTD 1 list.
- `make_cols()` excludes all `*_fold1/2/3` and `foldch_level*` names; includes `num_peaks_filt`, `reads_in_peaks_filt`, `frip_filt`.
- `frip_filt` in `PCT_COLS`; `num_peaks_filt`, `num_peaks_bl` in `INT_COLS`.
- `build_row` derives `frip_filt` from `reads_in_peaks_filt`/`mapped_reads` (e.g. 9000/90000 → `10.0`).
- `run_csv` writes a header with `frip_filt` and without `frip_fold1`.
- HTML header: `write_html` output contains the experiment date and gDNA batch strings when provided; contains neither label when both are `None`.
- HTML header: filter note renders the provided fold-change threshold value.

### U4. Thread experiment metadata into HTML render path (`render_report_html.py`, `stats.smk`)

**Goal:** Supply `experiment_date`/`gdna_batch` (and the selected fold value) to the HTML renderer from both the pipeline and the standalone CLI.
**Requirements:** R5, R8
**Dependencies:** U3 (`write_html` signature)
**Files:**
- `workflow/scripts/render_report_html.py`
- `workflow/rules/stats.smk` (`report_html` rule params; `report` rule if any param is needed there)

**Approach:**
- `render()` gains optional `experiment_date`, `gdna_batch`, `filter_foldch` params and forwards them to `write_html`.
- `render_report_html.main()` (Snakemake entry) reads them from `sm.params` (added to the `report_html` rule from `config.get("experiment_date")`, `config.get("gdna_batch")`, and the selected fold value `FOLD_LEVELS[MEME_FOLD_IDX-1]`).
- `render_report_html.cli_main()` parses `experiment_date`/`gdna_batch` from `{out_dir}/config_used.yaml` when that file exists (best-effort: ignore if missing or unparseable), and computes the fold value from the same snapshot's `macs3.foldch_levels`/`meme_foldch_level` when available.
- Add the three params to the `report_html` rule in `stats.smk`.
**Patterns to follow:** the `report`/`report_html` rules' existing `params:` blocks; the `FOLD_LEVELS`/`MEME_FOLD_IDX` globals from `common.smk`; the `config_used.yaml` snapshot written by the `onstart` hook in `Snakefile`.
**Test scenarios:**
- `render()` passes provided metadata through to the written HTML (header contains the values).
- CLI path with a `config_used.yaml` present picks up the values; CLI path without one omits them and still writes valid HTML.
- `Covers R8.` Regenerating HTML from a completed out dir succeeds whether or not `config_used.yaml` exists.

### U5. Update database schema and writes (`update_db.py`, `db.smk`)

**Goal:** Store experiment metadata per run; swap fold columns for the single `_filt` set; preserve backward compatibility.
**Requirements:** R6, R7
**Dependencies:** U2 (report CSV column names read by `read_report`)
**Files:**
- `workflow/scripts/update_db.py`
- `workflow/rules/db.smk`
- `tests/test_collect_stats_and_report.py` or a new `tests/test_update_db.py` (see U6)

**Approach:**
- `COLS`: remove `reads_in_peaks_fold1/2/3` and `num_peaks_fold1/2/3`; add `num_peaks_filt`, `reads_in_peaks_filt`. In `main()`'s per-sample loop, replace the `*_fold1/2/3` `stats.get(...)` reads with `num_peaks_filt`/`reads_in_peaks_filt`.
- `META_COLS`: add `experiment_date`, `gdna_batch` (next to `author`/`run_date`). Populate them in the `meta` dict from new `sm.params.experiment_date`/`sm.params.gdna_batch` (empty string when unset, like `author`).
- `db.smk`: add `experiment_date = config.get("experiment_date") or ""` and `gdna_batch = config.get("gdna_batch") or ""` to the `update_db` rule params.
- Rely on `_ensure_columns` for additive migration; no destructive DDL.
**Patterns to follow:** the `author`/`gene_annotation` param+meta flow already in `update_db.py`/`db.smk`; the `_ensure_columns` additive-migration helper.
**Test scenarios:**
- `update_db.COLS` excludes `reads_in_peaks_fold1`, `num_peaks_fold3`, etc.; includes `num_peaks_filt`, `reads_in_peaks_filt`.
- `update_db.META_COLS` includes `experiment_date`, `gdna_batch`.
- `Covers R7.` Writing rows to a temp SQLite DB that already has the old `*_fold1/2/3` columns succeeds; a second write with the new schema adds `num_peaks_filt`/`experiment_date`/`gdna_batch` columns and inserts without error.
- `experiment_date`/`gdna_batch` values from params land in the `run_metadata` row; empty string when unset.
- Re-running the same `output_dir` replaces its rows (idempotency preserved) in both tables.

### U6. Update and extend tests

**Goal:** Realign `tests/test_collect_stats_and_report.py` with the simplified scheme and cover the new fields.
**Requirements:** R3, R5, R6 (verification)
**Dependencies:** U2, U3, U5
**Files:**
- `tests/test_collect_stats_and_report.py`
- optionally `tests/test_update_db.py` (new, if DB tests are cleaner separated)

**Approach:** Rewrite the four currently-failing tests (`test_make_cols_order`, `test_renamed_columns_present`, `test_run_csv_excludes_mapping_rate_and_includes_renamed_frip`, `test_num_peaks_bl_after_num_peaks_filt_in_make_cols`) to assert the KTD 1 column list and `_filt` naming. Drop references to `reads_in_peaks_5fold`. Update the `_write_stats_csv` helper's `extras` usage to the new column names. Add tests for: HTML header rendering of experiment metadata (U3), and DB `COLS`/`META_COLS` membership + backward-compatible migration (U5).
**Execution note:** These tests currently fail against the diverged code — treat them as the definition of "done" for the column scheme; update them to the agreed list rather than to whatever the code happens to emit.
**Test scenarios:** (the tests themselves — enumerated per unit above under U2/U3/U5). Ensure the full file passes: `pixi run pytest tests/test_collect_stats_and_report.py`.
**Verification:** `pixi run pytest tests/` is green.

### U7. Update README documentation

**Goal:** Document the two new config fields and the simplified report columns.
**Requirements:** R9
**Dependencies:** U1, U3 (final field names + column set)
**Files:**
- `README.md`

**Approach:** Add `experiment_date` and `gdna_batch` to the config documentation near `author` (Required/Notable fields, marked optional). Update any report-columns description to the simplified `base + _filt` set and note that `_filt` reflects the `meme_foldch_level`-selected fold peaks. Scope the edit to these sections; do not attempt the broader README cleanup (deferred).
**Patterns to follow:** the existing `author: "Your Name"  # optional` documentation and the "Notable options" table format (`name=default` per the repo's param-doc convention).
**Test scenarios:** Test expectation: none — documentation only.
**Verification:** README reflects the shipped config keys and column names; no lingering references to the removed fold columns in the sections touched.

---

## Risks & Dependencies

- **Validation guard is a hard gate (R1).** If the new keys are not added to `_KNOWN_CONFIG_KEYS`, *every* run with the new config fails at DAG build. U1 must land before any config carrying the fields is used. Mitigation: U1 is the first unit and has no dependencies.
- **Column-name drift across four files.** `collect_stats.py`, `report.py`, `update_db.py`, and the tests must agree on the exact column strings. Mitigation: KTD 1 is the single source of truth; U6 asserts the list in each module.
- **`_filt` semantics could confuse users expecting a fixed threshold.** The report still uses one of the three configured fold levels. Mitigation: the HTML header states the actual fold-change value `_filt` represents (KTD 7).
- **CLI HTML regeneration without a snapshot.** Older output dirs may lack `config_used.yaml`. Mitigation: best-effort parse; omit header entries when absent (R8).
- **DB `_ensure_columns` runs DDL outside a transaction** (documented in `update_db.py`). No change to that mechanism here; new columns follow the existing safe path (KTD 6).

---

## Sources & Research

- Current 3-level implementation: `workflow/rules/filter.smk` (`filter_peaks`), `workflow/rules/common.smk` (`FOLD_LEVELS`, `MEME_FOLD_IDX`, fold validation), `workflow/scripts/collect_stats.py`, `workflow/scripts/report.py`, `workflow/scripts/update_db.py`, `workflow/rules/stats.smk`, `workflow/rules/db.smk`.
- Diverged tests (4 failing): `tests/test_collect_stats_and_report.py` (`test_make_cols_order` et al.) — encode the single-filtered-set target.
- README single-threshold docs: `README.md:100` (`macs3.min_foldch`), `README.md:311` (`min_foldch=2.0`, `*_peaks_filt.narrowPeak`).
- Origin of the 3-level design: `docs/plans/2026-06-24-001-feat-three-level-fold-filter-plan.md`.
- Config snapshot mechanism for the CLI path: `workflow/Snakefile` `onstart` hook (`config_used.yaml`), documented in `docs/plans/2026-07-08-001-feat-file-path-validation-and-config-snapshot-plan.md`.
- No relevant `docs/solutions/` learnings (directory empty); no matching brainstorm requirements doc.
