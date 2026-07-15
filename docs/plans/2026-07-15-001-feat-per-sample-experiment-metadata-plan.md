---
title: "feat: Make experiment_date and gdna_batch per-sample instead of experiment-wide"
status: completed
created: 2026-07-15
type: feat
---

## Summary

Move the `experiment_date` and `gdna_batch` provenance fields from **experiment-wide** (one value per run) to **per-sample** (one value per sample). This reverses the reporting-surface decision made in `docs/plans/2026-07-14-001-feat-simplify-report-experiment-metadata-plan.md`, which placed both as single top-level config keys shown once in the HTML report header and repeated across the DB rows.

Concretely:

1. **Config** — remove the two top-level keys; each sample under `samples:` gains its own optional `experiment_date` and `gdna_batch` (next to `r1`/`r2`). Field names are unchanged.
2. **Report** — the two values become **per-sample columns** in `report.csv`/`report.html` (positioned right after `sample`), and the single HTML header metadata line is removed. The `_filt` fold-change header note stays.
3. **Database** — each `run_metadata` row is populated from its own sample's config instead of the shared run-wide value. The `run_metadata` schema (`META_COLS`) already carries both columns, so no schema change is needed — only the value source moves.

This is a clean **replace** (per user decision): the experiment-wide keys are gone, not kept as a fallback default.

---

## Problem Frame

The just-shipped design treats `experiment_date`/`gdna_batch` as a single value for the whole run. In practice, samples in one pipeline run can come from different bench dates and different gDNA batches, so a single experiment-wide value loses that per-sample provenance. The DB already writes **one `run_metadata` row per sample** (`update_db.py` loops over `samples_cfg`), so the per-sample storage shape is already in place — today every row just gets the same repeated value. The report, by contrast, shows the metadata only once (a header line), which cannot represent per-sample variation at all.

**Success criteria:**
1. Each sample can set its own `experiment_date` and `gdna_batch` under `samples: <name>:`, alongside `r1`/`r2`; both remain optional.
2. The top-level `experiment_date`/`gdna_batch` keys no longer exist in either config file and are rejected by the validation guard with a clear "moved to per-sample" hint.
3. `report.csv` and `report.html` show `experiment_date` and `gdna_batch` as per-sample columns; two samples with different values render differently. The old single-line header metadata is gone.
4. Each `run_metadata` row stores its own sample's `experiment_date`/`gdna_batch` (empty string when unset), sourced from that sample's config block.
5. `pixi run pytest tests/` is green, with tests updated to the per-sample scheme and new coverage for per-sample rendering/storage.
6. Standalone HTML regeneration (`render_report_html.py <out_dir>`) still works and shows the per-sample columns (they come from `report.csv`, so no config snapshot dependency for these fields).

---

## Requirements

- **R1** — Remove top-level `experiment_date`/`gdna_batch` from `config.yaml` and `config/config.yaml`; add both as optional per-sample keys under each sample block. Remove both from `_KNOWN_CONFIG_KEYS` in `workflow/rules/common.smk`, and add a migration hint to the unrecognized-key error when either appears at top level.
- **R2** — Both fields stay optional and per-sample. When unset for a sample, the report cell shows `NA` (consistent with other unset cells) and the DB stores an empty string (consistent with `author`).
- **R3** — Add `experiment_date` and `gdna_batch` to `collect_stats.py`'s `COLUMNS` (immediately after `sample`), populated from new per-sample `sample_stats` rule params.
- **R4** — Add `experiment_date` and `gdna_batch` to `report.py`'s `make_cols()` (same position). Render them as text (sortable-as-text) columns in the HTML table.
- **R5** — Remove the experiment-metadata header line: drop the `experiment_date`/`gdna_batch` parameters from `_report_header_html`, `write_html`, `run_html` (report.py) and `render` (render_report_html.py), and simplify `read_run_metadata` to return only the `_filt` fold value. The filter-summary header note is retained.
- **R6** — Populate each `run_metadata` row's `experiment_date`/`gdna_batch` from its sample's config (`samples_cfg[sample]`) in `update_db.py`; remove the run-wide `experiment_date`/`gdna_batch` params from the `update_db` rule and the module-level shared reads.
- **R7** — No `run_metadata` / `pipeline_runs` schema change: `META_COLS` already lists both columns; existing DBs continue to open and accept new rows via the existing `_ensure_columns` path.
- **R8** — Update `tests/test_collect_stats_and_report.py` and `tests/test_update_db.py` to the per-sample scheme and add coverage for per-sample rendering and per-sample DB values.
- **R9** — Update `README.md` to document `experiment_date`/`gdna_batch` as optional per-sample fields (shown as report columns and stored per-sample in the DB), removing the experiment-wide description.

---

## Key Technical Decisions

**1. Per-sample values ride through the per-sample stats CSV — one source of truth.**
`collect_stats.py` writes one `{sample}.stats.csv`; that file already is the per-sample source that both `report.csv` and the HTML renderer consume. Baking `experiment_date`/`gdna_batch` into it (via new `sample_stats` rule params reading `config["samples"][wc.sample]`) means `report.py`, `report.csv`, `report.html`, and standalone CLI regeneration all pick them up with no extra plumbing and **no config-snapshot dependency**. This is cleaner than passing a separate date/batch map into the `report` rule, which would force the CLI regeneration path to re-parse `config_used.yaml` for these fields.

**2. The database sources per-sample values directly from `samples_cfg`, not from `report.csv`.**
`update_db.py:main()` already iterates `samples_cfg` and has each sample's config dict (`scfg`) in hand. Read `scfg.get("experiment_date")`/`scfg.get("gdna_batch")` (with `or ""`) in the loop. This keeps the DB write independent of report-column presence and mirrors the existing `r1`/`r2` per-sample sourcing right beside it.

**3. Column position: immediately after `sample`.**
`sample, experiment_date, gdna_batch, total_reads, …`. Provenance reads naturally next to the sample name and stays out of the numeric metric block. Applied identically in `collect_stats.COLUMNS` and `report.make_cols()` (the two must agree).

**4. Report cells: `NA` when unset; DB: empty string when unset.**
`collect_stats.py` emits the value or `NA` (so the report table is visually consistent with every other unset cell and the existing `.na` styling). The DB independently stores `""` for unset, matching how `author` behaves. Each convention is internally consistent; the divergence is intentional and localized.

**5. Text columns for HTML sorting.**
`experiment_date`/`gdna_batch` are free-form strings, so their `<th>` cells get `data-col-type="text"` (like `sample`), not the default `numeric`. Introduce a small `TEXT_COLS` set in `write_html` rather than special-casing `sample` alone. `_fmt_html` already routes non-INT/non-PCT columns through `html.escape`, so no formatting change is needed.

**6. Header metadata line is removed, filter note kept.**
Per-sample values now live in table columns, so the run-wide header entry is meaningless. `_report_header_html` loses its `experiment_date`/`gdna_batch` arguments and keeps only the `filter_foldch` note. `read_run_metadata` narrows from `(filter_foldch, experiment_date, gdna_batch)` to just the fold value; all callers (`report.cli_main`, `render_report_html.cli_main`) are updated.

**7. Replace, don't augment — with a helpful migration error.**
The top-level keys are removed from both config files and from `_KNOWN_CONFIG_KEYS`, so an old config carrying them at top level will fail validation. Add a hint to the `common.smk` unrecognized-key error (alongside the existing `input_control`/`slurm_*` hints): "`experiment_date`/`gdna_batch` are now set per-sample under `samples: <name>:`." This turns a silent break into an actionable message.

**8. No schema migration.**
`run_metadata`'s `META_COLS` already includes `experiment_date`/`gdna_batch`. Only the value *source* changes (run-wide param → per-sample config). `_ensure_columns` and the additive-migration guarantees are untouched.

---

## Scope Boundaries

**In scope:**
- Move `experiment_date`/`gdna_batch` from top-level to per-sample in both config files + validation allow-list and error hint.
- Per-sample columns in `collect_stats.py` / `report.py` and the `sample_stats` rule params.
- Removal of the HTML header metadata line and the now-unused params/return values through the report + render path.
- Per-sample DB sourcing in `update_db.py` + `db.smk` param removal.
- Test updates + README documentation for the above.

**Deferred to Follow-Up Work:**
- Per-sample-key validation (rejecting unknown keys *inside* a sample block, e.g. a typo like `expriment_date`). Currently only top-level keys are validated; adding nested validation is a separate hardening task.
- Backfilling the per-sample values into `pipeline_runs` (they live in `run_metadata`; no request to duplicate them into the metrics table).
- Broader README staleness unrelated to these two fields.

**Non-goals:**
- No change to peak calling, fold filtering, MEME/FIMO, or any QC metric computation.
- No change to the `run_metadata`/`pipeline_runs` schemas beyond the value source.
- No fallback/default behavior — per-sample fully replaces experiment-wide.

---

## Implementation Units

### U1. Move config fields to per-sample and update validation

**Goal:** Relocate `experiment_date`/`gdna_batch` from top-level to per-sample keys and make the validation guard reject the old top-level form with a clear hint.
**Requirements:** R1
**Dependencies:** none
**Files:**
- `config.yaml` (user template)
- `config/config.yaml` (pipeline defaults)
- `workflow/rules/common.smk`

**Approach:**
- `config.yaml`: delete the top-level `experiment_date`/`gdna_batch` block (currently ~lines 7-10). Under the `samples:` block, add both keys (set to `null`) to each sample next to `r1`/`r2`, with a one-line comment noting they are optional per-sample provenance (e.g. `experiment_date: null   # optional, e.g. "2026-07-14"`).
- `config/config.yaml`: delete the top-level keys (~lines 19-21); add both per-sample under each sample in the defaults `samples:` block.
- `common.smk`: remove `"experiment_date"` and `"gdna_batch"` from `_KNOWN_CONFIG_KEYS`. In the `_unexpected_keys` error builder, add a hint when either key is present at top level: point the user to the per-sample location.

**Patterns to follow:** the existing `input_control` and `slurm_partition`/`slurm_account` hint branches in `common.smk`; the existing per-sample `r1`/`r2` layout in both configs.
**Test scenarios:**
- Config with per-sample `experiment_date`/`gdna_batch` set loads without a validation error (exercised indirectly via U6 stats/DB tests using a `samples_cfg` dict carrying these keys).
- A top-level `experiment_date` triggers the unrecognized-key `ValueError`, and the message mentions the per-sample location. (Optional unit test if `common.smk` logic is extracted; otherwise verification is manual — see Verification.)
- Static config files: `Test expectation: none` — behavior is exercised through U6.
**Verification:** A dry run (`pixi run snakemake -n --configfile config.yaml`) with per-sample fields set builds the DAG; a dry run with a stray top-level `experiment_date` fails with the hinted message.

### U2. Emit per-sample metadata from `collect_stats.py`

**Goal:** Add `experiment_date`/`gdna_batch` to the per-sample stats CSV, populated from new `sample_stats` params.
**Requirements:** R2, R3
**Dependencies:** none (independent of U1; reads params, not top-level config)
**Files:**
- `workflow/scripts/collect_stats.py`
- `workflow/rules/stats.smk` (`sample_stats` rule params)

**Approach:**
- `COLUMNS`: insert `"experiment_date"` and `"gdna_batch"` immediately after `"sample"`.
- `main()`: set `row["experiment_date"] = str(sm.params.experiment_date) if sm.params.experiment_date else NA` and likewise for `gdna_batch`.
- `stats.smk` `sample_stats` rule `params:`: add `experiment_date = lambda wc: config["samples"][wc.sample].get("experiment_date")` and `gdna_batch = lambda wc: config["samples"][wc.sample].get("gdna_batch")`.

**Patterns to follow:** the existing per-sample lambda params in `sample_stats` (`is_pe = lambda wc: wc.sample in PE_SAMPLES`); the `row[...] = ... if ... else NA` idiom already in the file.
**Test scenarios:**
- `cs.COLUMNS` equals the new list exactly (order-sensitive), with `experiment_date`/`gdna_batch` right after `sample`.
- `cs.COLUMNS` still contains none of the removed `foldch_level*`/`*_fold1/2/3` names (regression guard preserved).
- Existing helper tests (`_count_lines`, `_safe_pct`, `select_meme_peak_file`, etc.) continue to pass unchanged.

### U3. Add per-sample report columns and remove header metadata (`report.py`)

**Goal:** Surface the two columns in `report.csv`/`report.html`, render them as text, and strip the run-wide header metadata line.
**Requirements:** R2, R4, R5
**Dependencies:** U2 (column names must agree)
**Files:**
- `workflow/scripts/report.py`
- `tests/test_collect_stats_and_report.py` (see U6)

**Approach:**
- `make_cols()`: insert `"experiment_date"`, `"gdna_batch"` immediately after `"sample"`.
- `build_row()`: no change — it copies `dict(data)` from the stats CSV, so the new columns flow through automatically; `run_csv`'s `{c: row.get(c, "NA") …}` fills `NA` when absent.
- `write_html()`: introduce `TEXT_COLS = {"sample", "experiment_date", "gdna_batch"}` and use it for the `data-col-type="text"` header-cell branch (replacing the `col == "sample"` special case). No change to `_fmt_html` (non-INT/non-PCT → `html.escape`).
- Remove `experiment_date`/`gdna_batch` params from `_report_header_html`, `write_html`, and `run_html`; keep `filter_foldch`.
- `read_run_metadata()`: return only the fold value (drop the two metadata elements). Update `cli_main()`'s unpacking and its `run_html(...)` call accordingly.

**Technical design (directional):** `write_html(rows, cols, …, out_path, filter_foldch=None, control_samples=None)` — metadata kwargs removed. `read_run_metadata(out_dir) -> filter_foldch | None`.
**Patterns to follow:** existing `INT_COLS`/`PCT_COLS` set-membership checks; the current `_report_header_html` structure (keep the `filter_foldch` `<p>` block).
**Test scenarios:**
- `make_cols()` equals the new list; `experiment_date`/`gdna_batch` appear right after `sample`.
- `run_csv` writes a header containing `experiment_date` and `gdna_batch`.
- Two samples with different `experiment_date`/`gdna_batch` values in their stats CSVs produce two report rows with the respective values (per-sample, not shared).
- `write_html` output contains each row's `experiment_date`/`gdna_batch` cell values; the corresponding `<th>` carries `data-col-type="text"`.
- `write_html` no longer emits an "Experiment date"/"gDNA batch" **header** line (only per-row cells); the `filter_foldch` note still renders when provided.
- `read_run_metadata` returns the fold value alone from a snapshot; returns `None` when the snapshot is missing.

### U4. Update HTML render path (`render_report_html.py`, `stats.smk`)

**Goal:** Align the renderer and the `report_html` rule with the reduced `write_html`/`render`/`read_run_metadata` signatures.
**Requirements:** R5
**Dependencies:** U3 (`write_html` / `read_run_metadata` signatures)
**Files:**
- `workflow/scripts/render_report_html.py`
- `workflow/rules/stats.smk` (`report_html` rule params)

**Approach:**
- `render()`: drop `experiment_date`/`gdna_batch` params; forward only `filter_foldch` and `control_samples` to `write_html`.
- `main()` (Snakemake entry): stop reading `sm.params.experiment_date`/`sm.params.gdna_batch`.
- `cli_main()`: update the `read_run_metadata` unpacking to the single return value and the `render(...)` call.
- `stats.smk` `report_html` rule: remove the `experiment_date` and `gdna_batch` params; keep `filter_foldch`, `control_samples`, dirs, etc.

**Patterns to follow:** the existing `report_html` `params:` block and the `render()`/`write_html` call wiring.
**Test scenarios:**
- `Covers R6/CLI.` Standalone regeneration (`render_report_html.py <out_dir>`) against a completed out dir writes valid HTML whose table includes the per-sample `experiment_date`/`gdna_batch` columns (values sourced from `report.csv`), with or without `config_used.yaml` present.
- `render()` no longer accepts (or requires) metadata kwargs — call succeeds with the reduced signature.
**Verification:** Regenerate HTML from a completed output directory; confirm the two columns show per-sample values and the header no longer has a metadata line.

### U5. Source per-sample metadata in the database (`update_db.py`, `db.smk`)

**Goal:** Populate each `run_metadata` row from its own sample's config; remove the run-wide params.
**Requirements:** R2, R6, R7
**Dependencies:** none (reads `samples_cfg`, which is already a param)
**Files:**
- `workflow/scripts/update_db.py`
- `workflow/rules/db.smk`

**Approach:**
- `update_db.py`: remove the module-level `experiment_date`/`gdna_batch = sm.params...` reads. In the per-sample loop, set `meta["experiment_date"] = scfg.get("experiment_date") or ""` and `meta["gdna_batch"] = scfg.get("gdna_batch") or ""` (right next to `r1`/`r2` sourcing).
- `db.smk`: remove `experiment_date = config.get("experiment_date") or ""` and `gdna_batch = config.get("gdna_batch") or ""` from the `update_db` rule params.
- `META_COLS`/`COLS` unchanged; `_ensure_columns` path unchanged.

**Patterns to follow:** the adjacent `r1 = get_r1(scfg)` / `row["r2"] = get_r2(scfg) or ""` per-sample sourcing in the same loop.
**Test scenarios:**
- `update_db.META_COLS` still includes `experiment_date`, `gdna_batch` (unchanged).
- Writing meta rows where each sample carries a **different** `experiment_date`/`gdna_batch` stores the distinct values per row (extend the existing `test_experiment_metadata_lands_in_run_metadata`).
- `Covers R7.` The existing legacy-schema migration test still passes (no schema change).
- Re-running the same `output_dir` replaces its `run_metadata` rows (idempotency preserved).
**Verification:** `main()`'s per-sample sourcing is exercised by a pipeline smoke run (not unit-tested, since `main()` needs the Snakemake object); confirm two samples with different config values land distinct rows in `run_metadata`.

### U6. Update and extend tests

**Goal:** Realign existing tests to the per-sample scheme and add per-sample coverage.
**Requirements:** R8
**Dependencies:** U2, U3, U5
**Files:**
- `tests/test_collect_stats_and_report.py`
- `tests/test_update_db.py`

**Approach:**
- `test_make_cols_order` and `test_collect_stats_columns_simplified`: update both expected lists to include `experiment_date`/`gdna_batch` after `sample`.
- Replace the header-metadata tests (`test_html_header_shows_experiment_metadata`, `test_html_header_omits_metadata_when_unset`) with: (a) per-row rendering of the two columns for two samples with different values, and (b) absence of the old header metadata line. Keep `test_html_header_shows_filter_foldch`.
- Update `test_read_run_metadata_from_snapshot` / `test_read_run_metadata_missing_snapshot` to the new single-value return shape.
- `test_update_db.py`: extend `test_experiment_metadata_lands_in_run_metadata` (or add a sibling) to assert distinct per-sample values across two rows.
**Execution note:** Treat the failing/updated column-order tests as the definition of "done" for the column list — update them to the agreed order rather than to whatever the code emits.
**Test scenarios:** (enumerated per unit under U2/U3/U5.) Ensure the full suite passes.
**Verification:** `pixi run pytest tests/` is green.

### U7. Update README documentation

**Goal:** Document the two fields as per-sample.
**Requirements:** R9
**Dependencies:** U1, U3 (final field placement + column set)
**Files:**
- `README.md`

**Approach:** Remove the top-level `experiment_date`/`gdna_batch` example (~lines 62-63) and the two experiment-wide table rows (~lines 89-90). Document both as **optional per-sample** fields in the samples section (next to `r1`/`r2`), noting they appear as report columns and are stored per-sample in the results database. Follow the repo's `name=default` param-doc convention.
**Patterns to follow:** the existing samples/`r1`/`r2` documentation and the "Notable options" table format.
**Test scenarios:** `Test expectation: none` — documentation only.
**Verification:** README shows the fields under samples; no lingering top-level/experiment-wide references.

---

## Risks & Dependencies

- **Old configs break at validation (intended).** Any existing config with top-level `experiment_date`/`gdna_batch` will fail after U1. Mitigation: the added `common.smk` hint gives an actionable "moved to per-sample" message (KTD 7); README updated (U7).
- **Column-name/order agreement across `collect_stats.py`, `report.py`, and tests.** All must list the two new columns in the same position. Mitigation: KTD 3 fixes the position; U6 asserts the list in both modules.
- **`main()` per-sample DB sourcing is not unit-tested.** `update_db.main()` requires the Snakemake object. Mitigation: unit tests cover `write_meta_rows` with distinct per-sample values; the sourcing wiring is verified by a smoke run (U5 Verification).
- **Mixed-provenance DBs.** Existing `run_metadata` rows hold the old repeated value; new runs hold per-sample values in the same column. This is expected and needs no migration (KTD 8).

---

## Sources & Research

- Prior (experiment-wide) design being reversed: `docs/plans/2026-07-14-001-feat-simplify-report-experiment-metadata-plan.md`.
- Config + validation: `config.yaml`, `config/config.yaml`, `workflow/rules/common.smk` (`_KNOWN_CONFIG_KEYS`, unrecognized-key hints).
- Stats/report/render path: `workflow/scripts/collect_stats.py` (`COLUMNS`, `main`), `workflow/scripts/report.py` (`make_cols`, `write_html`, `_report_header_html`, `read_run_metadata`, `cli_main`), `workflow/scripts/render_report_html.py` (`render`, `cli_main`), `workflow/rules/stats.smk` (`sample_stats`, `report`, `report_html`).
- Database path: `workflow/scripts/update_db.py` (`META_COLS`, `main` per-sample loop, `_ensure_columns`), `workflow/rules/db.smk` (`update_db` params).
- Tests to update: `tests/test_collect_stats_and_report.py`, `tests/test_update_db.py`.
- No matching brainstorm requirements doc (the one under `docs/brainstorms/` covers multi-user concurrency, unrelated); `docs/solutions/` has no relevant learnings.
