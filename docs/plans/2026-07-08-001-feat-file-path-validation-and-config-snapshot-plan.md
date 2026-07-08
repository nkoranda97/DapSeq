---
title: "feat: Add file-path validation and config snapshot"
status: active
created: 2026-07-08
type: feat
---

## Summary

Two early-startup improvements: (1) a batched file-path existence checker in `common.smk` that reports every bad path in a single clear error before any jobs run, and (2) an `onstart` hook in `Snakefile` that writes the full merged effective config to `{output_dir}/config_used.yaml` at the start of every real run.

---

## Problem Frame

Users currently get opaque Snakemake `MissingInputException` errors or tool-level failures mid-pipeline when a file path in their config is wrong — often after SLURM jobs have already been submitted. There is also no persistent record of the exact parameters used for a completed run; the results directory alone cannot reproduce the run.

**Success criteria:**
1. A config with any non-existent path fails immediately at DAG-build time, before any jobs are submitted, with one message naming every bad path and its config key.
2. Every real pipeline invocation writes `{output_dir}/config_used.yaml` containing the full merged effective configuration (user values + all pipeline defaults) with a header identifying the run timestamp and source config file(s).
3. Path validation runs on `--dryrun`; the config snapshot does not write on `--dryrun` (Snakemake's `onstart` hook is intentionally skipped on dry runs — correct behavior).

---

## Requirements

- **R1** — Validate existence of all file-path config fields: `genome_ref`, per-sample `r1`/`r2`, optional filter files when their feature is enabled, `gene_annotation` (if set), and `factorbook` paths (if explicitly set).
- **R2** — Accumulate all path errors and raise a single `ValueError` listing every bad path; do not fail on the first error encountered.
- **R3** — Write the full merged `config` dict as YAML to `{output_dir}/config_used.yaml` via the `onstart` hook on every real run (not dry-runs).
- **R4** — Include a `#`-prefixed header block in `config_used.yaml` with the run timestamp and absolute path(s) to the source config file(s).
- **R5** — Create `{output_dir}` with `exist_ok=True` before writing the snapshot (Snakemake only creates rule output directories when a rule fires, not before).
- **R6** — All existing validation in `common.smk` must remain intact and unchanged.

---

## Key Technical Decisions

**1. Path validation lives in `common.smk` module-level code — not a rule.**
Consistent with all existing validation (`_missing_fields`, `_unexpected_keys`, fold-change guards). Module-level Python runs at DAG-build time, which happens on both real runs and `--dryrun`. No new Snakemake rule is added.

**2. Batch all errors before raising.**
Use a `_path_errors: list` accumulator, then raise one `ValueError` with all entries. This is the established pattern in `common.smk:40–52` (`_missing_fields` block).

**3. `onstart` hook for config archival.**
Fires after DAG build (so `OUT`, `SAMPLES`, `CONTROL`, etc. are resolved and path validation has already passed) but before any jobs run. Confirmed to NOT fire on `--dryrun` in Snakemake 9 (`workflow.py:1451–1452`).

**4. Dump the merged `config` dict via `yaml.dump`.**
User confirmed: merged effective config (captures every pipeline default, not just the user's overrides). PyYAML is available in the pixi environment (Snakemake itself depends on it). Use `sort_keys=True, default_flow_style=False` for readable YAML.

**5. Identify the user's config file via `workflow.configfiles`.**
At runtime, `workflow.configfiles` is a list of `Path` objects for CLI `--configfile` args followed by the string `"config/config.yaml"` appended by the `configfile:` directive. Filter out entries matching `"config/config.yaml"` to surface the user-supplied path(s) for the header comment.

**6. `bbduk.adapters` is a special case.**
When null/unset, the pipeline defaults to the string `"adapters"` (bbduk's built-in reference, not a filesystem path). Only report a path error if `adapters` is explicitly set to a truthy value.

**7. Placement in `common.smk`.**
Insert the path-validation block after `PEAKS_FILTER_SUFFIX` is defined (~line 108), so `_BL_ENABLED` and `_RMSK_ENABLED` are already available.

---

## Scope Boundaries

**In scope:**
- File existence checks for all user-supplied path fields
- `onstart`-based merged config YAML archival to `{output_dir}/config_used.yaml`

**Deferred to Follow-Up Work:**
- Adding a `user_config_path` column to the SQLite DB in `update_db.py` (schema migration is easy via the existing `_ensure_columns` helper, but this is a separate concern)
- File content validation (FASTA format, BED validity, BAM header correctness)
- Absolute-vs-relative path detection or normalization
- Validating numeric config fields (mapq ranges, thread counts, fold-change business logic beyond the existing checks)

---

## Implementation Units

### U1. File path existence validation in common.smk

**Goal:** Report every non-existent file path from the config as a single batched `ValueError` before any jobs run.

**Requirements:** R1, R2, R6

**Dependencies:** none

**Files:**
- `workflow/rules/common.smk`

**Approach:**

Add a `_path_errors = []` block in `common.smk` immediately after `PEAKS_FILTER_SUFFIX` is assigned (~line 108). For each file-path field:

- `genome_ref` — always check (`os.path.exists`); it is guaranteed non-null at this point by the earlier required-fields check.
- Sample r1 — iterate `SAMPLES`, call `get_r1(s)` (normalizes to list), check each path.
- Sample r2 — iterate `SAMPLES`, read `config["samples"][s].get("r2")`; if set, normalize to list (handle both `str` and `list[str]`), check each path.
- `bbduk.adapters` — check only if the value is truthy (explicitly set by user; skip the bare-string `"adapters"` fallback).
- `blacklist_filter.bed` — check only when `_BL_ENABLED` is true; also error if enabled but no bed path is set.
- `rmsk_filter.txt` — check only when `_RMSK_ENABLED` is true; also error if enabled but no txt path is set.
- `gene_annotation` — check only when the value is non-null.
- `factorbook.tsv`, `factorbook.meme` — check only when each is explicitly set to a non-null string.

If `_path_errors` is non-empty, raise a single `ValueError` with all entries formatted as indented lines plus a closing hint: `"Verify that paths are absolute and the files are accessible from this machine."`

**Patterns to follow:**
- `common.smk:40–52` — `_missing_fields` accumulator and single-raise style
- Error entry format: `"  <config_key>: '<path>' does not exist"` (two leading spaces, key first, value in quotes)

**Test scenarios:**

- `genome_ref` set to a non-existent file → error message names `genome_ref` and the bad path
- A sample's `r1` is a non-existent single file → error names `samples.<name>.r1` and path
- A sample's `r1` is a list with one missing file → that specific path is named in the error
- A PE sample's `r2` is a non-existent file → error names `samples.<name>.r2` and path
- `bbduk.adapters: null` (default) → no error raised
- `bbduk.adapters: /explicit/adapters.fa` where file is missing → error raised
- `blacklist_filter.enabled: true` with non-existent bed path → error
- `blacklist_filter.enabled: true` with valid bed path → no error
- `blacklist_filter.enabled: false` with bad bed path → no error (path not checked when disabled)
- `rmsk_filter` scenarios mirror the blacklist scenarios above
- `gene_annotation: null` → no error
- `gene_annotation: /missing.gtf` → error
- `factorbook.tsv: null` → no error (uses repo-relative default)
- `factorbook.tsv: /explicit/path.tsv` where file is missing → error
- Multiple bad paths across different fields → all reported together in one `ValueError`
- All paths valid → no error; pipeline proceeds normally

**Verification:** Run `pixi run snakemake --dryrun --configfile config.yaml` with a config pointing to a non-existent `genome_ref` and two bad sample `r1` paths; confirm all three appear in one error. Run with a fully valid config and confirm no error.

---

### U2. Merged config snapshot via onstart hook

**Goal:** Write `{output_dir}/config_used.yaml` containing the full merged effective configuration at the start of every real pipeline run.

**Requirements:** R3, R4, R5

**Dependencies:** U1 (path validation passes first — guaranteed by Snakemake's execution order: module-level Python runs before `onstart`)

**Files:**
- `workflow/Snakefile`

**Approach:**

Add `import yaml` and `from datetime import datetime` at the top of `Snakefile` alongside the existing `from snakemake.utils import min_version`. (`os` is already available in the Snakefile namespace via `common.smk`.)

After the final `include:` directive, add an `onstart:` block that:

1. Calls `os.makedirs(OUT, exist_ok=True)` to ensure the output directory exists.
2. Computes the timestamp: `datetime.now().isoformat(timespec="seconds")`.
3. Filters `workflow.configfiles` to identify the user-supplied path(s): `[str(p) for p in workflow.configfiles if str(p) != "config/config.yaml"]`.
4. Opens `os.path.join(OUT, "config_used.yaml")` for writing (overwrites on each run).
5. Writes `#`-prefixed header lines: pipeline label, `Run timestamp:`, `Source config:`.
6. Calls `yaml.dump(dict(config), f, default_flow_style=False, sort_keys=True)`.

Use private variable names prefixed with `_` (e.g., `_snapshot_path`, `_ts`, `_user_cfgs`) to avoid polluting the Snakemake globals namespace.

**Patterns to follow:**
- `workflow/Snakefile` top-of-file import style (one import per line, stdlib before third-party)
- `onstart:` indented body syntax (no `def`, no argument)

**Test scenarios:**

- Fresh run with valid config → `{output_dir}/config_used.yaml` is created; contains all merged keys (`samples`, `genome_ref`, `threads`, `macs3`, `meme`, etc.)
- Re-run with modified config (e.g., changed `threads`) → file is overwritten with new timestamp and new value
- `--dryrun` invocation → `config_used.yaml` is NOT created
- `output_dir` does not exist yet → `os.makedirs` creates it and file is written successfully
- Header comment in output file → contains `Run timestamp:` in ISO 8601 format and `Source config:` with absolute path to the user's config
- Pipeline defaults (e.g., `bamcoverage`, `fimo` sections) are present in the snapshot even when the user's config omits them (merged from `config/config.yaml`)

**Verification:** Run a full pipeline invocation; open `{output_dir}/config_used.yaml`; confirm all top-level keys are present and the header comment contains the timestamp and config path. Confirm `--dryrun` does not create the file.

---

## Sources & Research

- `workflow/rules/common.smk:5–52` — existing validation patterns (`_KNOWN_CONFIG_KEYS`, `_missing_fields`)
- `workflow/rules/common.smk:105–108` — `_BL_ENABLED`, `_RMSK_ENABLED`, `PEAKS_FILTER_SUFFIX`
- `workflow/Snakefile:1–7` — current import and directive structure
- Snakemake 9.20.0 `workflow.py:1451–1452` — `onstart` skips dry-runs
- Snakemake 9.20.0 `workflow.py:1683–1687` — `onstart` hook registration
- Snakemake 9.20.0 `workflow.py:209` — `workflow.configfiles` initialization (CLI paths first, then `configfile:` directive appends pipeline default path)
