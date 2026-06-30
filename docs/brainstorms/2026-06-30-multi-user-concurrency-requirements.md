# Multi-User Concurrency Requirements

**Date:** 2026-06-30
**Status:** Ready for planning

## Goal

Allow multiple users on a shared HPC cluster to run the pipeline simultaneously from the same shared repo directory, each with their own config and output directory, without stepping on each other.

## Context

Users already invoke the pipeline from the repo root with `--configfile /path/to/their/config.yaml`. Each user points `output_dir` to their own location. Snakemake's locks are keyed per output directory, so concurrent runs with different `output_dir` values do not conflict at the scheduling layer. Two gaps remain: silent failure when `--configfile` is omitted, and unprotected concurrent writes to the shared SQLite database.

## Requirements

### 1. Early validation for missing user config

- If any of the following are still null/unset after config merging, the pipeline must fail immediately with a clear, actionable error message (not a silent empty run or a cryptic Snakemake error):
  - `output_dir`
  - `genome_ref`
  - At least one sample with `r1` set
- The error message should tell the user to pass `--configfile /path/to/your/config.yaml`.
- Validation lives in `workflow/rules/common.smk`, alongside existing validation logic.

### 2. SQLite WAL mode for the shared database

- The shared `pipeline_db.db` (defaults to repo root) must be opened in WAL (Write-Ahead Logging) mode so concurrent writes from simultaneous `update_db` rule executions do not corrupt the database or raise errors.
- WAL mode is set once at connection open time in `workflow/scripts/update_db.py`; no schema changes required.
- No retry logic is needed beyond what WAL provides by default.

### 3. README multi-user documentation

- Add a short section to `README.md` (under Running or a new "Multi-user / shared installation" heading) that explains:
  - Multiple users can run from the shared repo root simultaneously.
  - Each user must pass `--configfile` pointing to their own config file with a unique `output_dir`.
  - The shared database at the repo root accumulates results from all users; `db_path` in a user's config can redirect to a private DB if needed.

## Non-goals

- Per-user working directories (`--directory` isolation) — not needed; Snakemake locks are per-output-dir.
- Per-user databases — the shared DB is intentional.
- A wrapper/launcher script — the existing README invocation is sufficient.
- Fixing `configfile: "config/config.yaml"` to an absolute path — only needed if running from a non-repo-root directory, which is not the target use case.

## Success criteria

- Running `pixi run snakemake --profile profiles/slurm` (no `--configfile`) produces a clear error message within seconds, before any jobs are submitted.
- Two users running simultaneously with different configs and `output_dir` values complete successfully with no DB corruption.
- README makes the multi-user pattern clear to a new user without further explanation.
