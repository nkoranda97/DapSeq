---
title: "feat: Standalone HTML report generation"
date: 2026-06-02
status: active
type: feat
---

# feat: Standalone HTML report generation

## Summary

Split the combined `report` Snakemake rule into two independent rules so that `report.html` can be regenerated from already-computed data (per-sample stats CSVs + logo PNGs) without re-running any heavy computation. Also wire `report.html` into the `qc_done` aggregate target and fix `motifs_done` to include the newly-added `logo_rc1.png` outputs.

---

## Problem Frame

The `report` rule currently produces both `report.csv` and `report.html` as outputs of a single job. Because Snakemake treats all outputs of a rule as one atomic unit, regenerating the HTML (e.g., after a styling change or a logo re-render) forces the rule to re-run from scratch and re-read all per-sample stats CSVs. Additionally, `report.html` is not listed in `qc_done`, so it is never requested by the normal pipeline run — it is a dead output.

---

## Requirements

- R1: `report.html` can be generated or regenerated independently, given that `report.csv` and all logo PNGs already exist.
- R2: The normal full pipeline run (`qc_done`) produces `report.html` without any manual intervention.
- R3: The CLI entry point (`report.py`) continues to regenerate both CSV and HTML for convenience.
- R4: A standalone CLI entry point exists that regenerates only the HTML from an existing output directory.
- R5: `motifs_done` correctly declares `logo_rc1.png` as part of its outputs.

---

## Key Technical Decisions

**Split at the rule boundary, not the script boundary.** `report.py` keeps `run()` (does both CSV + HTML) for the CLI case. A new `render_report_html.py` script handles Snakemake's HTML-only mode, importing `write_html` and `logo_to_base64` from `report.py`. This avoids duplicating HTML rendering logic while giving Snakemake explicit, separate dependency edges.

**`report` rule keeps its name.** The rule that produces `report.csv` keeps the existing `report` name so any `snakemake report` invocations by users continue to work. The new rule is `report_html`.

**`report.py` `main()` becomes CSV-only.** The Snakemake-facing `main()` in `report.py` drops HTML outputs and logo inputs; `cli_main()` continues to do both (R3).

---

## Implementation Units

### U1. Make `report.py` CSV-only for Snakemake, keep CLI doing both

**Goal:** Strip logo inputs and HTML output from `report.py`'s Snakemake entry point while preserving the combined CLI behavior.

**Requirements:** R1, R3

**Dependencies:** none

**Files:**
- `workflow/scripts/report.py`

**Approach:** In `main()`, remove `sm.output.html`, `sm.params.meme_dir`, and `sm.params.factorbook_dir`. Call only `run_csv()` (or rename `run()` → `run_csv()` to clarify intent). `cli_main()` continues to call the combined `run()` which generates both outputs from an existing out directory.

Refactor `run()` into two functions:
- `run_csv(treatment_samples, stats_dir, csv_out)` — builds rows, writes CSV
- `run_html(treatment_samples, stats_dir, meme_dir, factorbook_dir, csv_out, html_out)` — calls `run_csv` then builds logo maps and calls `write_html`

`main()` calls `run_csv`. `cli_main()` calls `run_html` (which does both).

**Patterns to follow:** Existing `main()` / `cli_main()` split in `collect_stats.py` and `report.py`.

**Test scenarios:**
- Happy path: `main()` is called with Snakemake context having only `sm.output.csv` and `sm.params.stats_dir` + `sm.params.treatment_samples`; `report.csv` is written, no HTML is attempted.
- Happy path: `cli_main()` with a populated out directory writes both `report.csv` and `report.html`.
- Edge case: `cli_main()` with no `*.stats.csv` files exits with a clear error message.

**Verification:** `report.py`'s `main()` references no logo paths and no `sm.output.html`; `cli_main()` still produces an HTML file alongside the CSV.

---

### U2. Add `render_report_html.py`

**Goal:** New script that reads `report.csv` and logo directories and writes `report.html`. Used by the `report_html` Snakemake rule and callable standalone from the CLI.

**Requirements:** R1, R4

**Dependencies:** U1

**Files:**
- `workflow/scripts/render_report_html.py` *(new)*

**Approach:** Import `write_html`, `logo_to_base64`, and `make_cols` from `report`. The Snakemake `main()` reads `sm.input.report_csv`, `sm.params.meme_dir`, `sm.params.factorbook_dir`, `sm.params.treatment_samples`, and writes to `sm.output.html`. The `cli_main()` takes `out_dir` as a positional argument, discovers treatment samples by reading `report.csv` (first column), and writes `report.html` alongside it.

**Patterns to follow:** `render_meme_logo.py` for the single-output Snakemake script pattern; `report.py`'s `cli_main()` for the out-directory CLI pattern.

**Test scenarios:**
- Happy path: called with a populated meme dir and factorbook dir; produces `report.html` with embedded base64 logo images.
- Edge case: logo file is empty or missing for a sample; `logo_to_base64` returns `None`; HTML renders "NA" cell for that sample without raising.
- Happy path (CLI): `cli_main()` with an out directory that has `report.csv` but no logos; produces valid HTML with all logo cells as "NA".
- Edge case (CLI): `report.csv` absent; exits with a clear error.

**Verification:** Running `pixi run python workflow/scripts/render_report_html.py <out_dir>` on a directory with `report.csv` and logo PNGs produces a valid, openable `report.html` without errors.

---

### U3. Split `report` rule; add `report_html` rule; fix `qc_done` and `motifs_done`

**Goal:** Wire the new rule split into Snakemake so `report.html` is a first-class target and `qc_done` requests it automatically.

**Requirements:** R1, R2, R5

**Dependencies:** U1, U2

**Files:**
- `workflow/rules/stats.smk`
- `workflow/rules/common.smk`

**Approach:**

In `stats.smk`:
- `report` rule: inputs are `stats_csvs` only; output is `csv` only; params drop `meme_dir` and `factorbook_dir`; script remains `report.py`.
- `report_html` rule (new): inputs are `report_csv = OUT + "/stats/report.csv"`, `meme_logos`, `meme_logos_rc`, `factorbook_logos`; output is `html = OUT + "/stats/report.html"`; params carry `meme_dir`, `factorbook_dir`, `treatment_samples`; script is `render_report_html.py`. Inherits `report` resource profile.

In `common.smk`:
- `qc_done`: add `OUT + "/stats/report.html"` to its input list.
- `motifs_done`: add `expand(OUT + "/meme/{sample}/summits/logo_rc1.png", sample=TREATMENT_SAMPLES)` and `expand(OUT + "/meme/{sample}/peaks/logo_rc1.png", sample=TREATMENT_SAMPLES)` to its input list (these are now outputs of `meme_logo_summits` and `meme_logo_peaks`).

**Patterns to follow:** Existing `report` rule structure; `motifs_done` expand pattern.

**Test scenarios:**
- Happy path: `snakemake <out>/stats/report.html` with all logos and `report.csv` present runs `report_html` rule only (no upstream stats recomputation).
- Happy path: fresh pipeline run hits `qc_done` target and produces both `report.csv` and `report.html`.
- Integration: `snakemake motifs_done` succeeds after adding `logo_rc1.png` to its inputs (previously this would pass even with `logo_rc1.png` missing).

**Verification:** `snakemake --dry-run <out>/stats/report.html` when `report.csv` and logos exist shows only the `report_html` job in the execution plan with no upstream jobs.

---

## Scope Boundaries

### In scope
- Splitting the combined rule into `report` (CSV) + `report_html` (HTML)
- New `render_report_html.py` script with Snakemake and CLI entry points
- `qc_done` including `report.html`
- `motifs_done` including `logo_rc1.png`

### Out of scope
- Changing the HTML content or styling
- Adding other standalone re-render targets (e.g., re-running MEME logos)
- MultiQC or any other reporting outputs

---

## Risks & Dependencies

- `render_report_html.py` imports from `report.py` using a relative `sys.path.insert`. Snakemake resolves `script:` paths relative to the Snakefile, so the import path must match how `report.py` currently does it (same directory). Low risk — same pattern as `render_meme_logo.py` importing `logo_utils`.
- If `report.csv` doesn't exist when `report_html` runs, Snakemake will refuse to schedule it. This is correct behavior — `report_html` declares `report.csv` as input, so Snakemake will run `report` first.
