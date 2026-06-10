---
title: "fix: mapping_pct always returns NA"
type: fix
status: active
created: 2026-06-05
---

# fix: mapping\_pct always returns NA

## Summary

`mapping_pct` is never computed in `collect_stats.py` — the column is absent from `COLUMNS` and no calculation exists. Every downstream consumer reads a key that was never written, so `update_db.py` always stores `"NA"`. The fix adds the computation to `collect_stats.py` and adds the column to `report.py` so the value flows through the full pipeline.

---

## Problem Frame

**Expected:** `mapping_pct` = `mapped_reads` (post-aligner, post-MAPQ filter) / `trimmed_reads` (post-bbduk) × 100, rounded to 2 decimal places.

**Actual:** `mapping_pct` is always `"NA"` in the SQLite database.

**Data flow:**
```
collect_stats.py → {sample}.stats.csv
                → report.py → report.csv
                → update_db.py → SQLite
```

`update_db.py` reads from `report.csv` (not directly from `{sample}.stats.csv`). `report.py` writes only the columns in `make_cols()`. Because neither `collect_stats.py` nor `report.py` define `mapping_pct`, the value is never present at any stage.

---

## Requirements

- R1: `mapping_pct` is computed as `mapped_reads / trimmed_reads * 100` (rounded to 2 dp).
- R2: `"NA"` is returned when either input is `"NA"` or the denominator is zero.
- R3: The computed value flows through `report.csv` so `update_db.py` receives it.
- R4: `mapping_pct` is visible in the HTML report alongside `mapping_rate`.

---

## Key Technical Decisions

**KTD-1: Compute in `collect_stats.py`, not `report.py`.**
`collect_stats.py` already holds both `mapped_reads` and `trimmed_reads` at the point of computation. Computing there keeps `{sample}.stats.csv` as the single source of truth for all per-sample metrics and avoids re-deriving values in a downstream aggregator.

**KTD-2: Use a `_safe_pct` helper consistent with `_safe_frip` in `report.py`.**
`report.py` already has `_safe_frip` for guarded division. A parallel private `_safe_pct` in `collect_stats.py` (or a brief inline guard) follows the same defensive pattern the codebase already uses.

**KTD-3: `mapping_pct` is a percentage column in the HTML report.**
`PCT_COLS` in `report.py` controls `%`-suffix formatting. Adding `mapping_pct` there is consistent with `frip` / `frip_top_n_fold` formatting.

---

## Scope Boundaries

**In scope:** Add `mapping_pct` computation; expose it through `report.py` and the HTML report.

**Out of scope:**
- `alignment_rate` / `mapping_rate` key-name mismatch between `collect_stats.py` and `update_db.py` — separate, unrelated bug.
- Backfilling existing SQLite rows — only future pipeline runs will have the value.

### Deferred to Follow-Up Work
- Fix the `alignment_rate` vs `mapping_rate` column-name mismatch in `update_db.py` (also always `"NA"`).

---

## Implementation Units

### U1. Compute `mapping_pct` in `collect_stats.py`

**Goal:** Add `mapping_pct` to the per-sample stats CSV.

**Requirements:** R1, R2

**Dependencies:** none

**Files:**
- `workflow/scripts/collect_stats.py`

**Approach:**
- Add `"mapping_pct"` to the `COLUMNS` list (after `mapping_rate`).
- Add a `_safe_pct(numerator, denominator)` helper that mirrors `_safe_frip` in `report.py`: return `"NA"` if either arg is `"NA"`, guard against zero denominator, otherwise return `str(round(int(numerator) / int(denominator) * 100, 2))`.
- After `row["mapped_reads"]` and `row["trimmed_reads"]` are both set, assign `row["mapping_pct"] = _safe_pct(row["mapped_reads"], row["trimmed_reads"])`.

**Patterns to follow:** `_parse_bbduk_log` / `_samtools_count` pattern; `_safe_frip` in `report.py` for the guarded division shape.

**Test scenarios:**
- Happy path: both `mapped_reads` and `trimmed_reads` are valid integers → returns correctly rounded percentage string.
- `trimmed_reads` is `"NA"` → returns `"NA"`.
- `mapped_reads` is `"NA"` → returns `"NA"`.
- `trimmed_reads` is `"0"` → returns `"NA"` (zero-division guard).
- `mapped_reads` equals `trimmed_reads` → returns `"100.0"`.
- Result written to `{sample}.stats.csv` as the `mapping_pct` column.

**Verification:** A sample stats CSV contains a numeric `mapping_pct` value (not `"NA"`) when both bbduk and aligner logs are present.

---

### U2. Surface `mapping_pct` in `report.py`

**Goal:** Include `mapping_pct` in `report.csv` and format it correctly in the HTML report.

**Requirements:** R3, R4

**Dependencies:** U1 (value must exist in `{sample}.stats.csv`)

**Files:**
- `workflow/scripts/report.py`

**Approach:**
- Add `"mapping_pct"` to the list returned by `make_cols()`, after `"mapping_rate"` (keeps alignment metrics grouped).
- Add `"mapping_pct"` to `PCT_COLS` so `_fmt_html` appends `%` in the HTML table.
- No other changes needed — `build_row` passes all keys from the stats CSV through to the row dict, and `run_csv` uses `extrasaction="ignore"` against `make_cols()`, so it will automatically include the new column.

**Patterns to follow:** How `frip` and `frip_top_n_fold` are added to `make_cols()` and `PCT_COLS`.

**Test scenarios:**
- `report.csv` header contains `mapping_pct` column.
- A treatment sample row has a numeric `mapping_pct` value when the underlying stats CSV has one.
- `mapping_pct` missing from an individual stats CSV → row contains `"NA"` for that column.
- HTML report renders `mapping_pct` with a `%` suffix.

**Verification:** `report.csv` written by `run_csv` includes `mapping_pct` with a non-`"NA"` value; HTML report displays it with `%` formatting.

---

## Open Questions

None — root cause is confirmed, fix is fully determined.

---

## Sources & Research

- Traced data flow through `collect_stats.py`, `report.py`, and `update_db.py`.
- `_safe_frip` in `report.py:50–57` is the reference pattern for guarded percentage division.
- `PCT_COLS` at `report.py:22` is the reference for `%`-suffix formatting.
