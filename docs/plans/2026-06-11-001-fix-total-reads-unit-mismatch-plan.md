# fix: correct total_reads unit mismatch for paired-end samples

## Summary

For paired-end samples, `report.csv` shows `total_reads` smaller than `trimmed_reads` and `mapped_reads`, which is nonsensical (trimming and alignment can only remove reads, never add them). The cause is a unit mismatch introduced in `collect_stats.py`: `total_reads` is halved (treated as read-pair count) while `trimmed_reads` and `mapped_reads` count individual reads (R1+R2 combined). This plan fixes `_parse_subsample_log` so `total_reads` is reported in the same units as the other read-count columns.

---

## Problem Frame

`workflow/scripts/collect_stats.py:_parse_subsample_log` parses a `reformat.sh` subsample log line `Input: N reads` to populate `total_reads`. For PE samples it divides `N` by 2, on the theory that `reformat.sh` counts R1+R2 separately and the column should represent fragment/pair counts.

However:
- `trimmed_reads` comes from `_parse_bbduk_log`, which parses bbduk's `Result: N reads` line. For PE input, bbduk also counts R1+R2 combined and is **not** halved.
- `mapped_reads` comes from `_samtools_count` (`samtools view -c -F 2308`), which counts individual alignment records (R1+R2 combined) and is **not** halved.

Result: for PE samples, `total_reads` is reported in "read pairs" while `trimmed_reads`/`mapped_reads` are reported in "individual reads" (~2x larger). This produces the visibly broken ordering `total_reads < trimmed_reads` and `total_reads < mapped_reads` seen in `results/human/bwa/stats/report.html` (e.g. CXXC01: total_reads=62,835,988, trimmed_reads=107,076,682, mapped_reads=90,704,660 — the unhalved input count is 125,671,976, which correctly sits above trimmed and mapped).

`subsampled_frags` is unaffected: it is compared against `config["bbduk"]["max_frags"]`, which `trim.smk` computes as a fragment/pair count (`TOTAL_LINES / 4` per R1 file), so halving for PE is correct there.

---

## Requirements

- `total_reads` must be expressed in the same units as `trimmed_reads` and `mapped_reads` (individual reads, R1+R2 combined for PE) so the natural ordering `total_reads >= trimmed_reads >= mapped_reads` holds.
- `subsampled_frags` must continue to represent fragment/pair counts (no change to its computation).
- No other column's computation (`mapping_pct`, `frip`, `frip_filt`, etc.) should change, since none of them currently consume `total_reads`.

---

## Key Technical Decisions

### KTD1: Stop halving `total` in `_parse_subsample_log`; keep halving `subsampled`

**Decision:** In `_parse_subsample_log`, the `Input: N reads` line should populate `total` as `str(n)` for both SE and PE (no `// 2` for PE). The `Output: N reads` line continues to populate `subsampled` as `str(n // 2 if is_pe else n)`, unchanged.

**Rationale:** This is the minimal change that aligns `total_reads` with `trimmed_reads`/`mapped_reads` units (verified against real subsample/trim logs: PE `Input: 125671976 reads` should map to `total_reads=125671976`, matching the unhalved trim log `Result: 107076682 reads`). `subsampled_frags` is a distinct concept (gated on `max_frags`, a fragment-count config value) and must stay in fragment units — changing it would require also changing `trim.smk`'s subsampling math, which is out of scope and not broken.

**Alternatives considered:**
- *Halve `trimmed_reads` and `mapped_reads` for PE instead* — rejected. This would require threading `is_pe` into `_parse_bbduk_log` and `_samtools_count`, touches two more code paths, and produces "fragment" semantics for columns explicitly named `*_reads`. The `total_reads` fix is smaller and matches the column-naming convention (`*_reads` = individual reads, `*_frags` = pairs).

---

## Scope Boundaries

**In scope:**
- `workflow/scripts/collect_stats.py` — `_parse_subsample_log` fix and docstring update.
- `tests/test_collect_stats_and_report.py` — new unit tests for `_parse_subsample_log` covering SE and PE total/subsampled behavior.

**Out of scope:**
- Regenerating existing `*.stats.csv` / `report.csv` / `report.html` files under `results/` — these reflect the old (buggy) units and will only be corrected when the relevant Snakemake rules are re-run. Not a code change; left as an operational follow-up for the user.
- Any change to `subsampled_frags`, `trim.smk` subsampling math, or `mapping_pct`/`frip` calculations — none of these consume `total_reads` and none are affected by this bug.

---

## Implementation Units

### U1. Fix `_parse_subsample_log` total_reads unit

**Goal:** Make `total_reads` reflect individual-read counts (matching `trimmed_reads`/`mapped_reads`) rather than fragment/pair counts for PE samples.

**Requirements:** Addresses the unit-mismatch root cause described in Problem Frame.

**Dependencies:** None.

**Files:**
- `workflow/scripts/collect_stats.py` (modify `_parse_subsample_log`, lines ~21-40)

**Approach:**
- In the `Input:` branch, change `total = str(n // 2 if is_pe else n)` to `total = str(n)`.
- Leave the `Output:` branch (`subsampled`) unchanged.
- Update the function's docstring to clarify that `total` is always the raw individual-read count from the `Input:` line, while `subsampled` (fragment count) is halved for PE.

**Patterns to follow:** Existing style in `collect_stats.py` — small pure functions with regex-based log parsing (see `_parse_bbduk_log` for the analogous unhalved pattern).

**Test scenarios:**
- Happy path (PE): given a log with `Input: 125671976 reads` and `Output: 125671976 reads (100.00%)`, `is_pe=True` → `total == "125671976"`, `subsampled == "62835988"`.
- Happy path (SE): given a log with `Input: 3536279 reads` and `Output: 3536279 reads (100.00%)`, `is_pe=False` → `total == "3536279"`, `subsampled == "3536279"` (unchanged from current behavior).
- Subsampling applied (PE): `Input: 125671976 reads` / `Output: 10000000 reads (7.96%)`, `is_pe=True` → `total == "125671976"`, `subsampled == "5000000"`.
- Regression check: `total_reads >= trimmed_reads >= mapped_reads` ordering is restorable for the CXXC01 fixture values (total=125671976 unhalved vs. trimmed=107076682 vs. mapped=90704660) — can be expressed as an assertion in the new test rather than a full pipeline run.

**Verification:** New tests for `_parse_subsample_log` pass; existing tests in `tests/test_collect_stats_and_report.py` continue to pass unchanged (no other column computation depends on `total_reads`).

---

### U2. Add regression tests for `_parse_subsample_log`

**Goal:** Cover `_parse_subsample_log` with direct unit tests (currently untested), using fixture log content derived from real `reformat.sh` output, so the total/subsampled unit relationship is locked in.

**Requirements:** Provides test coverage for KTD1; prevents regression of the unit mismatch.

**Dependencies:** U1 (tests should be written against the corrected implementation; can be written test-first per the execution note below).

**Files:**
- `tests/test_collect_stats_and_report.py` (add new test functions; follow existing `tmp_path`-based fixture-writing helpers such as `_write_stats_csv`)

**Execution note:** Write these tests first (they should fail against the current `// 2` behavior), then apply the U1 fix to make them pass.

**Approach:**
- Add a small helper that writes a minimal `reformat.sh`-style log to a temp file, using the real captured log formats in the Appendix below as ground truth (these are not checked into the repo — they were captured from a completed PE and SE run during planning).
- Add test functions for: PE no-subsampling, SE no-subsampling, PE with subsampling applied (`Output:` < `Input:`).

**Test scenarios:**
- `test_parse_subsample_log_pe_no_subsampling`: PE log with `Input: 125671976 reads`, `Output: 125671976 reads (100.00%)` → `total == "125671976"`, `subsampled == "62835988"`.
- `test_parse_subsample_log_se_no_subsampling`: SE log with `Input: 3536279 reads`, `Output: 3536279 reads (100.00%)` → `total == "3536279"`, `subsampled == "3536279"`.
- `test_parse_subsample_log_pe_with_subsampling`: PE log with `Input: 125671976 reads`, `Output: 10000000 reads (7.96%)` → `total == "125671976"`, `subsampled == "5000000"`.
- `test_parse_subsample_log_total_not_halved_for_pe`: explicit assertion that `total` from a PE log with `Input: 125671976 reads` equals `"125671976"`, not `"62835988"` — directly encodes the bug being fixed.

**Verification:** `pixi run pytest tests/test_collect_stats_and_report.py -k subsample_log` passes; full test suite (`pixi run pytest`) passes.

---

## Open Questions / Follow-Up

- Existing `results/**/stats/*.stats.csv`, `report.csv`, and `report.html` files (across all completed runs, not just `results/human/bwa`) were generated with the old halved `total_reads` for PE samples and will remain inconsistent until the relevant `sample_stats_*` and `report`/`report_html` rules are re-run. This plan does not trigger that re-run — flagging it so the user can decide when/whether to regenerate reports for affected runs.
- `workflow/scripts/update_db.py` persists `total_reads` into the tracking database. Any existing PE rows written before this fix will retain the halved value, the same staleness as the CSV/HTML reports above. No code change is needed (the writer just stores whatever `collect_stats.py` produces), but a database backfill/correction may be worth considering alongside the report regeneration follow-up.

---

## Appendix: Verified Log Formats

Captured during planning from a completed PE run (`results/human/bwa/logs/bbduk/CXXC01.*`) and an existing SE run (`results/a_thaliana_7_v2/logs/bbduk/DEL1.subsample.log`). These are not checked into the repo — use them as the ground-truth format for U2's test fixtures.

**PE subsample log (`reformat.sh`, `int=t`, two input files) — relevant lines:**
```
Set INTERLEAVED to true
Reset INTERLEAVED to false because paired input files were specified.
Input is being processed as paired
Writing interleaved.
Input:                  	125671976 reads          	19856172208 bases
Output:                 	125671976 reads (100.00%) 	19856172208 bases (100.00%)
```

**PE trim log (`bbduk.BBDukS`) — relevant line:**
```
Result:                 	107076682 reads (85.20%) 	14322340110 bases (72.13%)
```

Cross-check against `CXXC01.stats.csv` (current, buggy output): `total_reads=62835988` (= 125671976 / 2, the bug), `trimmed_reads=107076682` (unhalved). After the fix, `total_reads` should be `125671976`.

**SE subsample log (`reformat.sh`, single input file, `int=f`) — relevant lines:**
```
Set INTERLEAVED to false
Input is being processed as unpaired
Input:                  	3536279 reads          	353627900 bases
Output:                 	3536279 reads (100.00%) 	353627900 bases (100.00%)
```

For SE, `total` and `subsampled` are already equal and unaffected by this fix (`is_pe=False` → no `// 2` either way).
