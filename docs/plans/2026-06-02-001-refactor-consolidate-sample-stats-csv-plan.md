---
title: "refactor: Consolidate per-sample statistics into one stats CSV per sample"
date: 2026-06-02
status: completed
type: refactor
---

## Summary

Replace the ~10 individual statistics `.txt` files written per sample with a single `{sample}.stats.csv` that accumulates all tracked metrics. A new `collect_stats.py` script and matching Snakemake rules write these CSVs after the pipeline's main computation steps complete. `report.py` is simplified to aggregate them, and `report.csv` gains four new QC columns (`subsampled_frags`, `median_frag_size`, `num_peaks`, `num_peaks_filt`). `update_db.py` gains a schema migration helper to add those columns to existing databases.

---

## Problem Frame

The pipeline currently writes ~10 individual `.txt` stat files per sample across the `stats/` and `MACS/` directories (`total_frags.txt`, `subsampled_frags.txt`, `align_rate.txt`, `median_frag_size.txt`, `filtered_stats.txt`, `frip_macs.txt`, `frip_macs_filt.txt`, `frip_macs_5fold.txt`, `numpeaks.txt`, `numpeaks_filt.txt`). The storage server bills by file count, so these accumulate cost quickly at scale. All of these values are statistics computed by the pipeline — not raw tool outputs — and belong together in a single structured file per sample.

`idxstats.txt` is raw `samtools idxstats` output and is kept as-is.

---

## Requirements

- **R1** — All per-sample stat `.txt` files are eliminated; each sample produces exactly one `{sample}.stats.csv`
- **R2** — Both treatment and control samples receive a stats CSV
- **R3** — `report.csv` gains `subsampled_frags`, `median_frag_size`, `num_peaks`, `num_peaks_filt` as new columns appended after existing columns
- **R4** — `update_db.py` migrates existing SQLite databases to include the new columns without data loss
- **R5** — `idxstats.txt` is retained unchanged
- **R6** — All existing `report.csv` columns and `update_db.py` pipeline_runs columns remain at the same positions and formats

---

## Key Technical Decisions

**KTD1: One stats CSV per sample, not one global CSV**
A single global CSV requires concurrent writes from parallel Snakemake jobs — a race condition. Per-sample CSVs are Snakemake-native: one rule owns exactly one output file. Aggregation into `report.csv` continues as the final step.

**KTD2: Two Snakemake rules — `sample_stats_treatment` and `sample_stats_control`**
Treatment and control samples have different available inputs (narrowPeak files and FIMO TSVs exist only for treatment). Snakemake's `input:` declarations must be concrete at DAG-build time; two rules with wildcard constraints (`TREATMENT_SAMPLES`, `CONTROL_SAMPLES`) are cleaner than optional inputs with lambdas.

**KTD3: Python script (`collect_stats.py`) runs all samtools/bedtools computations**
Mixing subprocess calls with structured CSV output is cleaner in Python than in a shell rule. Consistent with the existing `report.py` pattern.

**KTD4: `total_reads` and `subsampled_frags` extracted from the existing subsample log**
`reformat.sh` logs "Input: N reads" and "Output: N reads". The trim rules already use `grep "^Output:"` on this log to derive `subsampled_frags.txt`, confirming the format. Adding `grep "^Input:"` gives `total_frags`. For PE, reformat.sh counts individual reads, so Input/2 and Output/2 give fragment-pair counts. This eliminates `total_frags.txt` and `subsampled_frags.txt` without data loss.

**KTD5: `align_rate` extracted from the existing bowtie2 log**
`align_rate.txt` is already just `grep "overall alignment rate" {log.align} | awk …`. The bowtie2 log is kept as a rule log file. Parsing it in `collect_stats.py` eliminates the intermediate file.

**KTD6: `update_db.py` uses `ALTER TABLE ... ADD COLUMN` for schema migration**
`CREATE TABLE IF NOT EXISTS` does not add columns to an existing table. A helper checks `PRAGMA table_info` and runs `ALTER TABLE pipeline_runs ADD COLUMN "col" TEXT` for any column in `COLS` that is missing. This is idempotent — safe against both new and previously-run databases.

---

## Scope Boundaries

**In scope:**
- All per-sample stat `.txt` files produced by trim, align, peaks, and stats rules
- Treatment and control samples
- New columns in `report.csv` and matching DB schema migration in `update_db.py`

**Deferred to Follow-Up Work:**
- Including control sample stats in `report.csv` (controls get their own stats CSV but are excluded from the report by design)
- Log file count reduction (bbduk, bowtie2, bamcoverage logs — a separate concern)

**Out of scope:**
- Any changes to trimming, alignment, or peak-calling logic
- Changes to `idxstats.txt`
- Changing or reordering the 11 existing `report.csv` columns

---

## High-Level Technical Design

```
trim rules          align rules         peaks/stats/motifs rules
     │                   │                        │
     │ (log kept)        │ (log kept)             │ (BAM, peaks, FIMO kept)
     └──────────────┬────┘                        │
                    │                             │
           sample_stats_treatment ◄───────────────┘
           sample_stats_control ◄── (BAM + logs only)
                    │
            {sample}.stats.csv
                    │
              report rule
                    │
             report.csv  (15 cols)
                    │
            update_db rule
                    │
          pipeline_db.db (29 cols)
```

The four existing stats rules (`filtered_stats`, `frip_macs`, `frip_macs_filt`, `frip_macs_5fold`) are replaced by the two `sample_stats_*` rules. The `filtered_stats` rule is split: its `idxstats.txt` output moves to a standalone `idxstats` rule; the `filtered_stats.txt` output is absorbed into the stats CSV.

---

## Implementation Units

### U1. Remove stat file outputs from trim rules

**Goal:** Stop writing `total_frags.txt` and `subsampled_frags.txt`; the subsample log already contains equivalent data.

**Requirements:** R1, R2

**Dependencies:** none

**Files:**
- `workflow/rules/trim.smk`

**Approach:**
- In both `trim_se` and `trim_pe`: remove `total_frags` and `subsampled_frags` from the `output:` block
- Remove the corresponding shell lines: `echo "$TOTAL_FRAGS" > {output.total_frags}` and the `grep "^Output:"` pipeline that wrote subsampled_frags
- Keep the `$TOTAL_FRAGS` shell variable — it is still required for the subsampling fraction computation within the rule
- The subsample log (`{log.subsample}`) is already retained and will be read by `collect_stats.py`

**Test scenarios:**
- `trim_se` with an SE sample: rule completes without error; `{sample}.total_frags.txt` is not created; `{log.subsample}` exists and contains "Input:" and "Output:" lines
- `trim_pe` with a PE sample: same; `{sample}.subsampled_frags.txt` is not created
- Subsampling still produces the correct fragment count: trimmed output read count matches expected when `max_frags` is set (subsampling logic uses the `$TOTAL_FRAGS` variable, not the file)

**Verification:** `snakemake --dry-run` succeeds with no DAG errors; no references to `total_frags.txt` or `subsampled_frags.txt` remain in trim.smk outputs

---

### U2. Remove stat file outputs from align rules

**Goal:** Stop writing `align_rate.txt` (SE and PE) and `median_frag_size.txt` (PE only); the bowtie2 log already contains the alignment rate.

**Requirements:** R1, R2

**Dependencies:** none

**Files:**
- `workflow/rules/align.smk`

**Approach:**
- In `align_se`: remove `align_rate` from `output:` block; remove the `grep "overall alignment rate"` shell command
- In `align_pe`: same as SE, plus remove `median_frag_size` from `output:` and remove the `bamPEFragmentSize` command
- The bowtie2 log (`{log.align}`) is already kept as a rule log file and will be parsed by `collect_stats.py`
- BAM, BAI, and BigWig outputs are unchanged

**Test scenarios:**
- `align_se`: rule completes; `{sample}.align_rate.txt` is not created; `{log.align}` contains "overall alignment rate"
- `align_pe`: rule completes; `{sample}.align_rate.txt` and `{sample}.median_frag_size.txt` are not created
- BAM, BAI, and BigWig are still produced correctly for both SE and PE

**Verification:** `snakemake --dry-run` succeeds; align.smk `output:` blocks no longer reference `align_rate` or `median_frag_size`

---

### U3. Remove numpeaks stat file outputs from filter_peaks

**Goal:** Stop writing `numpeaks.txt` and `numpeaks_filt.txt`; these are simple line counts of files that remain in the DAG.

**Requirements:** R1

**Dependencies:** none

**Files:**
- `workflow/rules/peaks.smk`

**Approach:**
- Remove `numpeaks` and `numpeaks_filt` from the `filter_peaks` `output:` block
- Remove the two `wc -l` commands that wrote those files
- The narrowPeak files (`_peaks.narrowPeak`, `_peaks_filt.narrowPeak`, `_peaks_5fold.narrowPeak`) remain as outputs and are consumed by downstream rules

**Test scenarios:**
- `filter_peaks`: rule completes; `{sample}_numpeaks.txt` and `{sample}_numpeaks_filt.txt` are not created
- `{sample}_peaks_filt.narrowPeak` and `{sample}_peaks_5fold.narrowPeak` are still produced with correct content

**Verification:** `snakemake --dry-run` succeeds; filter_peaks output block references only peak files

---

### U4. Create collect_stats.py and sample_stats rules

**Goal:** New Python script + two Snakemake rules that compute all per-sample statistics and write `{sample}.stats.csv`.

**Requirements:** R1, R2, R3

**Dependencies:** U1, U2, U3

**Files:**
- `workflow/scripts/collect_stats.py` (new)
- `workflow/rules/stats.smk`

**Approach — collect_stats.py:**

Receives all paths via `snakemake.input`, `snakemake.output`, and `snakemake.params.is_pe` / `snakemake.params.is_treatment`.

Log parsing (no subprocesses needed):
- `total_reads`: `grep "^Input:" subsample_log`, value field; divide by 2 if `is_pe`
- `subsampled_frags`: `grep "^Output:" subsample_log`, value field; divide by 2 if `is_pe`
- `trimmed_reads`: existing `parse_bbduk_log` regex — `"Result:\s+(\d+)\s+reads"` from trim log
- `alignment_rate`: `grep "overall alignment rate" align_log`, strip `% overall...` suffix

Subprocess calls (samtools / bedtools — treatment and control):
- `mapped_reads`: `samtools view -c -F 2308 {bam}`

Subprocess calls (treatment only):
- `reads_in_peaks`: `bedtools intersect -abam {bam} -b {peaks} -u | samtools view -c`
- `reads_5fold`: same with `{peaks_5fold}`
- `reads_nfold`: same with `{peaks_filt}`

Subprocess call (PE treatment only):
- `median_frag_size`: `bamPEFragmentSize -b {bam}`, parse first "Median:" line

File reads (treatment only):
- `num_peaks`: line count of raw narrowPeak
- `num_peaks_filt`: line count of filtered narrowPeak
- `max_peak_score`: max value of column 7 (0-indexed) of filtered narrowPeak — same logic as current `parse_narrowpeak_max`
- `motif_peaks`: unique `sequence_name` values in FIMO TSV — same logic as current `parse_fimo`

Output: two-row CSV (header + data row) at `snakemake.output.stats_csv`. Non-applicable columns (frip fields for controls, `median_frag_size` for SE) are written as `"NA"`.

**Approach — Snakemake rules (added to stats.smk):**

`sample_stats_treatment` — wildcard constrained to `TREATMENT_SAMPLES`:
```
input:  bam, bai, peaks (raw), peaks_filt, peaks_5fold,
        fimo (peaks/fimo.tsv), subsample_log, trim_log, align_log
output: stats_csv = OUT + "/stats/{sample}.stats.csv"
params: is_pe = lambda wc: wc.sample in PE_SAMPLES, is_treatment = True
script: collect_stats.py
```

`sample_stats_control` — wildcard constrained to `CONTROL_SAMPLES` (if empty, no jobs fire):
```
input:  bam, bai, subsample_log, trim_log, align_log
output: stats_csv = OUT + "/stats/{sample}.stats.csv"
params: is_pe = lambda wc: wc.sample in PE_SAMPLES, is_treatment = False
script: collect_stats.py
```

**Test scenarios:**
- Treatment SE sample: `{sample}.stats.csv` produced; all columns present; `median_frag_size` is `"NA"`; frip and peaks columns are numeric
- Treatment PE sample: `{sample}.stats.csv` produced; `median_frag_size` is an integer; all other columns numeric
- Control sample: `{sample}.stats.csv` produced; `reads_in_peaks`, `reads_5fold`, `reads_nfold`, `num_peaks`, `num_peaks_filt`, `max_peak_score`, `motif_peaks`, `median_frag_size` are all `"NA"`; `total_reads`, `mapped_reads` are numeric
- Sample with empty FIMO TSV (no motif hits): `motif_peaks` is `"NA"`, not `0`
- No subsampling (`max_frags = None`): `subsampled_frags` equals `total_reads`; subsampling fraction = 1.0 does not cause divide-by-zero
- Pipeline without control samples: `sample_stats_control` rule fires for zero wildcards; DAG still resolves cleanly

**Verification:** `{sample}.stats.csv` exists for every sample; CSV header row matches the expected column list; data row values match what the removed rules would have produced

---

### U5. Refactor old stats rules in stats.smk

**Goal:** Replace `filtered_stats`, `frip_macs`, `frip_macs_filt`, `frip_macs_5fold` with a standalone `idxstats` rule; the other outputs are now in the stats CSV.

**Requirements:** R1, R5

**Dependencies:** U4

**Files:**
- `workflow/rules/stats.smk`
- `config/config.yaml`

**Approach:**
- The current `filtered_stats` rule produces two outputs: `{sample}.filtered_stats.txt` and `{sample}.idxstats.txt`. `idxstats.txt` is kept (R5); `filtered_stats.txt` is moved to the stats CSV (U4).
- Rename `filtered_stats` → `idxstats`: keep only `samtools idxstats {input.bam} > {output.idxstats}`; remove the `samtools view -c` and `printf` lines; remove `filtered_stats.txt` from output
- Delete `frip_macs`, `frip_macs_filt`, `frip_macs_5fold` rules in their entirety
- Remove the `filtered_stats` and `frip_macs` resource entries from `config/config.yaml`'s `resources:` block (these keys become unreferenced once the rules are deleted)

**Test scenarios:**
- `{sample}.idxstats.txt` still produced for all samples with correct per-chromosome counts
- `{sample}.filtered_stats.txt` is not produced by any rule
- `{sample}.frip_macs.txt`, `{sample}.frip_macs_filt.txt`, `{sample}.frip_macs_5fold.txt` are not produced

**Verification:** `snakemake --dry-run` succeeds; `idxstats.txt` appears in the DAG; no frip txt files in expected outputs

---

### U6. Update report rule and report.py

**Goal:** Simplify `report.py` to read from per-sample stats CSVs; add four new columns to `report.csv`; update the Snakemake `report` rule inputs accordingly.

**Requirements:** R1, R3, R6

**Dependencies:** U4, U5

**Files:**
- `workflow/rules/stats.smk` (report rule)
- `workflow/scripts/report.py`

**Approach — report rule in stats.smk:**
- Replace the long `input:` list with: `stats_csvs = expand(OUT + "/stats/{sample}.stats.csv", sample=TREATMENT_SAMPLES)`
- Remove: `trim_logs`, `total_frags`, `filt_stats`, `frip`, `frip_filt`, `frip_5fold`, `align_rates`, `narrowpeaks`, `fimo_peaks` (all now consumed by collect_stats.py)
- Simplify `params:` to only `treatment_samples` and `stats_dir`

**Approach — report.py:**
- `build_row(sample, stats_dir)`: opens `{stats_dir}/{sample}.stats.csv` as a `csv.DictReader`, reads the single data row, computes `mapping_pct` from the row's `reads_in_peaks` and `mapped_reads` values (not stored in stats CSV to avoid duplication), returns a dict
- `make_cols()`: add `subsampled_frags`, `median_frag_size`, `num_peaks`, `num_peaks_filt` at the end — after the existing 11 columns
- Remove all standalone parsing functions (`parse_kv_file`, `parse_bbduk_log`, `parse_align_rate`, `parse_narrowpeak_max`, `parse_fimo`) — their logic now lives in `collect_stats.py`
- `run()` and the Snakemake entrypoint: update params to use `stats_dir` only
- `cli_main()`: discover samples from `*.stats.csv` glob instead of `*.frip_macs.txt`

**New report.csv column order:**
```
sample, total_reads, trimmed_reads, mapped_reads, reads_in_peaks,
reads_5fold, reads_nfold, alignment_rate, mapping_pct,
max_peak_score, motif_peaks,
subsampled_frags, median_frag_size, num_peaks, num_peaks_filt
```

**Test scenarios:**
- `report.csv` produced with 15 columns; header matches expected order
- PE sample row: `median_frag_size` is numeric; SE sample row: `median_frag_size` is `"NA"`
- `mapping_pct` computed correctly from values in stats CSV (not carried through from stats CSV directly)
- `cli_main()` discovers sample names from `*.stats.csv` files and produces report without error
- Missing stats CSV for a sample: report row for that sample has `"NA"` for all columns (graceful degradation)

**Verification:** `report.csv` schema matches expected 15-column header; values in previously-existing columns are numerically unchanged from the pre-refactor pipeline

---

### U7. Update update_db.py for new report.csv columns

**Goal:** Add four new columns to `pipeline_runs` and handle schema migration on existing databases.

**Requirements:** R4

**Dependencies:** U6

**Files:**
- `workflow/scripts/update_db.py`

**Approach:**
- Add `"subsampled_frags"`, `"median_frag_size"`, `"num_peaks"`, `"num_peaks_filt"` to the `COLS` list (append at end to preserve existing column positions)
- Add corresponding `row["col"] = stats.get("col", "NA")` lines in `main()`
- Add a `_ensure_columns(con, table, cols)` helper that runs `PRAGMA table_info({table})`, finds any columns in `cols` missing from the schema, and runs `ALTER TABLE {table} ADD COLUMN "{col}" TEXT` for each. Call this in `write_rows()` after `CREATE TABLE IF NOT EXISTS`, before the `DELETE`/`INSERT`.
- `pipeline_runs` table goes from 25 to 29 QC metric columns (total row width unchanged: 29 QC cols + `id` auto-increment)

**Test scenarios:**
- Fresh DB: `pipeline_runs` table created with all 29 QC columns; INSERT succeeds
- Existing DB with 25 columns: `_ensure_columns` adds the 4 missing columns without error; existing rows retain their data; new rows insert successfully
- Re-running `_ensure_columns` on a DB that already has all 29 columns: no error (columns already present, ALTER skipped)
- Control sample rows in DB: `num_peaks`, `reads_in_peaks`, etc. are `"NA"` (controls are not in report.csv; `stats.get` defaults to `"NA"`)

**Verification:** DB schema has 29 QC columns + `id`; new pipeline run inserts correctly; existing run data is intact after migration

---

## Open Questions / Implementation Notes

- **Subsample log "Input:" format**: The trim rules use `grep "^Output:"` on the reformat.sh log, confirming that format. The parallel `"^Input:"` line is assumed to have the same structure (`Input:\t<N> reads`). Verify against an actual reformat.sh log before finalizing the parsing regex in `collect_stats.py`. If the format differs, only the parsing logic in U4 needs adjustment.
- **CONTROL_SAMPLES could be empty**: When no control samples are configured, `sample_stats_control` fires for zero wildcards. Snakemake handles this correctly — no action needed.
- **`mapping_pct` stays in report.py**: It is derived from other stats CSV columns (`reads_in_peaks / mapped_reads`) and should not be stored in the stats CSV to avoid duplication. `report.py` continues to compute it.

---

## System-Wide Impact

- **File count**: Treatment samples: ~11 stat files → `{sample}.stats.csv` + `idxstats.txt` = 2 files. Control samples: ~5 stat files → 2 files.
- **`report.csv`**: Gains 4 columns appended at the end. Consumers using positional column access will break — they should use column-name access (DictReader). `update_db.py` already uses `stats.get(col)` so it is safe.
- **`pipeline_db.db`**: Schema migrated in-place when `update_db.py` runs. One-way change (no downgrade path). Existing rows gain four `NULL`/empty columns for the new metrics.
- **No behavioral changes**: Trimming, alignment, peak calling, and MEME/FIMO logic are untouched. Only the stat file outputs and their consumers change.
