# DapSeq Pipeline Audit

## Codebase Changes Since Last Audit

These architectural changes affect line-number references throughout this document:

1. **`update_db.py` completely rewritten** — now writes two tables (`pipeline_runs`, `run_metadata`), reads QC stats from `report.csv` via `DictReader` instead of raw stat files, and has an expanded column set (`reads_nfold_*`, `mapping_pct_*`, `motif_peaks_*`, etc.).
2. **`report.py` completely rewritten** — now uses named helper functions (`parse_kv_file`, `parse_bbduk_log`, `parse_narrowpeak_max`, `parse_fimo`, `_safe_frip`, `build_row`), handles both `full` and `chr` bam types, and has a standalone CLI mode.
3. **`filter_chr` rule relocated** — lives in `align.smk` (lines 124–163), not in `filter.smk`.
4. **Logomaker replaces ghostscript** — `meme_summits` and `factorbook_logo` render PNG logos via `logo_utils._render_logo_with_logomaker`; no more ghostscript PDF step.
5. **`{bam_type}` wildcard added throughout** — peaks, motifs, stats, and report rules now use a `full|chr` wildcard, doubling rule instances per sample.

---

## Critical Bugs

### 1. `align.smk:150-153` — `samtools -e` expression filter silently broken in container
**File:** `workflow/rules/align.smk:150-153` (rule `filter_chr`)

The shell command uses `samtools view -e 'rname =~ "{params.pattern}"'`. The `dapseq.sif` samtools silently fails with `-e 'rname =~ ...'` — reads are not filtered and an unfiltered BAM is produced with no error.

**Fix:** Replace the `-e` filter with header-parsing + explicit region args:
```bash
REGIONS=$(samtools view -H {input.bam} \
  | awk '/^@SQ/{{split($2,a,":");if(a[2]~"{params.pattern}")printf "%s ", a[2]}}')
samtools view -b -@ {threads} -o {output.bam} {input.bam} $REGIONS 2>{log.filter}
```

---

## High Priority Bugs

### 2. `report.py:48` — `fields[6]` accessed without bounds check in `parse_narrowpeak_max`
**File:** `workflow/scripts/report.py:47-50`

If any narrowPeak line is malformed (fewer than 7 columns), `float(fields[6])` raises an `IndexError` with no informative message.

**Fix:** Add `if len(fields) < 10: continue` (narrowPeak requires 10 columns).

---

### 3. `report.py:112` — Unclosed file handle
**File:** `workflow/scripts/report.py:112`

`open(total_frags_path).read().strip()` in a conditional expression — no `with` statement, no explicit close. All other reads in this file correctly use `with open(...)`.

**Fix:** Replace with:
```python
with open(total_frags_path) as fh:
    row["reads_original"] = fh.read().strip()
```

---

## Medium Priority — Anti-patterns & Snakemake Guidelines

### 4. `ref.smk:25` and `filter.smk:25` — Append (`2>>`) instead of overwrite (`2>`) for log files
**Files:** `workflow/rules/ref.smk:25`, `workflow/rules/filter.smk:25`

`samtools faidx` (ref.smk) and `samtools index` inside `merge_control` (filter.smk) both use `2>>` (append). On reruns these logs grow indefinitely and mix output from different runs. All other rules correctly use `2>`.

**Fix:** Change both to `2>{log}`.

---

### 5. `motifs.smk` — `meme_summits` and `factorbook_logo` use `run:` with embedded `shell()` calls
**Files:** `workflow/rules/motifs.smk:159-199` (`meme_summits`), `workflow/rules/motifs.smk:18-89` (`factorbook_logo`)

Both rules use a `run:` block that mixes Python logic with `shell()` calls. Snakemake guidelines recommend `run:` only for pure Python logic; compute-heavy rules should use `shell:` or `script:`. A `run:` block cannot be benchmarked and does not integrate cleanly with the SLURM executor's resource tracking. The companion rule `meme_peaks` correctly uses `shell:`.

**Fix:** Extract each rule's logic into a script in `workflow/scripts/` (`run_meme.py`, `run_factorbook_logo.py`) and replace each `run:` block with `script: "../scripts/<name>.py"`.

---

### 6. `align.smk`, `peaks.smk`, `filter.smk`, `motifs.smk` — `{config[...]}` direct access in shell blocks
**Files:**
- `align.smk:51,116`: `{config[genome_size]}` (in `align_se` and `align_pe` bamCoverage calls)
- `align.smk:160`: `{config[genome_size]}` (in `filter_chr` bamCoverage call)
- `peaks.smk:28-29`: `{config[macs3][format]}`, `{config[genome_size]}`
- `filter.smk:53-56`: `{config[bamcompare][bin_size]}`, `{config[bamcompare][operation]}`, `{config[bamcompare][scale_factors_method]}`, `{config[bamcompare][n]}`
- `motifs.smk:260,288`: `{config[fimo][thresh]}`

Snakemake guidelines require config values used in `shell:` blocks to be declared in `params:` — this makes them visible in the DAG, validates types at parse time, and avoids silent failures from missing keys.

**Fix:** Move each `{config[...]}` reference used in `shell:` into a corresponding `params:` entry and reference `{params.xxx}` in the shell command.

---

### 7. No `set -euo pipefail` in multi-command shell blocks
**Files:** `workflow/rules/peaks.smk:22-31`, `workflow/rules/filter.smk:22-26`, `workflow/rules/ref.smk:22-27`, and others.

If any command in a chained shell block fails, subsequent commands run with incomplete/corrupt data and the rule appears to succeed. Snakemake guidelines recommend `set -euo pipefail` at the top of all non-trivial shell blocks.

**Fix:** Add `set -euo pipefail` as the first line of every shell block that contains more than one command. (The `exec 2>{log}` pattern already used in `stats.smk` is a good model — combine with `set -euo pipefail`.)

---

### 8. `update_db.py:184-188` — Stored BAM/bigWig paths don't match actual output paths
**File:** `workflow/scripts/update_db.py:184-188`

`_build_meta_paths` stores `bam: f"{o}/bam/{sample}.bam"` and `bigwig: f"{o}/bigWig/{sample}.bw"`. Actual outputs are `{sample}.full.bam` / `{sample}.chr.bam` and `{sample}.full.bw` / `{sample}.chr.bw`. The DB will never contain a valid path to either file.

**Fix:** Store both paths as `bam_full` / `bam_chr` (and `bigwig_full` / `bigwig_chr`), updating `META_COLS` accordingly, or pick the primary (`full`) and fix the path strings to include the bam_type suffix.

---

### 9. `ref.smk:1-27` — `bowtie2_index` has no `threads:` directive
**File:** `workflow/rules/ref.smk`

`bowtie2-build` supports `--threads`. Without a `threads:` directive the rule runs single-threaded and the SLURM scheduler allocates only 1 core, wasting wall-time on large genomes.

**Fix:** Add `threads: config["threads"]` and pass `--threads {threads}` to the `bowtie2-build` call.

---

### 10. Duplicate rules that should be unified
**File:** `workflow/rules/motifs.smk`

- `narrow_peak_to_fasta_summits` and `narrow_peak_to_fasta_peaks` differ only in `params.extend_bp` — they call the same script.
- `fimo_summits` and `fimo_peaks` are identical except for the `summits` vs `peaks` path segment.

Both pairs can be collapsed into a single rule with a `{peak_type}` wildcard (`summits|peaks`), eliminating duplication and simplifying future maintenance.

**Fix:** Introduce `{peak_type}` wildcard constrained to `"summits|peaks"` and merge each pair into one rule.

---

## Low Priority

### 11. Missing `benchmark:` directives on compute-heavy rules
**Files:** `align_se`, `align_pe`, `meme_summits`, `meme_peaks`, `bowtie2_index`

Without `benchmark:`, there is no per-rule runtime/memory profiling. Snakemake guidelines recommend benchmarks on any rule expected to run >1 min.

**Fix:** Add `benchmark: OUT + "/benchmarks/{rule_name}/{sample}.{bam_type}.tsv"` to the five heaviest rules.

---

### 12. `motifs.smk:45,89,164,199` — Bare `open(...).close()` instead of `Path.touch()`
**File:** `workflow/rules/motifs.smk:45,89,164,199`

`open(str(f), "w").close()` for creating empty sentinel files is non-idiomatic and appears four times (`factorbook_logo` ×2, `meme_summits` ×2).

**Fix:** Replace with `from pathlib import Path; Path(str(f)).touch()`.

---

## Verification

After applying fixes:

1. **Dry-run:** `pixi run snakemake -n --rerun-triggers mtime` — confirms DAG is valid.
2. **filter_chr fix:** `samtools view -H <output.bam> | grep @SQ` — should only contain chromosomes matching the pattern.
3. **DB path fix:** `sqlite3 pipeline_db.db 'SELECT sample, bam_full, bigwig_full FROM run_metadata'` — should contain real paths including `.full.bam` / `.full.bw` suffixes.
4. **Report script fix:** Mock a narrowPeak with a truncated line; confirm `report.py` does not crash.
5. **Regression:** `pixi run pytest tests/` — existing update_db tests should still pass.
