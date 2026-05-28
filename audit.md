o# DapSeq Pipeline Audit

## Critical Bugs

### 1. `update_db.py:262` — DB column always written as "NA" (field name mismatch)
**File:** `workflow/scripts/update_db.py:262`

`report.py` writes a TSV column named `f"peaks_{n}fold#"` (e.g., `"peaks_2fold#"` when `min_foldch=2`). `update_db.py` looks it up with the hardcoded key `"min5fold peak#"` — wrong name and hardcoded value. The DB column `min5fold_peak#` will always be `"NA"`.

**Fix:** Derive `n` from `sm.params.macs3_min_foldch` in `main()` (same formula as `_fmt_n` in report.py: `str(int(x)) if x == int(x) else str(x)`), then use `stats.get(f"peaks_{n}fold#", "NA")`.

---

### 2. `align.smk:150-153` — `samtools -e` expression filter silently broken in container
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

### 3. `report.py:138,145` — ZeroDivisionError when `filtered_reads == "0"`
**File:** `workflow/scripts/report.py:137-146`

Both FRiP calculations check `!= "NA"` but not `> 0`. A sample with zero filtered reads crashes the entire report rule.

**Fix:** Change both conditions to also guard `int(filtered_reads) > 0`.

---

### 4. `report.py:116,122` — Unclosed file handles
**File:** `workflow/scripts/report.py:116,122`

`open(path).read().strip()` used in conditional expressions — no `with` statement, no explicit close. File handles leak.

**Fix:** Wrap each in a `with open(path) as fh: ...` block or a small helper function.

---

### 5. `report.py:54` — `fields[6]` accessed without bounds check in `parse_narrowpeak`
**File:** `workflow/scripts/report.py:54`

If any narrowPeak line is malformed (fewer than 7 columns), an `IndexError` crashes the script with no informative message.

**Fix:** Add `if len(fields) < 10: continue` (narrowPeak requires 10 columns).

---

## Medium Priority — Anti-patterns & Snakemake Guidelines

### 6. `ref.smk:25` and `filter.smk:25` — Append (`2>>`) instead of overwrite (`2>`) for log files
**Files:** `workflow/rules/ref.smk:25`, `workflow/rules/filter.smk:25`

`samtools faidx` and `samtools index` both use `2>>` (append). On reruns these logs grow indefinitely and mix output from different runs. All other rules correctly use `2>`.

**Fix:** Change both to `2>{log}`.

---

### 7. `motifs.smk:68-106` — `meme_summits` uses `run:` with embedded `shell()` calls
**File:** `workflow/rules/motifs.smk:68-106`

`meme_summits` uses a `run:` block that mixes Python logic with `shell()` calls. Snakemake guidelines recommend `run:` only for pure Python logic; compute-heavy rules should use `shell:` or `script:`. A `run:` block also cannot be benchmarked and does not integrate cleanly with the SLURM executor's resource tracking. The companion rule `meme_peaks` correctly uses `shell:`.

**Fix:** Extract the custom-alphabet + meme + ghostscript logic into `workflow/scripts/run_meme.py` and replace the `run:` block with `script: "../scripts/run_meme.py"`.

---

### 8. `align.smk`, `peaks.smk`, `filter.smk`, `motifs.smk` — `{config[...]}` direct access in shell blocks
**Files:**
- `align.smk:51`: `{config[genome_size]}`
- `peaks.smk:28-29`: `{config[macs3][format]}`, `{config[genome_size]}`
- `filter.smk:54-56`: `{config[bamcompare][bin_size]}`, `{config[bamcompare][operation]}`, etc.
- `motifs.smk:167,195`: `{config[fimo][thresh]}`

Snakemake guidelines require config values used in `shell:` blocks to be declared in `params:` — this makes them visible in the DAG, validates types at parse time, and avoids silent failures from missing keys.

**Fix:** Move each `{config[...]}` reference used in `shell:` into a corresponding `params:` entry and reference `{params.xxx}` in the shell command.

---

### 9. No `set -euo pipefail` in multi-command shell blocks
**Files:** `workflow/rules/peaks.smk:22-31`, `workflow/rules/filter.smk:22-26`, `workflow/rules/ref.smk:22-27`, and others.

If any command in a chained shell block fails, subsequent commands run with incomplete/corrupt data and the rule appears to succeed. Snakemake guidelines recommend `set -euo pipefail` at the top of all non-trivial shell blocks.

**Fix:** Add `set -euo pipefail` as the first line of every shell block that contains more than one command. (The `exec 2>{log}` pattern already used in `stats.smk` is a good model — combine with `set -euo pipefail`.)

---

### 10. `update_db.py:181-183` — Stored BAM/bigWig paths don't match actual output paths
**File:** `workflow/scripts/update_db.py:181-183`

`_build_meta_paths` stores `bam: f"{o}/bam/{sample}.bam"` and `bigwig: f"{o}/bigWig/{sample}.bw"`. Actual outputs are `{sample}.full.bam` / `{sample}.chr.bam` and `{sample}.full.bw` / `{sample}.chr.bw`. The DB will never contain a valid path to either file.

**Fix:** Store both paths as `bam_full` / `bam_chr` (and similarly for bigwig), or pick the primary (`full`) and fix the path strings to include the bam_type suffix.

---

### 11. `ref.smk:1-27` — `bowtie2_index` has no `threads:` directive
**File:** `workflow/rules/ref.smk`

`bowtie2-build` supports `--threads`. Without a `threads:` directive the rule runs single-threaded and the SLURM scheduler allocates only 1 core, wasting wall-time on large genomes.

**Fix:** Add `threads: config["threads"]` and pass `--threads {threads}` to the `bowtie2-build` call.

---

### 12. Duplicate rules that should be unified
**File:** `workflow/rules/motifs.smk`

- `narrow_peak_to_fasta_summits` and `narrow_peak_to_fasta_peaks` differ only in `params.extend_bp` — they call the same script.
- `fimo_summits` and `fimo_peaks` are identical except for the `summits` vs `peaks` path segment.

Both pairs can be collapsed into a single rule with a `{peak_type}` wildcard (`summits|peaks`), eliminating duplication and simplifying future maintenance.

**Fix:** Introduce `{peak_type}` wildcard constrained to `"summits|peaks"` and merge each pair into one rule.

---

## Low Priority

### 13. Missing `benchmark:` directives on compute-heavy rules
**Files:** `align_se`, `align_pe`, `meme_summits`, `meme_peaks`, `bowtie2_index`

Without `benchmark:`, there is no per-rule runtime/memory profiling. Snakemake guidelines recommend benchmarks on any rule expected to run >1 min.

**Fix:** Add `benchmark: OUT + "/benchmarks/{rule_name}/{sample}.{bam_type}.tsv"` to the five heaviest rules.

---

### 14. `motifs.smk:73,106` — Bare `open(...).close()` instead of `Path.touch()`
**File:** `workflow/rules/motifs.smk:73,106`

`open(str(f), "w").close()` for creating empty sentinel files is non-idiomatic.

**Fix:** Replace with `from pathlib import Path; Path(str(f)).touch()`.

---

## Verification

After applying fixes:

1. **Dry-run:** `pixi run snakemake -n --rerun-triggers mtime` — confirms DAG is valid.
2. **filter_chr fix:** `samtools view -H <output.bam> | grep @SQ` — should only contain chromosomes matching the pattern.
3. **DB column fix:** `sqlite3 pipeline_db.db 'SELECT sample, "min5fold_peak#" FROM pipeline_runs'` — should contain real integers, not "NA".
4. **Report script fix:** Mock a narrowPeak with a truncated line and a zero-read stats file; confirm `report.py` does not crash.
5. **Regression:** `pixi run pytest tests/` — existing update_db tests should still pass.
