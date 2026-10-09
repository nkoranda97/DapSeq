# Project audit — 2026-10-09

Read-only audit of the workflow rules, Python scripts, tests, config, docs, and environment.
Findings were spot-checked against the code; items marked *(unverified)* depend on tool behaviour not tested here.

## A. Correctness bugs (wrong results or crashes)

1. **Multi-lane paired-end trimming is broken** — `trim.smk` `trim_pe`: `R1="${R1_FILES[0]}"` keeps only the first R1 lane, and `R2={input.r2}` is unquoted, so two R2 lanes become a shell command. Lists are clearly intended (`get_r1`, `common.smk`).
2. **Single-end samples are peak-called as BAMPE** — `peaks.smk` `macs3` / `macs3_control` take `-f` from global `macs3.format` (default `BAMPE`) rather than per-sample SE/PE. SE or mixed runs fail at MACS *(MACS3 rejection unverified)*.
3. **Defaults leak into user configs via recursive merge** — `config/config.yaml` defines samples `TF1`/`TF2`/`input_DNA` and enables hg38 `blacklist_filter` + `rmsk_filter`. A user sample named `TF1` without `control:` inherits `input_DNA`; a user config omitting the filter sections gets hg38 filtering on any genome, and `blacklist/rmsk.txt.gz` is gitignored, so fresh clones fail path validation.
4. **"_filt" stats and DB paths describe the wrong peak set** — `collect_stats.py` and `update_db.py` use `peaks_fold{N}.narrowPeak` (before blacklist/rmsk), while MEME/FIMO use `get_final_filtered_peaks` (`..._bl_rmsk`). The report header says the opposite. HOMER annotation also uses the unfiltered set.
5. **`genome_size` handling** — comment says `"1e5"` notation is accepted, but the same value goes to `bamCoverage --effectiveGenomeSize` (int) *(deeptools parsing unverified)*. A null `genome_size` isn't caught by `_missing_fields`.
6. **MEME failures masked** — the MEME shell blocks in `motifs.smk` lack `set -e` and end with `rm -f`, so a MEME error exits 0 and shows up as MissingOutput.
7. **`ref.smk` bowtie2 branch** — `samtools faidx ... 2>{log}` overwrites the bowtie2-build log (should be `2>>`).

## B. Results database

8. **Idempotency key is the raw `output_dir` string** — `update_db.py` deletes `WHERE output_dir = ?` with an unnormalised path; `out/` vs `out` vs absolute leaves duplicates, and two projects using relative `results` clobber each other in the shared DB.
9. **Schema-migration race** — `PRAGMA table_info` then `ALTER TABLE ADD COLUMN` outside a transaction; two concurrent runs against an old DB → `duplicate column name`. `pipeline_runs` and `run_metadata` commit separately; connections aren't closed on exception.
10. **Missing/always-NA columns** — `alignment_rate` is always NA (not produced by `report.make_cols()`); `frip`, `num_peaks_bl`, `num_peaks_rmsk`, per-sample `experiment_date` not stored. HOMER paths recorded even though HOMER doesn't run by default.

## C. Container & environment

11. **logomaker missing from the container** — imported by `logo_utils.py`, used by logo rules that run in the container, but only in root `pixi.toml`. Logo rules likely fail under Apptainer *(built .sif not inspected)*.
12. **Container build ignores its lockfile** — `dapseq.def` copies only `pixi.toml`; tools pinned `"*"`, base image `pixi:latest` → non-reproducible builds.
13. **SLURM profile** — bind paths/account are lab-specific and committed; `--bind a, --bind b` has a stray comma; `jobscript:` is likely ignored by the Snakemake 8+ SLURM executor, so `module load apptainer` may never run *(unverified)*.
14. Minor drift: Python 3.13 (root) vs 3.12 (container); `min_version("8.0")` vs `snakemake >=9.19`.

## D. Tests

15. **Default test command fails** — pytest recurses into `apptainer_build/dapseq_sandbox/` (17 collection errors); no `testpaths`, no pixi `test` task.
16. **10 stale tests** in `tests/test_narrow_peak_to_fasta.py` reference the removed entropy filter (`_triplet_entropy`, `complexity_filter_enabled`). `pytest tests` → 84 passed, 10 failed.
17. **Coverage gaps** — no tests for `logo_utils`, `render_meme_logo`, `run_factorbook_logo`, `render_report_html`, `update_db.main()`, most `collect_stats` parsers, or HTML escaping in the report.

## E. Documentation

18. **README breaks runs** — required-fields example includes `slurm_partition`/`slurm_account`, which `common.smk` now rejects with `ValueError`.
19. **README documents removed options** — `macs3.min_foldch` and `*_peaks_filt.narrowPeak`; real keys are `macs3.foldch_levels` + `meme_foldch_level`, and outputs are `peaks_fold{N}[_bl][_rmsk]`. Setting `min_foldch` is silently ignored (the unknown-key guard doesn't check inside `macs3:`).
20. README lacks an outputs section (layout, `config_used.yaml`, HTML report), test instructions, and `blacklist_filter`/`rmsk_filter` in the options table; `pixi run apptainer build` claims apptainer is in pixi (it isn't).
21. Root template `config.yaml` has drifted from defaults: dead `resources.merge_control`, missing `bwa_mem2_*`/`sample_stats`, opposite filter defaults, different `threads`/`maxpeaks`. `samtools.extra_merge` is dead in both.
22. Stale planning docs: several implemented plans still `status: active`; the multi-user brainstorm still says "Ready for planning" and specifies WAL mode, but the code deliberately uses `journal_mode=DELETE`; claim that Snakemake locks are per-output-dir is inaccurate (shared `.snakemake/` in repo root). `todo.md` lists done items.

## F. Robustness & hygiene (low)

23. Sample names joined into regexes without `re.escape` (`common.smk`); a sample named `X_control` where `X` is a control → ambiguous MACS output.
24. Unquoted paths in shell blocks (spaces break); no `temp()` on trimmed FASTQs; no `benchmark:` directives; `trim_se` reads inputs three times.
25. Logo rules borrow `meme.mem_mb` (16–32 GB) with hardcoded runtimes not settable from config.
26. Index rules write next to `genome_ref` — shared-genome users race to build the index and need write access.
27. `homer_annotate` not in `rule all`, so setting `gene_annotation` does nothing by default.
28. Factorbook matching uses exact `sample.upper() == target`, so `CTCF_rep1` or plant TFs get no logo; unhandled `ValueError` if TSV columns are missing.
29. Repo: `.DS_Store` tracked; empty `yamls/`; no LICENSE / CI / lint config; `.gitattributes` misses `apptainer_build/pixi.lock`; uncommitted `.compound-engineering` template swap (new `config.yaml` and `config.example.yaml` are identical).

## Suggested grouping for follow-up work

| Group | Items | Why together |
|---|---|---|
| 1. SE/PE & multi-lane correctness | 1, 2, 5 | Runs fail or silently drop data |
| 2. Config defaults & template | 3, 21, 18, 19 | Same root cause: merged defaults and drifted template/docs |
| 3. Filtered-peak consistency | 4, 10, 27 | Stats, DB, HOMER should all describe the MEME peak set |
| 4. DB robustness | 8, 9 | Shared-DB integrity |
| 5. Container reproducibility | 11, 12, 13, 14 | Builds and cluster runs |
| 6. Test suite repair & coverage | 15, 16, 17 | Gate for every other change |
| 7. Docs & hygiene cleanup | 20, 22, 29 | Low-risk, independent |
