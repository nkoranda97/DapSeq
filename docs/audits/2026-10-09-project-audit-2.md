# Project audit 2 — 2026-10-09

**Status (2026-10-09):**
- Group 1 was narrowed to `motif_peaks` (items 1, 2 and 32, the version marker from item 18, and the `num_peaks*` wording from item 7) and merged in #14. See `docs/plans/2026-10-09-1415-fix-motif-peaks-count-plan.md`.
- Item 3 was reviewed: the filtered-peak stats are correct as they are, so it is not planned.
- The bowtie2-build log overwrite in item 26 (first-audit item 7) is fixed.
- Groups 2–6 are otherwise not started.

Read-only second audit, run after every group from the [first audit](2026-10-09-project-audit.md) merged. It covers the workflow rules, Python scripts, tests, config, profiles, container and docs.
Findings were checked against the code. Items 1, 2, 5, 7, 10 and 12 were also reproduced with the container's own tools (MEME/FIMO 5.5.9, MACS3 3.0.4, BBMap 39.81 from `apptainer_build/.pixi`) on simulated paired-end data. Items marked *(unverified)* depend on the built `.sif` or the cluster, neither of which was available.

Baseline: `pixi run test` passes (172 tests). Dry runs build the full DAG for PE (multi-lane, with and without a control, blacklist and rmsk on), SE, and mixed configs.

Severity: **high** means wrong results or a crash in a normal run. **Medium** means a supported path breaks, results change silently, or following the docs gives a wrong run. Low means hygiene.

## A. Wrong numbers in the report and database

1. **High: `motif_peaks` counts chromosomes, not peaks.**
   - FASTA headers are always `chr:start-end` (`motifs.smk:37`), and FIMO runs with `--parse-genomic-coord` (`motifs.smk:218`), so FIMO writes only the chromosome in `sequence_name`.
   - `_fimo_motif_peaks` (`collect_stats.py:174`) counts distinct `sequence_name` values. On a simulated run, 17 of 20 peaks had a hit and `motif_peaks` came out as 2.
   - Every report.csv, report.html and DB row carries this value.
   - The test fixture uses `peak1`-style names, which FIMO never writes with this flag.
   - Pre-existing.
2. **Medium: zero FIMO hits still report `NA`.** When nothing matches, FIMO 5.5.9 writes no header row, only a blank line and `#` comments, so the 930143c fix never fires. The README's "What's changed" line saying `motif_peaks` is now 0 is false.
3. **Medium: `chrom_filter` only reaches the MEME FASTA.**
   - `motifs.smk:38` is the only place it is read.
   - `reads_in_peaks_filt`, `frip_filt`, `max_peak_score`, HOMER and the DB's filtered-peak path still include the excluded chromosomes. For example, ChrC/ChrM peaks in an Arabidopsis run inflate `frip_filt`.
   - The settled "_filt describes the MEME input" decision holds only when `chrom_filter` is empty.
4. **Medium: the template ships a human `genome_size`.**
   - `config.yaml:70` sets `"3000000000"`, so the "must be set" check never fires for a template user (the defaults file has `null`).
   - The value is MACS3's `-g` and bamCoverage's `--effectiveGenomeSize` for RPGC. For Arabidopsis both are about 20x too large, with no error.
   - The README never says this is the effective genome size.
5. **Medium: MEME searches at most 100,000 letters, whatever `meme.maxsize` says.**
   - The README describes `maxsize` as the search-size limit. It is actually the dataset-size cap (`-maxsize`).
   - MEME's `-searchsize` default of 100,000 still samples the input. At the template's `maxpeaks: 500`, peaks-mode input likely exceeds that, so part of it is silently left out of motif discovery.
   - README:191 and README:405 contradict each other.
6. **Medium: the template overrides the documented defaults.**
   - `meme.maxpeaks` is 500 instead of 100, and `threads` is 16 instead of 8.
   - A full `resources:` block sits below the HPC-raised defaults from 279fbfa: `trim_align` 32 GB/480 min instead of 64 GB/120 min, `bowtie2_index` 32 GB instead of 64 GB, `fastqc` 2 GB instead of 16 GB.
   - This is left over from first-audit item 21.
7. Low: `num_peaks*` count summit rows, not peak regions. With `--call-summits`, MACS3 writes one row per summit and repeats the region (20 rows for 10 regions in the simulation). The "peaks" FASTA de-duplicates regions, so the counts and the MEME input use different units.
8. Low: on the bwa_mem2 path, `-F 1804` keeps supplementary alignments (`align.smk:150,205`). The FRiP numerator includes them; `mapped_reads` (`-F 2308`) does not. bowtie2 is not affected.
9. Low: falsy per-sample metadata such as `gdna_batch: 0` becomes NA or empty (`collect_stats.py:236`, `update_db.py:356`), because the code tests truthiness instead of `is not None`.

## B. Runs that halt or crash

10. **Medium: MACS3 exits 1 on weak treatments called as BAM.**
    - Only control calls get `--nomodel` (`peaks.smk:37`).
    - An SE treatment, or a PE treatment with an SE control, that has fewer than 100 paired peaks cannot build MACS3's shift model and exits 1 (`callpeak_cmd.py:202` in the MACS3 3.0.4 source).
    - Neither profile sets `keep-going`, so one failed TF stops the run.
    - The PE-with-SE-control case is new since 93907e5, which moved that pair from BAMPE (no model) to BAM. The settled per-call format still stands; the gap is that the fallback call can abort the run.
11. **Medium: `bbduk.max_frags` subsampling needs `bc`, which is not in the container** (`trim.smk:64,184`).
    - `bc` is in neither `apptainer_build/pixi.lock` nor the Ubuntu 24.04 base *(built .sif not inspected)*.
    - The same path calls `reformat.sh` without the `threads=` cap the rule gives bbduk to avoid its deadlock at 12 or more threads.
12. **Medium: setting `meme.base_colors` breaks both MEME rules** (`motifs.smk:87,162`). The alphabet file has unquoted names and an `END ALPHABET` line. MEME 5.5.9 rejects it, writes an empty `meme.txt` and exits non-zero. Quoted names with no END line work.
13. Low: config shape is not checked.
    - A numeric sample name gives a raw `TypeError` (`sample_names.py:24`).
    - A `factorbook:` header with all children commented out gives an `AttributeError` at `motifs.smk:2`.
    - Per-sample key typos are ignored, so `contol: input_A` silently runs the sample with no control.
    - A `macs3.format` typo fails only when MACS3 runs, and a forced `BAMPE` on SE samples suppresses the fallback warning.
14. Low: sample names with shell characters pass the startup checks and then fail in an unquoted shell rule. One example is `TF(2)`, which a test treats as valid. Allowing only `[A-Za-z0-9_+-]` would turn this into a startup error.
15. Low: a user-supplied Factorbook TSV with a blank line raises `IndexError`, and CRLF line endings make every lookup miss silently (`run_factorbook_logo.py:46`). The bundled TSV is clean.
16. Low: `rmsk_filter` silently removes nothing when the file is not in UCSC rmsk layout (`peaks.smk:126`), for example a plant TE BED or GFF.

## C. Results database

17. **Medium: only the first lane of each mate is recorded.**
    - `update_db.get_r1`/`get_r2` keep `r1[0]` in both tables (`update_db.py:154`).
    - All lanes have been trimmed since #4, so the DB's record of input files is incomplete.
    - `test_update_db.py:464` locks this in, with a comment that is wrong.
    - The read-layout plan deferred this, and no group picked it up.
18. Low: rows carry no pipeline version or filter flags.
    - 8a1728b changed what `frip_filt` and `max_peak_score` mean when filters are on, and fixing item 1 will change `motif_peaks`.
    - Nothing in a row says which definition it used. Add a marker before item 1's fix lands.
19. Low: `homer_annotations` is recorded even when HOMER does not run (`update_db.py:263`).
20. Low: the docstring's claim that fcntl locks "work on Lustre/GPFS via NLM" is wrong (`update_db.py:12`). GPFS honours them; Lustre does only with the `flock` mount option *(cluster filesystem not checked)*.

## D. Cluster and container

21. **Medium: the SLURM profile still hardcodes the lab's binds and account** (`profiles/slurm/config.yaml:16,23`). This is first-audit item 13.
    - Plan 2026-06-30 is marked `completed`, but its placeholder-bind step was not done.
    - The local profile binds only `$HOME`.
    - The README covers neither.
22. **Medium: the HPC instructions don't say how to get the SIF.**
    - Every job needs `apptainer_build/dapseq.sif`, which is gitignored.
    - The only build steps are under "Local workstation", and the new "Rebuild the container" step points there.
    - Not covered: whether the cluster build needs `--fakeroot` or a build node, and what a rebuild does to other users' running jobs.
23. Low: MEME jobs reserve `threads` CPUs (8, or 16 from the template) and pass `-p`, but the container's MEME is built without parallel support and runs serially *(seen in the binary, not at runtime)*.
24. Low: MEME writes its own `logo1.png`/`logo_rc1.png` into the directory the logo rules write to. If a run is interrupted between `meme_*` and `meme_logo_*`, the next run keeps MEME's image.
25. Low: the SLURM profile has no `latency-wait`, so it gets Snakemake's 5 s default *(not checked on the cluster filesystem)*. The local profile's `jobs: 4` means 4 cores, not the 4 parallel jobs the README describes.
26. Low: `ref.smk:27` still overwrites the bowtie2-build log (`samtools faidx … 2>{log}`). This is first-audit item 7, never fixed. Python is still 3.13 on the host and 3.12 in the container (item 14).

## E. Documentation

27. **Medium: "Default Pipeline Parameters" contradicts the code.**
    - It lists `call-summits=f`, but every call passes `--call-summits`.
    - It lists `nomodel=f`, but control calls pass `--nomodel`.
    - It lists `tbo=f` and `tpe=f`, but PE trimming passes `tbo tpe`.
    - It lists `minavgquality=0` and `maqb=10`, but the pipeline passes `maq=10`.
    - bowtie2's `--no-mixed --no-discordant` is not listed.
28. **Medium: keys the code reads are undocumented.**
    - The missing keys are `bbduk.max_frags` (which subsamples every sample), `bbduk.adapters`, `threads`, `container`, `resources.*`, `factorbook.*`, `bamcompare.*`, `bwa_mem2.extra*` and `bowtie2.extra_build`.
    - "Any tool accepts an `extra` field" is false for `samtools`.
    - bamCompare is build-on-request by design (plan 2026-07-20). The README presents `-b2` as part of the run and never says how to request `bigWig/<sample>.peaks.bw`.
29. **Medium: "What's changed" (7aa95f1, merged in #13) has one false line and misses two upgrade effects.**
    - The line saying `motif_peaks` is `0` is false (item 2).
    - Rule code changed, so re-running an existing `output_dir` recomputes every sample from trimming (Snakemake's `code` rerun trigger). `--rerun-triggers mtime` or a new `output_dir` avoids this.
    - DB rows recorded under a relative `output_dir` stay as stale duplicates.
30. **Medium: the first audit's status line (a22923e, merged in #13) overstates what closed.**
    - Items 7, 13, 14 and part of 29 are still open, as is the DB's first-lane deferral.
    - Item 6 (masked MEME failures) is moot rather than fixed: Snakemake prepends `set -euo pipefail` to every bash shell block.
31. Low: smaller doc and hygiene gaps.
    - Three implemented plans have no `status:` (2026-06-10, 2026-06-11, and the 2026-10-09 read-layout plan).
    - The multi-user brainstorm still says locks are per output directory.
    - The README does not warn that `snakemake --unlock` clears every user's locks in the shared `.snakemake/` *(unverified)*.
    - `.gitattributes` still misses `apptainer_build/pixi.lock`.
    - `.compound-engineering/config.yaml` is identical to `config.example.yaml`.

## F. Tests

32. **Medium: the FIMO fixtures don't match real FIMO output.** This is how items 1 and 2 pass.
    - Several tests don't check what they are named for. `test_html_header_shows_filter_foldch` passes because `"5"` appears in the CSS colour `#2c3e50`. The blank-cell test checks only that `TF1` appears. The concurrency test bypasses `write_run`.
    - `collect_stats.main()`, which holds the filtered-peak wiring, has no tests, and neither do the CLI entry points.

## Checked and sound

- **First-audit fixes:**
  - 93907e5: per-call MACS3 format, multi-lane trimming, and the lane-count and `genome_size` checks.
  - 1c8332f: defaults no longer leak into user configs.
  - 003ff01: filtered-peak wiring, apart from item 3.
  - 8db7e41: the container lock is in sync, and the profile parses under Snakemake 9.20.
  - afdfc0b: name escaping and name checks.
  - e33055c: Factorbook lookup.
- **DB:** keys are normalised, both tables are written in one transaction, the migration race is tolerated, and `journal_mode=DELETE` is kept.
- **Parsers and report:** read vs pair units match real BBMap 39.81 logs. Also sound: `bamPEFragmentSize` parsing, `narrow_peak_to_fasta` coordinates and empty input, MEME matrix parsing, and report HTML escaping and blank cells.
- **Template:** as shipped it stops with the required-fields message. Filled in, it dry-runs (74 jobs). A pre-audit template also passes, and `macs3.format: BAMPE` forces BAMPE as the upgrade note says.

## Suggested grouping

| Group | Items | Why together |
|---|---|---|
| 1. Motif and peak statistics | 1, 2, 3, 7, 18, 32 | Report and DB numbers are wrong; tests need real FIMO fixtures; rows need a version marker before definitions change |
| 2. Genome size and template defaults | 4, 5, 6 | Following the template silently changes results |
| 3. Run-halting failures | 10, 11, 12, 13, 14, 15, 16 | Supported options or weak samples stop the run |
| 4. DB provenance | 17, 19, 20 | Shared-DB records |
| 5. Cluster and container | 21, 22, 23, 24, 25, 26 | HPC setup and wasted allocation |
| 6. Docs | 27, 28, 29, 30, 31 | Low risk; 29 and 30 are now on `main` (#13) |
| — | 8, 9 | Small; fold into whichever group touches the file |
