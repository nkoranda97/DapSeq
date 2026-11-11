# Project audit 3 — 2026-11-11

**Status:**
- 2026-11-11, #17: second-audit item 5 and items 10–16, with N1, N3 and N4, merged. Item 15's CRLF claim was wrong (text-mode `open` already handles CRLF); only blank lines were a bug.
- 2026-11-11, branch `fix/db-provenance`: items 9, 17, 19 and 20 (group 4 plus item 9). Multi-lane `r1`/`r2` are stored joined with `;`.
- Still open: items 8, 21–25 and 27–32, and N2.

Read-only third audit, run after #14 (`motif_peaks`), #15 (bowtie2-build log) and #16 (template defaults, `genome_size` warning) merged. It covers the workflow rules, Python scripts, config, profiles and tests, and re-checks the open items of the [second audit](2026-10-09-project-audit-2.md).

Baseline: `pixi run test` passes (205 tests). Probes used the container's own tools from `apptainer_build/.pixi` (FIMO 5.5.9, bedtools, samtools, Biopython, logomaker). HOMER could not be run locally (see "Unverified").

No new high-severity bugs. The #14, #15 and #16 changes hold up.

Severity follows the second audit: **high** means wrong results or a crash in a normal run, **medium** means a supported path breaks or results change silently, low means hygiene.

## New findings

- **N1. Low (costly when hit): a gzipped `genome_ref` fails late.**
  - The new genome-size check reads `.gz` references (`layout_utils.py:91`), but the rest of the pipeline cannot use them.
  - With plain gzip, `samtools faidx` in the index rule fails.
  - With bgzip, indexing, alignment and MACS3 all finish, then `narrow_peak_to_fasta.py:38` crashes with `UnicodeDecodeError` (reproduced). HOMER probably fails the same way *(unverified)*.
  - The startup checks only test that the file exists. A startup error for `.gz` references would catch this before any compute runs.
- **N2. Low: duplicate sequences can reach MEME in peaks mode.**
  - The step that removes duplicate regions keeps every row that ties for the highest fold-change (`narrow_peak_to_fasta.py:70-71`).
  - When two summits of one region have the same fold-change, the region goes into the FASTA twice (reproduced: 2 entries for 1 region). MEME gets a duplicate site and a `maxpeaks` slot is used up.
- **N3. Low: trimming detects gzip with `gzip -t`.**
  - `gzip -t` decompresses the whole file just to check it, so in SE mode each input is decompressed 4 times (`trim.smk:97,108`). In PE mode a single-lane input is decompressed once extra.
  - A `.gz` with trailing bytes makes `gzip -t` return 2 (confirmed), so the compressed bytes are piped raw into bbduk. A truncated download behaves the same way. The user sees a bbduk parse error instead of "corrupt gzip".
  - The PE lane-joining code already checks the gzip magic bytes (`trim.smk:205`), which avoids both problems.
- **N4. Low (add to second-audit item 12): `meme.base_colors` values made only of digits are read as numbers.**
  - Unquoted, `000000` loads as `0`, so it is silently replaced by the default colour.
  - `123456` loads as an int, and `.lstrip` at `motifs.smk:64` raises `AttributeError` while the workflow loads.
  - Unquoted `#008000` is read as a YAML comment (null).
- **N5. Hygiene (add to second-audit item 32): `write_rows` and `write_meta_rows` (`update_db.py:236,245`) are only called by tests.** They skip the `motif_peaks` blanking that the real `write_run` applies, so those tests check a path the pipeline never uses.

## Unverified

- HOMER on a control's self-call peaks, which are often empty. This runs for every control when `gene_annotation` is set. The local `apptainer_build/.pixi` HOMER has its install path hard-coded to the repo's old location (`~/Projects/DapSeq`) in `annotatePeaks.pl` and `HomerConfig.pm`, so it could not be run. This affects only the local env, not the `.sif`.

## Checked and sound

- **#14:** real FIMO 5.5.9 with `--no-pgc` names hits by peak. `motif_peaks` counted 2 for 2 peaks with hits, and zero hits gives `0`, not `NA`.
- **Empty peak file:** FRiP gives 0 (bedtools intersect with an empty `-b` returns 0).
- **Logos:** logomaker renders motif rows with zero probabilities without warnings.
- **Scripts that import each other:** `update_db` imports `collect_stats`, and `render_report_html` imports `report`. The importing script never triggers the other's Snakemake entry point.
- **DB writes:** blanking, delete and insert happen in one transaction.
- **#16:** `genome_size` parsing handles scientific notation, and the oversized-size check only warns.

## Still open from the second audit

No commits have touched these since 2026-10-09. Items 5, 8, 9, 11, 12, 17, 20 and 25 were re-checked in the code today. The others are carried over from the second audit without being re-checked one by one.

- **Item 5 (medium), MEME sampling:** the rules pass `-maxsize`, not `-searchsize`, so MEME still samples its input down to 100,000 letters. With the new `maxpeaks: 100` default, peaks mode is near the cap. Anyone who raises `maxpeaks` loses part of their input without being told.
- **Items 10–16 (run-stopping):**
  - Item 10: a weak SE treatment, or a PE treatment with an SE control, can't build MACS3's shift model and exits 1. No profile sets `keep-going`.
  - Item 11: `max_frags` subsampling needs `bc`, which is not in `apptainer_build/pixi.toml`.
  - Item 12: setting `base_colors` breaks both MEME rules.
  - Items 13–16: config-shape checks, sample names with shell characters, the Factorbook TSV, and rmsk layout.
- **Items 8, 9:** bwa-mem2 keeps supplementary alignments, and falsy per-sample metadata becomes NA (`collect_stats.py:252,255`).
- **Items 17, 19, 20 (DB):**
  - Item 17: only the first lane of each mate is recorded (`update_db.py:165`).
  - Item 19: `homer_annotations` is recorded even when HOMER does not run.
  - Item 20: the docstring's Lustre claim is wrong (`update_db.py:12`).
- **Items 21–25 (cluster):**
  - Item 21: the SLURM profile hard-codes the lab's binds and account.
  - Item 22: the README does not say how to get the `.sif`.
  - Item 23: MEME is given `-p` but runs serially.
  - Item 24: MEME's own logo files share the logo rules' directory.
  - Item 25: there is no `latency-wait`.
- **Items 27–31:** docs.
- **Item 32:** tests.

## Suggested next group

Item 5, then the run-stopping group (items 10–16), with N1 and N4 folded in. Both turn a late crash or a silent fallback into a clear startup error. N2 and N3 are small and can go with whichever change touches their file.
