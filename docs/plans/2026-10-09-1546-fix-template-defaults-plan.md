---
title: Template Defaults - Plan
type: fix
date: 2026-10-09
topic: template-defaults
artifact_contract: ce-unified-plan/v1
product_contract_source: ce-brainstorm
execution: code
---

# Template Defaults - Plan

## Goal Capsule

- **Objective:** A lab member who starts a config from the template gets the documented default for every setting except the ones they fill in. They are told when their `genome_size` cannot be right for their genome.
- **Means:** Shrink the template to the per-run fields (KTD1), and add a startup warning that compares `genome_size` with the reference's length (KTD2, KTD3).
- **Product authority:** This plan covers audit group 2 items 4 and 6 in `docs/audits/2026-10-09-project-audit-2.md`: the template's human `genome_size`, and template values that override the documented defaults. Item 5 (MEME's search size) and audit groups 3–6 are not active scope. The Product Contract governs behavior; the Planning Contract governs how.
- **Execution profile:** U1–U3 in one pull request from `origin/main`.
- **Stop conditions:** Stop and ask in either of these cases:
  - Snakemake turns out not to merge a user config over `config/config.yaml` key by key.
  - Reading the reference's length at startup turns out to be impractical for a large genome.
- **Open blockers:** None.

---

## Product Contract

Product Contract preservation, restructured with no scope change:
- **Problem Frame:** the description of the template's `resources:` block is corrected. It differs from the defaults, some limits higher and some lower.
- **Outstanding Questions:** the two Deferred-to-Planning questions are resolved by KTD2 and KTD4 and removed.

### Summary

The template config holds only the per-run fields, with `genome_size` left blank, so every other setting comes from the defaults file. The template and README explain that `genome_size` is the effective genome size and give typical values. Each run warns, without stopping, when `genome_size` is larger than the reference genome's total length. The README tells users how to update configs copied from the old template, and how to find past runs that may have used a wrong size.

### Problem Frame

Each experiment config is merged on top of the defaults file (`config/config.yaml`), so any value in a config overrides the documented default. The root template (`config.yaml`) is what users copy, and it sets values that differ from those defaults:
- `genome_size` is `"3000000000"`, a human value.
- `meme.maxpeaks` is 500 instead of 100.
- `threads` is 16 instead of 8.
- Its own `resources:` block differs from the cluster-tuned defaults. Some limits are lower: `trim_align` 32 GB instead of 64 GB, `bowtie2_index` 32 GB instead of 64 GB, `fastqc` 2 GB instead of 16 GB. Others are higher: `meme` and `fimo` 32 GB instead of 16 GB.

`genome_size` is passed to MACS3 as `-g` and to bamCoverage as `--effectiveGenomeSize` for RPGC scaling. The defaults file leaves it null so the startup check forces each config to set it, but the template's pre-filled value means that check never fires.

The lab runs a mix of species with sizes set inconsistently. For an Arabidopsis run that kept the template value, both uses are about 20 times too large, and nothing reports it. Configs already copied from the template carry the same values, so fixing the template alone does not reach them.

### Key Decisions

- **The template holds only per-run fields.** Everything else comes from the defaults file, so the two cannot disagree again. Governs R1, R2. (session-settled: user-approved — chosen over keeping a full copy pinned to the defaults by a test, and over making the template's values the new defaults)
- **An impossible `genome_size` warns but does not stop the run.** Governs R6. (session-settled: user-directed — chosen over stopping the run with a message, and over no check)
- **This plan covers the template only; MEME's search size is separate work.** (session-settled: user-directed — chosen over planning MEME search size (audit item 5) first)

### Requirements

**Template**

- R1. The template sets only the per-run fields: `author`, `samples`, `output_dir`, `genome_ref`, `genome_size`, `gene_annotation` and `aligner`. Other keys appear only as commented examples.
- R2. A config copied from the template, with only the per-run fields filled in, runs with the defaults file's value for every other setting.
- R3. The template's `genome_size` is blank, so a config copied from it stops at startup with the required-field message until a size is set.

**genome_size**

- R4. No config the project ships carries a `genome_size` value: not the template, the defaults file, or the README's example config.
- R5. The template and README say that `genome_size` is the effective genome size (the part of the genome reads can map to), used by MACS3 and by bamCoverage's RPGC scaling. They give typical values for common genomes and cite where those values come from.
- R6. When `genome_size` is larger than the total length of the reference FASTA, the run log shows a warning naming both numbers, and the run continues with the configured value.

**Upgrade notes**

- R7. The README "What's changed" section tells users with configs copied from the old template two things:
  - Delete its `meme.maxpeaks`, `threads` and `resources` lines to get the defaults.
  - Set `genome_size` to their genome's effective size.
- R8. The README explains how to find past runs that may have used a wrong `genome_size`. The results database records `genome_size` and `genome_ref` for every run.

### Acceptance Examples

- AE1. **Covers R3.** **Given** a config copied from the template with `genome_size` left blank, **when** the run starts, **then** it stops with the required-field message naming `genome_size`.
- AE2. **Covers R6.** **Given** an Arabidopsis reference (about 135 Mb) and `genome_size: 3000000000`, **when** the run starts, **then** the log warns with both numbers. MACS3 still receives 3000000000.
- AE3. **Covers R6.** **Given** the same reference and `genome_size: 1.19e8`, **then** no warning appears.
- AE4. **Covers R2.** **Given** a copied template with only the per-run fields filled in, **then** the run uses `meme.maxpeaks` 100, `threads` 8 and the defaults file's `resources`.

### Scope Boundaries

- MEME's search size and the misdescribed `meme.maxsize` (audit item 5).
- Audit groups 3–6.
- Computing the effective genome size automatically.
- Catching a plausible but wrong size, such as a human value for a mouse genome. The warning covers only sizes larger than the reference.
- Changing any value in the defaults file.
- Rewriting lab members' existing configs.

#### Considered and not built

- **Creating the FASTA index at startup so the check can always read it.** The index is a declared output of the index rule (`workflow/rules/ref.smk`). Writing it outside that rule would bypass Snakemake's bookkeeping and could race another user's index build. Revisit only if streaming a FASTA once at startup proves too slow.

### Dependencies / Assumptions

- Each user config is merged on top of `config/config.yaml`, so a key left out of a config takes the default. This is verified in `workflow/Snakefile` (`configfile: "config/config.yaml"`) and in how runs pass `--configfile`.
- The unknown-key guard in `workflow/rules/common.smk` is a hard-coded key list, not one read from the template, so shrinking the template cannot start rejecting valid configs.
- PR #15 (open) also edits the README's "What's changed" section. Whichever merges second rebases onto the other.

### Sources / Research

- `docs/audits/2026-10-09-project-audit-2.md`, items 4 and 6, and first-audit item 21.
- Template versus defaults:
  - `config.yaml:70` (`genome_size: "3000000000"`), `:79` (`threads: 16`), `:144` (`maxpeaks: 500`), `:178` (`resources:`).
  - `config/config.yaml:21` (`genome_size: null`), `:26` (`threads: 8`), `:101` (`maxpeaks: 100`).
- Where `genome_size` is used and checked:
  - `workflow/rules/peaks.smk:51,85` (MACS3 `-g`).
  - `workflow/rules/align.smk:48,105,158,213` (`--effectiveGenomeSize`).
  - `workflow/rules/common.smk:78-79` (the required-field check).
- `README.md:145`: the example config block also carries `"3000000000"`.
- Grounding dossier with the full template-versus-defaults table and history: `/tmp/compound-engineering-1000/ce-brainstorm/audit2-20261009/grounding-template.md` (session scratch). The template values originate in commit 2e287c9.

---

## Planning Contract

### Key Technical Decisions

- KTD1. **The template keeps only the R1 fields, plus commented examples of commonly tuned keys.**
  - The examples are `macs3.foldch_levels`/`meme_foldch_level`, `meme.maxpeaks`/`nmotifs`, `fimo.thresh`, `chrom_filter`, `blacklist_filter`/`rmsk_filter`, `bbduk.max_frags`, `threads`, `db_path` and a `resources:` override. Each example is shown with the defaults file's current value and a pointer to the README options table.
  - The commented values are illustrations only, so they cannot override anything.

  Implements R1, R2. (session-settled: user-approved — chosen over keeping a full copy pinned to the defaults by a test: confirmed at brainstorm)
- KTD2. **The size check runs once per run in the Snakefile's `onstart`, beside the existing paired-end/single-end fallback warning, and logs through the same `logger.warning`.**
  - `onstart` runs only in the main process, before any job, so the warning shows once in the run log, not in every SLURM job's re-parse.
  - Dry runs skip `onstart`, as they already do for the fallback warning.

  Implements R6. (session-settled: user-directed — chosen over stopping the run: confirmed at brainstorm)
- KTD3. **The reference's total length comes from `<genome_ref>.fai` (the sum of its length column) when that index exists. Otherwise the FASTA is read once and its sequence letters are counted.**
  - Gzipped references are read through gzip.
  - Any error while reading or parsing either file makes the helper report the length as unknown, and the check is skipped. Snakemake does not catch exceptions raised in `onstart`, so a raised error would stop the run and break the warn-only decision. Such errors include an unreadable path, a `.fai` that another user's index job is still writing, and a truncated `.gz`. A truly missing reference still fails elsewhere, with its own error.

  The helper lives next to `parse_genome_size` in `workflow/scripts/layout_utils.py`, where the genome-size logic already sits and is tested.
- KTD4. **Typical values come from the container's own MACS3.** The README table takes MACS3's built-in `-g` presets (human, mouse, worm, fly) from `macs3 callpeak --help` in `apptainer_build/.pixi/envs/default/bin`, so the numbers match the tool that consumes them.
  - Genomes without a preset get the rule deepTools' effective-genome-size documentation gives: use the assembly's non-N length.
  - Arabidopsis appears as a worked example only when a value with a citable source is found. Otherwise the README gives the rule alone.
- KTD5. **A test pins the template's contract against the defaults file.** It merges a filled-in copy of the template over `config/config.yaml` with Snakemake's own `snakemake.utils.update_config` (the merge Snakemake applies to `--configfile` layers). Every non-per-run key must equal the default, and the template's top-level keys must be exactly the R1 set. Implements R1, R2, R3.

---

## Implementation Units

### U1. Shrink the template to the per-run fields

- **Goal:** A config copied from the template runs on the documented defaults and must set its own `genome_size`.
- **Requirements:** R1, R2, R3, R4 (template part), R5 (template part); KTD1, KTD5.
- **Dependencies:** None.
- **Files:**
  - `config.yaml`
  - `tests/test_template_config.py` (new)
- **Approach:**
  1. Keep the header comment, `author`, the `samples` skeleton (with its `control`, `experiment_date` and `gdna_batch` comments), `output_dir`, `genome_ref`, `gene_annotation` and `aligner`.
  2. Set `genome_size` to `null`. Its comment says it is the effective genome size, used by MACS3 `-g` and bamCoverage RPGC, must not exceed the reference length, and points to the README's typical values.
  3. Replace every other section with the commented examples KTD1 lists.
- **Patterns to follow:**
  - The existing template's comment style.
  - Module-level YAML loading in tests via the repo root, as `tests/conftest.py` locates `workflow/scripts`.
- **Test scenarios:**
  - The template's top-level keys are exactly `author`, `samples`, `output_dir`, `genome_ref`, `genome_size`, `gene_annotation`, `aligner`.
  - The template's `genome_size` is `None`.
  - Covers AE4. The template is filled in (sample paths, `output_dir`, `genome_ref`, a `genome_size`) and merged over `config/config.yaml` with `snakemake.utils.update_config`. Every key outside the per-run set then equals the defaults file's value, including `meme.maxpeaks` 100, `threads` 8 and the whole `resources` block.
  - The defaults file's `genome_size` is `None`, which pins R4 for that file.
- **Verification:**
  - Tests pass.
  - Covers AE1. A dry run with the template as shipped (sample paths filled, `genome_size` blank) stops with the required-field message naming `genome_size`.
  - The same template with a `genome_size` set dry-runs without errors.

### U2. Warn when `genome_size` exceeds the reference length

- **Goal:** A run whose `genome_size` cannot be right for its reference says so in the run log, then proceeds.
- **Requirements:** R6; KTD2, KTD3.
- **Dependencies:** None.
- **Files:**
  - `workflow/scripts/layout_utils.py` (new helpers next to `parse_genome_size`)
  - `workflow/Snakefile` (`onstart`)
  - `tests/test_layout_utils.py`
- **Approach:**
  1. Add a helper that returns the reference's total length per KTD3, or `None` when it cannot be read.
  2. Add a helper that, given the parsed genome size and that length, returns the warning text when the size is larger and `None` otherwise.
  3. The warning names both numbers, says `genome_size` is the effective genome size and cannot exceed the reference length, and points to the README's "Choosing genome_size" section (U3).
  4. In `onstart`, call both and `logger.warning` the message when there is one. Pass `GENOME_SIZE` and `config["genome_ref"]`.
- **Patterns to follow:**
  - The `fallback_pairs` call and its `logger.warning` in `workflow/Snakefile` `onstart`.
  - The `parse_genome_size` tests in `tests/test_layout_utils.py`.
- **Test scenarios:**
  - Covers AE2. With a `.fai` whose lengths sum to 135,000,000 and a size of 3,000,000,000, the warning names both numbers.
  - Covers AE3. With the same `.fai` and 119,000,000, there is no warning.
  - A size equal to the reference length gives no warning.
  - With no `.fai`, the length comes from the FASTA's sequence letters: headers and newlines are not counted, and lowercase and `N` are.
  - A gzipped FASTA with no `.fai` gives the same length as the plain one.
  - An unreadable path gives `None`, and therefore no warning.
  - A `.fai` with a malformed line (for example, a partially written last line) and a truncated `.fa.gz` both give `None`, never an exception.
- **Verification:**
  - Tests pass.
  - The Verification Contract's startup-warning check prints the warning, and the run completes.

### U3. Document `genome_size` and the template change

- **Goal:** Readers know what `genome_size` means, what values to use, and how to update old configs and check past runs.
- **Requirements:** R4 (README part), R5 (README part), R7, R8; KTD4.
- **Dependencies:** U1, U2.
- **Files:** `README.md`
- **Approach:**
  1. In "What's changed", add the two R7 steps to "Update your setup before the next run". The old template's `resources` block differs from the cluster-tuned defaults in both directions, so the note says "delete it to use the defaults", not that its values are smaller.
  2. Add a `#### Choosing genome_size` subsection directly after the Required fields example block (around `README.md:128-148`). The options table has no `genome_size` row to sit beside. The subsection covers:
     - What the effective genome size is.
     - The KTD4 table and rule.
     - That a value larger than the reference now triggers a startup warning.
  3. Set the example config block's `genome_size` (`README.md:145`) to blank, with a pointer to that subsection.
  4. In the MACS3 `gsize=hs` parameter entry (around `README.md:384`), replace its "estimate by removing Ns and simple repeats" advice with a pointer to the subsection, so the README gives one rule. Point the bamCoverage `effectiveGenomeSize=` entry (around `README.md:365`) there too.
  5. Add the R8 note: the results database's `pipeline_runs` table holds `genome_size` and `genome_ref` per row, so a lab member can list runs whose size does not fit their reference and re-run them.
- **Test scenarios:** Test expectation: none -- documentation only; U1 and U2 tests pin the behavior it describes.
- **Verification:** The README's `genome_size` text, example block and "What's changed" agree with the template and the warning. Every listed preset matches the container's `macs3 callpeak --help`.

---

## Risks & Dependencies

- **Startup cost on a genome's first run:** reading a 3 GB FASTA once to total its length takes tens of seconds before the index exists. Later runs read the `.fai`.
- **PR #15:** it adds a "What's changed" bullet in the same list as R7. A trivial rebase whichever merges second.

---

## Verification Contract

| Check | How | Proves |
|---|---|---|
| Unit tests | `pixi run test` | U1, U2 |
| Required-field stop | Dry run (`pixi run snakemake -n --configfile <filled template, genome_size blank>`) stops naming `genome_size` | AE1, R3 |
| Template validates | The same template with a `genome_size` set dry-runs cleanly | R2 |
| Startup warning | A real run that targets only the database step's output. Make a scratch output folder holding a hand-written `stats/report.csv`, and a scratch config with an oversize `genome_size`, a small scratch reference and `db_path` in scratch. Then run `pixi run snakemake <scratch_out>/stats/db_updated.flag --configfile <scratch config>`; `--allowed-rules` would swallow the target path. The run log shows the warning, and the run completes | AE2, R6 |
| Preset values | The README's values equal those in `apptainer_build/.pixi/envs/default/bin/macs3 callpeak --help` | R5, KTD4 |

---

## Definition of Done

- R1–R8 hold, and AE1–AE4 are each covered by a test or a Verification Contract check.
- `pixi run test` passes.
- The defaults file is unchanged.
- One pull request from `origin/main`, carrying this plan.
- No leftover code from abandoned approaches is in the diff.
