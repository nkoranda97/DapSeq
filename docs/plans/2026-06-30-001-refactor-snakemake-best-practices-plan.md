---
title: "refactor: Audit and remediate Snakemake best-practice violations"
date: 2026-06-30
status: completed
type: refactor
---

# refactor: Audit and remediate Snakemake best-practice violations

## Summary

Comprehensive audit of the DapSeq Snakemake workflow (11 `.smk` files, 2 profile configs). Findings are grouped into five categories: shell safety, SLURM boilerplate duplication, profile portability, Python/Snakemake idiom issues, and pipeline connectivity gaps. Two items — output directory casing (`Fastqc`, `MACS`) and single-global-threads — are noted but deferred.

---

## Problem Frame

The pipeline works but accumulates small deviations from Snakemake conventions across every rule file. The most consequential issues are (1) missing `set -euo pipefail` in complex `trim.smk` shell blocks, which can silently corrupt outputs, and (2) `slurm_partition`/`slurm_account` repeated verbatim in every rule's `resources:` block — ~38 lines of boilerplate that belong in the SLURM profile instead. Other issues are lower severity but harm maintainability and portability.

---

## Issues Catalog

### A — Shell safety

| # | Location | Issue |
|---|----------|-------|
| A1 | `workflow/rules/trim.smk` `trim_se` and `trim_pe` shell blocks | Multi-line bash with arrays, loops, `bc` arithmetic, and subshells — no `set -euo pipefail`. Silent failure at any step would leave a truncated or empty output file and Snakemake would consider the rule successful. |
| A2 | `workflow/rules/peaks.smk` `filter_peaks` shell | Three `awk` commands; only the first redirects stderr to `{log}`. Errors from the second and third awk invocations are silently discarded. |
| A3 | `workflow/rules/stats.smk` `idxstats` shell | Uses `exec 2>{log}` (unconventional) but no `set -euo pipefail`, leaving the `samtools idxstats` command unguarded against failure. |

### B — SLURM boilerplate repetition

| # | Location | Issue |
|---|----------|-------|
| B1 | All 9 rule files, every `resources:` block | `slurm_partition = config["slurm_partition"]` and `slurm_account = config["slurm_account"]` appear ~19 times. Snakemake's `default-resources:` in the profile is the correct home for these cluster-wide defaults. |
| B2 | `config.yaml` and `config/config.yaml` | `slurm_partition` and `slurm_account` keys exist in both config files and in `common.smk`'s `_KNOWN_CONFIG_KEYS`. Once moved to the profile, these become dead config entries. |

### C — Profile portability

| # | Location | Issue |
|---|----------|-------|
| C1 | `profiles/local/config.yaml:11` | `apptainer-args: "--bind /home/nkoranda"` — hardcoded username. Any other user gets a bind path that doesn't exist. |
| C2 | `profiles/slurm/config.yaml:12` | `apptainer-args: "--bind /project/zhuzhuzhang,/project/gzy8899"` — hardcoded HPC project directories that are specific to one lab/user. Collaborators must edit this manually with no guidance. |

### D — Python / Snakemake idiom issues

| # | Location | Issue |
|---|----------|-------|
| D1 | `workflow/rules/filter.smk:4-5` `merge_control` input | `lambda wc: expand(...)` wrappers on both `bams` and `bais` that never reference `wc`. Plain `expand()` calls suffice. |
| D2 | `workflow/rules/qc.smk:17` `fastqc` shell | `--outdir={OUT}/Fastqc/` uses the global Python variable `{OUT}` directly in the shell string. Directory paths should route through `params:` so they appear in the rule's provenance. |
| D3 | `workflow/rules/qc.smk:36` `multiqc` shell | Same issue: `{OUT}/Fastqc/` and `-o {OUT}/` hardcoded in the shell string. |
| D4 | `workflow/rules/align.smk:1` and `workflow/rules/ref.smk:1` | Both define `_aligner = config.get("aligner", "bowtie2")` identically. Should be defined once in `common.smk` and reused. |
| D5 | `workflow/rules/motifs.smk` `meme_summits`/`meme_peaks` params | Several params (`has_custom`, `a_color`, `c_color`, `g_color`, `t_color`) use `lambda wc:` but never reference `wc`. These should be plain values. |

### E — Stale resource config keys

| # | Location | Issue |
|---|----------|-------|
| E1 | `config.yaml` `resources:` | Keys `bowtie2`, `samtools_filter_sort`, `bamcoverage`, `filtered_stats`, `frip_macs` are defined but not referenced by any current rule. Relics from prior rule structures. |
| E2 | `config/config.yaml` `resources:` | Keys `bowtie2`, `samtools_filter_sort`, `bamcoverage`, `qc_summary` are stale for the same reason. |
| E3 | `workflow/rules/motifs.smk` `factorbook_logo`, `meme_logo_summits`, `meme_logo_peaks` | These rules borrow `config["resources"]["report"]` for their memory/runtime. They are not report generation; this semantic mismatch makes resource tuning confusing. They should have their own resource keys (e.g., `logo`) or use the `meme` resource entry. |

### F — Pipeline connectivity and path portability

| # | Location | Issue |
|---|----------|-------|
| F1 | `workflow/rules/stats.smk` `idxstats` rule | Produces `{OUT}/stats/{sample}.idxstats.txt` but this output is not requested by any aggregate rule (`qc_done`, `mapped`, `peaked`, `rule all:`). The rule is an orphan — it will never run unless explicitly targeted on the command line. |
| F2 | `workflow/rules/motifs.smk:3-4` `factorbook_logo` input | Hard-coded relative paths `"factorbook/factorbook_chip_seq_meme_motifs.tsv"` and `"factorbook/factorbook_chip_seq_meme_motif_catalog.meme"` depend on the current working directory. All other external file inputs are config-parameterized or use `workflow.basedir`. |
| F3 | `workflow/rules/db.smk:7` `update_db` params | `db_path = os.path.join(os.path.dirname(workflow.basedir), "pipeline_db.db")` computes the database path relative to the **parent** of the workflow directory. This assumes a specific installation layout and breaks if the pipeline is moved or symlinked. |

### G — Noted / deferred (not in scope)

| # | Issue | Reason deferred |
|---|-------|----------------|
| G1 | Output directory casing: `{OUT}/Fastqc/` (mixed), `{OUT}/MACS/` (all caps) vs lowercase elsewhere | Renaming would invalidate all existing run outputs — high blast radius for a cosmetic fix |
| G2 | Single global `config["threads"]` for all rules — alignment benefits from many threads; `awk`/filter rules don't | Requires per-rule tuning and performance profiling |

---

## Key Technical Decisions

**SLURM partition/account: profile-owned vs config-owned.** The current design routes `slurm_partition` and `slurm_account` through `config.yaml` so users fill them in once alongside sample/genome settings. Moving them to `profiles/slurm/config.yaml` `default-resources:` means users must edit the profile instead. The profile is the more semantically correct home — it is already cluster-specific (`use-apptainer`, `jobs`, bind paths), whereas `config.yaml` is experiment-specific. Post-fix, the slurm profile becomes the single place for cluster identity. This is the Snakemake-idiomatic approach for Snakemake ≥ 8.

**`idxstats` rule: wire into target vs remove.** The `idxstats` rule produces useful per-chromosome alignment counts. It should be wired into `qc_done` rather than removed, so it runs as part of the normal QC pass.

**`factorbook_logo` input paths: config key vs `ancient()`.** The factorbook files are large, rarely changing reference databases. The fix is to add `factorbook_tsv` and `factorbook_meme` config keys in `config.yaml` (with sensible defaults pointing to the existing `factorbook/` directory), matching the pattern used for `genome_ref`, `gene_annotation`, and blacklist files.

---

## Scope Boundaries

### In scope
- All `.smk` files under `workflow/rules/`
- `workflow/Snakefile`
- `profiles/local/config.yaml` and `profiles/slurm/config.yaml`
- `config.yaml` and `config/config.yaml` (resource key cleanup and SLURM key removal)
- `workflow/rules/common.smk` (`_KNOWN_CONFIG_KEYS` update)

### Deferred to Follow-Up Work
- Directory renames (`Fastqc` → `fastqc`, `MACS` → `macs`) — separate breaking change
- Per-rule thread tuning — requires benchmarking
- Adding `benchmark:` directives to all rules — significant boilerplate, low priority

---

## Implementation Units

### U1. Fix shell safety in trim.smk, filter_peaks, and idxstats

**Goal:** Ensure all multi-line shell blocks fail loudly on error rather than silently continuing.

**Requirements:** Addresses A1, A2, A3.

**Dependencies:** None.

**Files:**
- `workflow/rules/trim.smk`
- `workflow/rules/peaks.smk`
- `workflow/rules/stats.smk`

**Approach:**
- `trim_se` and `trim_pe`: Add `set -euo pipefail` as the first line of each shell block. Both blocks already use bash-specific syntax (arrays, `((...))`, `bc`) so bash mode is in effect; `set -euo pipefail` is idiomatic there.
- `filter_peaks`: Redirect stderr from all three `awk` invocations to `{log}`. The standard pattern for multiple commands sharing one log is to open the log once with `exec 2>{log}` at the top, or redirect each line. Given `set -euo pipefail` is already present, appending `2>>{log}` on commands 2 and 3 (append, not overwrite) is simplest.
- `idxstats`: Replace the `exec 2>{log}` idiom with the standard `set -euo pipefail` + `... > {output} 2>{log}` pattern to be consistent with the rest of the pipeline.

**Patterns to follow:** `workflow/rules/align.smk` shell blocks (have `set -euo pipefail`).

**Test scenarios:**
- Manually verify `trim_se` shell is syntactically valid bash by running `bash -n` on a temp file containing the block (with dummy substitutions)
- Confirm `filter_peaks` log captures stderr from all three awk invocations by inspecting the log after a dry run with an intentionally bad input
- Confirm `idxstats` log format matches other rules

**Verification:** After changes, `grep -n 'set -euo' workflow/rules/trim.smk` should match the start of both `trim_se` and `trim_pe` shell blocks.

---

### U2. Remove per-rule SLURM boilerplate; move to profile default-resources

**Goal:** Eliminate the ~38-line duplication of `slurm_partition`/`slurm_account` across all rule `resources:` blocks by routing them through `profiles/slurm/config.yaml` `default-resources:`.

**Requirements:** Addresses B1, B2.

**Dependencies:** None (can be applied independently).

**Files:**
- `workflow/rules/align.smk`
- `workflow/rules/trim.smk`
- `workflow/rules/ref.smk`
- `workflow/rules/filter.smk`
- `workflow/rules/peaks.smk`
- `workflow/rules/motifs.smk`
- `workflow/rules/stats.smk`
- `workflow/rules/qc.smk`
- `workflow/rules/annotation.smk`
- `workflow/rules/db.smk`
- `profiles/slurm/config.yaml`
- `config.yaml`
- `config/config.yaml`
- `workflow/rules/common.smk`

**Approach:**
1. In every rule's `resources:` block, delete the two lines:
   ```
   slurm_partition = config["slurm_partition"],
   slurm_account   = config["slurm_account"],
   ```
2. In `profiles/slurm/config.yaml`, add to `default-resources:`:
   ```yaml
   default-resources:
     mem_mb: 8000
     runtime: 120
     slurm_partition: "FILL_IN"    # e.g. caslake
     slurm_account: "FILL_IN"      # e.g. pi-yourlab
   ```
   Add an inline comment instructing users to fill these in when installing.
3. Remove `slurm_partition:` and `slurm_account:` keys from `config.yaml` and `config/config.yaml`.
4. In `common.smk`, remove `"slurm_partition"` and `"slurm_account"` from `_KNOWN_CONFIG_KEYS`.
5. Update the README's "SLURM" section to tell users to set partition/account in the slurm profile rather than in `config.yaml`.

**Patterns to follow:** `profiles/local/config.yaml` already shows the `default-resources:` block format. `profiles/slurm/config.yaml` existing structure.

**Test scenarios:**
- `grep -rn 'slurm_partition' workflow/rules/` returns zero results after the change
- `grep -rn 'slurm_partition' profiles/slurm/config.yaml` returns the `default-resources:` entry
- `grep -rn 'slurm_partition\|slurm_account' config.yaml config/config.yaml` returns zero results
- A Snakemake dry-run (`--dry-run`) with `--profile profiles/slurm` succeeds without KeyError

**Verification:** Dry-run passes cleanly with both profiles after the change.

---

### U3. Fix hardcoded paths in profile configs

**Goal:** Make both profile configs usable by any user without manual path editing.

**Requirements:** Addresses C1, C2.

**Dependencies:** None.

**Files:**
- `profiles/local/config.yaml`
- `profiles/slurm/config.yaml`
- `README.md` (update setup instructions)

**Approach:**

**Local profile (`profiles/local/config.yaml:11`):**
Replace the hardcoded user home directory with a more portable approach. Apptainer supports environment variable expansion in bind arguments. Change:
```yaml
apptainer-args: "--bind /home/nkoranda"
```
to:
```yaml
apptainer-args: "--bind $HOME"
```
Apptainer expands `$HOME` at bind time. Add a comment explaining why this is needed (Apptainer's `--home` only mounts the pipeline directory; this bind makes data outside it accessible).

**SLURM profile (`profiles/slurm/config.yaml:12`):**
The HPC project directories (`/project/zhuzhuzhang`, `/project/gzy8899`) are user-specific. Replace the hardcoded line with a placeholder and a comment instructing users to fill in their lab's project directory:
```yaml
# Bind your HPC project directory so data paths outside /home are visible.
# Replace with your actual project path(s), e.g. /project/pi-yourlab
apptainer-args: "--bind /project/FILL_IN"
```

**Test scenarios:**
- Local profile with `$HOME` expanded: verify `apptainer --bind $HOME` is accepted by running `apptainer exec --bind $HOME dapseq.sif echo ok` manually
- SLURM profile placeholder: confirm placeholder text is clearly instructional and won't be used as-is (the bind path doesn't exist, so Apptainer would fail loudly rather than silently)

**Verification:** Any developer who clones the repo and runs `pixi run snakemake --profile profiles/local` without editing profile files gets a working bind (not a path error).

---

### U4. Python / Snakemake idiom cleanup

**Goal:** Remove unnecessary lambda wrappers, eliminate direct Python variable use in shell commands, and deduplicate the `_aligner` variable.

**Requirements:** Addresses D1–D5.

**Dependencies:** None.

**Files:**
- `workflow/rules/filter.smk`
- `workflow/rules/qc.smk`
- `workflow/rules/align.smk`
- `workflow/rules/ref.smk`
- `workflow/rules/common.smk`
- `workflow/rules/motifs.smk`
- `config.yaml`
- `config/config.yaml`

**Approach:**

**D1 — `merge_control` unnecessary lambdas (`filter.smk:4-5`):**
```python
# Before
bams = lambda wc: expand(OUT + "/bam/{s}.bam", s=CONTROL_SAMPLES),
bais = lambda wc: expand(OUT + "/bam/{s}.bam.bai", s=CONTROL_SAMPLES),
# After
bams = expand(OUT + "/bam/{s}.bam", s=CONTROL_SAMPLES),
bais = expand(OUT + "/bam/{s}.bam.bai", s=CONTROL_SAMPLES),
```
Snakemake evaluates plain `expand()` calls at DAG-build time — correct here since `CONTROL_SAMPLES` is a module-level constant.

**D2/D3 — `{OUT}` in qc.smk shell commands:**
Add `outdir` to `params:` in both `fastqc` and `multiqc`, then reference `{params.outdir}` in the shell string:
```python
# fastqc params
outdir = OUT + "/Fastqc",
# fastqc shell
"fastqc {params.extra} {input} --outdir={params.outdir}/ 2>{log}"

# multiqc params
indir  = OUT + "/Fastqc",
outdir = OUT,
# multiqc shell
"multiqc {params.extra} {params.indir}/ -o {params.outdir}/ 2>{log}"
```

**D4 — Duplicate `_aligner` in `align.smk` and `ref.smk`:**
Remove the `_aligner = ...` definition from both files. Add it once to `common.smk` (after the config validation block):
```python
ALIGNER = config.get("aligner", "bowtie2")
```
Use `ALIGNER` (uppercase, module-level constant) in the `if/elif` branches of both files.

**D5 — Unnecessary lambdas in `motifs.smk` params:**
The color params (`has_custom`, `a_color`, etc.) are constant per config load. Replace the lambda wrappers with plain computed values:
```python
# Before
has_custom = lambda wc: "1" if config["meme"].get("base_colors") else "",
a_color    = lambda wc: ((config["meme"].get("base_colors") or {}).get("A") or "008000").lstrip("#"),
# After
has_custom = "1" if config["meme"].get("base_colors") else "",
a_color    = ((config["meme"].get("base_colors") or {}).get("A") or "008000").lstrip("#"),
```
These params appear in both `meme_summits` and `meme_peaks` — apply in both places.

**Stale resource keys (E1, E2):**
Remove from `config.yaml` `resources:` section: `bowtie2`, `samtools_filter_sort`, `bamcoverage`, `filtered_stats`, `frip_macs`.
Remove from `config/config.yaml` `resources:` section: `bowtie2`, `samtools_filter_sort`, `bamcoverage`, `qc_summary`.

**E3 — `motifs.smk` borrow of `report` resource:**
Rename the resource lookup in `factorbook_logo`, `meme_logo_summits`, and `meme_logo_peaks` to use the `meme` resource key (they are already gated by the same runtime/memory profile), or add a dedicated `logo` resource key to both config files if different sizing is ever needed. Using `meme` is the minimal correct fix.

**Patterns to follow:** Other plain-value params throughout the codebase (e.g., `peaks.smk` `_MACS3_KEEP_DUP`). The `ALIGNER` naming convention follows `CONTROL`, `OUT`, `SAMPLES`, etc.

**Test scenarios:**
- `grep -n 'lambda wc' workflow/rules/filter.smk` returns zero results after D1
- `grep -n '{OUT}' workflow/rules/qc.smk` returns zero results after D2/D3
- `grep -n '_aligner\|ALIGNER' workflow/rules/align.smk workflow/rules/ref.smk` shows only the usage (not the definition) in both files, and a definition in `common.smk`
- Stale keys absent from `resources:` subsection: `grep -n '  bowtie2:\|  samtools_filter_sort:\|  bamcoverage:\|  filtered_stats:\|  frip_macs:\|  qc_summary:' config.yaml config/config.yaml` returns zero results. Note: use leading-space prefix to distinguish `resources.bamcoverage` (stale) from the top-level `bamcoverage:` tool-config section (live — must remain). Confirm `bamcoverage` is still present as a top-level key in `_KNOWN_CONFIG_KEYS`.
- Snakemake dry-run passes after all changes

**Verification:** `snakemake --dry-run --profile profiles/local --configfile config.yaml` completes without error.

---

### U5. Pipeline connectivity: wire idxstats, fix factorbook and db paths

**Goal:** Ensure `idxstats` runs as part of QC, give `factorbook_logo` config-parameterized inputs, and make the database path portable.

**Requirements:** Addresses F1, F2, F3.

**Dependencies:** None.

**Files:**
- `workflow/rules/common.smk` (wire `idxstats` into `qc_done`)
- `workflow/rules/motifs.smk` (factorbook inputs)
- `workflow/rules/db.smk` (db_path)
- `config.yaml` (add factorbook config keys)
- `config/config.yaml` (add factorbook config keys)
- `workflow/rules/common.smk` (add factorbook keys to `_KNOWN_CONFIG_KEYS`)

**Approach:**

**F1 — Wire `idxstats` into `qc_done`:**
In `common.smk`, add idxstats outputs to the `rule qc_done:` input:
```python
rule qc_done:
    input:
        OUT + "/stats/report.csv",
        OUT + "/stats/report.html",
        OUT + "/multiqc_report.html",
        expand(OUT + "/stats/{sample}.idxstats.txt", sample=SAMPLES),
```
`idxstats` runs on all samples (treatment + control), which matches `sample=SAMPLES`.

**F2 — Factorbook config parameterization:**
Add two optional config keys to `config.yaml` and `config/config.yaml`:
```yaml
factorbook:
  tsv: "factorbook/factorbook_chip_seq_meme_motifs.tsv"
  meme: "factorbook/factorbook_chip_seq_meme_motif_catalog.meme"
```
In `motifs.smk:factorbook_logo`, change the `input:` to:
```python
input:
    tsv  = config.get("factorbook", {}).get("tsv", "factorbook/factorbook_chip_seq_meme_motifs.tsv"),
    meme = config.get("factorbook", {}).get("meme", "factorbook/factorbook_chip_seq_meme_motif_catalog.meme"),
```
The default values preserve current behavior. Users who store the factorbook files elsewhere can override in their `config.yaml`. Add `"factorbook"` to `_KNOWN_CONFIG_KEYS` in `common.smk`.

**F3 — Portable db_path:**
The current expression `os.path.dirname(workflow.basedir)` resolves to the parent of `workflow/` — i.e., the repository root — and places `pipeline_db.db` there. This is layout-dependent. Replace with an explicit config key with a sensible default:
```python
# In db.smk params
db_path = config.get("db_path", os.path.join(os.path.dirname(workflow.basedir), "pipeline_db.db")),
```
Add `db_path: null` (commented out) to `config.yaml` as a user-overridable option, and add `"db_path"` to `_KNOWN_CONFIG_KEYS`.

**Patterns to follow:** `config.get("gene_annotation")` pattern for optional config keys with defaults (`motifs.smk`, `annotation.smk`). `common.smk:66` shows the `workflow.basedir` idiom (`SCRIPTS = os.path.join(workflow.basedir, "scripts")`); F3 wraps the equivalent path with `config.get("db_path", os.path.join(os.path.dirname(workflow.basedir), "pipeline_db.db"))` so users can override it.

**Test scenarios:**
- After wiring idxstats: `snakemake --dry-run qc_done` shows `idxstats` in the planned jobs
- Factorbook: `snakemake --dry-run` with `factorbook.tsv` pointing to a nonexistent path fails with a clear missing-input error (not a KeyError in the Snakefile)
- db_path: set `db_path: /tmp/test.db` in config.yaml and verify `update_db.params.db_path` in the dry-run output reflects the override

**Verification:** `snakemake --dry-run --profile profiles/local --configfile config.yaml` shows `idxstats` in the job list and the database path resolves correctly.

---

## Risks & Dependencies

| Risk | Mitigation |
|------|-----------|
| U2 (SLURM DRY): users who have `slurm_partition` in their config.yaml will get an "Unrecognized config key" error after `_KNOWN_CONFIG_KEYS` update | Add a clear migration note in the README and a helpful error hint in `common.smk` similar to the existing `input_control` hint |
| U3 (profile paths): `$HOME` expansion in apptainer-args is shell-expanded by the calling process (Snakemake runs apptainer via `shell=True`), so it is reliable on all POSIX systems — no Apptainer version dependency | Add a comment in the local profile: `# $HOME is expanded by the shell before apptainer sees it` |
| U5 (idxstats wire): adds N more jobs to every run; may increase runtime if SLURM job queuing is a bottleneck | `idxstats` is fast (seconds per sample); the impact is negligible |
| U4 (stale key removal): users who had custom values for now-removed keys (e.g., custom `bamcoverage` resource) won't get a warning — the key will simply be ignored | The `_KNOWN_CONFIG_KEYS` validator in `common.smk` will catch any removed key that users still have in their config.yaml, producing a clear error |

---

## Sources & Research

Snakemake 8 best practices referenced:
- [Snakemake docs: Profiles](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles) — `default-resources` as the home for cluster-wide resource defaults
- [Snakemake docs: Rules](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html) — `params:` should not use lambdas when wildcards are not needed
- [Snakemake style guide](https://snakemake.readthedocs.io/en/stable/snakefiles/best_practices.html) — `set -euo pipefail` recommended for all multi-command shell blocks
