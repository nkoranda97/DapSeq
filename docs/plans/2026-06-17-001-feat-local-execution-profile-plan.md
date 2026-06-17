---
title: "feat: Add local (no-SLURM) execution profile"
date: 2026-06-17
status: completed
type: feat
---

# feat: Add local (no-SLURM) execution profile

## Summary

Add a `profiles/local/` Snakemake profile so the pipeline can run on a local workstation without submitting jobs to SLURM. The local executor uses the same Apptainer container as the HPC run, keeping tool behavior identical. Only two files need to be created or modified: the new profile and the user-facing `config.yaml` template.

---

## Problem Frame

The pipeline's only profile (`profiles/default/`) hard-codes `executor: slurm` and an HPC-specific jobscript that loads the Apptainer module via `module load`. Local machines don't have SLURM or the `module` system, so there is currently no supported way to run the pipeline outside of an HPC cluster.

The `slurm_partition` and `slurm_account` resource keys in every rule are already null-safe: `config/config.yaml` defaults both to `null`, and the Snakemake local executor ignores unknown resource keys. No rule-level changes are required.

---

## Requirements

- R1: A local profile exists that runs all pipeline rules without submitting SLURM jobs.
- R2: Tool execution is identical to the HPC run (same Apptainer container, same tool versions).
- R3: HPC-specific container bind paths are not applied when running locally.
- R4: The user config template clearly signals that `slurm_partition`/`slurm_account` are optional for local runs.
- R5: README documents how to build the SIF locally and invoke the local profile.

---

## Key Technical Decisions

| Decision | Rationale |
|---|---|
| Omit `executor:` from local profile | Snakemake 8/9 defaults to local execution when no executor plugin is named; avoids activating `snakemake-executor-plugin-slurm`. |
| Keep `use-apptainer: true` in local profile | Bioinformatics tools (bowtie2, bbduk, macs3, meme, fimo, homer, bamCoverage) are only packaged in the Apptainer SIF, not in `pixi.toml`. Keeps tool versions identical to HPC. |
| Omit `apptainer-args` from local profile | The HPC bind paths (`/project/zhuzhuzhang`, `/project/gzy8899`) don't exist on local machines; omitting prevents Apptainer from erroring on missing directories. |
| No rule-level changes | `config/config.yaml` already defaults `slurm_partition: null` and `slurm_account: null`; the local executor silently ignores those resource keys. |
| `jobs: 4` default in local profile | A conservative parallel-job limit suited for a workstation; user can override with `--jobs N`. |

---

## Implementation Units

### U1. Create `profiles/local/config.yaml`

**Goal:** A minimal Snakemake profile for local execution that uses Apptainer but skips all SLURM settings.

**Requirements:** R1, R2, R3

**Dependencies:** none

**Files:**
- `profiles/local/config.yaml` *(create)*

**Approach:** Omit `executor:` (defaults to local), set `use-apptainer: true`, set `jobs: 4`. Do **not** include `jobscript:` or `apptainer-args:`. Copy the `default-resources` mem/runtime hints from the default profile so per-rule resource declarations still resolve.

**Patterns to follow:** `profiles/default/config.yaml` structure (subset).

**Test scenarios:**
- `pixi run snakemake --profile profiles/local --configfile config.yaml --dry-run` completes without error and shows no SLURM-related output.
- Running with a test dataset and the local profile produces the same output files as the default profile on HPC (validates container invocation is working).
- Running when `slurm_partition` and `slurm_account` are absent from the user config does not raise a KeyError.

**Verification:** Dry-run succeeds; a sample run reaches peak-calling stage without binding errors from Apptainer.

---

### U2. Update user config template (`config.yaml`)

**Goal:** Make the SLURM fields clearly optional so local users don't feel required to fill them in.

**Requirements:** R4

**Dependencies:** none

**Files:**
- `config.yaml` *(modify)*

**Approach:** Move `slurm_partition` and `slurm_account` under a `# ── SLURM (HPC only) ─` comment block and set their defaults to `null` with an inline note that they are ignored when using the local profile. No structural changes — just clarifying comments and default values.

**Test scenarios:**
- A config.yaml with both fields set to `null` runs successfully with `--profile profiles/local`.

**Verification:** No change in existing HPC behavior; local dry-run does not KeyError on the slurm fields.

---

### U3. Add local execution section to README

**Goal:** Document the one-time SIF build step and the local run command so users can discover and use the new profile.

**Requirements:** R5

**Dependencies:** U1

**Files:**
- `README.md` *(modify)*

**Approach:** Add a `## Running locally` section (or `### Local workstation` sub-section under the existing Running section) with:
1. One-time SIF build: `pixi run apptainer build apptainer_build/dapseq.sif apptainer_build/dapseq.def`
2. Run command: `pixi run snakemake --profile profiles/local --configfile /path/to/config.yaml`
3. Note that `slurm_partition`/`slurm_account` can be omitted from the local config.

**Test scenarios:**
- Test expectation: none — prose-only change.

**Verification:** README renders correctly; the build and run commands are copy-pasteable.

---

## Scope Boundaries

**In scope:**
- New `profiles/local/config.yaml`
- Config template comment clarifications
- README local-run documentation

**Deferred to Follow-Up Work:**
- Adding all bioinformatics tools to `pixi.toml` as an alternative to Apptainer (would enable `use-apptainer: false` local runs)
- A `pixi run` task shortcut for the local profile invocation
- Automated test that exercises the local profile end-to-end in CI

**Out of scope:**
- macOS or Windows support (Apptainer is Linux-only; `pixi.toml` already restricts to `linux-64`)
- Changing any rule logic or resource declarations
