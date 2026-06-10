---
title: "refactor: Rewrite README to reflect current pipeline state"
date: 2026-06-05
status: active
type: refactor
---

# refactor: Rewrite README to reflect current pipeline state

## Summary

The current README documents tool parameters in exhaustive detail, references only Bowtie2 for alignment, uses a stale `mapq` default (10 vs. actual 30), and is missing `chrom_filter` and `base_colors` (MEME). The goal is a minimal, accurate document: setup, running, required config fields, and notable options only.

## Problem Frame

Users reading the README get incorrect defaults, no guidance on the aligner choice, and a wall of parameter docs that duplicates upstream tool documentation. The README should orient a new user and point to `config/config.yaml` for defaults — not reproduce every flag.

## Requirements

- R1: Reflect both supported aligners (`bowtie2`, `bwa_mem2`) in setup and index build sections
- R2: Remove verbose per-parameter documentation for all tools
- R3: Fix `mapq` default (currently shown as 10; actual default is 30)
- R4: Document `chrom_filter` config option
- R5: Keep the multi-user `--directory` warning (it is non-obvious and load-bearing)
- R6: Keep content minimal — link to tool docs rather than reproducing them

## Key Technical Decisions

| Decision | Choice | Rationale |
|---|---|---|
| Config depth | Required + notable fields only | User said concise; full defaults live in `config/config.yaml` |
| Per-tool parameter tables | Remove entirely | Duplicate upstream docs; user said remove irrelevant content |
| `base_colors` (MEME) | Include as commented example | Documented in config but not README; low noise |
| HOMER section | One line — link only | Existing README had just a link stub anyway |

## Implementation Units

### U1. Rewrite README.md

**Goal:** Replace the README with a minimal, accurate document that reflects the current pipeline.

**Requirements:** R1, R2, R3, R4, R5, R6

**Files:**
- `README.md`

**Approach:**

Structure:

```
# DAP-seq Pipeline
[1-2 sentence description]

## Setup
[pixi install command]

## Running
[snakemake command with --directory, --configfile, --snakefile]
[note about --directory matching output_dir]
[note about multi-user shared installation]

## Configuration
### Required fields
[table or bullet list: samples, input_control, output_dir, genome_ref, genome_size, gene_annotation, aligner, slurm_partition, slurm_account]

### Notable options
- aligner: bowtie2 (default) | bwa_mem2
- samtools.mapq: 30 (default) — MAPQ filter applied after alignment
- chrom_filter: [] — list of chromosomes to exclude before MEME analysis (e.g. chrEBV)
- macs3.min_foldch: 2.0 — peak fold-change filter
- meme.base_colors — optional hex color overrides for sequence logos

### Extra arguments
[one-liner showing extra: "--param value" pattern]

## Index build (admin, one-time per genome)
### Bowtie2
[command]
### bwa-mem2
[command]
[note about racing on first run]
```

Remove entirely:
- Detailed `bbduk` parameter table
- Detailed `samtools view` parameter table  
- Detailed `bamcoverage` / `bamcompare` parameter tables
- Detailed `macs3` parameter table
- Detailed `meme` / `fimo` parameter tables
- The HOMER section (already a stub — just link in a tools list if mentioned)

**Test scenarios:** None — documentation change with no behavioral impact.

**Verification:** Read the resulting README cold; a new user should be able to run the pipeline without consulting any other file.

## Scope Boundaries

### Deferred to Follow-Up Work
- Adding a pipeline overview diagram or DAG description
- Documenting the output directory structure
