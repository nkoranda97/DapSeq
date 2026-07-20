---
title: "feat: Per-sample control assignment for multi-experiment runs"
type: feat
status: completed
date: 2026-07-20
depth: standard
---

# feat: Per-sample control assignment for multi-experiment runs

## Summary

Today the pipeline has a **single, run-wide control**. The top-level `control:` field
names one sample (or a list that gets merged into `merged_control`), and **every**
treatment sample uses that same control as its MACS3 `-c` background and its bamCompare
`-b2`. That makes it impossible to run two independent experiments — each with its own
input/control — in the same Snakemake invocation.

This plan moves control assignment to **per-sample**: each sample carries an optional
`control:` key naming the control sample it should be called against. A "control sample"
is any sample referenced by another sample's `control:`. Treatments with different
controls coexist in one run, so several experiments execute together. The constraint
"exactly one control per sample" is preserved — `control:` is a scalar, not a list.

The top-level `control:` field is **removed** (breaking change). Old configs raise a
clear migration error. Because there is no longer a run-wide control and one control per
sample, the multi-control **merge** feature (`merged_control`) is dropped.

Control samples keep their current behavior: they are peak-called against themselves and
appear as flagged rows in the report (see origin: `docs/plans/2026-07-14-002-feat-control-self-macs-report-plan.md`).

---

## Problem Frame

Control is currently a **run-level scalar**, resolved once in
[workflow/rules/common.smk](workflow/rules/common.smk#L59-L95) into a single `CONTROL`
name that all treatments share:

- [workflow/rules/peaks.smk:13-16](workflow/rules/peaks.smk#L13-L16) — `rule macs3` resolves
  `control_bam` to `{CONTROL}.bam` for every treatment.
- [workflow/rules/filter.smk:26-32](workflow/rules/filter.smk#L26-L32) — `rule bamcompare`
  is gated on the global `if CONTROL:` and uses the fixed `{CONTROL}.bam` as `-b2`.
- [workflow/rules/filter.smk:1-23](workflow/rules/filter.smk#L1-L23) — `rule merge_control`
  merges a list of controls into one `merged_control.bam`.
- [workflow/rules/db.smk:11](workflow/rules/db.smk#L11) and
  [workflow/scripts/update_db.py:240](workflow/scripts/update_db.py#L240) — `control` is a
  single run-level value written identically onto every DB row.
- [workflow/scripts/report.py:325-334](workflow/scripts/report.py#L325-L334) —
  `read_control` reads the top-level `control` field from the config snapshot.

The single lever that couples everything is the run-wide `CONTROL`. Replacing it with a
per-sample map (`SAMPLE_CONTROL`) and a resolver (`control_for`) is the core change; every
consumer above then reads the map instead of the scalar.

---

## Requirements

- **R1** — A sample may declare exactly one control via a per-sample `control:` key naming
  another sample.
- **R2** — Multiple distinct controls can be used in a single run; treatments grouped under
  different controls all execute together (multi-experiment runs).
- **R3** — A sample with no `control:` is peak-called with no `-c` background (current
  no-control behavior preserved).
- **R4** — Control samples (any sample referenced as a control) are peak-called against
  themselves and appear as flagged report rows, exactly as today.
- **R5** — The top-level `control:` field is removed; configs that still set it fail fast
  with a migration message pointing to the per-sample key.
- **R6** — Each DB row records the control used for **that** sample (per-sample), not a
  single run-wide value.
- **R7** — Invalid control references (unknown sample, control with no `r1`, a sample used
  as both a treatment-with-its-own-control and a control for others) fail fast with a clear
  error.

---

## Key Technical Decisions

- **Per-sample `control:` scalar key** (chosen over top-level experiment blocks). It fits
  the existing `samples:` dict, keeps sample definitions in one place, and mirrors the
  recent per-sample `experiment_date` / `gdna_batch` change
  (commit `dfe0e40`). One control per sample is enforced by the value being a scalar.
- **Remove the top-level `control:` field** (chosen over keeping it as a run-wide default).
  Single source of truth; no ambiguity between a run default and per-sample overrides. Cost:
  breaking change (mitigated by a migration error) and loss of the `merged_control` merge
  feature, which has no meaning under one-control-per-sample.
- **`SAMPLE_CONTROL` map + `control_for(wc)` resolver** as the single new abstraction in
  `common.smk`. All rules read through it, so the peaks/bamcompare/db/report call sites each
  shrink to a lookup.
- **Control set is derived, not declared.** `CONTROL_SAMPLES = {sample referenced by any
  other sample's control}`. No separate list to keep in sync.
- **`is_treatment` / report-flagging semantics unchanged.** A control is any referenced
  sample; treatments are the rest. This keeps `REPORT_SAMPLES` ordering (controls last) and
  the control-row badge behavior intact.

---

## High-Level Technical Design

Control resolution today vs. after this change:

```
BEFORE (run-wide)                        AFTER (per-sample)
------------------                       ------------------
config.control ─┐                        samples:
                │                          TF1: {control: inA}   ┐
   ┌────────────▼───────────┐             TF2: {control: inA}    │ experiment A
   │ CONTROL (single name)   │             TF3: {control: inB}   ┘ experiment B
   └────────────┬───────────┘             inA: {}   inB: {}
                │ used by ALL treatments
                ▼                          SAMPLE_CONTROL = {TF1:inA, TF2:inA, TF3:inB}
   macs3 -c {CONTROL}                      CONTROL_SAMPLES = {inA, inB}
   bamcompare -b2 {CONTROL}                control_for(wc) → per-sample bam
   merged_control (if list)                (merge feature removed)
```

Directional only — not implementation specification. Control samples (`inA`, `inB`) are
still peak-called against themselves and flagged in the report.

---

## Implementation Units

### U1. Config schema: per-sample `control:`, remove top-level field

**Goal:** Express control assignment per sample and drop the run-wide field, with example
configs that demonstrate a two-experiment run.

**Requirements:** R1, R2, R3, R5

**Dependencies:** none

**Files:**
- [config.yaml](config.yaml) — user template
- [config/config.yaml](config/config.yaml) — pipeline defaults

**Approach:** Under each entry in `samples:`, add an optional `control:` key (scalar sample
name, `null` to omit) alongside `r1`/`r2`/`experiment_date`/`gdna_batch`. Remove the
top-level `control:` block and its explanatory comment. Update the template to show two
experiments sharing distinct controls (e.g., `TF1`/`TF2` → `input_A`, `TF3` → `input_B`).
Preserve the two-config-file split (see project memory: root template vs. `config/config.yaml`
defaults).

**Patterns to follow:** the per-sample `experiment_date`/`gdna_batch` layout already present
in both files.

**Test scenarios:** `Test expectation: none -- config/example data only; behavior is
exercised by U2 validation tests.`

**Verification:** Both YAML files parse; no top-level `control:` remains; each sample block
shows the optional `control:` key with a comment.

---

### U2. Control derivation and validation in `common.smk`

**Goal:** Replace the run-wide `CONTROL` scalar with a per-sample map and resolver, derive
the control set, and validate references. Add the top-level-`control` migration error.

**Requirements:** R1, R2, R3, R4, R5, R7

**Dependencies:** U1

**Files:**
- [workflow/rules/common.smk](workflow/rules/common.smk)

**Approach:**
- Remove `"control"` from `_KNOWN_CONFIG_KEYS` (line 8) — it is no longer a top-level key.
- In the unexpected-keys block (lines 14-30), add a hint when `control` appears at top level:
  point users to the per-sample `samples: <name>: control:` key (mirror the existing
  `experiment_date`/`gdna_batch` hint at lines 25-29).
- Replace the control-resolution block (lines 59-95): build
  `SAMPLE_CONTROL = {sample: control_name}` from each sample's `control:` value (skip
  `null`); derive `CONTROL_SAMPLES = sorted(set(SAMPLE_CONTROL.values()))`;
  `TREATMENT_SAMPLES = [s for s in SAMPLES if s not in CONTROL_SAMPLES]`.
- Keep `REPORT_SAMPLES = TREATMENT_SAMPLES + CONTROL_SAMPLES` (controls last) and the
  `CONTROL_SAMPLE_CONSTRAINT` / `TREATMENT_SAMPLE_CONSTRAINT` / `REPORT_SAMPLE_CONSTRAINT`
  derivations unchanged in shape.
- Delete the `CONTROL` scalar and any `merged_control` reference.
- Add `control_for(wc)` returning `SAMPLE_CONTROL.get(wc.sample)` (or `None`).
- Add a `CONTROL_TREATMENT_SAMPLES` constraint helper (samples that *have* a control) for the
  bamcompare wildcard constraint in U4.
- Validation (raise `ValueError` with a clear message): each referenced control must be a key
  in `SAMPLES` with a valid `r1`; a sample that is referenced as a control must not itself
  declare a `control:` (a control can't be a treatment); a sample may not name itself as its
  own control.

**Patterns to follow:** existing `ValueError` raising style and the `get_r1` helper in the
same file; wildcard-constraint `"(?!)"` empty-set idiom (lines 93-95).

**Test scenarios:**
- Covers R1. A single sample with `control: inA` → `SAMPLE_CONTROL == {"TF1": "inA"}`,
  `CONTROL_SAMPLES == ["inA"]`, `TREATMENT_SAMPLES == ["TF1"]`.
- Covers R2. Two samples with different controls → both controls appear in `CONTROL_SAMPLES`;
  both treatments in `TREATMENT_SAMPLES`; `REPORT_SAMPLES` orders controls last.
- Covers R3. A sample with `control: null` → not in `SAMPLE_CONTROL`; `control_for` returns
  `None`.
- Covers R7. `control:` naming an unknown sample → `ValueError` naming the offending sample.
- Covers R7. `control:` naming a sample whose `r1` is unset → `ValueError`.
- Covers R7. A control sample that itself declares a `control:` → `ValueError`.
- Covers R7. A sample naming itself → `ValueError`.
- Covers R5. Top-level `control:` present → `ValueError` with the per-sample migration hint.
- `control_for` returns the correct control name for a treatment and `None` for a control
  sample.

**Note on testability:** `common.smk` runs inside Snakemake's namespace. If the derivation
logic is not directly importable, extract the pure map/validation logic into a small helper
function (in `common.smk` or a tiny importable module under `workflow/scripts/`) so the
scenarios above can be unit-tested without a full Snakemake run. Decide the exact seam during
implementation.

**Verification:** A dry run (`snakemake -n`) with a two-experiment config resolves without
error; malformed configs raise the specified messages.

---

### U3. Per-sample control in MACS peak calling

**Goal:** `rule macs3` uses each sample's assigned control as the `-c` background.

**Requirements:** R1, R3, R4

**Dependencies:** U2

**Files:**
- [workflow/rules/peaks.smk](workflow/rules/peaks.smk)

**Approach:** Replace the `control_bam` input lambda (lines 13-16): for a control sample,
keep self (`{wc.sample}.bam`); for a treatment, resolve via `control_for(wc)` →
`{control}.bam`, or `[]` when the sample has no control. `params.ctrl` and the `--nomodel`
gate (control-against-itself) are unchanged in logic. `rule macs3_control` is unchanged.

**Patterns to follow:** the existing input-lambda / `params.ctrl` idiom in the same rule.

**Test scenarios:** `Test expectation: none -- pure input-wiring change; correctness covered
by the U2 resolver tests and the integration dry run in U8.`

**Verification:** Dry run shows each treatment's `macs3` job depends on the correct control
BAM; a no-control sample emits no `-c`.

---

### U4. Per-sample control in bamCompare; remove `merge_control`

**Goal:** `rule bamcompare` runs per sample against that sample's control; delete the
now-meaningless control-merge rule.

**Requirements:** R1, R2, R6-adjacent

**Dependencies:** U2

**Files:**
- [workflow/rules/filter.smk](workflow/rules/filter.smk)

**Approach:**
- Delete `rule merge_control` (lines 1-23) — no run-wide merge under one-control-per-sample.
- Replace the global `if CONTROL:` gate: `rule bamcompare` now applies to the set of samples
  that **have** a control. Add a `wildcard_constraints: sample = CONTROL_TREATMENT_SAMPLES`
  (from U2) so it only fires for samples with an assigned control. Resolve `control_bam` /
  `control_bai` inputs via `control_for(wc)`.
- Remove the `merge_control` entry from `resources:` in [config/config.yaml](config/config.yaml).
- Also remove `"merged_control"` handling anywhere else it lingers.

**Note:** the `.peaks.bw` output of `rule bamcompare` is not consumed by any target rule
today (no `all`/`mapped`/`peaked` target requires it). This change keeps the rule correct but
does not add it to a default target — call out in the PR that `.peaks.bw` remains build-on-request.

**Patterns to follow:** wildcard-constraint gating used by `rule macs3_control`
([peaks.smk:60-61](workflow/rules/peaks.smk#L60-L61)).

**Test scenarios:** `Test expectation: none -- rule-wiring change; exercised by the U8 dry
run. Verify no rule references merge_control or merged_control.`

**Verification:** `grep -rn "merge_control\|merged_control" workflow/` returns nothing; dry
run builds `{sample}.peaks.bw` only for samples with a control.

---

### U5. Per-sample control in the results database

**Goal:** Record the control used for each sample on its own DB row instead of a single
run-wide value.

**Requirements:** R6

**Dependencies:** U2

**Files:**
- [workflow/rules/db.smk](workflow/rules/db.smk)
- [workflow/scripts/update_db.py](workflow/scripts/update_db.py)

**Approach:**
- In `db.smk`, replace `control = config.get("control")` (line 11) with the per-sample map,
  e.g. `sample_control = SAMPLE_CONTROL`.
- In `update_db.py`, remove `"control"` from the run-level `shared` dict
  ([line 240](workflow/scripts/update_db.py#L240)) and instead set it per-sample when building
  each row: a treatment records its control name; a control sample records `""`. Keep
  `"control"` in `COLS` (line 30) — its position and column name are unchanged, only its
  source moves from run-level to per-row.

**Patterns to follow:** how per-sample fields (`sample`, `r1`, `r2`, `is_treatment`) are
already set per-row in `main()` vs. the `shared` dict; mirrors the per-sample
`experiment_date`/`gdna_batch` handling.

**Test scenarios:**
- Covers R6. Build rows for a two-experiment config → each treatment row's `control` equals
  its assigned control; each control row's `control` is `""`.
- A no-control treatment row has `control == ""`.
- `COLS` still contains `control` exactly once and the column ordering is unchanged (guards
  the existing DB schema).

**Verification:** `tests/test_update_db.py` passes; inspecting the DB after a run shows
distinct `control` values across experiments.

---

### U6. Report control derivation from per-sample config

**Goal:** `report.read_control` derives control samples from per-sample `control:` keys in
the config snapshot, since the top-level field is gone.

**Requirements:** R4, R5

**Dependencies:** U1

**Files:**
- [workflow/scripts/report.py](workflow/scripts/report.py)

**Approach:** Rewrite `read_control` (lines 325-334): load the snapshot, iterate
`samples:` entries, collect each non-null `control:` value into a de-duplicated, ordered
list. Return `[]` when the snapshot is missing/unparseable or no sample sets a control. The
Snakemake path (`render_report_html.py` receiving `control_samples=CONTROL_SAMPLES`) is
unchanged; only the snapshot-fallback derivation changes.

**Patterns to follow:** existing `_load_config_snapshot` best-effort loader and the
list-normalization style in the current `read_control`.

**Test scenarios:**
- Covers R4. Snapshot with `samples: {TF1: {control: inA}, TF2: {control: inB}, inA: {}, inB: {}}`
  → `read_control` returns `["inA", "inB"]` (order stable, de-duplicated when shared).
- Two treatments sharing one control → returns that control once.
- Snapshot with no `control:` on any sample → returns `[]`.
- Missing/unparseable snapshot → returns `[]`.

**Verification:** Report flags control rows correctly for a multi-experiment run; the
existing `test_write_html_flags_control_row` still passes with derived control names.

---

### U7. Documentation

**Goal:** README reflects per-sample control, the removed top-level field, and multi-experiment
usage.

**Requirements:** R1, R2, R3, R5

**Dependencies:** U1-U6

**Files:**
- [README.md](README.md)

**Approach:** Update the sample-config example (around line 68) to show per-sample `control:`
and a two-experiment layout. Replace the top-level `control` row in the config parameter table
(line 87) with a per-sample `control` entry documented as `name=default` per the project's
param-doc format (see feedback memory). Add a short "Running multiple experiments in one run"
note. Document the removal of the top-level field and the migration error. Note that
control-merging is no longer supported.

**Patterns to follow:** existing README config-table format and the `name=default` param
convention.

**Test scenarios:** `Test expectation: none -- documentation.`

**Verification:** README examples parse as valid YAML and match the U1 template; no stale
references to top-level `control:` or `merged_control`.

---

### U8. Test suite updates and integration dry run

**Goal:** Update fixtures and assertions for the per-sample model; confirm an end-to-end dry
run of a two-experiment config.

**Requirements:** R1-R7

**Dependencies:** U2-U6

**Files:**
- [tests/test_collect_stats_and_report.py](tests/test_collect_stats_and_report.py)
- [tests/test_update_db.py](tests/test_update_db.py)
- [tests/conftest.py](tests/conftest.py)

**Approach:**
- Replace the `read_control` snapshot tests (lines 542-552: `test_read_control_scalar`,
  `test_read_control_list`, `test_read_control_null`) with per-sample-derivation tests from
  U6. The old top-level `control:` snapshot format no longer exists.
- Keep the control-flagging tests (they take `control_samples` directly) but ensure fixture
  sample names reflect referenced-control naming.
- Update `test_update_db.py` per-sample-control assertions from U5.
- Add/adjust any fixture config in `conftest.py` to include per-sample `control:` and a
  two-experiment shape.
- Run `pixi run` Snakemake dry run (`-n`) against a two-experiment config to confirm the DAG
  resolves (per project memory: all execution via `pixi run`, never bare `python`).

**Test scenarios:** the enumerated scenarios in U2, U5, and U6 live here. Additionally:
- Covers R2 (integration). A two-experiment dry run resolves without error and schedules the
  correct per-treatment `macs3`/`bamcompare` jobs.
- Covers R5 (integration). A dry run with a top-level `control:` config fails with the
  migration message.

**Verification:** `pixi run pytest` passes; the two dry runs behave as specified.

---

## Scope Boundaries

**In scope:** per-sample control assignment; removal of the top-level `control:` field and
migration error; per-sample control in MACS, bamCompare, DB, and report; removal of the
`merge_control` rule; docs and tests.

**Out of scope / non-goals:**
- Multiple controls **per single sample** — explicitly not supported (one control per sample).
- Any change to how controls are peak-called against themselves or flagged in the report
  (behavior preserved from `docs/plans/2026-07-14-002-feat-control-self-macs-report-plan.md`).
- Experiment-level grouping metadata (naming/labeling experiments as first-class objects in
  the report or DB).

### Deferred to Follow-Up Work
- Adding `.peaks.bw` to a default target rule, if desired (currently build-on-request; this
  plan preserves that).
- An explicit per-experiment label column in the report/DB, if experiment grouping later needs
  to be surfaced to users.

---

## Risks & Dependencies

- **Breaking config change (R5).** Existing user configs with top-level `control:` will fail.
  Mitigation: fast, explicit migration error in `common.smk` (U2) pointing to the per-sample
  key; README migration note (U7).
- **`common.smk` testability.** The derivation logic runs in Snakemake's namespace; the U2
  note flags extracting a pure helper if direct unit testing is otherwise hard. Low risk —
  the logic is small.
- **DB schema stability.** U5 moves `control` from run-level to per-row but must keep the
  column name and position (guarded by a U5 test) so existing DB readers/queries keep working.
- **Loss of control-merging.** Anyone relying on `control: [a, b]` merge behavior loses it.
  Acceptable per the chosen design; called out in README.

---

## Sources & Research

- Local code: `workflow/rules/common.smk`, `peaks.smk`, `filter.smk`, `db.smk`,
  `workflow/scripts/update_db.py`, `report.py`.
- Prior plan (control-self behavior to preserve):
  `docs/plans/2026-07-14-002-feat-control-self-macs-report-plan.md`.
- Prior plan (per-sample field precedent):
  `docs/plans/2026-07-15-001-feat-per-sample-experiment-metadata-plan.md`.
- No external research required — fully internal Snakemake refactor with strong local patterns.
