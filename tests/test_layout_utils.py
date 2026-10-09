"""
Tests for workflow/scripts/layout_utils.py

Covers MACS3 format resolution per peak call, PE->BAM fallback detection,
r1/r2 lane-list normalisation and validation, and genome_size parsing.
"""

import pytest

import layout_utils as m


PE = {"TF_A", "input_A", "TF_B", "TF_D"}
CONTROLS = {"input_A", "input_B"}


# ── macs3_format ──────────────────────────────────────────────────────────────

def test_pe_treatment_with_pe_control_is_bampe():
    assert m.macs3_format("TF_A", "input_A", PE, CONTROLS, None) == "BAMPE"


def test_pe_treatment_with_se_control_falls_back_to_bam():
    assert m.macs3_format("TF_B", "input_B", PE, CONTROLS, None) == "BAM"


def test_se_treatment_without_control_is_bam():
    assert m.macs3_format("TF_C", None, PE, CONTROLS, None) == "BAM"


def test_pe_treatment_without_control_is_bampe():
    assert m.macs3_format("TF_D", None, PE, CONTROLS, None) == "BAMPE"


def test_control_self_call_follows_control_layout():
    assert m.macs3_format("input_A", None, PE, CONTROLS, None) == "BAMPE"
    assert m.macs3_format("input_B", None, PE, CONTROLS, None) == "BAM"


def test_override_wins_for_every_call():
    assert m.macs3_format("TF_C", None, PE, CONTROLS, "BAMPE") == "BAMPE"
    assert m.macs3_format("TF_A", "input_A", PE, CONTROLS, "BAM") == "BAM"


def test_empty_override_is_unset():
    assert m.macs3_format("TF_C", None, PE, CONTROLS, "") == "BAM"


# ── fallback_pairs ────────────────────────────────────────────────────────────

def test_fallback_pairs_lists_pe_treatments_with_se_controls():
    sample_control = {"TF_A": "input_A", "TF_B": "input_B"}
    assert m.fallback_pairs(sample_control, PE, None) == [("TF_B", "input_B")]


def test_fallback_pairs_ignores_se_treatments():
    assert m.fallback_pairs({"TF_C": "input_B"}, PE, None) == []


def test_fallback_pairs_empty_when_override_set():
    assert m.fallback_pairs({"TF_B": "input_B"}, PE, "BAMPE") == []


# ── lanes ─────────────────────────────────────────────────────────────────────

def test_as_list_normalises_null_path_and_list():
    assert m.as_list(None) == []
    assert m.as_list("a.fq.gz") == ["a.fq.gz"]
    assert m.as_list(["a.fq.gz", "b.fq.gz"]) == ["a.fq.gz", "b.fq.gz"]


def test_lane_count_mismatch_names_sample_and_counts():
    samples = {"TF_A": {"r1": ["l1_R1", "l2_R1"], "r2": "l1_R2"}}
    errors = m.lane_count_errors(samples, {"TF_A"})
    assert len(errors) == 1
    assert "TF_A" in errors[0]
    assert "2" in errors[0] and "1" in errors[0]


def test_matching_lane_counts_have_no_errors():
    samples = {
        "TF_A": {"r1": ["l1_R1", "l2_R1"], "r2": ["l1_R2", "l2_R2"]},
        "TF_B": {"r1": "R1", "r2": "R2"},
    }
    assert m.lane_count_errors(samples, {"TF_A", "TF_B"}) == []


def test_se_samples_are_never_flagged():
    samples = {"TF_C": {"r1": ["l1", "l2"], "r2": None}}
    assert m.lane_count_errors(samples, set()) == []


# ── parse_genome_size ─────────────────────────────────────────────────────────

@pytest.mark.parametrize("value, expected", [
    ("2.7e9", 2_700_000_000),
    (3_000_000_000, 3_000_000_000),
    ("3000000000", 3_000_000_000),
    ("1e5", 100_000),
    (" 1.2E8 ", 120_000_000),
])
def test_parse_genome_size_accepts_int_and_scientific(value, expected):
    result = m.parse_genome_size(value)
    assert result == expected
    assert isinstance(result, int)


@pytest.mark.parametrize("value", [None, "abc", "0", "-5", "1.5", True, ""])
def test_parse_genome_size_rejects_invalid(value):
    with pytest.raises(ValueError, match="genome_size"):
        m.parse_genome_size(value)
