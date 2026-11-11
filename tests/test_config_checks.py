"""
Tests for workflow/scripts/config_checks.py: startup checks on the shape of
the merged config.
"""

import gzip

import pytest
import yaml

import config_checks as m


# ── fill_empty_sections ───────────────────────────────────────────────────────

def test_section_with_all_keys_commented_out_gets_its_defaults():
    # what Snakemake's merge leaves for a user file containing "factorbook:\n#  tsv: x"
    config = {"factorbook": yaml.safe_load("factorbook:\n#  tsv: x\n")["factorbook"],
              "threads": 4}
    defaults = {"factorbook": {"tsv": None, "meme": None}, "threads": 8}
    assert m.fill_empty_sections(config, defaults) == ["factorbook"]
    assert config == {"factorbook": {"tsv": None, "meme": None}, "threads": 4}


def test_set_or_absent_sections_are_left_alone():
    config = {"bbduk": {"k": 23}, "gene_annotation": None}
    defaults = {"bbduk": {"k": 21}, "gene_annotation": None, "homer": {"extra": ""}}
    assert m.fill_empty_sections(config, defaults) == []
    assert config == {"bbduk": {"k": 23}, "gene_annotation": None}


def test_nested_section_with_all_keys_commented_out_gets_its_defaults():
    # "resources:\n  meme:\n#    mem_mb: 32000" merges to resources.meme = None
    config = {"resources": {"meme": None, "fimo": {"mem_mb": 1, "runtime": 2}}}
    defaults = {"resources": {"meme": {"mem_mb": 16000, "runtime": 120},
                              "fimo": {"mem_mb": 16000, "runtime": 60}}}
    assert m.fill_empty_sections(config, defaults) == ["resources.meme"]
    assert config["resources"]["meme"] == {"mem_mb": 16000, "runtime": 120}
    assert config["resources"]["fimo"] == {"mem_mb": 1, "runtime": 2}


# ── max_frags ─────────────────────────────────────────────────────────────────

@pytest.mark.parametrize("value", [None, "", 0])
def test_unset_max_frags_means_no_subsampling(value):
    assert m.max_frags(value) is None


@pytest.mark.parametrize("text, expected", [
    ("max_frags: 10000000", 10_000_000),
    ('max_frags: "1e7"', 10_000_000),
    ("max_frags: 1.0e+7", 10_000_000),
])
def test_max_frags_accepts_whole_numbers(text, expected):
    assert m.max_frags(yaml.safe_load(text)["max_frags"]) == expected


@pytest.mark.parametrize("text", ["max_frags: 10,000,000", "max_frags: 10M", "max_frags: -5",
                                  "max_frags: 2.5"])
def test_max_frags_typos_are_rejected(text):
    # mawk would read "10,000,000" and "10M" as 10 and keep ~10 fragments per sample
    with pytest.raises(ValueError, match="bbduk.max_frags"):
        m.max_frags(yaml.safe_load(text)["max_frags"])


# ── sample_config_errors ──────────────────────────────────────────────────────

def test_known_sample_keys_pass():
    cfg = {"TF_A": {"r1": "a", "r2": "b", "control": "in", "tf": "ARF5",
                    "experiment_date": "2026", "gdna_batch": 3}}
    assert m.sample_config_errors(cfg) == []


def test_misspelled_sample_key_is_named_with_a_suggestion():
    errors = m.sample_config_errors({"TF_A": {"r1": "a", "contol": "input_A"}})
    assert len(errors) == 1
    assert "'TF_A' sets unknown key 'contol' (did you mean 'control'?)" in errors[0]


def test_sample_with_no_settings_is_reported():
    errors = m.sample_config_errors(yaml.safe_load("TF_A:\nTF_B: {r1: x}\n"))
    assert len(errors) == 1 and "sample 'TF_A' has no settings" in errors[0]


# ── macs3_format_override ─────────────────────────────────────────────────────

@pytest.mark.parametrize("value", [None, ""])
def test_unset_format_chooses_per_call(value):
    assert m.macs3_format_override(value, {"SE1"}) is None


def test_format_is_case_insensitive():
    assert m.macs3_format_override("bam", {"SE1"}) == "BAM"
    assert m.macs3_format_override("Bampe", set()) == "BAMPE"


def test_unknown_format_is_rejected():
    with pytest.raises(ValueError, match="must be null .* BAM or BAMPE, got 'BAMP'"):
        m.macs3_format_override("BAMP", set())


def test_bampe_with_single_end_samples_is_rejected():
    with pytest.raises(ValueError, match="single-end .*: SE1, SE2"):
        m.macs3_format_override("BAMPE", {"SE2", "SE1"})


# ── compressed_reference_error ────────────────────────────────────────────────

def test_gzipped_reference_is_rejected(tmp_path):
    path = tmp_path / "genome.fa.gz"
    with gzip.open(path, "wt") as fh:
        fh.write(">chr1\nACGT\n")
    error = m.compressed_reference_error(str(path))
    assert error and "gzip-compressed" in error and str(path) in error


def test_plain_reference_passes_even_with_a_gz_name(tmp_path):
    path = tmp_path / "genome.fa.gz"
    path.write_text(">chr1\nACGT\n")
    assert m.compressed_reference_error(str(path)) is None


def test_missing_reference_is_left_to_the_path_check(tmp_path):
    assert m.compressed_reference_error(str(tmp_path / "missing.fa")) is None
