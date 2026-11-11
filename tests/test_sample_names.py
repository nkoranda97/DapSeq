"""
Tests for workflow/scripts/sample_names.py: the regex alternations used in
wildcard constraints, and the sample names the pipeline must reject.
"""

import re

import sample_names as m


def test_sample_regex_matches_names_literally():
    rx = m.sample_regex(["TF+1", "TF(2)", "plain"])
    assert re.fullmatch(rx, "TF+1")
    assert re.fullmatch(rx, "TF(2)")
    assert not re.fullmatch(rx, "TFF1")      # '+' must not act as a quantifier


def test_sample_regex_for_no_samples_never_matches():
    assert not re.fullmatch(m.sample_regex([]), "")
    assert not re.fullmatch(m.sample_regex([]), "anything")


def test_names_with_dot_slash_or_whitespace_are_rejected():
    errors = m.invalid_sample_names(["ok_1", "rep.1", "a/b", "has space", "tab\tx"])
    named = " ".join(errors)
    for bad in ("rep.1", "a/b", "has space"):
        assert repr(bad) in named
    assert "ok_1" not in named
    assert len(errors) == 4


def test_only_letters_digits_underscore_plus_and_minus_are_allowed():
    ok = ["TF_A", "TF+1", "input-B", "ARF5_rep2", "_x"]
    assert m.invalid_sample_names(ok) == []
    # shell metacharacters pass Snakemake's matching but break unquoted shell rules
    bad = ["TF(2)", "a;b", "x$y", "q'r", "star*", "c:d"]
    errors = m.invalid_sample_names(bad)
    assert len(errors) == len(bad)
    assert all("letters, digits, '_', '+' and '-'" in e for e in errors)


def test_name_starting_with_minus_is_rejected():
    # MACS3 would read "-n -TF" as an option
    errors = m.invalid_sample_names(["-TF"])
    assert len(errors) == 1 and "may not start with '-'" in errors[0]


def test_numeric_sample_name_gets_a_message_not_a_type_error():
    # YAML reads an unquoted key such as 1234: as an int
    errors = m.invalid_sample_names([1234, "TF_A"])
    assert len(errors) == 1
    assert "1234" in errors[0] and "quote it" in errors[0]


def test_sample_named_like_a_controls_self_call_is_rejected():
    errors = m.control_name_collisions(["input_A", "input_A_control", "TF_control"], ["input_A"])
    assert len(errors) == 1 and "input_A_control" in errors[0]
