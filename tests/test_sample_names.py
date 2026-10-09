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


def test_sample_named_like_a_controls_self_call_is_rejected():
    errors = m.control_name_collisions(["input_A", "input_A_control", "TF_control"], ["input_A"])
    assert len(errors) == 1 and "input_A_control" in errors[0]
