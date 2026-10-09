"""
Tests for workflow/scripts/pipeline_version.py: the version label stored with
every results-DB row. git is replaced by a stub; no repository is needed.
"""

import subprocess
from types import SimpleNamespace

import pipeline_version as pv


def _git(monkeypatch, result=None, raises=None):
    calls = []

    def fake_run(cmd, **kwargs):
        calls.append(cmd)
        if raises is not None:
            raise raises
        return result

    monkeypatch.setattr(subprocess, "run", fake_run)
    return calls


def test_returns_git_describe_output(monkeypatch):
    calls = _git(monkeypatch, SimpleNamespace(returncode=0, stdout="7836373-dirty\n"))
    assert pv.pipeline_version("/repo") == "7836373-dirty"
    # Other lab members run from the shared checkout, which git would refuse
    # as "dubious ownership" without the safe.directory override.
    assert "safe.directory=*" in calls[0] and "/repo" in calls[0]


def test_git_missing_is_unknown(monkeypatch):
    _git(monkeypatch, raises=FileNotFoundError("git"))
    assert pv.pipeline_version("/repo") == "unknown"


def test_git_failure_is_unknown(monkeypatch):
    _git(monkeypatch, SimpleNamespace(returncode=128, stdout=""))
    assert pv.pipeline_version("/repo") == "unknown"


def test_git_timeout_is_unknown(monkeypatch):
    _git(monkeypatch, raises=subprocess.TimeoutExpired("git", 10))
    assert pv.pipeline_version("/repo") == "unknown"
