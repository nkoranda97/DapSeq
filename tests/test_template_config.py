"""
The root config.yaml is the template users copy. It must hold only the
per-run fields, so a copy runs on the defaults in config/config.yaml for
every other setting.
"""

import copy
from pathlib import Path

import yaml
from snakemake.utils import update_config

ROOT = Path(__file__).parent.parent
TEMPLATE = yaml.safe_load((ROOT / "config.yaml").read_text())
DEFAULTS = yaml.safe_load((ROOT / "config" / "config.yaml").read_text())

PER_RUN = {"author", "samples", "output_dir", "genome_ref", "genome_size",
           "gene_annotation", "aligner"}


def test_template_sets_only_the_per_run_fields():
    assert set(TEMPLATE) == PER_RUN


def test_template_genome_size_is_blank():
    assert TEMPLATE["genome_size"] is None


def test_defaults_genome_size_is_blank():
    assert DEFAULTS["genome_size"] is None


def test_filled_in_template_runs_on_the_defaults():
    filled = copy.deepcopy(TEMPLATE)
    filled.update(output_dir="/out", genome_ref="/ref/genome.fa", genome_size=119000000)
    for sample in filled["samples"].values():
        sample["r1"] = "/fq/R1.fq.gz"

    merged = copy.deepcopy(DEFAULTS)
    update_config(merged, filled)   # how Snakemake layers --configfile over the defaults

    for key in set(DEFAULTS) - PER_RUN:
        assert merged[key] == DEFAULTS[key], key
    assert merged["meme"]["maxpeaks"] == 100 and merged["threads"] == 8
