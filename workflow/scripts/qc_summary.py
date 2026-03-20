"""
Aggregate per-sample QC statistics into a single TSV summary table.

Called via Snakemake's script: directive; uses the snakemake object for I/O.
"""

import os


def parse_kv_file(path):
    """
    Read a key<TAB>value file, skipping blank lines and lines starting with '#'.
    Returns a dict of {key: value}.
    """
    result = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t", 1)
            if len(parts) == 2:
                result[parts[0]] = parts[1]
    return result


def main():
    sm = snakemake  # noqa: F821 — injected by Snakemake

    samples           = list(sm.params.samples)
    treatment_samples = set(sm.params.treatment_samples)
    use_macs          = sm.params.use_macs
    use_gem           = sm.params.use_gem
    out_dir           = sm.params.out_dir

    # Base columns present for every sample
    base_cols = [
        "sample",
        "total_reads",
        "mapped_reads",
        "unmapped_reads",
        "mapping_rate",
        "unique_mapped",
        "multimapped",
        "multimapped_rate",
        "soft_clipped_reads",
        "properly_paired",
        "properly_paired_rate",
        "filtered_reads",
        "reads_in_peaks",
        "frip",
    ]
    optional_cols = []
    if use_macs:
        optional_cols += ["reads_in_peaks_macs", "frip_macs"]
    if use_gem:
        optional_cols += ["reads_in_peaks_gem", "frip_gem"]

    all_cols = base_cols + optional_cols

    rows = []
    for sample in samples:
        raw_path      = os.path.join(out_dir, f"{sample}.raw_alignment_stats.txt")
        filtered_path = os.path.join(out_dir, f"{sample}.filtered_stats.txt")

        raw      = parse_kv_file(raw_path)
        filtered = parse_kv_file(filtered_path)

        row = {
            "sample":              sample,
            "total_reads":         raw.get("total_reads",         "NA"),
            "mapped_reads":        raw.get("mapped_reads",        "NA"),
            "unmapped_reads":      raw.get("unmapped_reads",      "NA"),
            "mapping_rate":        raw.get("mapping_rate",        "NA"),
            "unique_mapped":       raw.get("unique_mapped",       "NA"),
            "multimapped":         raw.get("multimapped",         "NA"),
            "multimapped_rate":    raw.get("multimapped_rate",    "NA"),
            "soft_clipped_reads":  raw.get("soft_clipped_reads",  "NA"),
            # PE-only fields: absent for SE → "NA"
            "properly_paired":     raw.get("properly_paired",     "NA"),
            "properly_paired_rate": raw.get("properly_paired_rate", "NA"),
            "filtered_reads":      filtered.get("filtered_reads", "NA"),
        }

        if sample in treatment_samples:
            frip_path = os.path.join(out_dir, f"{sample}.frip.txt")
            frip      = parse_kv_file(frip_path)
            row["reads_in_peaks"] = frip.get("reads_in_peaks", "NA")
            row["frip"]           = frip.get("frip",           "NA")

            if use_macs:
                frip_macs_path = os.path.join(out_dir, f"{sample}.frip_macs.txt")
                frip_macs      = parse_kv_file(frip_macs_path)
                row["reads_in_peaks_macs"] = frip_macs.get("reads_in_peaks_macs", "NA")
                row["frip_macs"]           = frip_macs.get("frip_macs",           "NA")

            if use_gem:
                frip_gem_path = os.path.join(out_dir, f"{sample}.frip_gem.txt")
                frip_gem      = parse_kv_file(frip_gem_path)
                row["reads_in_peaks_gem"] = frip_gem.get("reads_in_peaks_gem", "NA")
                row["frip_gem"]           = frip_gem.get("frip_gem",           "NA")
        else:
            # Control sample — no peak-based metrics
            row["reads_in_peaks"] = "NA"
            row["frip"]           = "NA"
            if use_macs:
                row["reads_in_peaks_macs"] = "NA"
                row["frip_macs"]           = "NA"
            if use_gem:
                row["reads_in_peaks_gem"] = "NA"
                row["frip_gem"]           = "NA"

        rows.append(row)

    with open(sm.output[0], "w") as out:
        out.write("\t".join(all_cols) + "\n")
        for row in rows:
            out.write("\t".join(str(row.get(c, "NA")) for c in all_cols) + "\n")


if "snakemake" in dir():
    main()
