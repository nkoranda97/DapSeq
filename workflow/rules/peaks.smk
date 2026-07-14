_MACS3_KEEP_DUP = config["macs3"].get("keep_dup", 1)
_MACS3_EXTRA    = config["macs3"].get("extra", "")


# Peak calling for every report sample. Treatment samples use CONTROL as the
# -c background; control samples use *themselves* (control-against-itself run,
# expected to yield few or no peaks) so the control flows through the same
# downstream filtering/motif/stats/report path as treatments. control_bam is
# resolved per-sample; params.ctrl emits -c only when a background bam exists.
rule macs3:
    input:
        sample_bam  = OUT + "/bam/{sample}.bam",
        control_bam = lambda wc: (
            OUT + f"/bam/{wc.sample}.bam" if wc.sample in CONTROL_SAMPLES
            else (OUT + f"/bam/{CONTROL}.bam" if CONTROL else [])
        ),
    output:
        summits    = OUT + "/MACS/{sample}_summits.bed",
        narrowpeak = OUT + "/MACS/{sample}_peaks.narrowPeak",
    wildcard_constraints:
        sample = REPORT_SAMPLE_CONSTRAINT,
    params:
        ctrl        = lambda wc, input: f"-c {input.control_bam}" if input.control_bam else "",
        outdir      = OUT + "/MACS",
        name        = lambda wc: wc.sample,
        keep_dup    = _MACS3_KEEP_DUP,
        extra       = _MACS3_EXTRA,
        tmpdir      = OUT + "/temp",
        macs3_format = config["macs3"]["format"],
        genome_size  = config["genome_size"],
    resources:
        mem_mb          = config["resources"]["macs3"]["mem_mb"],
        runtime         = config["resources"]["macs3"]["runtime"],
    log:
        OUT + "/logs/macs3/{sample}.log"
    shell:
        """
        set -euo pipefail
        mkdir -p {params.tmpdir}
        export TMPDIR={params.tmpdir}
        macs3 callpeak \
          -t {input.sample_bam} {params.ctrl} \
          -f {params.macs3_format} --outdir {params.outdir} \
          -g {params.genome_size} -n {params.name} \
          --call-summits --keep-dup {params.keep_dup} {params.extra} --verbose=0 2>{log}
        """


rule macs3_control:
    input:
        sample_bam = OUT + "/bam/{sample}.bam",
    output:
        summits    = OUT + "/MACS/{sample}_control_summits.bed",
        narrowpeak = OUT + "/MACS/{sample}_control_peaks.narrowPeak",
    wildcard_constraints:
        sample = CONTROL_SAMPLE_CONSTRAINT,
    params:
        outdir      = OUT + "/MACS",
        name        = lambda wc: f"{wc.sample}_control",
        keep_dup    = _MACS3_KEEP_DUP,
        extra       = _MACS3_EXTRA,
        tmpdir      = OUT + "/temp",
        macs3_format = config["macs3"]["format"],
        genome_size  = config["genome_size"],
    resources:
        mem_mb          = config["resources"]["macs3"]["mem_mb"],
        runtime         = config["resources"]["macs3"]["runtime"],
    log:
        OUT + "/logs/macs3_control/{sample}.log"
    shell:
        """
        set -euo pipefail
        mkdir -p {params.tmpdir}
        export TMPDIR={params.tmpdir}
        macs3 callpeak \
          -t {input.sample_bam} \
          -f {params.macs3_format} --outdir {params.outdir} \
          -g {params.genome_size} -n {params.name} \
          --call-summits --keep-dup {params.keep_dup} {params.extra} --verbose=0 2>{log}
        """


if config.get("blacklist_filter", {}).get("enabled", False):
    rule blacklist_filter_peaks:
        input:
            peaks     = OUT + f"/MACS/{{sample}}_peaks_fold{MEME_FOLD_IDX}.narrowPeak",
            blacklist = config["blacklist_filter"]["bed"],
        output:
            OUT + f"/MACS/{{sample}}_peaks_fold{MEME_FOLD_IDX}_bl.narrowPeak",
        resources:
            mem_mb          = config["resources"]["blacklist_filter_peaks"]["mem_mb"],
            runtime         = config["resources"]["blacklist_filter_peaks"]["runtime"],
        log:
            OUT + "/logs/blacklist_filter/{sample}.log"
        shell:
            """
            set -euo pipefail
            bedtools intersect -v -a {input.peaks} -b {input.blacklist} \
              > {output} 2>{log}
            """


if config.get("rmsk_filter", {}).get("enabled", False):
    _RMSK_IN_SUFFIX = ("_bl" if config.get("blacklist_filter", {}).get("enabled", False) else "")
    rule rmsk_filter_peaks:
        input:
            peaks = get_pre_rmsk_peaks,
            rmsk  = config["rmsk_filter"]["txt"],
        output:
            OUT + f"/MACS/{{sample}}_peaks_fold{MEME_FOLD_IDX}{_RMSK_IN_SUFFIX}_rmsk.narrowPeak",
        resources:
            mem_mb          = config["resources"]["rmsk_filter_peaks"]["mem_mb"],
            runtime         = config["resources"]["rmsk_filter_peaks"]["runtime"],
        log:
            OUT + "/logs/rmsk_filter/{sample}.log"
        shell:
            """
            set -euo pipefail
            zcat {input.rmsk} \
              | awk '($7+0)==$7 {{OFS="\\t"; print $6, $7, $8}}' \
              | bedtools intersect -v -a {input.peaks} -b stdin \
              > {output} 2>{log}
            """


rule filter_peaks:
    input:
        OUT + "/MACS/{sample}_peaks.narrowPeak",
    output:
        fold1 = OUT + "/MACS/{sample}_peaks_fold1.narrowPeak",
        fold2 = OUT + "/MACS/{sample}_peaks_fold2.narrowPeak",
        fold3 = OUT + "/MACS/{sample}_peaks_fold3.narrowPeak",
    params:
        fch1 = FOLD_LEVELS[0],
        fch2 = FOLD_LEVELS[1],
        fch3 = FOLD_LEVELS[2],
    resources:
        mem_mb          = config["resources"]["filter_peaks"]["mem_mb"],
        runtime         = config["resources"]["filter_peaks"]["runtime"],
    log:
        OUT + "/logs/filter_peaks/{sample}.log"
    shell:
        """
        set -euo pipefail
        awk -v FCH={params.fch1} '$7 >= FCH' {input} > {output.fold1} 2>{log}
        awk -v FCH={params.fch2} '$7 >= FCH' {input} > {output.fold2} 2>>{log}
        awk -v FCH={params.fch3} '$7 >= FCH' {input} > {output.fold3} 2>>{log}
        """
