# Per-sample stats for every report sample (treatment samples and controls
# alike). Controls now have peaks from the against-itself MACS run (rule macs3
# with the control as its own -c background), so they get the same full
# peak/FRiP/motif stats as treatment samples. All report samples have peaks, so
# peak-stat computation runs unconditionally; the semantic treatment/control
# distinction lives in the report (flagged row) and the DB (is_treatment column).
rule sample_stats:
    input:
        bam          = OUT + "/bam/{sample}.bam",
        bai          = OUT + "/bam/{sample}.bam.bai",
        peaks        = OUT + "/MACS/{sample}_peaks.narrowPeak",
        peaks_fold1  = OUT + "/MACS/{sample}_peaks_fold1.narrowPeak",
        peaks_fold2  = OUT + "/MACS/{sample}_peaks_fold2.narrowPeak",
        peaks_fold3  = OUT + "/MACS/{sample}_peaks_fold3.narrowPeak",
        fimo         = OUT + "/fimo/{sample}/peaks/fimo.tsv",
        subsample_log= OUT + "/logs/bbduk/{sample}.subsample.log",
        trim_log     = OUT + "/logs/bbduk/{sample}.trim.log",
        peaks_bl     = lambda wc: (
            [OUT + f"/MACS/{wc.sample}_peaks_fold{MEME_FOLD_IDX}_bl.narrowPeak"]
            if config.get("blacklist_filter", {}).get("enabled", False)
            else []
        ),
        peaks_rmsk   = lambda wc: (
            [OUT + f"/MACS/{wc.sample}_peaks_fold{MEME_FOLD_IDX}{PEAKS_FILTER_SUFFIX}.narrowPeak"]
            if config.get("rmsk_filter", {}).get("enabled", False)
            else []
        ),
    output:
        stats_csv = OUT + "/stats/{sample}.stats.csv",
    wildcard_constraints:
        sample = REPORT_SAMPLE_CONSTRAINT,
    params:
        is_pe             = lambda wc: wc.sample in PE_SAMPLES,
        max_frags         = config["bbduk"].get("max_frags"),
        blacklist_enabled = config.get("blacklist_filter", {}).get("enabled", False),
        rmsk_enabled      = config.get("rmsk_filter",     {}).get("enabled", False),
        meme_foldch_level = MEME_FOLD_IDX,
        experiment_date   = lambda wc: config["samples"][wc.sample].get("experiment_date"),
        gdna_batch        = lambda wc: config["samples"][wc.sample].get("gdna_batch"),
    resources:
        mem_mb          = config["resources"]["sample_stats"]["mem_mb"],
        runtime         = config["resources"]["sample_stats"]["runtime"],
    log:
        OUT + "/logs/sample_stats/{sample}.log"
    script:
        "../scripts/collect_stats.py"


rule idxstats:
    input:
        bam = OUT + "/bam/{sample}.bam",
        bai = OUT + "/bam/{sample}.bam.bai",
    output:
        idxstats = OUT + "/stats/{sample}.idxstats.txt",
    resources:
        mem_mb          = config["resources"]["alignment_stats"]["mem_mb"],
        runtime         = config["resources"]["alignment_stats"]["runtime"],
    log:
        OUT + "/logs/idxstats/{sample}.log"
    shell:
        """
        set -euo pipefail
        samtools idxstats {input.bam} > {output.idxstats} 2>{log}
        """


rule report:
    input:
        stats_csvs = expand(OUT + "/stats/{sample}.stats.csv", sample=REPORT_SAMPLES),
    output:
        csv = OUT + "/stats/report.csv",
    params:
        report_samples = REPORT_SAMPLES,
        stats_dir      = OUT + "/stats",
    resources:
        mem_mb          = config["resources"]["report"]["mem_mb"],
        runtime         = config["resources"]["report"]["runtime"],
    log:
        OUT + "/logs/report.log"
    script:
        "../scripts/report.py"


rule report_html:
    input:
        report_csv      = OUT + "/stats/report.csv",
        meme_logos      = expand(OUT + "/meme/{sample}/summits/logo1.png",    sample=REPORT_SAMPLES),
        meme_logos_rc   = expand(OUT + "/meme/{sample}/summits/logo_rc1.png", sample=REPORT_SAMPLES),
        factorbook_logos= expand(OUT + "/factorbook/{sample}.logo.png",       sample=REPORT_SAMPLES),
    output:
        html = OUT + "/stats/report.html",
    params:
        report_samples    = REPORT_SAMPLES,
        control_samples   = CONTROL_SAMPLES,
        meme_dir          = OUT + "/meme",
        factorbook_dir    = OUT + "/factorbook",
        filter_foldch     = FOLD_LEVELS[MEME_FOLD_IDX - 1],
    resources:
        mem_mb          = config["resources"]["report"]["mem_mb"],
        runtime         = config["resources"]["report"]["runtime"],
    log:
        OUT + "/logs/report_html.log"
    script:
        "../scripts/render_report_html.py"
