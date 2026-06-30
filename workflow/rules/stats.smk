rule sample_stats_treatment:
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
        sample = TREATMENT_SAMPLE_CONSTRAINT,
    params:
        is_pe             = lambda wc: wc.sample in PE_SAMPLES,
        is_treatment      = True,
        max_frags         = config["bbduk"].get("max_frags"),
        blacklist_enabled = config.get("blacklist_filter", {}).get("enabled", False),
        rmsk_enabled      = config.get("rmsk_filter",     {}).get("enabled", False),
        foldch_levels     = FOLD_LEVELS,
        meme_foldch_level = MEME_FOLD_IDX,
    resources:
        mem_mb          = config["resources"]["sample_stats"]["mem_mb"],
        runtime         = config["resources"]["sample_stats"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
    log:
        OUT + "/logs/sample_stats/{sample}.log"
    script:
        "../scripts/collect_stats.py"


rule sample_stats_control:
    input:
        bam          = OUT + "/bam/{sample}.bam",
        bai          = OUT + "/bam/{sample}.bam.bai",
        subsample_log= OUT + "/logs/bbduk/{sample}.subsample.log",
        trim_log     = OUT + "/logs/bbduk/{sample}.trim.log",
    output:
        stats_csv = OUT + "/stats/{sample}.stats.csv",
    wildcard_constraints:
        sample = CONTROL_SAMPLE_CONSTRAINT,
    params:
        is_pe        = lambda wc: wc.sample in PE_SAMPLES,
        is_treatment = False,
        max_frags    = config["bbduk"].get("max_frags"),
    resources:
        mem_mb          = config["resources"]["sample_stats"]["mem_mb"],
        runtime         = config["resources"]["sample_stats"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
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
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
    log:
        OUT + "/logs/idxstats/{sample}.log"
    shell:
        """
        set -euo pipefail
        samtools idxstats {input.bam} > {output.idxstats} 2>{log}
        """


rule report:
    input:
        stats_csvs = expand(OUT + "/stats/{sample}.stats.csv", sample=TREATMENT_SAMPLES),
    output:
        csv = OUT + "/stats/report.csv",
    params:
        treatment_samples = TREATMENT_SAMPLES,
        stats_dir         = OUT + "/stats",
        foldch_levels     = FOLD_LEVELS,
        meme_foldch_level = MEME_FOLD_IDX,
    resources:
        mem_mb          = config["resources"]["report"]["mem_mb"],
        runtime         = config["resources"]["report"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
    log:
        OUT + "/logs/report.log"
    script:
        "../scripts/report.py"


rule report_html:
    input:
        report_csv      = OUT + "/stats/report.csv",
        meme_logos      = expand(OUT + "/meme/{sample}/summits/logo1.png",    sample=TREATMENT_SAMPLES),
        meme_logos_rc   = expand(OUT + "/meme/{sample}/summits/logo_rc1.png", sample=TREATMENT_SAMPLES),
        factorbook_logos= expand(OUT + "/factorbook/{sample}.logo.png",       sample=TREATMENT_SAMPLES),
    output:
        html = OUT + "/stats/report.html",
    params:
        treatment_samples = TREATMENT_SAMPLES,
        meme_dir          = OUT + "/meme",
        factorbook_dir    = OUT + "/factorbook",
    resources:
        mem_mb          = config["resources"]["report"]["mem_mb"],
        runtime         = config["resources"]["report"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
    log:
        OUT + "/logs/report_html.log"
    script:
        "../scripts/render_report_html.py"
