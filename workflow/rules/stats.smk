_align_log_dir = "bwa_mem2" if config.get("aligner", "bowtie2") == "bwa_mem2" else "bowtie2"
_is_bwa_mem2   = config.get("aligner", "bowtie2") == "bwa_mem2"


rule sample_stats_treatment:
    input:
        bam          = OUT + "/bam/{sample}.bam",
        bai          = OUT + "/bam/{sample}.bam.bai",
        peaks        = OUT + "/MACS/{sample}_peaks.narrowPeak",
        peaks_filt   = OUT + "/MACS/{sample}_peaks_filt.narrowPeak",
        peaks_5fold  = OUT + "/MACS/{sample}_peaks_5fold.narrowPeak",
        fimo         = OUT + "/fimo/{sample}/peaks/fimo.tsv",
        subsample_log= OUT + "/logs/bbduk/{sample}.subsample.log",
        trim_log     = OUT + "/logs/bbduk/{sample}.trim.log",
        align_log    = OUT + f"/logs/{_align_log_dir}/{{sample}}.log",
        **({} if not _is_bwa_mem2 else
           {"flagstat_log": OUT + "/logs/bwa_mem2/{sample}.flagstat.log"}),
    output:
        stats_csv = OUT + "/stats/{sample}.stats.csv",
    wildcard_constraints:
        sample = "|".join(TREATMENT_SAMPLES) if TREATMENT_SAMPLES else "(?!)",
    params:
        is_pe        = lambda wc: wc.sample in PE_SAMPLES,
        is_treatment = True,
        aligner      = config.get("aligner", "bowtie2"),
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
        align_log    = OUT + f"/logs/{_align_log_dir}/{{sample}}.log",
        **({} if not _is_bwa_mem2 else
           {"flagstat_log": OUT + "/logs/bwa_mem2/{sample}.flagstat.log"}),
    output:
        stats_csv = OUT + "/stats/{sample}.stats.csv",
    wildcard_constraints:
        sample = "|".join(CONTROL_SAMPLES) if CONTROL_SAMPLES else "(?!)",
    params:
        is_pe        = lambda wc: wc.sample in PE_SAMPLES,
        is_treatment = False,
        aligner      = config.get("aligner", "bowtie2"),
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
        exec 2>{log}
        samtools idxstats {input.bam} > {output.idxstats}
        """


rule report:
    input:
        stats_csvs = expand(OUT + "/stats/{sample}.stats.csv", sample=TREATMENT_SAMPLES),
    output:
        csv = OUT + "/stats/report.csv",
    params:
        treatment_samples = TREATMENT_SAMPLES,
        stats_dir         = OUT + "/stats",
    resources:
        mem_mb          = config["resources"]["report"]["mem_mb"],
        runtime         = config["resources"]["report"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
    log:
        OUT + "/logs/report.log"
    script:
        "../scripts/report.py"
