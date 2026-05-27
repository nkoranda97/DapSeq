rule filtered_stats:
    input:
        bam = OUT + "/bam/{sample}.{bam_type}.bam",
        bai = OUT + "/bam/{sample}.{bam_type}.bam.bai",
    output:
        stats    = OUT + "/stats/{sample}.{bam_type}.filtered_stats.txt",
        idxstats = OUT + "/stats/{sample}.{bam_type}.idxstats.txt",
    resources:
        mem_mb          = config["resources"]["filtered_stats"]["mem_mb"],
        runtime         = config["resources"]["filtered_stats"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
    log:
        OUT + "/logs/filtered_stats/{sample}.{bam_type}.log"
    shell:
        """
        exec 2>{log}
        FILTERED=$(samtools view -c {input.bam})
        printf "filtered_reads\t%d\n" "$FILTERED" >  {output.stats}
        printf "\\n### samtools flagstat ###\\n"   >> {output.stats}
        samtools flagstat {input.bam}             >> {output.stats}
        samtools idxstats {input.bam}              > {output.idxstats}
        """


rule frip_macs:
    input:
        bam   = OUT + "/bam/{sample}.{bam_type}.bam",
        bai   = OUT + "/bam/{sample}.{bam_type}.bam.bai",
        peaks = OUT + "/MACS/{sample}.{bam_type}_peaks.narrowPeak",
    output:
        OUT + "/stats/{sample}.{bam_type}.frip_macs.txt"
    wildcard_constraints:
        sample = "|".join(TREATMENT_SAMPLES) if TREATMENT_SAMPLES else "(?!)",
    resources:
        mem_mb          = config["resources"]["frip_macs"]["mem_mb"],
        runtime         = config["resources"]["frip_macs"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
    log:
        OUT + "/logs/frip_macs/{sample}.{bam_type}.log"
    shell:
        """
        exec 2>{log}
        # -F 2308 excludes: unmapped (4), not primary alignment (256), supplementary (2048)
        TOTAL=$(samtools view -c -F 2308 {input.bam})
        IN_PEAKS=$(bedtools intersect -abam {input.bam} -b {input.peaks} -u | samtools view -c)
        if [ "$TOTAL" -gt 0 ]; then
          FRIP=$(awk "BEGIN {{printf \\"%.4f\\", $IN_PEAKS / $TOTAL}}")
        else
          FRIP=0.0000
        fi
        printf "reads_in_peaks_macs\t%d\n" "$IN_PEAKS" >  {output}
        printf "frip_macs\t%s\n"           "$FRIP"      >> {output}
        """


rule frip_macs_filt:
    input:
        bam   = OUT + "/bam/{sample}.{bam_type}.bam",
        bai   = OUT + "/bam/{sample}.{bam_type}.bam.bai",
        peaks = OUT + "/MACS/{sample}.{bam_type}_peaks_filt.narrowPeak",
    output:
        OUT + "/stats/{sample}.{bam_type}.frip_macs_filt.txt"
    wildcard_constraints:
        sample = "|".join(TREATMENT_SAMPLES) if TREATMENT_SAMPLES else "(?!)",
    resources:
        mem_mb          = config["resources"]["frip_macs"]["mem_mb"],
        runtime         = config["resources"]["frip_macs"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
    log:
        OUT + "/logs/frip_macs_filt/{sample}.{bam_type}.log"
    shell:
        """
        exec 2>{log}
        TOTAL=$(samtools view -c -F 2308 {input.bam})
        IN_PEAKS=$(bedtools intersect -abam {input.bam} -b {input.peaks} -u | samtools view -c)
        if [ "$TOTAL" -gt 0 ]; then
          FRIP=$(awk "BEGIN {{printf \\"%.4f\\", $IN_PEAKS / $TOTAL}}")
        else
          FRIP=0.0000
        fi
        printf "reads_in_peaks_macs_filt\t%d\n" "$IN_PEAKS" >  {output}
        printf "frip_macs_filt\t%s\n"            "$FRIP"     >> {output}
        """


rule report:
    input:
        trim_logs    = expand(OUT + "/logs/bbduk/{sample}.trim.log",               sample=TREATMENT_SAMPLES),
        total_frags  = expand(OUT + "/stats/{sample}.total_frags.txt",             sample=TREATMENT_SAMPLES),
        align_rates  = expand(OUT + "/stats/{sample}.align_rate.txt",              sample=TREATMENT_SAMPLES),
        filt_stats   = expand(OUT + "/stats/{sample}.{bam_type}.filtered_stats.txt",
                              sample=TREATMENT_SAMPLES, bam_type=BAM_TYPES),
        frip_macs      = expand(OUT + "/stats/{sample}.{bam_type}.frip_macs.txt",
                                sample=TREATMENT_SAMPLES, bam_type=BAM_TYPES),
        frip_macs_filt = expand(OUT + "/stats/{sample}.{bam_type}.frip_macs_filt.txt",
                                sample=TREATMENT_SAMPLES, bam_type=BAM_TYPES),
        narrowpeaks    = expand(OUT + "/MACS/{sample}.{bam_type}_peaks.narrowPeak",
                                sample=TREATMENT_SAMPLES, bam_type=BAM_TYPES),
        logo_pngs    = expand(OUT + "/meme/{sample}.{bam_type}/summits/logo1.png",
                              sample=TREATMENT_SAMPLES, bam_type=BAM_TYPES),
    output:
        tsv  = OUT + "/stats/report.tsv",
        html = OUT + "/stats/report.html",
    params:
        treatment_samples = TREATMENT_SAMPLES,
        bam_types         = BAM_TYPES,
        trim_log_dir      = OUT + "/logs/bbduk",
        macs_dir          = OUT + "/MACS",
        meme_dir          = OUT + "/meme",
        stats_dir         = OUT + "/stats",
        min_foldch        = config["macs3"]["min_foldch"],
    resources:
        mem_mb          = config["resources"]["report"]["mem_mb"],
        runtime         = config["resources"]["report"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
    log:
        OUT + "/logs/report.log"
    script:
        "../scripts/report.py"
