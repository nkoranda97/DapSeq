rule filtered_stats:
    input:
        bam = OUT + "/bam/{sample}.bam",
        bai = OUT + "/bam/{sample}.bam.bai",
    output:
        stats    = OUT + "/stats/{sample}.filtered_stats.txt",
        idxstats = OUT + "/stats/{sample}.idxstats.txt",
    resources:
        mem_mb          = config["resources"]["filtered_stats"]["mem_mb"],
        runtime         = config["resources"]["filtered_stats"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
    log:
        OUT + "/logs/filtered_stats/{sample}.log"
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
        bam   = OUT + "/bam/{sample}.bam",
        bai   = OUT + "/bam/{sample}.bam.bai",
        peaks = OUT + "/MACS/{sample}_peaks.narrowPeak",
    output:
        OUT + "/stats/{sample}.frip_macs.txt"
    wildcard_constraints:
        sample = "|".join(TREATMENT_SAMPLES) if TREATMENT_SAMPLES else "(?!)",
    resources:
        mem_mb          = config["resources"]["frip_macs"]["mem_mb"],
        runtime         = config["resources"]["frip_macs"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
    log:
        OUT + "/logs/frip_macs/{sample}.log"
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
        bam   = OUT + "/bam/{sample}.bam",
        bai   = OUT + "/bam/{sample}.bam.bai",
        peaks = OUT + "/MACS/{sample}_peaks_filt.narrowPeak",
    output:
        OUT + "/stats/{sample}.frip_macs_filt.txt"
    wildcard_constraints:
        sample = "|".join(TREATMENT_SAMPLES) if TREATMENT_SAMPLES else "(?!)",
    resources:
        mem_mb          = config["resources"]["frip_macs"]["mem_mb"],
        runtime         = config["resources"]["frip_macs"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
    log:
        OUT + "/logs/frip_macs_filt/{sample}.log"
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


rule frip_macs_5fold:
    input:
        bam   = OUT + "/bam/{sample}.bam",
        bai   = OUT + "/bam/{sample}.bam.bai",
        peaks = OUT + "/MACS/{sample}_peaks_5fold.narrowPeak",
    output:
        OUT + "/stats/{sample}.frip_macs_5fold.txt"
    wildcard_constraints:
        sample = "|".join(TREATMENT_SAMPLES) if TREATMENT_SAMPLES else "(?!)",
    resources:
        mem_mb          = config["resources"]["frip_macs"]["mem_mb"],
        runtime         = config["resources"]["frip_macs"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
    log:
        OUT + "/logs/frip_macs_5fold/{sample}.log"
    shell:
        """
        exec 2>{log}
        IN_PEAKS=$(bedtools intersect -abam {input.bam} -b {input.peaks} -u | samtools view -c)
        printf "reads_in_peaks_macs_5fold\t%d\n" "$IN_PEAKS" > {output}
        """


rule report:
    input:
        trim_logs  = expand(OUT + "/logs/bbduk/{sample}.trim.log",        sample=TREATMENT_SAMPLES),
        total_frags= expand(OUT + "/stats/{sample}.total_frags.txt",       sample=TREATMENT_SAMPLES),
        filt_stats = expand(OUT + "/stats/{sample}.filtered_stats.txt",    sample=TREATMENT_SAMPLES),
        frip       = expand(OUT + "/stats/{sample}.frip_macs.txt",         sample=TREATMENT_SAMPLES),
        frip_filt  = expand(OUT + "/stats/{sample}.frip_macs_filt.txt",    sample=TREATMENT_SAMPLES),
        frip_5fold = expand(OUT + "/stats/{sample}.frip_macs_5fold.txt",   sample=TREATMENT_SAMPLES),
        align_rates= expand(OUT + "/stats/{sample}.align_rate.txt",        sample=TREATMENT_SAMPLES),
        narrowpeaks= expand(OUT + "/MACS/{sample}_peaks_filt.narrowPeak",  sample=TREATMENT_SAMPLES),
        fimo_peaks = expand(OUT + "/fimo/{sample}/peaks/fimo.tsv",         sample=TREATMENT_SAMPLES),
    output:
        csv = OUT + "/stats/report.csv",
    params:
        treatment_samples = TREATMENT_SAMPLES,
        trim_log_dir      = OUT + "/logs/bbduk",
        macs_dir          = OUT + "/MACS",
        stats_dir         = OUT + "/stats",
        fimo_dir          = OUT + "/fimo",
    resources:
        mem_mb          = config["resources"]["report"]["mem_mb"],
        runtime         = config["resources"]["report"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
    log:
        OUT + "/logs/report.log"
    script:
        "../scripts/report.py"
