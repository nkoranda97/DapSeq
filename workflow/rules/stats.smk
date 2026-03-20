rule alignment_stats:
    input:
        OUT + "/sam/{sample}.sam"
    output:
        OUT + "/stats/{sample}.raw_alignment_stats.txt"
    params:
        mapq  = config["samtools"]["mapq"],
        is_pe = lambda wc: "true" if wc.sample in PE_SAMPLES else "false",
    resources:
        mem_mb          = config["resources"]["alignment_stats"]["mem_mb"],
        runtime         = config["resources"]["alignment_stats"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
    log:
        OUT + "/logs/alignment_stats/{sample}.log"
    shell:
        """
        exec 2>{log}
        TOTAL=$(samtools view -c -F 2308 {input})
        MAPPED=$(samtools view -c -F 2312 {input})
        UNMAPPED=$(samtools view -c -f 4 -F 2304 {input})
        if [ "$TOTAL" -gt 0 ]; then
          MAPPING_RATE=$(awk "BEGIN {{printf \\"%.2f\\", 100.0 * $MAPPED / $TOTAL}}")
        else
          MAPPING_RATE=0.00
        fi
        UNIQUE=$(samtools view -c -F 2312 -q {params.mapq} {input})
        MULTI=$((MAPPED - UNIQUE))
        if [ "$TOTAL" -gt 0 ]; then
          MULTI_RATE=$(awk "BEGIN {{printf \\"%.2f\\", 100.0 * $MULTI / $TOTAL}}")
        else
          MULTI_RATE=0.00
        fi
        SOFT=$(samtools view -F 2312 {input} | awk '$6 ~ /S/' | wc -l)

        printf "total_reads\t%d\n"          "$TOTAL"        >  {output}
        printf "mapped_reads\t%d\n"         "$MAPPED"       >> {output}
        printf "unmapped_reads\t%d\n"       "$UNMAPPED"     >> {output}
        printf "mapping_rate\t%s\n"         "$MAPPING_RATE" >> {output}
        printf "unique_mapped\t%d\n"        "$UNIQUE"       >> {output}
        printf "multimapped\t%d\n"          "$MULTI"        >> {output}
        printf "multimapped_rate\t%s\n"     "$MULTI_RATE"   >> {output}
        printf "soft_clipped_reads\t%d\n"   "$SOFT"         >> {output}

        if [ "{params.is_pe}" = "true" ]; then
          PP=$(samtools view -c -f 2 -F 2308 {input})
          if [ "$TOTAL" -gt 0 ]; then
            PP_RATE=$(awk "BEGIN {{printf \\"%.2f\\", 100.0 * $PP / $TOTAL}}")
          else
            PP_RATE=0.00
          fi
          printf "properly_paired\t%d\n"      "$PP"      >> {output}
          printf "properly_paired_rate\t%s\n" "$PP_RATE" >> {output}
        fi

        printf "\\n### samtools flagstat ###\\n" >> {output}
        samtools flagstat {input}               >> {output}
        """


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


rule frip:
    input:
        bam   = OUT + "/bam/{sample}.bam",
        bai   = OUT + "/bam/{sample}.bam.bai",
        peaks = OUT + "/combined_bed/{sample}.combined.bed",
    output:
        OUT + "/stats/{sample}.frip.txt"
    wildcard_constraints:
        sample = "|".join(TREATMENT_SAMPLES) if TREATMENT_SAMPLES else "(?!)",
    resources:
        mem_mb          = config["resources"]["frip"]["mem_mb"],
        runtime         = config["resources"]["frip"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
    log:
        OUT + "/logs/frip/{sample}.log"
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
        printf "reads_in_peaks\t%d\n" "$IN_PEAKS" >  {output}
        printf "frip\t%s\n"           "$FRIP"      >> {output}
        """

if USE_MACS:
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

if USE_GEM:
    rule frip_gem:
        input:
            bam   = OUT + "/bam/{sample}.bam",
            bai   = OUT + "/bam/{sample}.bam.bai",
            peaks = OUT + "/compare_bed/{sample}.GEM.bed",
        output:
            OUT + "/stats/{sample}.frip_gem.txt"
        wildcard_constraints:
            sample = "|".join(TREATMENT_SAMPLES) if TREATMENT_SAMPLES else "(?!)",
        resources:
            mem_mb          = config["resources"]["frip_gem"]["mem_mb"],
            runtime         = config["resources"]["frip_gem"]["runtime"],
            slurm_partition = config["slurm_partition"],
            slurm_account   = config["slurm_account"],
        log:
            OUT + "/logs/frip_gem/{sample}.log"
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
            printf "reads_in_peaks_gem\t%d\n" "$IN_PEAKS" >  {output}
            printf "frip_gem\t%s\n"           "$FRIP"      >> {output}
            """


rule qc_summary:
    input:
        raw_stats      = expand(OUT + "/stats/{sample}.raw_alignment_stats.txt", sample=SAMPLES),
        filtered_stats = expand(OUT + "/stats/{sample}.filtered_stats.txt",      sample=SAMPLES),
        frip           = expand(OUT + "/stats/{sample}.frip.txt",                sample=TREATMENT_SAMPLES),
        frip_macs      = (expand(OUT + "/stats/{sample}.frip_macs.txt", sample=TREATMENT_SAMPLES) if USE_MACS else []),
        frip_gem       = (expand(OUT + "/stats/{sample}.frip_gem.txt",  sample=TREATMENT_SAMPLES) if USE_GEM  else []),
    output:
        OUT + "/stats/qc_summary.tsv"
    params:
        samples           = SAMPLES,
        treatment_samples = TREATMENT_SAMPLES,
        use_macs          = USE_MACS,
        use_gem           = USE_GEM,
        out_dir           = OUT + "/stats",
    resources:
        mem_mb          = config["resources"]["qc_summary"]["mem_mb"],
        runtime         = config["resources"]["qc_summary"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
    log:
        OUT + "/logs/qc_summary.log"
    script:
        "../scripts/qc_summary.py"
