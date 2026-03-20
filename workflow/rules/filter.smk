rule samtools_filter_sort_dedup:
    input:
        OUT + "/sam/{sample}.sam"
    output:
        bam = OUT + "/bam/{sample}.bam",
        bai = OUT + "/bam/{sample}.bam.bai",
    params:
        prefix = OUT + "/bam/{sample}.tmp.bam",
        is_se  = lambda wc: "true" if wc.sample in SE_SAMPLES else "false",
    threads:
        config["threads"]
    resources:
        mem_mb          = config["resources"]["samtools_filter_sort_dedup"]["mem_mb"],
        runtime         = config["resources"]["samtools_filter_sort_dedup"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
    log:
        OUT + "/logs/samtools/{sample}.log"
    shell:
        """
        printf "primary_mapped %d\n" "$(samtools view -c -F 2308 {input})" >> {log}
        printf "mapq_passed %d\n" "$(samtools view -c -F 2308 -q {config[samtools][mapq]} {input})" >> {log}
        if [ "{params.is_se}" = "true" ]; then
          samtools view -@ {threads} -h -F 4 -q {config[samtools][mapq]} -u {input} \
            | samtools sort -@ {threads} -o {params.prefix} - 2>>{log}
        else
          samtools view -@ {threads} -h -F 4 -q {config[samtools][mapq]} -u {input} \
            | samtools fixmate -m -@ {threads} - - \
            | samtools sort -@ {threads} -o {params.prefix} - 2>>{log}
        fi
        samtools markdup -r -@ {threads} {params.prefix} {output.bam} 2>>{log}
        samtools index {output.bam} 2>>{log}
        rm -f {params.prefix}
        """


rule bamcoverage:
    input:
        bam = OUT + "/bam/{sample}.bam",
        bai = OUT + "/bam/{sample}.bam.bai",
    output:
        OUT + "/bigWig/{sample}.bw"
    threads:
        config["threads"]
    resources:
        mem_mb          = config["resources"]["bamcoverage"]["mem_mb"],
        runtime         = config["resources"]["bamcoverage"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
    log:
        OUT + "/logs/bamcoverage/{sample}.log"
    shell:
        "bamCoverage -b {input.bam} --normalizeUsing {config[bamcoverage][normalize_using]} --extendReads {config[bamcoverage][extend_reads]} -p {threads} -o {output} 2>{log}"


if CONTROL:
    rule bamcompare:
        input:
            sample_bam  = OUT + "/bam/{sample}.bam",
            sample_bai  = OUT + "/bam/{sample}.bam.bai",
            control_bam = OUT + f"/bam/{CONTROL}.bam",
            control_bai = OUT + f"/bam/{CONTROL}.bam.bai",
        output:
            OUT + "/bigWig/{sample}.peaks.bw"
        threads:
            config["threads"]
        resources:
            mem_mb          = config["resources"]["bamcompare"]["mem_mb"],
            runtime         = config["resources"]["bamcompare"]["runtime"],
            slurm_partition = config["slurm_partition"],
            slurm_account   = config["slurm_account"],
        log:
            OUT + "/logs/bamcompare/{sample}.log"
        shell:
            """
            bamCompare \
              -b1 {input.sample_bam} -b2 {input.control_bam} \
              -o {output} \
              --binSize {config[bamcompare][bin_size]} --operation {config[bamcompare][operation]} \
              --scaleFactorsMethod {config[bamcompare][scale_factors_method]} -n {config[bamcompare][n]} \
              -p {threads} 2>{log}
            """
