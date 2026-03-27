if len(CONTROL_SAMPLES) > 1:
    rule merge_control:
        input:
            bams = expand(OUT + "/bam/{s}.bam",     s=CONTROL_SAMPLES),
            bais = expand(OUT + "/bam/{s}.bam.bai", s=CONTROL_SAMPLES),
        output:
            bam = OUT + "/bam/merged_control.bam",
            bai = OUT + "/bam/merged_control.bam.bai",
        params:
            extra_merge = config["samtools"].get("extra_merge", ""),
        threads:
            config["threads"]
        resources:
            mem_mb          = config["resources"]["merge_control"]["mem_mb"],
            runtime         = config["resources"]["merge_control"]["runtime"],
            slurm_partition = config["slurm_partition"],
            slurm_account   = config["slurm_account"],
        log:
            OUT + "/logs/samtools/merge_control.log"
        shell:
            """
            samtools merge -@ {threads} -f {params.extra_merge} {output.bam} {input.bams} 2>{log}
            samtools index {output.bam} 2>>{log}
            """


if CONTROL:
    rule bamcompare:
        input:
            sample_bam  = OUT + "/bam/{sample}.bam",
            sample_bai  = OUT + "/bam/{sample}.bam.bai",
            control_bam = OUT + f"/bam/{CONTROL}.bam",
            control_bai = OUT + f"/bam/{CONTROL}.bam.bai",
        output:
            OUT + "/bigWig/{sample}.peaks.bw"
        params:
            extra = config["bamcompare"].get("extra", ""),
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
              -p {threads} {params.extra} 2>{log}
            """
