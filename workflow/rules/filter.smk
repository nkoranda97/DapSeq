if len(CONTROL_SAMPLES) > 1:
    rule merge_control:
        input:
            bams = lambda wc: expand(OUT + "/bam/{s}.{bam_type}.bam",
                                     s=CONTROL_SAMPLES, bam_type=wc.bam_type),
            bais = lambda wc: expand(OUT + "/bam/{s}.{bam_type}.bam.bai",
                                     s=CONTROL_SAMPLES, bam_type=wc.bam_type),
        output:
            bam = OUT + "/bam/merged_control.{bam_type}.bam",
            bai = OUT + "/bam/merged_control.{bam_type}.bam.bai",
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
            OUT + "/logs/samtools/merge_control.{bam_type}.log"
        shell:
            """
            samtools merge -@ {threads} -f {params.extra_merge} {output.bam} {input.bams} 2>{log}
            samtools index {output.bam} 2>>{log}
            """


if CONTROL:
    rule bamcompare:
        input:
            sample_bam  = OUT + "/bam/{sample}.{bam_type}.bam",
            sample_bai  = OUT + "/bam/{sample}.{bam_type}.bam.bai",
            control_bam = OUT + f"/bam/{CONTROL}.{{bam_type}}.bam",
            control_bai = OUT + f"/bam/{CONTROL}.{{bam_type}}.bam.bai",
        output:
            OUT + "/bigWig/{sample}.{bam_type}.peaks.bw"
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
            OUT + "/logs/bamcompare/{sample}.{bam_type}.log"
        shell:
            """
            bamCompare \
              -b1 {input.sample_bam} -b2 {input.control_bam} \
              -o {output} \
              --binSize {config[bamcompare][bin_size]} --operation {config[bamcompare][operation]} \
              --scaleFactorsMethod {config[bamcompare][scale_factors_method]} -n {config[bamcompare][n]} \
              -p {threads} {params.extra} 2>{log}
            """
