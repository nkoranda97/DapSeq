rule fastqc:
    input:
        OUT + "/trimmed/{sample}.{read}.fastq.gz"
    output:
        html = OUT + "/Fastqc/{sample}.{read}_fastqc.html",
        zip  = OUT + "/Fastqc/{sample}.{read}_fastqc.zip",
    params:
        extra = config.get("fastqc", {}).get("extra", ""),
    resources:
        mem_mb          = config["resources"]["fastqc"]["mem_mb"],
        runtime         = config["resources"]["fastqc"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
    log:
        OUT + "/logs/fastqc/{sample}.{read}.log"
    shell:
        "fastqc {params.extra} {input} --outdir={OUT}/Fastqc/ 2>{log}"


rule multiqc:
    input:
        expand(OUT + "/Fastqc/{sample}.R1_fastqc.zip", sample=SAMPLES),
        expand(OUT + "/Fastqc/{sample}.R2_fastqc.zip", sample=PE_SAMPLES),
    output:
        OUT + "/multiqc_report.html"
    params:
        extra = config.get("multiqc", {}).get("extra", ""),
    resources:
        mem_mb          = config["resources"]["multiqc"]["mem_mb"],
        runtime         = config["resources"]["multiqc"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
    log:
        OUT + "/logs/multiqc.log"
    shell:
        "multiqc {params.extra} {OUT}/Fastqc/ -o {OUT}/ 2>{log}"
