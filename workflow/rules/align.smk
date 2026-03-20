if ALIGNER == "bowtie2":
    rule bowtie2:
        input:
            r1  = OUT + "/trimmed/{sample}.R1.fastq.gz",
            r2  = get_r2,
            idx = config["genome_ref"] + ".1.bt2",
        output:
            temp(OUT + "/sam/{sample}.sam")
        params:
            idx     = config["genome_ref"],
            r2_arg  = lambda wc: "" if wc.sample in SE_SAMPLES else "-2 " + OUT + f"/trimmed/{wc.sample}.R2.fastq.gz",
            pe_flag = lambda wc: "-U" if wc.sample in SE_SAMPLES else "-1",
            extra   = config["bowtie2"].get("extra", ""),
        threads:
            config["threads"]
        resources:
            mem_mb          = config["resources"]["bowtie2"]["mem_mb"],
            runtime         = config["resources"]["bowtie2"]["runtime"],
            slurm_partition = config["slurm_partition"],
            slurm_account   = config["slurm_account"],
        log:
            OUT + "/logs/bowtie2/{sample}.log"
        shell:
            "bowtie2 -p {threads} {params.extra} -x {params.idx} {params.pe_flag} {input.r1} {params.r2_arg} -S {output} 2>{log}"

elif ALIGNER == "bwa":
    rule bwa_mem:
        input:
            r1  = OUT + "/trimmed/{sample}.R1.fastq.gz",
            r2  = get_r2,
            idx = config["genome_ref"] + ".bwt",
        output:
            temp(OUT + "/sam/{sample}.sam")
        params:
            r2_arg = lambda wc: "" if wc.sample in SE_SAMPLES else OUT + f"/trimmed/{wc.sample}.R2.fastq.gz",
            k      = config["bwa"]["k"],
            B      = config["bwa"]["B"],
            O      = config["bwa"]["O"],
        threads:
            config["threads"]
        resources:
            mem_mb          = config["resources"]["bwa_mem"]["mem_mb"],
            runtime         = config["resources"]["bwa_mem"]["runtime"],
            slurm_partition = config["slurm_partition"],
            slurm_account   = config["slurm_account"],
        log:
            OUT + "/logs/bwa/{sample}.log"
        shell:
            "bwa mem -t {threads} -k {params.k} -B {params.B} -O {params.O} -v 3 {config[genome_ref]} {input.r1} {params.r2_arg} > {output} 2>{log}"
