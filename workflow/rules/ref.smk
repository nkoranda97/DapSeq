rule bowtie2_index:
    input:
        config["genome_ref"]
    output:
        bt2_1  = config["genome_ref"] + ".1.bt2",
        bt2_2  = config["genome_ref"] + ".2.bt2",
        bt2_3  = config["genome_ref"] + ".3.bt2",
        bt2_4  = config["genome_ref"] + ".4.bt2",
        bt2_r1 = config["genome_ref"] + ".rev.1.bt2",
        bt2_r2 = config["genome_ref"] + ".rev.2.bt2",
        sizes  = config["genome_ref"] + ".sizes",
        fai    = config["genome_ref"] + ".fai",
    params:
        extra_build = config["bowtie2"].get("extra_build", ""),
    resources:
        mem_mb          = config["resources"]["bowtie2_index"]["mem_mb"],
        runtime         = config["resources"]["bowtie2_index"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
    log:
        OUT + "/logs/bowtie2_index.log"
    shell:
        """
        bowtie2-build {params.extra_build} {input} {input} 2>{log}
        samtools faidx {input} 2>>{log}
        cut -f1,2 {input}.fai > {output.sizes}
        """
