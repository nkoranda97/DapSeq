if ALIGNER == "bowtie2":
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
        resources:
            mem_mb          = config["resources"]["bowtie2_index"]["mem_mb"],
            runtime         = config["resources"]["bowtie2_index"]["runtime"],
            slurm_partition = config["slurm_partition"],
            slurm_account   = config["slurm_account"],
        log:
            OUT + "/logs/bowtie2_index.log"
        shell:
            """
            bowtie2-build {input} {input} 2>{log}
            samtools faidx {input} 2>>{log}
            cut -f1,2 {input}.fai > {output.sizes}
            """

elif ALIGNER == "bwa":
    rule bwa_index:
        input:
            config["genome_ref"]
        output:
            bwt   = config["genome_ref"] + ".bwt",
            pac   = config["genome_ref"] + ".pac",
            ann   = config["genome_ref"] + ".ann",
            amb   = config["genome_ref"] + ".amb",
            sa    = config["genome_ref"] + ".sa",
            sizes = config["genome_ref"] + ".sizes",
            fai   = config["genome_ref"] + ".fai",
        resources:
            mem_mb          = config["resources"]["bwa_index"]["mem_mb"],
            runtime         = config["resources"]["bwa_index"]["runtime"],
            slurm_partition = config["slurm_partition"],
            slurm_account   = config["slurm_account"],
        log:
            OUT + "/logs/bwa_index.log"
        shell:
            """
            bwa index -a bwtsw {input} 2>{log}
            samtools faidx {input} 2>>{log}
            cut -f1,2 {input}.fai > {output.sizes}
            """


rule split_genome:
    input:
        fa   = config["genome_ref"],
        fai  = config["genome_ref"] + ".fai",
    output:
        directory(GENOME_SPLIT_DIR)
    resources:
        mem_mb          = config["resources"]["split_genome"]["mem_mb"],
        runtime         = config["resources"]["split_genome"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
    log:
        OUT + "/logs/split_genome.log"
    shell:
        """
        mkdir -p {output}
        cut -f1 {input.fai} | while read chr; do
            samtools faidx {input.fa} "$chr" \
              | awk '/^>/ {{print ">" substr($1,2); next}} {{print}}' \
              > {output}/${{chr}}.fa
        done 2>{log}
        """
