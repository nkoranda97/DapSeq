_aligner = config.get("aligner", "bowtie2")

if _aligner == "bowtie2":
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
        threads:
            config["threads"]
        resources:
            mem_mb          = config["resources"]["bowtie2_index"]["mem_mb"],
            runtime         = config["resources"]["bowtie2_index"]["runtime"],
        log:
            OUT + "/logs/bowtie2_index.log"
        shell:
            """
            set -euo pipefail
            bowtie2-build --threads {threads} {params.extra_build} {input} {input} 2>{log}
            samtools faidx {input} 2>{log}
            cut -f1,2 {input}.fai > {output.sizes}
            """

elif _aligner == "bwa_mem2":
    rule bwa_mem2_index:
        input:
            config["genome_ref"]
        output:
            idx    = config["genome_ref"] + ".0123",
            amb    = config["genome_ref"] + ".amb",
            ann    = config["genome_ref"] + ".ann",
            bwt    = config["genome_ref"] + ".bwt.2bit.64",
            pac    = config["genome_ref"] + ".pac",
            sizes  = config["genome_ref"] + ".sizes",
            fai    = config["genome_ref"] + ".fai",
        params:
            extra_index = config["bwa_mem2"].get("extra_index", ""),
        threads:
            config["threads"]
        resources:
            mem_mb          = config["resources"]["bwa_mem2_index"]["mem_mb"],
            runtime         = config["resources"]["bwa_mem2_index"]["runtime"],
        log:
            OUT + "/logs/bwa_mem2_index.log"
        shell:
            """
            set -euo pipefail
            bwa-mem2 index {params.extra_index} {input} 2>{log}
            samtools faidx {input} 2>>{log}
            cut -f1,2 {input}.fai > {output.sizes}
            """
