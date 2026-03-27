rule align_se:
    wildcard_constraints:
        sample = "|".join(sorted(SE_SAMPLES)) if SE_SAMPLES else "(?!)",
    input:
        r1  = OUT + "/trimmed/{sample}.R1.fastq.gz",
        idx = config["genome_ref"] + ".1.bt2",
    output:
        bam        = OUT + "/bam/{sample}.bam",
        bai        = OUT + "/bam/{sample}.bam.bai",
        bw         = OUT + "/bigWig/{sample}.bw",
        align_rate = OUT + "/stats/{sample}.align_rate.txt",
    params:
        idx   = config["genome_ref"],
        mapq  = config["samtools"]["mapq"],
        extra = config["bowtie2"].get("extra", ""),
    threads:
        config["threads"]
    resources:
        mem_mb          = config["resources"]["trim_align"]["mem_mb"],
        runtime         = config["resources"]["trim_align"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
    log:
        align       = OUT + "/logs/bowtie2/{sample}.log",
        bamcoverage = OUT + "/logs/bamcoverage/{sample}.log",
    shell:
        """
        bowtie2 \
          --reorder \
          --threads {threads} \
          -x {params.idx} \
          -U {input.r1} \
          {params.extra} \
          2>{log.align} \
        | samtools view -h -F 4 -q {params.mapq} -b - \
        | samtools sort -@ {threads} -o {output.bam} -
        samtools index {output.bam}

        grep "overall alignment rate" {log.align} \
          | awk -v FS="% " '{{print $1}}' > {output.align_rate}

        bamCoverage -b {output.bam} -o {output.bw} -p {threads} \
          --normalizeUsing RPGC \
          --effectiveGenomeSize {config[genome_size]} \
          --binSize 1 --ignoreDuplicates \
          2>{log.bamcoverage}
        """


rule align_pe:
    wildcard_constraints:
        sample = "|".join(sorted(PE_SAMPLES)) if PE_SAMPLES else "(?!)",
    input:
        r1  = OUT + "/trimmed/{sample}.R1.fastq.gz",
        r2  = OUT + "/trimmed/{sample}.R2.fastq.gz",
        idx = config["genome_ref"] + ".1.bt2",
    output:
        bam              = OUT + "/bam/{sample}.bam",
        bai              = OUT + "/bam/{sample}.bam.bai",
        bw               = OUT + "/bigWig/{sample}.bw",
        align_rate       = OUT + "/stats/{sample}.align_rate.txt",
        median_frag_size = OUT + "/stats/{sample}.median_frag_size.txt",
    params:
        idx   = config["genome_ref"],
        mapq  = config["samtools"]["mapq"],
        extra = config["bowtie2"].get("extra", ""),
    threads:
        config["threads"]
    resources:
        mem_mb          = config["resources"]["trim_align"]["mem_mb"],
        runtime         = config["resources"]["trim_align"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
    log:
        align       = OUT + "/logs/bowtie2/{sample}.log",
        bamcoverage = OUT + "/logs/bamcoverage/{sample}.log",
    shell:
        """
        bowtie2 \
          --reorder \
          --threads {threads} \
          -x {params.idx} \
          -1 {input.r1} -2 {input.r2} \
          --no-mixed --no-discordant \
          {params.extra} \
          2>{log.align} \
        | samtools view -h -F 4 -q {params.mapq} -b - \
        | samtools sort -@ {threads} -o {output.bam} -
        samtools index {output.bam}

        grep "overall alignment rate" {log.align} \
          | awk -v FS="% " '{{print $1}}' > {output.align_rate}

        bamPEFragmentSize -b {output.bam} \
          | grep Median: | awk 'NR==1 {{print $2}}' > {output.median_frag_size}

        bamCoverage -b {output.bam} -o {output.bw} -p {threads} \
          --normalizeUsing RPGC \
          --effectiveGenomeSize {config[genome_size]} \
          --binSize 1 --ignoreDuplicates --maxFragmentLength 600 --extendReads \
          2>{log.bamcoverage}
        """
