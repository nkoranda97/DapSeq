rule align_se:
    wildcard_constraints:
        sample = "|".join(sorted(SE_SAMPLES)) if SE_SAMPLES else "(?!)",
    input:
        r1  = OUT + "/trimmed/{sample}.R1.fastq.gz",
        idx = config["genome_ref"] + ".1.bt2",
    output:
        bam        = OUT + "/bam/{sample}.full.bam",
        bai        = OUT + "/bam/{sample}.full.bam.bai",
        bw         = OUT + "/bigWig/{sample}.full.bw",
        align_rate = OUT + "/stats/{sample}.align_rate.txt",
    params:
        idx             = config["genome_ref"],
        mapq            = config["samtools"]["mapq"],
        extra           = config["bowtie2"].get("extra", ""),
        bw_normalize    = config["bamcoverage"].get("normalize_using", "RPGC"),
        bw_bin_size     = config["bamcoverage"].get("bin_size", 1),
        bw_ignore_dups  = "--ignoreDuplicates"q if config["bamcoverage"].get("ignore_duplicates", True) else "",
        bw_extra        = config["bamcoverage"].get("extra", ""),
        bw_tempdir      = OUT + "/temp",
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

        mkdir -p {params.bw_tempdir}
        export TMPDIR={params.bw_tempdir}
        bamCoverage -b {output.bam} -o {output.bw} -p {threads} \
          --normalizeUsing {params.bw_normalize} \
          --effectiveGenomeSize {config[genome_size]} \
          --binSize {params.bw_bin_size} {params.bw_ignore_dups} \
          {params.bw_extra} \
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
        bam              = OUT + "/bam/{sample}.full.bam",
        bai              = OUT + "/bam/{sample}.full.bam.bai",
        bw               = OUT + "/bigWig/{sample}.full.bw",
        align_rate       = OUT + "/stats/{sample}.align_rate.txt",
        median_frag_size = OUT + "/stats/{sample}.median_frag_size.txt",
    params:
        idx             = config["genome_ref"],
        mapq            = config["samtools"]["mapq"],
        extra           = config["bowtie2"].get("extra", ""),
        bw_normalize    = config["bamcoverage"].get("normalize_using", "RPGC"),
        bw_bin_size     = config["bamcoverage"].get("bin_size", 1),
        bw_ignore_dups  = "--ignoreDuplicates" if config["bamcoverage"].get("ignore_duplicates", True) else "",
        bw_max_frag     = config["bamcoverage"].get("max_fragment_length", 600),
        bw_extend_reads = "--extendReads" if config["bamcoverage"].get("extend_reads", True) else "",
        bw_extra        = config["bamcoverage"].get("extra", ""),
        bw_tempdir      = OUT + "/temp",
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

        mkdir -p {params.bw_tempdir}
        export TMPDIR={params.bw_tempdir}
        bamCoverage -b {output.bam} -o {output.bw} -p {threads} \
          --normalizeUsing {params.bw_normalize} \
          --effectiveGenomeSize {config[genome_size]} \
          --binSize {params.bw_bin_size} {params.bw_ignore_dups} \
          --maxFragmentLength {params.bw_max_frag} {params.bw_extend_reads} \
          {params.bw_extra} \
          2>{log.bamcoverage}
        """


rule filter_chr:
    input:
        bam = OUT + "/bam/{sample}.full.bam",
        bai = OUT + "/bam/{sample}.full.bam.bai",
    output:
        bam = OUT + "/bam/{sample}.chr.bam",
        bai = OUT + "/bam/{sample}.chr.bam.bai",
        bw  = OUT + "/bigWig/{sample}.chr.bw",
    params:
        pattern        = config["contig_filter"]["pattern"],
        bw_normalize   = config["bamcoverage"].get("normalize_using", "RPGC"),
        bw_bin_size    = config["bamcoverage"].get("bin_size", 1),
        bw_ignore_dups = "--ignoreDuplicates" if config["bamcoverage"].get("ignore_duplicates", True) else "",
        bw_extra       = config["bamcoverage"].get("extra", ""),
        bw_tempdir     = OUT + "/temp",
    threads:
        config["threads"]
    resources:
        mem_mb          = config["resources"]["filter_chr"]["mem_mb"],
        runtime         = config["resources"]["filter_chr"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
    log:
        filter      = OUT + "/logs/filter_chr/{sample}.log",
        bamcoverage = OUT + "/logs/bamcoverage/{sample}.chr.log",
    shell:
        """
        samtools view -b -@ {threads} \
          -e 'rname =~ "{params.pattern}"' \
          -o {output.bam} {input.bam} 2>{log.filter}
        samtools index {output.bam}

        mkdir -p {params.bw_tempdir}
        export TMPDIR={params.bw_tempdir}
        bamCoverage -b {output.bam} -o {output.bw} -p {threads} \
          --normalizeUsing {params.bw_normalize} \
          --effectiveGenomeSize {config[genome_size]} \
          --binSize {params.bw_bin_size} {params.bw_ignore_dups} \
          {params.bw_extra} 2>{log.bamcoverage}
        """
