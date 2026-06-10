_aligner = config.get("aligner", "bowtie2")

if _aligner == "bowtie2":
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
        params:
            idx             = config["genome_ref"],
            mapq            = config["samtools"]["mapq"],
            extra           = config["bowtie2"].get("extra", ""),
            bw_normalize    = config["bamcoverage"].get("normalize_using", "RPGC"),
            bw_bin_size     = config["bamcoverage"].get("bin_size", 1),
            bw_ignore_dups  = "--ignoreDuplicates" if config["bamcoverage"].get("ignore_duplicates", True) else "",
            bw_extra        = config["bamcoverage"].get("extra", ""),
            bw_tempdir      = OUT + "/temp",
            genome_size     = config["genome_size"],
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
            set -euo pipefail
            bowtie2 \
              --reorder \
              --threads {threads} \
              -x {params.idx} \
              -U {input.r1} \
              {params.extra} \
              2>{log.align} \
            | samtools view -h -F 1804 -q {params.mapq} -b - \
            | samtools sort -@ {threads} -o {output.bam} -
            samtools index {output.bam}

            mkdir -p {params.bw_tempdir}
            export TMPDIR={params.bw_tempdir}
            bamCoverage -b {output.bam} -o {output.bw} -p {threads} \
              --normalizeUsing {params.bw_normalize} \
              --effectiveGenomeSize {params.genome_size} \
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
            bam              = OUT + "/bam/{sample}.bam",
            bai              = OUT + "/bam/{sample}.bam.bai",
            bw               = OUT + "/bigWig/{sample}.bw",
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
            genome_size     = config["genome_size"],
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
            set -euo pipefail
            bowtie2 \
              --reorder \
              --threads {threads} \
              -x {params.idx} \
              -1 {input.r1} -2 {input.r2} \
              --no-mixed --no-discordant \
              {params.extra} \
              2>{log.align} \
            | samtools view -h -F 1804 -q {params.mapq} -b - \
            | samtools sort -@ {threads} -o {output.bam} -
            samtools index {output.bam}

            mkdir -p {params.bw_tempdir}
            export TMPDIR={params.bw_tempdir}
            bamCoverage -b {output.bam} -o {output.bw} -p {threads} \
              --normalizeUsing {params.bw_normalize} \
              --effectiveGenomeSize {params.genome_size} \
              --binSize {params.bw_bin_size} {params.bw_ignore_dups} \
              --maxFragmentLength {params.bw_max_frag} {params.bw_extend_reads} \
              {params.bw_extra} \
              2>{log.bamcoverage}
            """

elif _aligner == "bwa_mem2":
    rule align_se:
        wildcard_constraints:
            sample = "|".join(sorted(SE_SAMPLES)) if SE_SAMPLES else "(?!)",
        input:
            r1  = OUT + "/trimmed/{sample}.R1.fastq.gz",
            idx = config["genome_ref"] + ".0123",
        output:
            bam        = OUT + "/bam/{sample}.bam",
            bai        = OUT + "/bam/{sample}.bam.bai",
            bw         = OUT + "/bigWig/{sample}.bw",
        params:
            idx             = config["genome_ref"],
            mapq            = config["samtools"]["mapq"],
            extra           = config["bwa_mem2"].get("extra", ""),
            bw_normalize    = config["bamcoverage"].get("normalize_using", "RPGC"),
            bw_bin_size     = config["bamcoverage"].get("bin_size", 1),
            bw_ignore_dups  = "--ignoreDuplicates" if config["bamcoverage"].get("ignore_duplicates", True) else "",
            bw_extra        = config["bamcoverage"].get("extra", ""),
            bw_tempdir      = OUT + "/temp",
            genome_size     = config["genome_size"],
        threads:
            config["threads"]
        resources:
            mem_mb          = config["resources"]["bwa_mem2_align"]["mem_mb"],
            runtime         = config["resources"]["bwa_mem2_align"]["runtime"],
            slurm_partition = config["slurm_partition"],
            slurm_account   = config["slurm_account"],
        log:
            align       = OUT + "/logs/bwa_mem2/{sample}.log",
            bamcoverage = OUT + "/logs/bamcoverage/{sample}.log",
        shell:
            """
            set -euo pipefail
            bwa-mem2 mem \
              -t {threads} \
              {params.idx} \
              {input.r1} \
              {params.extra} \
              2>{log.align} \
            | samtools view -h -F 1804 -q {params.mapq} -b - \
            | samtools sort -@ {threads} -o {output.bam} -
            samtools index {output.bam}

            mkdir -p {params.bw_tempdir}
            export TMPDIR={params.bw_tempdir}
            bamCoverage -b {output.bam} -o {output.bw} -p {threads} \
              --normalizeUsing {params.bw_normalize} \
              --effectiveGenomeSize {params.genome_size} \
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
            idx = config["genome_ref"] + ".0123",
        output:
            bam              = OUT + "/bam/{sample}.bam",
            bai              = OUT + "/bam/{sample}.bam.bai",
            bw               = OUT + "/bigWig/{sample}.bw",
        params:
            idx             = config["genome_ref"],
            mapq            = config["samtools"]["mapq"],
            extra           = config["bwa_mem2"].get("extra", ""),
            bw_normalize    = config["bamcoverage"].get("normalize_using", "RPGC"),
            bw_bin_size     = config["bamcoverage"].get("bin_size", 1),
            bw_ignore_dups  = "--ignoreDuplicates" if config["bamcoverage"].get("ignore_duplicates", True) else "",
            bw_max_frag     = config["bamcoverage"].get("max_fragment_length", 600),
            bw_extend_reads = "--extendReads" if config["bamcoverage"].get("extend_reads", True) else "",
            bw_extra        = config["bamcoverage"].get("extra", ""),
            bw_tempdir      = OUT + "/temp",
            genome_size     = config["genome_size"],
        threads:
            config["threads"]
        resources:
            mem_mb          = config["resources"]["bwa_mem2_align"]["mem_mb"],
            runtime         = config["resources"]["bwa_mem2_align"]["runtime"],
            slurm_partition = config["slurm_partition"],
            slurm_account   = config["slurm_account"],
        log:
            align       = OUT + "/logs/bwa_mem2/{sample}.log",
            bamcoverage = OUT + "/logs/bamcoverage/{sample}.log",
        shell:
            """
            set -euo pipefail
            bwa-mem2 mem \
              -t {threads} \
              {params.idx} \
              {input.r1} {input.r2} \
              {params.extra} \
              2>{log.align} \
            | samtools view -h -F 1804 -q {params.mapq} -b - \
            | samtools sort -@ {threads} -o {output.bam} -
            samtools index {output.bam}

            mkdir -p {params.bw_tempdir}
            export TMPDIR={params.bw_tempdir}
            bamCoverage -b {output.bam} -o {output.bw} -p {threads} \
              --normalizeUsing {params.bw_normalize} \
              --effectiveGenomeSize {params.genome_size} \
              --binSize {params.bw_bin_size} {params.bw_ignore_dups} \
              --maxFragmentLength {params.bw_max_frag} {params.bw_extend_reads} \
              {params.bw_extra} \
              2>{log.bamcoverage}
            """

else:
    raise ValueError(
        f"Unknown aligner: '{_aligner}'. Must be 'bowtie2' or 'bwa_mem2'."
    )
