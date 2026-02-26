# ─────────────────────────────────────────────────────────────────────────────
# DAP-seq Snakemake pipeline
# ─────────────────────────────────────────────────────────────────────────────

configfile: "config/config.yaml"

container: "docker://nkoranda/dapseq:latest"

# Container-internal paths for tools bundled in the image
GEM_JAR       = "/opt/gem/gem.jar"
GEM_READ_DIST = "/opt/gem/Read_Distribution_default.txt"

# ─────────────────────────────────────────────────────────────────────────────
# Sample setup
# CONTROL is None when input_control is absent or null — all rules adapt.
# ─────────────────────────────────────────────────────────────────────────────
SAMPLES           = list(config["samples"].keys())
CONTROL           = config.get("input_control") or None
TREATMENT_SAMPLES = [s for s in SAMPLES if s != CONTROL]

wildcard_constraints:
    sample = "[^/]+",
    read   = "R[12]",

# ─────────────────────────────────────────────────────────────────────────────
# rule all
# ─────────────────────────────────────────────────────────────────────────────
rule all:
    input:
        config["genome_ref"] + ".sizes",
        expand("Trimmed_output/Fastqc/{sample}.{read}_fastqc.html",
               sample=SAMPLES, read=["R1", "R2"]),
        "multiqc_report.html",
        expand("bigWig/{sample}.bw",                          sample=SAMPLES),
        *(expand("bigWig/{sample}.peaks.bw",                  sample=TREATMENT_SAMPLES) if CONTROL else []),
        expand("compare_bed/{sample}.MACS_peak.bed",          sample=TREATMENT_SAMPLES),
        expand("compare_bed/{sample}.GEM_peak.bed",           sample=TREATMENT_SAMPLES),
        expand("compare_bed/{sample}.compare_peak.bed",       sample=TREATMENT_SAMPLES),


# ─────────────────────────────────────────────────────────────────────────────
# Step 0: Index reference genome
# ─────────────────────────────────────────────────────────────────────────────
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
    resources:
        mem_mb=8000, runtime=60,
        slurm_partition=config["slurm_partition"],
        slurm_account=config["slurm_account"],
    log:
        "logs/bwa_index.log"
    shell:
        """
        bwa index -a bwtsw {input} 2>{log}
        samtools faidx {input} 2>>{log}
        cut -f1,2 {input}.fai > {output.sizes}
        """


# ─────────────────────────────────────────────────────────────────────────────
# Step 1: Trim paired-end reads (ALL samples)
# ─────────────────────────────────────────────────────────────────────────────
rule trimmomatic:
    input:
        r1 = lambda wc: config["samples"][wc.sample]["r1"],
        r2 = lambda wc: config["samples"][wc.sample]["r2"],
    output:
        r1          = "Trimmed_output/{sample}.R1.fastq.gz",
        r2          = "Trimmed_output/{sample}.R2.fastq.gz",
        r1_unpaired = temp("Trimmed_output/{sample}.R1.unpaired.fastq.gz"),
        r2_unpaired = temp("Trimmed_output/{sample}.R2.unpaired.fastq.gz"),
    resources:
        mem_mb=4000, runtime=60,
        slurm_partition=config["slurm_partition"],
        slurm_account=config["slurm_account"],
    log:
        "logs/trimmomatic/{sample}.log"
    shell:
        """
        ADAPTERS=$(find $CONDA_PREFIX/share -name "TruSeq3-PE.fa" | head -1)
        trimmomatic PE -phred33 \
          {input.r1} {input.r2} \
          {output.r1} {output.r1_unpaired} \
          {output.r2} {output.r2_unpaired} \
          ILLUMINACLIP:$ADAPTERS:2:30:10 \
          LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 2>{log}
        """


# ─────────────────────────────────────────────────────────────────────────────
# Step 2: FastQC on trimmed R1 and R2 (ALL samples)
# ─────────────────────────────────────────────────────────────────────────────
rule fastqc:
    input:
        "Trimmed_output/{sample}.{read}.fastq.gz"
    output:
        html = "Trimmed_output/Fastqc/{sample}.{read}_fastqc.html",
        zip  = "Trimmed_output/Fastqc/{sample}.{read}_fastqc.zip",
    resources:
        mem_mb=2000, runtime=30,
        slurm_partition=config["slurm_partition"],
        slurm_account=config["slurm_account"],
    log:
        "logs/fastqc/{sample}.{read}.log"
    shell:
        "fastqc {input} --outdir=Trimmed_output/Fastqc/ 2>{log}"


# ─────────────────────────────────────────────────────────────────────────────
# Step 2.1: MultiQC
# ─────────────────────────────────────────────────────────────────────────────
rule multiqc:
    input:
        expand("Trimmed_output/Fastqc/{sample}.{read}_fastqc.zip",
               sample=SAMPLES, read=["R1", "R2"])
    output:
        "multiqc_report.html"
    resources:
        mem_mb=2000, runtime=30,
        slurm_partition=config["slurm_partition"],
        slurm_account=config["slurm_account"],
    log:
        "logs/multiqc.log"
    shell:
        "multiqc Trimmed_output/Fastqc/ -o . 2>{log}"


# ─────────────────────────────────────────────────────────────────────────────
# Step 3: Align with BWA MEM, paired-end (ALL samples)
# ─────────────────────────────────────────────────────────────────────────────
rule bwa_mem:
    input:
        r1  = "Trimmed_output/{sample}.R1.fastq.gz",
        r2  = "Trimmed_output/{sample}.R2.fastq.gz",
        idx = config["genome_ref"] + ".sizes",
    output:
        temp("Sam_File/{sample}.sam")
    threads:
        config["threads"]
    resources:
        mem_mb=16000, runtime=120,
        slurm_partition=config["slurm_partition"],
        slurm_account=config["slurm_account"],
    log:
        "logs/bwa_mem/{sample}.log"
    shell:
        "bwa mem -t {threads} -k 60 -B 7 -O 6 -v 0 {config[genome_ref]} {input.r1} {input.r2} > {output} 2>{log}"


# ─────────────────────────────────────────────────────────────────────────────
# Step 4: Filter, sort, dedup, index BAM (ALL samples)
# ─────────────────────────────────────────────────────────────────────────────
rule samtools_filter_sort_dedup:
    input:
        "Sam_File/{sample}.sam"
    output:
        bam = "Sam_sorted/{sample}._mapped_sorted.bam",
        bai = "Sam_sorted/{sample}._mapped_sorted.bam.bai",
    params:
        prefix = "Sam_sorted/{sample}._temp",
    threads:
        config["threads"]
    resources:
        mem_mb=8000, runtime=60,
        slurm_partition=config["slurm_partition"],
        slurm_account=config["slurm_account"],
    log:
        "logs/samtools/{sample}.log"
    shell:
        """
        samtools view -@ {threads} -h -F 4 -q 30 -u -S {input} \
          | samtools sort -@ {threads} - {params.prefix} 2>>{log}
        samtools rmdup -s {params.prefix}.bam {output.bam} 2>>{log}
        samtools index {output.bam} 2>>{log}
        rm -f {params.prefix}.bam
        """


# ─────────────────────────────────────────────────────────────────────────────
# Step 4.1: Coverage bigWig (ALL samples)
# ─────────────────────────────────────────────────────────────────────────────
rule bamcoverage:
    input:
        bam = "Sam_sorted/{sample}._mapped_sorted.bam",
        bai = "Sam_sorted/{sample}._mapped_sorted.bam.bai",
    output:
        "bigWig/{sample}.bw"
    threads:
        config["threads"]
    resources:
        mem_mb=8000, runtime=60,
        slurm_partition=config["slurm_partition"],
        slurm_account=config["slurm_account"],
    log:
        "logs/bamcoverage/{sample}.log"
    shell:
        "bamCoverage -b {input.bam} --normalizeUsing BPM --extendReads 300 -p {threads} -o {output} 2>{log}"


# ─────────────────────────────────────────────────────────────────────────────
# Step 4.2: Sample vs control ratio bigWig — only created when CONTROL is set
# ─────────────────────────────────────────────────────────────────────────────
if CONTROL:
    rule bamcompare:
        input:
            sample_bam  = "Sam_sorted/{sample}._mapped_sorted.bam",
            sample_bai  = "Sam_sorted/{sample}._mapped_sorted.bam.bai",
            control_bam = f"Sam_sorted/{CONTROL}._mapped_sorted.bam",
            control_bai = f"Sam_sorted/{CONTROL}._mapped_sorted.bam.bai",
        output:
            "bigWig/{sample}.peaks.bw"
        threads:
            config["threads"]
        resources:
            mem_mb=8000, runtime=60,
            slurm_partition=config["slurm_partition"],
            slurm_account=config["slurm_account"],
        log:
            "logs/bamcompare/{sample}.log"
        shell:
            """
            bamCompare \
              -b1 {input.sample_bam} -b2 {input.control_bam} \
              -o {output} \
              --binSize 200 --operation ratio \
              --scaleFactorsMethod SES -n 1000 \
              -p {threads} 2>{log}
            """


# ─────────────────────────────────────────────────────────────────────────────
# Step 5: MACS2 peak calling (TREATMENT samples)
# Runs with or without a control — -c flag omitted when CONTROL is null.
# ─────────────────────────────────────────────────────────────────────────────
rule macs2:
    input:
        sample_bam  = "Sam_sorted/{sample}._mapped_sorted.bam",
        control_bam = (f"Sam_sorted/{CONTROL}._mapped_sorted.bam" if CONTROL else []),
    output:
        summits      = "MACS/{sample}_summits.bed",
        narrowpeak   = "MACS/{sample}_peaks.narrowPeak",
        treat_pileup = "MACS/{sample}_treat_pileup.bdg",
        ctrl_lambda  = "MACS/{sample}_control_lambda.bdg",
    params:
        ctrl = lambda wc, input: f"-c {input.control_bam}" if CONTROL else "",
    resources:
        mem_mb=8000, runtime=60,
        slurm_partition=config["slurm_partition"],
        slurm_account=config["slurm_account"],
    log:
        "logs/macs2/{sample}.log"
    shell:
        """
        macs2 callpeak \
          -t {input.sample_bam} {params.ctrl} \
          -f BAM --outdir MACS \
          -g {config[genome_size]} -n {wildcards.sample} \
          -B -q 0.01 -m 2 50 --verbose=0 2>{log}
        """


# ─────────────────────────────────────────────────────────────────────────────
# Step 6: GEM peak calling (TREATMENT samples)
# Runs with or without a control — --ctrl flag omitted when CONTROL is null.
# ─────────────────────────────────────────────────────────────────────────────
rule gem:
    input:
        sample_bam  = "Sam_sorted/{sample}._mapped_sorted.bam",
        control_bam = (f"Sam_sorted/{CONTROL}._mapped_sorted.bam" if CONTROL else []),
        chrom_sizes = config["genome_ref"] + ".sizes",
    output:
        gem_events = "GEM_BED/{sample}/{sample}.GEM_events.txt",
        gps_events = "GEM_BED/{sample}/{sample}.GPS_events.txt",
    params:
        ctrl = lambda wc, input: f"--ctrl {input.control_bam}" if CONTROL else "",
    threads:
        config["threads"]
    resources:
        mem_mb=32000, runtime=240,
        slurm_partition=config["slurm_partition"],
        slurm_account=config["slurm_account"],
    log:
        "logs/gem/{sample}.log"
    shell:
        """
        java -Xmx16G -jar {GEM_JAR} \
          --t {threads} --f BAM \
          --d {GEM_READ_DIST} \
          --g {input.chrom_sizes} \
          --genome {config[genome_dir]} \
          --expt {input.sample_bam} {params.ctrl} \
          --outBED \
          --out GEM_BED/{wildcards.sample}/{wildcards.sample} \
          --k_min 6 --kmax 20 --k_seqs 600 --k_neg_dinu_shuffle 2>{log}
        """


# ─────────────────────────────────────────────────────────────────────────────
# Step 7: Combine MACS2 + GEM peaks
# ─────────────────────────────────────────────────────────────────────────────
rule combine_peaks:
    input:
        macs_summits = "MACS/{sample}_summits.bed",
        gem_events   = "GEM_BED/{sample}/{sample}.GEM_events.txt",
    output:
        combined_bed = "Combined_BED/{sample}.combined.bed",
        macs_bed     = "compare_bed/{sample}.MACS.bed",
        gem_bed      = "compare_bed/{sample}.GEM.bed",
    params:
        window_size = 200,
        min_score   = 1,
    resources:
        mem_mb=4000, runtime=30,
        slurm_partition=config["slurm_partition"],
        slurm_account=config["slurm_account"],
    log:
        "logs/combine_peaks/{sample}.log"
    script:
        "scripts/combine_peaks.py"


# ─────────────────────────────────────────────────────────────────────────────
# Step 7.2a: Extract FASTA sequences for peak regions
# ─────────────────────────────────────────────────────────────────────────────
rule getfasta:
    input:
        bed    = "Combined_BED/{sample}.combined.bed",
        genome = config["genome_ref"],
    output:
        "fasta/{sample}.fasta"
    resources:
        mem_mb=4000, runtime=30,
        slurm_partition=config["slurm_partition"],
        slurm_account=config["slurm_account"],
    log:
        "logs/getfasta/{sample}.log"
    shell:
        "bedtools getfasta -fo {output} -fi {input.genome} -bed {input.bed} 2>{log}"


# ─────────────────────────────────────────────────────────────────────────────
# Step 7.2b: Remove duplicate FASTA records
# ─────────────────────────────────────────────────────────────────────────────
rule dedup_fasta:
    input:
        "fasta/{sample}.fasta"
    output:
        "fasta/{sample}.fasta.nodup"
    resources:
        mem_mb=4000, runtime=30,
        slurm_partition=config["slurm_partition"],
        slurm_account=config["slurm_account"],
    log:
        "logs/dedup_fasta/{sample}.log"
    script:
        "scripts/dedup_fasta.py"


# ─────────────────────────────────────────────────────────────────────────────
# Step 7.3: Filter tandem repeats
# ─────────────────────────────────────────────────────────────────────────────
rule tandem_filter:
    input:
        "fasta/{sample}.fasta.nodup"
    output:
        "fasta/{sample}.fasta.filtered.fasta"
    params:
        k     = 6,
        k_max = 4,
    resources:
        mem_mb=4000, runtime=30,
        slurm_partition=config["slurm_partition"],
        slurm_account=config["slurm_account"],
    log:
        "logs/tandem_filter/{sample}.log"
    script:
        "scripts/tandem_filter.py"


# ─────────────────────────────────────────────────────────────────────────────
# Step 8: MEME motif discovery
# ─────────────────────────────────────────────────────────────────────────────
rule meme:
    input:
        "fasta/{sample}.fasta.filtered.fasta"
    output:
        txt  = "GEM_BED/{sample}-meme/meme.txt",
        logo = "GEM_BED/{sample}-meme/logo1.png",
    threads:
        config["threads"]
    resources:
        mem_mb=8000, runtime=120,
        slurm_partition=config["slurm_partition"],
        slurm_account=config["slurm_account"],
    log:
        "logs/meme/{sample}.log"
    shell:
        """
        meme {input} \
          -nmotifs 1 -minw 4 -maxw 12 \
          -dna -mod oops -nostatus \
          -p {threads} -oc GEM_BED/{wildcards.sample}-meme 2>{log}
        """


# ─────────────────────────────────────────────────────────────────────────────
# Step 9: FIMO motif scanning
# ─────────────────────────────────────────────────────────────────────────────
rule fimo:
    input:
        meme_txt = "GEM_BED/{sample}-meme/meme.txt",
        genome   = config["genome_ref"],
    output:
        bed = "GEM_BED/{sample}-fimo/fimo.bed",
        gff = "GEM_BED/{sample}-fimo/fimo.gff",
    resources:
        mem_mb=8000, runtime=60,
        slurm_partition=config["slurm_partition"],
        slurm_account=config["slurm_account"],
    log:
        "logs/fimo/{sample}.log"
    shell:
        """
        fimo -verbosity 1 --thresh 1e-5 \
          -oc GEM_BED/{wildcards.sample}-fimo \
          {input.meme_txt} {input.genome} 2>{log}
        """


# ─────────────────────────────────────────────────────────────────────────────
# Step 9.2: Intersect motif sites with peak beds
# ─────────────────────────────────────────────────────────────────────────────
rule bedtools_intersect:
    input:
        fimo_bed = "GEM_BED/{sample}-fimo/fimo.bed",
        macs_bed = "compare_bed/{sample}.MACS.bed",
        gem_bed  = "compare_bed/{sample}.GEM.bed",
    output:
        macs_peak    = "compare_bed/{sample}.MACS_peak.bed",
        gem_peak     = "compare_bed/{sample}.GEM_peak.bed",
        compare_peak = "compare_bed/{sample}.compare_peak.bed",
    resources:
        mem_mb=4000, runtime=30,
        slurm_partition=config["slurm_partition"],
        slurm_account=config["slurm_account"],
    log:
        "logs/bedtools_intersect/{sample}.log"
    shell:
        """
        bedtools intersect -wa -a {input.fimo_bed} -b {input.macs_bed} > {output.macs_peak} 2>{log}
        bedtools intersect -wa -a {input.fimo_bed} -b {input.gem_bed}  > {output.gem_peak}  2>>{log}
        bedtools intersect -wa -a {input.macs_bed} -b {input.gem_bed}  > {output.compare_peak} 2>>{log}
        """
