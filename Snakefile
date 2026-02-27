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
SAMPLES           = [s for s in config["samples"] if config["samples"][s]["r1"] is not None]
CONTROL           = config.get("input_control") or None
TREATMENT_SAMPLES = [s for s in SAMPLES if s != CONTROL]
OUT               = config["output_dir"]
GENOME_SPLIT_DIR  = config["genome_dir"] + "/split"

wildcard_constraints:
    sample = "[^/]+",
    read   = "R[12]",

# ─────────────────────────────────────────────────────────────────────────────
# rule all
# ─────────────────────────────────────────────────────────────────────────────
rule all:
    input:
        config["genome_ref"] + ".sizes",
        expand(OUT + "/Fastqc/{sample}.{read}_fastqc.html",
               sample=SAMPLES, read=["R1", "R2"]),
        OUT + "/multiqc_report.html",
        expand(OUT + "/bigWig/{sample}.bw",                      sample=SAMPLES),
        *(expand(OUT + "/bigWig/{sample}.peaks.bw",              sample=TREATMENT_SAMPLES) if CONTROL else []),
        expand(OUT + "/compare_bed/{sample}.MACS_peak.bed",      sample=TREATMENT_SAMPLES),
        expand(OUT + "/compare_bed/{sample}.GEM_peak.bed",       sample=TREATMENT_SAMPLES),
        expand(OUT + "/compare_bed/{sample}.compare_peak.bed",   sample=TREATMENT_SAMPLES),
        expand(OUT + "/meme/{sample}-intersection/meme.txt",      sample=TREATMENT_SAMPLES),


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
        OUT + "/logs/bwa_index.log"
    shell:
        """
        bwa index -a bwtsw {input} 2>{log}
        samtools faidx {input} 2>>{log}
        cut -f1,2 {input}.fai > {output.sizes}
        """


# ─────────────────────────────────────────────────────────────────────────────
# Step 0.1: Split genome into per-chromosome FASTAs for GEM
# ─────────────────────────────────────────────────────────────────────────────
rule split_genome:
    input:
        config["genome_ref"]
    output:
        directory(GENOME_SPLIT_DIR)
    resources:
        mem_mb=8000, runtime=60,
        slurm_partition=config["slurm_partition"],
        slurm_account=config["slurm_account"],
    log:
        OUT + "/logs/split_genome.log"
    shell:
        """
        mkdir -p {output}
        awk '/^>/ {{ filename="{output}/" substr($1,2) ".fa"; print > filename; next }} {{ print >> filename }}' {input} 2>{log}
        """


# ─────────────────────────────────────────────────────────────────────────────
# Step 1: Trim paired-end reads (ALL samples)
# ─────────────────────────────────────────────────────────────────────────────
rule trimmomatic:
    input:
        r1 = lambda wc: config["samples"][wc.sample]["r1"],
        r2 = lambda wc: config["samples"][wc.sample]["r2"],
    output:
        r1          = OUT + "/trimmed/{sample}.R1.fastq.gz",
        r2          = OUT + "/trimmed/{sample}.R2.fastq.gz",
        r1_unpaired = temp(OUT + "/trimmed/{sample}.R1.unpaired.fastq.gz"),
        r2_unpaired = temp(OUT + "/trimmed/{sample}.R2.unpaired.fastq.gz"),
    resources:
        mem_mb=4000, runtime=60,
        slurm_partition=config["slurm_partition"],
        slurm_account=config["slurm_account"],
    log:
        OUT + "/logs/trimmomatic/{sample}.log"
    shell:
        """
        ADAPTERS=$(find /opt/conda/envs/dapseq/share -name "TruSeq3-PE.fa" | head -1)
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
        OUT + "/trimmed/{sample}.{read}.fastq.gz"
    output:
        html = OUT + "/Fastqc/{sample}.{read}_fastqc.html",
        zip  = OUT + "/Fastqc/{sample}.{read}_fastqc.zip",
    resources:
        mem_mb=2000, runtime=30,
        slurm_partition=config["slurm_partition"],
        slurm_account=config["slurm_account"],
    log:
        OUT + "/logs/fastqc/{sample}.{read}.log"
    shell:
        "fastqc {input} --outdir={OUT}/Fastqc/ 2>{log}"


# ─────────────────────────────────────────────────────────────────────────────
# Step 2.1: MultiQC
# ─────────────────────────────────────────────────────────────────────────────
rule multiqc:
    input:
        expand(OUT + "/Fastqc/{sample}.{read}_fastqc.zip",
               sample=SAMPLES, read=["R1", "R2"])
    output:
        OUT + "/multiqc_report.html"
    resources:
        mem_mb=2000, runtime=30,
        slurm_partition=config["slurm_partition"],
        slurm_account=config["slurm_account"],
    log:
        OUT + "/logs/multiqc.log"
    shell:
        "multiqc {OUT}/Fastqc/ -o {OUT}/ 2>{log}"


# ─────────────────────────────────────────────────────────────────────────────
# Step 3: Align with BWA MEM (ALL samples)
# ─────────────────────────────────────────────────────────────────────────────
rule bwa_mem:
    input:
        r1  = OUT + "/trimmed/{sample}.R1.fastq.gz",
        r2  = OUT + "/trimmed/{sample}.R2.fastq.gz",
        idx = config["genome_ref"] + ".sizes",
    output:
        temp(OUT + "/sam/{sample}.sam")
    threads:
        config["threads"]
    resources:
        mem_mb=16000, runtime=120,
        slurm_partition=config["slurm_partition"],
        slurm_account=config["slurm_account"],
    log:
        OUT + "/logs/bwa_mem/{sample}.log"
    shell:
        "bwa mem -t {threads} -k 60 -B 7 -O 6 -v 0 {config[genome_ref]} {input.r1} {input.r2} > {output} 2>{log}"


# ─────────────────────────────────────────────────────────────────────────────
# Step 4: Filter, sort, dedup, index BAM (ALL samples)
# ─────────────────────────────────────────────────────────────────────────────
rule samtools_filter_sort_dedup:
    input:
        OUT + "/sam/{sample}.sam"
    output:
        bam = OUT + "/bam/{sample}.bam",
        bai = OUT + "/bam/{sample}.bam.bai",
    params:
        prefix = OUT + "/bam/{sample}.tmp.bam",
    threads:
        config["threads"]
    resources:
        mem_mb=8000, runtime=60,
        slurm_partition=config["slurm_partition"],
        slurm_account=config["slurm_account"],
    log:
        OUT + "/logs/samtools/{sample}.log"
    shell:
        """
        samtools view -@ {threads} -h -F 4 -q 30 -u {input} \
          | samtools fixmate -m -@ {threads} - - \
          | samtools sort -@ {threads} -o {params.prefix} - 2>>{log}
        samtools markdup -r -@ {threads} {params.prefix} {output.bam} 2>>{log}
        samtools index {output.bam} 2>>{log}
        rm -f {params.prefix}
        """


# ─────────────────────────────────────────────────────────────────────────────
# Step 4.1: Coverage bigWig (ALL samples)
# ─────────────────────────────────────────────────────────────────────────────
rule bamcoverage:
    input:
        bam = OUT + "/bam/{sample}.bam",
        bai = OUT + "/bam/{sample}.bam.bai",
    output:
        OUT + "/bigWig/{sample}.bw"
    threads:
        config["threads"]
    resources:
        mem_mb=8000, runtime=60,
        slurm_partition=config["slurm_partition"],
        slurm_account=config["slurm_account"],
    log:
        OUT + "/logs/bamcoverage/{sample}.log"
    shell:
        "bamCoverage -b {input.bam} --normalizeUsing BPM --extendReads 300 -p {threads} -o {output} 2>{log}"


# ─────────────────────────────────────────────────────────────────────────────
# Step 4.2: Sample vs control ratio bigWig (only when CONTROL is set)
# ─────────────────────────────────────────────────────────────────────────────
if CONTROL:
    rule bamcompare:
        input:
            sample_bam  = OUT + "/bam/{sample}.bam",
            sample_bai  = OUT + "/bam/{sample}.bam.bai",
            control_bam = OUT + f"/bam/{CONTROL}.bam",
            control_bai = OUT + f"/bam/{CONTROL}.bam.bai",
        output:
            OUT + "/bigWig/{sample}.peaks.bw"
        threads:
            config["threads"]
        resources:
            mem_mb=8000, runtime=60,
            slurm_partition=config["slurm_partition"],
            slurm_account=config["slurm_account"],
        log:
            OUT + "/logs/bamcompare/{sample}.log"
        shell:
            """
            bamCompare \
              -b1 {input.sample_bam} -b2 {input.control_bam} \
              -o {output} \
              --binSize 80 --operation ratio \
              --scaleFactorsMethod SES -n 1000 \
              -p {threads} 2>{log}
            """


# ─────────────────────────────────────────────────────────────────────────────
# Step 5: MACS3 peak calling (TREATMENT samples)
# ─────────────────────────────────────────────────────────────────────────────
rule macs3:
    input:
        sample_bam  = OUT + "/bam/{sample}.bam",
        control_bam = (OUT + f"/bam/{CONTROL}.bam" if CONTROL else []),
    output:
        summits      = OUT + "/MACS/{sample}_summits.bed",
        narrowpeak   = OUT + "/MACS/{sample}_peaks.narrowPeak",
        treat_pileup = OUT + "/MACS/{sample}_treat_pileup.bdg",
        ctrl_lambda  = OUT + "/MACS/{sample}_control_lambda.bdg",
    params:
        ctrl   = lambda wc, input: f"-c {input.control_bam}" if CONTROL else "",
        outdir = OUT + "/MACS",
    resources:
        mem_mb=8000, runtime=60,
        slurm_partition=config["slurm_partition"],
        slurm_account=config["slurm_account"],
    log:
        OUT + "/logs/macs3/{sample}.log"
    shell:
        """
        macs3 callpeak \
          -t {input.sample_bam} {params.ctrl} \
          -f BAM --outdir {params.outdir} \
          -g {config[genome_size]} -n {wildcards.sample} \
          -B -q 0.01 -m 2 50 --verbose=0 2>{log}
        """


# ─────────────────────────────────────────────────────────────────────────────
# Step 6: GEM peak calling (TREATMENT samples)
# ─────────────────────────────────────────────────────────────────────────────
rule gem:
    input:
        sample_bam  = OUT + "/bam/{sample}.bam",
        control_bam = (OUT + f"/bam/{CONTROL}.bam" if CONTROL else []),
        chrom_sizes = config["genome_ref"] + ".sizes",
        genome_dir  = GENOME_SPLIT_DIR,
    output:
        gem_events = OUT + "/GEM/{sample}/{sample}.GEM_events.txt",
        gps_events = OUT + "/GEM/{sample}/{sample}.GPS_events.txt",
    params:
        ctrl       = lambda wc, input: f"--ctrl {input.control_bam}" if CONTROL else "",
        out_prefix = OUT + "/GEM/{sample}/{sample}",
    threads:
        config["threads"]
    resources:
        mem_mb=32000, runtime=240,
        slurm_partition=config["slurm_partition"],
        slurm_account=config["slurm_account"],
    log:
        OUT + "/logs/gem/{sample}.log"
    shell:
        """
        java -Xmx16G -jar {GEM_JAR} \
          --t {threads} --f BAM \
          --d {GEM_READ_DIST} \
          --g {input.chrom_sizes} \
          --genome {input.genome_dir} \
          --expt {input.sample_bam} {params.ctrl} \
          --outBED \
          --out {params.out_prefix} \
          --k_min 6 --k_max 20 --k_seqs 600 --k_neg_dinu_shuffle 2>{log}
        """


# ─────────────────────────────────────────────────────────────────────────────
# Step 7: Combine MACS3 + GEM peaks
# ─────────────────────────────────────────────────────────────────────────────
rule combine_peaks:
    input:
        macs_summits = OUT + "/MACS/{sample}_summits.bed",
        gem_events   = OUT + "/GEM/{sample}/{sample}.GEM_events.txt",
    output:
        combined_bed = OUT + "/combined_bed/{sample}.combined.bed",
        macs_bed     = OUT + "/compare_bed/{sample}.MACS.bed",
        gem_bed      = OUT + "/compare_bed/{sample}.GEM.bed",
    params:
        window_size = 80,
        min_score   = 1,
    resources:
        mem_mb=4000, runtime=30,
        slurm_partition=config["slurm_partition"],
        slurm_account=config["slurm_account"],
    log:
        OUT + "/logs/combine_peaks/{sample}.log"
    script:
        "scripts/combine_peaks.py"


# ─────────────────────────────────────────────────────────────────────────────
# Step 7.2a: Extract FASTA sequences for peak regions
# ─────────────────────────────────────────────────────────────────────────────
rule getfasta:
    input:
        bed    = OUT + "/combined_bed/{sample}.combined.bed",
        genome = config["genome_ref"],
    output:
        OUT + "/fasta/{sample}.fasta"
    resources:
        mem_mb=4000, runtime=30,
        slurm_partition=config["slurm_partition"],
        slurm_account=config["slurm_account"],
    log:
        OUT + "/logs/getfasta/{sample}.log"
    shell:
        "bedtools getfasta -fo {output} -fi {input.genome} -bed {input.bed} 2>{log}"


# ─────────────────────────────────────────────────────────────────────────────
# Step 7.2b: Remove duplicate FASTA records
# ─────────────────────────────────────────────────────────────────────────────
rule dedup_fasta:
    input:
        OUT + "/fasta/{sample}.fasta"
    output:
        OUT + "/fasta/{sample}.fasta.nodup"
    resources:
        mem_mb=4000, runtime=30,
        slurm_partition=config["slurm_partition"],
        slurm_account=config["slurm_account"],
    log:
        OUT + "/logs/dedup_fasta/{sample}.log"
    script:
        "scripts/dedup_fasta.py"


# ─────────────────────────────────────────────────────────────────────────────
# Step 7.3: Filter tandem repeats
# ─────────────────────────────────────────────────────────────────────────────
rule tandem_filter:
    input:
        OUT + "/fasta/{sample}.fasta.nodup"
    output:
        OUT + "/fasta/{sample}.fasta.filtered.fasta"
    params:
        k     = 6,
        k_max = 4,
    resources:
        mem_mb=4000, runtime=30,
        slurm_partition=config["slurm_partition"],
        slurm_account=config["slurm_account"],
    log:
        OUT + "/logs/tandem_filter/{sample}.log"
    script:
        "scripts/tandem_filter.py"


# ─────────────────────────────────────────────────────────────────────────────
# Step 8: MEME motif discovery
# ─────────────────────────────────────────────────────────────────────────────
rule meme:
    input:
        OUT + "/fasta/{sample}.fasta.filtered.fasta"
    output:
        txt  = OUT + "/meme/{sample}/meme.txt",
        logo = OUT + "/meme/{sample}/logo1.png",
    params:
        outdir = OUT + "/meme/{sample}",
    threads:
        config["threads"]
    resources:
        mem_mb=8000, runtime=120,
        slurm_partition=config["slurm_partition"],
        slurm_account=config["slurm_account"],
    log:
        OUT + "/logs/meme/{sample}.log"
    shell:
        """
        meme {input} \
          -nmotifs 1 -minw 4 -maxw 12 \
          -dna -mod oops -nostatus \
          -p {threads} -oc {params.outdir} 2>{log}
        """


# ─────────────────────────────────────────────────────────────────────────────
# Step 9: FIMO motif scanning
# ─────────────────────────────────────────────────────────────────────────────
rule fimo:
    input:
        meme_txt = OUT + "/meme/{sample}/meme.txt",
        genome   = config["genome_ref"],
    output:
        bed = OUT + "/fimo/{sample}/fimo.bed",
        gff = OUT + "/fimo/{sample}/fimo.gff",
    params:
        outdir = OUT + "/fimo/{sample}",
    resources:
        mem_mb=8000, runtime=60,
        slurm_partition=config["slurm_partition"],
        slurm_account=config["slurm_account"],
    log:
        OUT + "/logs/fimo/{sample}.log"
    shell:
        """
        fimo -verbosity 1 --thresh 1e-5 \
          -oc {params.outdir} \
          {input.meme_txt} {input.genome} 2>{log}
        """


# ─────────────────────────────────────────────────────────────────────────────
# Step 9.2: Intersect motif sites with peak beds
# ─────────────────────────────────────────────────────────────────────────────
rule bedtools_intersect:
    input:
        fimo_bed = OUT + "/fimo/{sample}/fimo.bed",
        macs_bed = OUT + "/compare_bed/{sample}.MACS.bed",
        gem_bed  = OUT + "/compare_bed/{sample}.GEM.bed",
    output:
        macs_peak    = OUT + "/compare_bed/{sample}.MACS_peak.bed",
        gem_peak     = OUT + "/compare_bed/{sample}.GEM_peak.bed",
        compare_peak = OUT + "/compare_bed/{sample}.compare_peak.bed",
    resources:
        mem_mb=4000, runtime=30,
        slurm_partition=config["slurm_partition"],
        slurm_account=config["slurm_account"],
    log:
        OUT + "/logs/bedtools_intersect/{sample}.log"
    shell:
        """
        bedtools intersect -wa -a {input.fimo_bed} -b {input.macs_bed} > {output.macs_peak} 2>{log}
        bedtools intersect -wa -a {input.fimo_bed} -b {input.gem_bed}  > {output.gem_peak}  2>>{log}
        bedtools intersect -wa -a {input.macs_bed} -b {input.gem_bed}  > {output.compare_peak} 2>>{log}
        """


# ─────────────────────────────────────────────────────────────────────────────
# Step 11: Union of MACS and GEM motif hits (intersection.bed)
# ─────────────────────────────────────────────────────────────────────────────
rule motif_intersect:
    input:
        fimo_bed = OUT + "/fimo/{sample}/fimo.bed",
        macs_bed = OUT + "/compare_bed/{sample}.MACS.bed",
        gem_bed  = OUT + "/compare_bed/{sample}.GEM.bed",
    output:
        OUT + "/compare_bed/{sample}.intersection.bed"
    resources:
        mem_mb=4000, runtime=30,
        slurm_partition=config["slurm_partition"],
        slurm_account=config["slurm_account"],
    log:
        OUT + "/logs/motif_intersect/{sample}.log"
    shell:
        """
        bedtools intersect -wa -a {input.fimo_bed} -b {input.macs_bed} {input.gem_bed} \
          | sort -u > {output} 2>{log}
        """


# ─────────────────────────────────────────────────────────────────────────────
# Step 11.1: Extract FASTA for intersection motif sites
# ─────────────────────────────────────────────────────────────────────────────
rule getfasta_intersection:
    input:
        bed    = OUT + "/compare_bed/{sample}.intersection.bed",
        genome = config["genome_ref"],
    output:
        OUT + "/fasta/{sample}.motif.fasta"
    resources:
        mem_mb=4000, runtime=30,
        slurm_partition=config["slurm_partition"],
        slurm_account=config["slurm_account"],
    log:
        OUT + "/logs/getfasta_intersection/{sample}.log"
    shell:
        "bedtools getfasta -fo {output} -fi {input.genome} -bed {input.bed} 2>{log}"


# ─────────────────────────────────────────────────────────────────────────────
# Step 12: MEME motif refinement on intersection sequences
# ─────────────────────────────────────────────────────────────────────────────
rule meme_intersection:
    input:
        OUT + "/fasta/{sample}.motif.fasta"
    output:
        txt  = OUT + "/meme/{sample}-intersection/meme.txt",
        logo = OUT + "/meme/{sample}-intersection/logo1.png",
    params:
        outdir = OUT + "/meme/{sample}-intersection",
    threads:
        config["threads"]
    resources:
        mem_mb=8000, runtime=120,
        slurm_partition=config["slurm_partition"],
        slurm_account=config["slurm_account"],
    log:
        OUT + "/logs/meme_intersection/{sample}.log"
    shell:
        """
        meme {input} \
          -nmotifs 1 -minw 4 -maxw 12 \
          -dna -mod oops -nostatus \
          -p {threads} -oc {params.outdir} 2>{log}
        """
