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
GENOME_SPLIT_DIR = OUT + "/genome_split"

PEAK_CALLER = config.get("peak_caller", "both")
USE_MACS    = PEAK_CALLER in ("both", "macs3")
USE_GEM     = PEAK_CALLER in ("both", "gem")
USE_BOTH    = PEAK_CALLER == "both"

ALIGNER = config.get("aligner", "bowtie2")

SE_SAMPLES = {s for s in SAMPLES if config["samples"][s].get("r2") is None}
PE_SAMPLES = {s for s in SAMPLES if config["samples"][s].get("r2") is not None}

wildcard_constraints:
    sample = "[^/.]+",
    read   = "R[12]",

# ─────────────────────────────────────────────────────────────────────────────
# rule all
# ─────────────────────────────────────────────────────────────────────────────
rule all:
    input:
        config["genome_ref"] + ".sizes",
        expand(OUT + "/Fastqc/{sample}.R1_fastqc.html", sample=SAMPLES),
        expand(OUT + "/Fastqc/{sample}.R2_fastqc.html", sample=PE_SAMPLES),
        OUT + "/multiqc_report.html",
        expand(OUT + "/bigWig/{sample}.bw",                      sample=SAMPLES),
        *(expand(OUT + "/bigWig/{sample}.peaks.bw",              sample=TREATMENT_SAMPLES) if CONTROL else []),
        *(expand(OUT + "/compare_bed/{sample}.MACS_peak.bed",    sample=TREATMENT_SAMPLES) if USE_MACS else []),
        *(expand(OUT + "/compare_bed/{sample}.GEM_peak.bed",     sample=TREATMENT_SAMPLES) if USE_GEM  else []),
        *(expand(OUT + "/compare_bed/{sample}.compare_peak.bed", sample=TREATMENT_SAMPLES) if USE_BOTH else []),
        expand(OUT + "/meme/{sample}-intersection/meme.txt",      sample=TREATMENT_SAMPLES),

# ─────────────────────────────────────────────────────────────────────────────
# Checkpoint targets — run snakemake <target> to stop at a natural stage
# ─────────────────────────────────────────────────────────────────────────────
rule mapped:
    input:
        expand(OUT + "/bam/{sample}.bam",     sample=SAMPLES),
        expand(OUT + "/bam/{sample}.bam.bai", sample=SAMPLES),

rule peaked:
    input:
        expand(OUT + "/combined_bed/{sample}.combined.bed", sample=TREATMENT_SAMPLES),

rule motifs:
    input:
        expand(OUT + "/meme/{sample}/meme.txt",  sample=TREATMENT_SAMPLES),
        expand(OUT + "/fimo/{sample}/fimo.tsv",  sample=TREATMENT_SAMPLES),


# ─────────────────────────────────────────────────────────────────────────────
# Step 0: Index reference genome
# ─────────────────────────────────────────────────────────────────────────────
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


# ─────────────────────────────────────────────────────────────────────────────
# Step 0.1: Split genome into per-chromosome FASTAs for GEM
# ─────────────────────────────────────────────────────────────────────────────
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


# ─────────────────────────────────────────────────────────────────────────────
# Step 1: Trim reads — SE and PE variants
# ─────────────────────────────────────────────────────────────────────────────
rule trimmomatic_se:
    wildcard_constraints:
        sample = "|".join(SE_SAMPLES) if SE_SAMPLES else "(?!)",
    input:
        r1 = lambda wc: config["samples"][wc.sample]["r1"],
    output:
        r1 = OUT + "/trimmed/{sample}.R1.fastq.gz",
    params:
        adapters = config["trimmomatic"].get("adapters") or "",
    resources:
        mem_mb          = config["resources"]["trimmomatic"]["mem_mb"],
        runtime         = config["resources"]["trimmomatic"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
    log:
        OUT + "/logs/trimmomatic/{sample}.log"
    shell:
        """
        ADAPTERS={params.adapters:q}
        if [ -z "$ADAPTERS" ]; then
          ADAPTERS=$(find /opt/conda/envs/dapseq/share -name "TruSeq3-SE.fa" | head -1)
        elif [[ "$ADAPTERS" != /* ]]; then
          ADAPTERS=$(find /opt/conda/envs/dapseq/share -name "$ADAPTERS" | head -1)
        fi
        trimmomatic SE -phred33 \
          {input.r1} {output.r1} \
          ILLUMINACLIP:$ADAPTERS:2:30:10 \
          LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:{config[trimmomatic][minlen]} 2>{log}
        """


rule trimmomatic_pe:
    wildcard_constraints:
        sample = "|".join(PE_SAMPLES) if PE_SAMPLES else "(?!)",
    input:
        r1 = lambda wc: config["samples"][wc.sample]["r1"],
        r2 = lambda wc: config["samples"][wc.sample]["r2"],
    output:
        r1          = OUT + "/trimmed/{sample}.R1.fastq.gz",
        r2          = OUT + "/trimmed/{sample}.R2.fastq.gz",
        r1_unpaired = temp(OUT + "/trimmed/{sample}.R1.unpaired.fastq.gz"),
        r2_unpaired = temp(OUT + "/trimmed/{sample}.R2.unpaired.fastq.gz"),
    params:
        adapters = config["trimmomatic"].get("adapters") or "",
    resources:
        mem_mb          = config["resources"]["trimmomatic"]["mem_mb"],
        runtime         = config["resources"]["trimmomatic"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
    log:
        OUT + "/logs/trimmomatic/{sample}.log"
    shell:
        """
        ADAPTERS={params.adapters:q}
        if [ -z "$ADAPTERS" ]; then
          ADAPTERS=$(find /opt/conda/envs/dapseq/share -name "TruSeq3-PE-2.fa" | head -1)
        elif [[ "$ADAPTERS" != /* ]]; then
          ADAPTERS=$(find /opt/conda/envs/dapseq/share -name "$ADAPTERS" | head -1)
        fi
        trimmomatic PE -phred33 \
          {input.r1} {input.r2} \
          {output.r1} {output.r1_unpaired} \
          {output.r2} {output.r2_unpaired} \
          ILLUMINACLIP:$ADAPTERS:2:30:10:1:true \
          LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:{config[trimmomatic][minlen]} 2>{log}
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
        mem_mb          = config["resources"]["fastqc"]["mem_mb"],
        runtime         = config["resources"]["fastqc"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
    log:
        OUT + "/logs/fastqc/{sample}.{read}.log"
    shell:
        "fastqc {input} --outdir={OUT}/Fastqc/ 2>{log}"


# ─────────────────────────────────────────────────────────────────────────────
# Step 2.1: MultiQC
# ─────────────────────────────────────────────────────────────────────────────
rule multiqc:
    input:
        expand(OUT + "/Fastqc/{sample}.R1_fastqc.zip", sample=SAMPLES),
        expand(OUT + "/Fastqc/{sample}.R2_fastqc.zip", sample=PE_SAMPLES),
    output:
        OUT + "/multiqc_report.html"
    resources:
        mem_mb          = config["resources"]["multiqc"]["mem_mb"],
        runtime         = config["resources"]["multiqc"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
    log:
        OUT + "/logs/multiqc.log"
    shell:
        "multiqc {OUT}/Fastqc/ -o {OUT}/ 2>{log}"


# ─────────────────────────────────────────────────────────────────────────────
# Step 3: Align reads (ALL samples) — aligner selected by config["aligner"]
# ─────────────────────────────────────────────────────────────────────────────
if ALIGNER == "bowtie2":
    rule bowtie2:
        input:
            r1  = OUT + "/trimmed/{sample}.R1.fastq.gz",
            r2  = lambda wc: [] if wc.sample in SE_SAMPLES else [OUT + f"/trimmed/{wc.sample}.R2.fastq.gz"],
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
            r2  = lambda wc: [] if wc.sample in SE_SAMPLES else [OUT + f"/trimmed/{wc.sample}.R2.fastq.gz"],
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
            "bwa mem -t {threads} -k {params.k} -B {params.B} -O {params.O} -v 2 {config[genome_ref]} {input.r1} {params.r2_arg} > {output} 2>{log}"


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
        is_se  = lambda wc: "true" if wc.sample in SE_SAMPLES else "false",
    threads:
        config["threads"]
    resources:
        mem_mb          = config["resources"]["samtools_filter_sort_dedup"]["mem_mb"],
        runtime         = config["resources"]["samtools_filter_sort_dedup"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
    log:
        OUT + "/logs/samtools/{sample}.log"
    shell:
        """
        if [ "{params.is_se}" = "true" ]; then
          samtools view -@ {threads} -h -F 4 -q {config[samtools][mapq]} -u {input} \
            | samtools sort -@ {threads} -o {params.prefix} - 2>>{log}
        else
          samtools view -@ {threads} -h -F 4 -q {config[samtools][mapq]} -u {input} \
            | samtools fixmate -m -@ {threads} - - \
            | samtools sort -@ {threads} -o {params.prefix} - 2>>{log}
        fi
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
        mem_mb          = config["resources"]["bamcoverage"]["mem_mb"],
        runtime         = config["resources"]["bamcoverage"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
    log:
        OUT + "/logs/bamcoverage/{sample}.log"
    shell:
        "bamCoverage -b {input.bam} --normalizeUsing {config[bamcoverage][normalize_using]} --extendReads {config[bamcoverage][extend_reads]} -p {threads} -o {output} 2>{log}"


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
            mem_mb          = config["resources"]["bamcompare"]["mem_mb"],
            runtime         = config["resources"]["bamcompare"]["runtime"],
            slurm_partition = config["slurm_partition"],
            slurm_account   = config["slurm_account"],
        log:
            OUT + "/logs/bamcompare/{sample}.log"
        shell:
            """
            bamCompare \
              -b1 {input.sample_bam} -b2 {input.control_bam} \
              -o {output} \
              --binSize {config[bamcompare][bin_size]} --operation {config[bamcompare][operation]} \
              --scaleFactorsMethod {config[bamcompare][scale_factors_method]} -n {config[bamcompare][n]} \
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
        mem_mb          = config["resources"]["macs3"]["mem_mb"],
        runtime         = config["resources"]["macs3"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
    log:
        OUT + "/logs/macs3/{sample}.log"
    shell:
        """
        macs3 callpeak \
          -t {input.sample_bam} {params.ctrl} \
          -f {config[macs3][format]} --outdir {params.outdir} \
          -g {config[genome_size]} -n {wildcards.sample} \
          -B -q {config[macs3][qvalue]} -m {config[macs3][min_fold]} {config[macs3][max_fold]} --verbose=0 2>{log}
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
        out_prefix = OUT + "/GEM/{sample}",
    threads:
        config["threads"]
    resources:
        mem_mb          = config["resources"]["gem"]["mem_mb"],
        runtime         = config["resources"]["gem"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
    log:
        OUT + "/logs/gem/{sample}.log"
    shell:
        """
        mkdir -p {params.out_prefix}
        java -Xmx{config[gem][xmx]} -jar {GEM_JAR} \
          --t {threads} --f BAM \
          --d {GEM_READ_DIST} \
          --g {input.chrom_sizes} \
          --genome {input.genome_dir} \
          --expt {input.sample_bam} {params.ctrl} \
          --outBED \
          --out {params.out_prefix} \
          --k_min {config[gem][k_min]} --k_max {config[gem][k_max]} --k_seqs {config[gem][k_seqs]} --k_neg_dinu_shuffle 2>{log}
        """


# ─────────────────────────────────────────────────────────────────────────────
# Step 7: Combine MACS3 + GEM peaks
# ─────────────────────────────────────────────────────────────────────────────
rule combine_peaks:
    input:
        macs_summits = (OUT + "/MACS/{sample}_summits.bed"            if USE_MACS else []),
        gem_events   = (OUT + "/GEM/{sample}/{sample}.GEM_events.txt" if USE_GEM  else []),
    output:
        combined_bed = OUT + "/combined_bed/{sample}.combined.bed",
        macs_bed     = (OUT + "/compare_bed/{sample}.MACS.bed" if USE_MACS else []),
        gem_bed      = (OUT + "/compare_bed/{sample}.GEM.bed"  if USE_GEM  else []),
    params:
        window_size = config["combine_peaks"]["window_size"],
        min_score   = config["combine_peaks"]["min_score"],
        peak_caller = PEAK_CALLER,
    resources:
        mem_mb          = config["resources"]["combine_peaks"]["mem_mb"],
        runtime         = config["resources"]["combine_peaks"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
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
        mem_mb          = config["resources"]["getfasta"]["mem_mb"],
        runtime         = config["resources"]["getfasta"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
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
        mem_mb          = config["resources"]["dedup_fasta"]["mem_mb"],
        runtime         = config["resources"]["dedup_fasta"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
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
        k     = config["tandem_filter"]["k"],
        k_max = config["tandem_filter"]["k_max"],
    resources:
        mem_mb          = config["resources"]["tandem_filter"]["mem_mb"],
        runtime         = config["resources"]["tandem_filter"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
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
    params:
        outdir  = OUT + "/meme/{sample}",
        nmotifs = config["meme"]["nmotifs"],
        minw    = config["meme"]["minw"],
        maxw    = config["meme"]["maxw"],
        mod     = config["meme"]["mod"],
    threads:
        config["threads"]
    resources:
        mem_mb          = config["resources"]["meme"]["mem_mb"],
        runtime         = config["resources"]["meme"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
    log:
        OUT + "/logs/meme/{sample}.log"
    shell:
        "bash scripts/run_meme.sh {input} {params.outdir} {threads} {log} 'filtered FASTA (step 8)' {params.nmotifs} {params.minw} {params.maxw} {params.mod}"


# ─────────────────────────────────────────────────────────────────────────────
# Step 9: FIMO motif scanning
# ─────────────────────────────────────────────────────────────────────────────
rule fimo:
    input:
        meme_txt = OUT + "/meme/{sample}/meme.txt",
        genome   = config["genome_ref"],
    output:
        tsv = OUT + "/fimo/{sample}/fimo.tsv",
        gff = OUT + "/fimo/{sample}/fimo.gff",
    params:
        outdir = OUT + "/fimo/{sample}",
    resources:
        mem_mb          = config["resources"]["fimo"]["mem_mb"],
        runtime         = config["resources"]["fimo"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
    log:
        OUT + "/logs/fimo/{sample}.log"
    shell:
        """
        fimo -verbosity 1 --thresh {config[fimo][thresh]} \
          -oc {params.outdir} \
          {input.meme_txt} {input.genome} 2>{log}
        """
# ─────────────────────────────────────────────────────────────────────────────
# Step 9.2: FIMO to Bed
# ─────────────────────────────────────────────────────────────────────────────
rule fimo_to_bed:
    input:
        OUT + "/fimo/{sample}/fimo.tsv"
    output:
        OUT + "/fimo/{sample}/fimo.bed"
    resources:
        mem_mb          = config["resources"]["fimo_to_bed"]["mem_mb"],
        runtime         = config["resources"]["fimo_to_bed"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
    script:
        "scripts/fimo_to_bed.py"

# ─────────────────────────────────────────────────────────────────────────────
# Step 9.3: Intersect motif sites with peak beds
# ─────────────────────────────────────────────────────────────────────────────
if USE_MACS:
    rule fimo_macs_intersect:
        input:
            fimo_bed = OUT + "/fimo/{sample}/fimo.bed",
            macs_bed = OUT + "/compare_bed/{sample}.MACS.bed",
        output:
            OUT + "/compare_bed/{sample}.MACS_peak.bed",
        resources:
            mem_mb          = config["resources"]["bedtools_intersect"]["mem_mb"],
            runtime         = config["resources"]["bedtools_intersect"]["runtime"],
            slurm_partition = config["slurm_partition"],
            slurm_account   = config["slurm_account"],
        log:
            OUT + "/logs/bedtools_intersect/{sample}.MACS.log"
        shell:
            "bedtools intersect -wa -a {input.fimo_bed} -b {input.macs_bed} > {output} 2>{log}"

if USE_GEM:
    rule fimo_gem_intersect:
        input:
            fimo_bed = OUT + "/fimo/{sample}/fimo.bed",
            gem_bed  = OUT + "/compare_bed/{sample}.GEM.bed",
        output:
            OUT + "/compare_bed/{sample}.GEM_peak.bed",
        resources:
            mem_mb          = config["resources"]["bedtools_intersect"]["mem_mb"],
            runtime         = config["resources"]["bedtools_intersect"]["runtime"],
            slurm_partition = config["slurm_partition"],
            slurm_account   = config["slurm_account"],
        log:
            OUT + "/logs/bedtools_intersect/{sample}.GEM.log"
        shell:
            "bedtools intersect -wa -a {input.fimo_bed} -b {input.gem_bed} > {output} 2>{log}"

if USE_BOTH:
    rule compare_callers_intersect:
        input:
            macs_bed = OUT + "/compare_bed/{sample}.MACS.bed",
            gem_bed  = OUT + "/compare_bed/{sample}.GEM.bed",
        output:
            OUT + "/compare_bed/{sample}.compare_peak.bed",
        resources:
            mem_mb          = config["resources"]["bedtools_intersect"]["mem_mb"],
            runtime         = config["resources"]["bedtools_intersect"]["runtime"],
            slurm_partition = config["slurm_partition"],
            slurm_account   = config["slurm_account"],
        log:
            OUT + "/logs/bedtools_intersect/{sample}.compare.log"
        shell:
            "bedtools intersect -wa -a {input.macs_bed} -b {input.gem_bed} > {output} 2>{log}"


# ─────────────────────────────────────────────────────────────────────────────
# Step 11: Union of MACS and GEM motif hits (intersection.bed)
# ─────────────────────────────────────────────────────────────────────────────
rule motif_intersect:
    input:
        fimo_bed  = OUT + "/fimo/{sample}/fimo.bed",
        peak_beds = (
            [OUT + "/compare_bed/{sample}.MACS.bed",
             OUT + "/compare_bed/{sample}.GEM.bed"] if USE_BOTH else
            [OUT + "/compare_bed/{sample}.MACS.bed"] if USE_MACS else
            [OUT + "/compare_bed/{sample}.GEM.bed"]
        ),
    output:
        OUT + "/compare_bed/{sample}.intersection.bed"
    resources:
        mem_mb          = config["resources"]["motif_intersect"]["mem_mb"],
        runtime         = config["resources"]["motif_intersect"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
    log:
        OUT + "/logs/motif_intersect/{sample}.log"
    shell:
        """
        bedtools intersect -wa -a {input.fimo_bed} -b {input.peak_beds} 2>{log} \
          | sort -u > {output}
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
        mem_mb          = config["resources"]["getfasta_intersection"]["mem_mb"],
        runtime         = config["resources"]["getfasta_intersection"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
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
    params:
        outdir  = OUT + "/meme/{sample}-intersection",
        nmotifs = config["meme"]["nmotifs"],
        minw    = config["meme"]["minw"],
        maxw    = config["meme"]["maxw"],
        mod     = config["meme"]["mod"],
    threads:
        config["threads"]
    resources:
        mem_mb          = config["resources"]["meme_intersection"]["mem_mb"],
        runtime         = config["resources"]["meme_intersection"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
    log:
        OUT + "/logs/meme_intersection/{sample}.log"
    shell:
        "bash scripts/run_meme.sh {input} {params.outdir} {threads} {log} 'intersection FASTA (step 12)' {params.nmotifs} {params.minw} {params.maxw} {params.mod}"
