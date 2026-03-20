if USE_MACS:
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

if USE_GEM:
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
        "../scripts/combine_peaks.py"
