rule factorbook_logo:
    input:
        tsv  = "factorbook/factorbook_chip_seq_meme_motifs.tsv",
        meme = "factorbook/factorbook_chip_seq_meme_motif_catalog.meme",
    output:
        OUT + "/factorbook/{sample}.logo.png",
    wildcard_constraints:
        sample = "|".join(TREATMENT_SAMPLES) if TREATMENT_SAMPLES else "(?!)",
    resources:
        mem_mb          = config["resources"]["report"]["mem_mb"],
        runtime         = 20,
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
    params:
        base_colors = config["meme"].get("base_colors") or {},
    log:
        OUT + "/logs/factorbook/{sample}.log"
    script:
        "../scripts/run_factorbook_logo.py"


rule narrow_peak_to_fasta:
    input:
        narrowpeak = OUT + "/MACS/{sample}.{bam_type}_peaks_filt.narrowPeak",
        genome     = config["genome_ref"],
    output:
        OUT + "/fasta/{sample}.{bam_type}.{peak_type}.fasta",
    params:
        maxpeaks   = config["meme"]["maxpeaks"],
        extend_bp  = lambda wc: config["meme"]["summit_extend"] if wc.peak_type == "summits" else "all",
        fimocoords = True,
    resources:
        mem_mb          = config["resources"]["narrow_peak_to_fasta"]["mem_mb"],
        runtime         = config["resources"]["narrow_peak_to_fasta"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
    log:
        OUT + "/logs/narrow_peak_to_fasta/{sample}.{bam_type}.{peak_type}.log"
    script:
        "../scripts/narrow_peak_to_fasta.py"


rule meme_summits:
    input:
        OUT + "/fasta/{sample}.{bam_type}.summits.fasta",
    output:
        txt  = OUT + "/meme/{sample}.{bam_type}/summits/meme.txt",
        xml  = OUT + "/meme/{sample}.{bam_type}/summits/meme.xml",
        logo = OUT + "/meme/{sample}.{bam_type}/summits/logo1.png",
    params:
        outdir      = OUT + "/meme/{sample}.{bam_type}/summits",
        nmotifs     = config["meme"]["nmotifs"],
        minw        = config["meme"]["minw"],
        maxw        = config["meme"]["maxw"],
        mod         = config["meme"]["mod"],
        maxsize     = config["meme"].get("maxsize", 10000000),
        extra       = config["meme"].get("extra", ""),
        base_colors = config["meme"].get("base_colors") or {},
    threads:
        config["threads"]
    resources:
        mem_mb          = config["resources"]["meme"]["mem_mb"],
        runtime         = config["resources"]["meme"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
    log:
        OUT + "/logs/meme/{sample}.{bam_type}.summits.log"
    script:
        "../scripts/run_meme.py"


rule meme_peaks:
    input:
        OUT + "/fasta/{sample}.{bam_type}.peaks.fasta",
    output:
        txt  = OUT + "/meme/{sample}.{bam_type}/peaks/meme.txt",
        xml  = OUT + "/meme/{sample}.{bam_type}/peaks/meme.xml",
        logo = OUT + "/meme/{sample}.{bam_type}/peaks/logo1.png",
    params:
        outdir      = OUT + "/meme/{sample}.{bam_type}/peaks",
        nmotifs     = config["meme"]["nmotifs"],
        minw        = config["meme"]["minw"],
        maxw        = config["meme"]["maxw"],
        mod         = config["meme"]["mod"],
        maxsize     = config["meme"].get("maxsize", 10000000),
        extra       = config["meme"].get("extra", ""),
        base_colors = config["meme"].get("base_colors") or {},
    threads:
        config["threads"]
    resources:
        mem_mb          = config["resources"]["meme"]["mem_mb"],
        runtime         = config["resources"]["meme"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
    log:
        OUT + "/logs/meme/{sample}.{bam_type}.peaks.log"
    script:
        "../scripts/run_meme.py"


rule fimo:
    input:
        meme_xml = OUT + "/meme/{sample}.{bam_type}/{peak_type}/meme.xml",
        fasta    = OUT + "/fasta/{sample}.{bam_type}.{peak_type}.fasta",
    output:
        tsv = OUT + "/fimo/{sample}.{bam_type}/{peak_type}/fimo.tsv",
    params:
        outdir      = OUT + "/fimo/{sample}.{bam_type}/{peak_type}",
        extra       = config["fimo"].get("extra", ""),
        fimo_thresh = config["fimo"]["thresh"],
    resources:
        mem_mb          = config["resources"]["fimo"]["mem_mb"],
        runtime         = config["resources"]["fimo"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
    log:
        OUT + "/logs/fimo/{sample}.{bam_type}.{peak_type}.log"
    shell:
        """
        if [ ! -s {input.meme_xml} ]; then
            mkdir -p {params.outdir} && touch {output.tsv}
        else
            fimo --parse-genomic-coord --thresh {params.fimo_thresh} \
              {params.extra} -oc {params.outdir} \
              {input.meme_xml} {input.fasta} 2>{log}
        fi
        """
