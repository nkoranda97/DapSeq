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
        narrowpeak = OUT + "/MACS/{sample}_peaks_filt.narrowPeak",
        genome     = config["genome_ref"],
    output:
        OUT + "/fasta/{sample}.{peak_type}.fasta",
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
        OUT + "/logs/narrow_peak_to_fasta/{sample}.{peak_type}.log"
    script:
        "../scripts/narrow_peak_to_fasta.py"


rule meme_summits:
    input:
        OUT + "/fasta/{sample}.summits.fasta",
    output:
        txt = OUT + "/meme/{sample}/summits/meme.txt",
        xml = OUT + "/meme/{sample}/summits/meme.xml",
    params:
        outdir     = OUT + "/meme/{sample}/summits",
        nmotifs    = config["meme"]["nmotifs"],
        minw       = config["meme"]["minw"],
        maxw       = config["meme"]["maxw"],
        mod        = config["meme"]["mod"],
        maxsize    = config["meme"].get("maxsize", 10000000),
        extra      = config["meme"].get("extra", ""),
        has_custom = lambda wc: "1" if config["meme"].get("base_colors") else "",
        a_color    = lambda wc: ((config["meme"].get("base_colors") or {}).get("A") or "008000").lstrip("#"),
        c_color    = lambda wc: ((config["meme"].get("base_colors") or {}).get("C") or "0000FF").lstrip("#"),
        g_color    = lambda wc: ((config["meme"].get("base_colors") or {}).get("G") or "FFB300").lstrip("#"),
        t_color    = lambda wc: ((config["meme"].get("base_colors") or {}).get("T") or "FF0000").lstrip("#"),
    threads:
        config["threads"]
    resources:
        mem_mb          = config["resources"]["meme"]["mem_mb"],
        runtime         = config["resources"]["meme"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
    log:
        OUT + "/logs/meme/{sample}.summits.log"
    shell:
        """
        if [ ! -s {input} ]; then
            touch {output.txt} {output.xml}
        else
            mkdir -p {params.outdir}
            if [ -n "{params.has_custom}" ]; then
                printf 'ALPHABET "DNA"\\nA ADENINE {params.a_color} ~ T THYMINE {params.t_color}\\nC CYTOSINE {params.c_color} ~ G GUANINE {params.g_color}\\nEND ALPHABET\\n' \
                    > {params.outdir}/custom_alph.txt
                ALPH="-alph {params.outdir}/custom_alph.txt"
            else
                ALPH="-dna"
            fi
            meme {input} -oc {params.outdir} $ALPH -revcomp \
                -mod {params.mod} -nmotifs {params.nmotifs} \
                -minw {params.minw} -maxw {params.maxw} \
                -maxsize {params.maxsize} -p {threads} -nostatus {params.extra} \
                2>{log}
        fi
        """


rule meme_logo_summits:
    input:
        OUT + "/meme/{sample}/summits/meme.txt",
    output:
        logo    = OUT + "/meme/{sample}/summits/logo1.png",
        logo_rc = OUT + "/meme/{sample}/summits/logo_rc1.png",
    params:
        base_colors = config["meme"].get("base_colors") or {},
    resources:
        mem_mb          = config["resources"]["report"]["mem_mb"],
        runtime         = 10,
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
    log:
        OUT + "/logs/meme/{sample}.summits.logo.log"
    script:
        "../scripts/render_meme_logo.py"


rule meme_peaks:
    input:
        OUT + "/fasta/{sample}.peaks.fasta",
    output:
        txt = OUT + "/meme/{sample}/peaks/meme.txt",
        xml = OUT + "/meme/{sample}/peaks/meme.xml",
    params:
        outdir     = OUT + "/meme/{sample}/peaks",
        nmotifs    = config["meme"]["nmotifs"],
        minw       = config["meme"]["minw"],
        maxw       = config["meme"]["maxw"],
        mod        = config["meme"]["mod"],
        maxsize    = config["meme"].get("maxsize", 10000000),
        extra      = config["meme"].get("extra", ""),
        has_custom = lambda wc: "1" if config["meme"].get("base_colors") else "",
        a_color    = lambda wc: ((config["meme"].get("base_colors") or {}).get("A") or "008000").lstrip("#"),
        c_color    = lambda wc: ((config["meme"].get("base_colors") or {}).get("C") or "0000FF").lstrip("#"),
        g_color    = lambda wc: ((config["meme"].get("base_colors") or {}).get("G") or "FFB300").lstrip("#"),
        t_color    = lambda wc: ((config["meme"].get("base_colors") or {}).get("T") or "FF0000").lstrip("#"),
    threads:
        config["threads"]
    resources:
        mem_mb          = config["resources"]["meme"]["mem_mb"],
        runtime         = config["resources"]["meme"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
    log:
        OUT + "/logs/meme/{sample}.peaks.log"
    shell:
        """
        if [ ! -s {input} ]; then
            touch {output.txt} {output.xml}
        else
            mkdir -p {params.outdir}
            if [ -n "{params.has_custom}" ]; then
                printf 'ALPHABET "DNA"\\nA ADENINE {params.a_color} ~ T THYMINE {params.t_color}\\nC CYTOSINE {params.c_color} ~ G GUANINE {params.g_color}\\nEND ALPHABET\\n' \
                    > {params.outdir}/custom_alph.txt
                ALPH="-alph {params.outdir}/custom_alph.txt"
            else
                ALPH="-dna"
            fi
            meme {input} -oc {params.outdir} $ALPH -revcomp \
                -mod {params.mod} -nmotifs {params.nmotifs} \
                -minw {params.minw} -maxw {params.maxw} \
                -maxsize {params.maxsize} -p {threads} -nostatus {params.extra} \
                2>{log}
        fi
        """


rule meme_logo_peaks:
    input:
        OUT + "/meme/{sample}/peaks/meme.txt",
    output:
        logo    = OUT + "/meme/{sample}/peaks/logo1.png",
        logo_rc = OUT + "/meme/{sample}/peaks/logo_rc1.png",
    params:
        base_colors = config["meme"].get("base_colors") or {},
    resources:
        mem_mb          = config["resources"]["report"]["mem_mb"],
        runtime         = 10,
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
    log:
        OUT + "/logs/meme/{sample}.peaks.logo.log"
    script:
        "../scripts/render_meme_logo.py"


rule fimo:
    input:
        meme_xml = OUT + "/meme/{sample}/{peak_type}/meme.xml",
        fasta    = OUT + "/fasta/{sample}.{peak_type}.fasta",
    output:
        tsv = OUT + "/fimo/{sample}/{peak_type}/fimo.tsv",
    params:
        outdir      = OUT + "/fimo/{sample}/{peak_type}",
        extra       = config["fimo"].get("extra", ""),
        fimo_thresh = config["fimo"]["thresh"],
    resources:
        mem_mb          = config["resources"]["fimo"]["mem_mb"],
        runtime         = config["resources"]["fimo"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
    log:
        OUT + "/logs/fimo/{sample}.{peak_type}.log"
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
