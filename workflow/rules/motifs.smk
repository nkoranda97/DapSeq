_FB_DIR  = os.path.join(workflow.basedir, "..", "factorbook")
_FB_TSV  = (config.get("factorbook", {}).get("tsv") or
            os.path.join(_FB_DIR, "factorbook_chip_seq_meme_motifs.tsv"))
_FB_MEME = (config.get("factorbook", {}).get("meme") or
            os.path.join(_FB_DIR, "factorbook_chip_seq_meme_motif_catalog.meme"))


rule factorbook_logo:
    input:
        tsv  = _FB_TSV,
        meme = _FB_MEME,
    output:
        OUT + "/factorbook/{sample}.logo.png",
    wildcard_constraints:
        sample = "|".join(TREATMENT_SAMPLES) if TREATMENT_SAMPLES else "(?!)",
    resources:
        mem_mb          = config["resources"]["meme"]["mem_mb"],
        runtime         = 20,
    params:
        base_colors = config["meme"].get("base_colors") or {},
    log:
        OUT + "/logs/factorbook/{sample}.log"
    script:
        "../scripts/run_factorbook_logo.py"


rule narrow_peak_to_fasta:
    input:
        narrowpeak = get_final_filtered_peaks,
        genome     = config["genome_ref"],
    output:
        OUT + "/fasta/{sample}.{peak_type}.fasta",
    params:
        maxpeaks              = config["meme"]["maxpeaks"],
        extend_bp             = lambda wc: config["meme"]["summit_extend"] if wc.peak_type == "summits" else "all",
        fimocoords            = True,
        filter_chroms         = config.get("chrom_filter", []),
    resources:
        mem_mb          = config["resources"]["narrow_peak_to_fasta"]["mem_mb"],
        runtime         = config["resources"]["narrow_peak_to_fasta"]["runtime"],
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
        markov_order = config["meme"].get("markov_order", 0),
        maxsize    = config["meme"].get("maxsize", 10000000),
        extra      = config["meme"].get("extra", ""),
        has_custom = "1" if config["meme"].get("base_colors") else "",
        a_color    = ((config["meme"].get("base_colors") or {}).get("A") or "008000").lstrip("#"),
        c_color    = ((config["meme"].get("base_colors") or {}).get("C") or "0000FF").lstrip("#"),
        g_color    = ((config["meme"].get("base_colors") or {}).get("G") or "FFB300").lstrip("#"),
        t_color    = ((config["meme"].get("base_colors") or {}).get("T") or "FF0000").lstrip("#"),
    threads:
        config["threads"]
    resources:
        mem_mb          = config["resources"]["meme"]["mem_mb"],
        runtime         = config["resources"]["meme"]["runtime"],
    log:
        OUT + "/logs/meme/{sample}.summits.log"
    shell:
        """
        if [ ! -s {input} ]; then
            touch {output.txt} {output.xml}
        else
            mkdir -p {params.outdir}
            MEME_IN=$(mktemp)
            awk -v MINW={params.minw} '/^>/ {{h=$0; next}} length($0) >= MINW {{print h; print}}' {input} > "$MEME_IN"
            if [ ! -s "$MEME_IN" ]; then
                touch {output.txt} {output.xml}
                rm -f "$MEME_IN"
            else
            if [ -n "{params.has_custom}" ]; then
                printf 'ALPHABET "DNA"\\nA ADENINE {params.a_color} ~ T THYMINE {params.t_color}\\nC CYTOSINE {params.c_color} ~ G GUANINE {params.g_color}\\nEND ALPHABET\\n' \
                    > {params.outdir}/custom_alph.txt
                ALPH="-alph {params.outdir}/custom_alph.txt"
            else
                ALPH="-dna"
            fi
            meme "$MEME_IN" -oc {params.outdir} $ALPH -revcomp \
                -mod {params.mod} -nmotifs {params.nmotifs} \
                -minw {params.minw} -maxw {params.maxw} \
                -markov_order {params.markov_order} \
                -maxsize {params.maxsize} -p {threads} -nostatus {params.extra} \
                2>{log}
            rm -f "$MEME_IN"
            fi
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
        mem_mb          = config["resources"]["meme"]["mem_mb"],
        runtime         = 10,
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
        markov_order = config["meme"].get("markov_order", 0),
        maxsize    = config["meme"].get("maxsize", 10000000),
        extra      = config["meme"].get("extra", ""),
        has_custom = "1" if config["meme"].get("base_colors") else "",
        a_color    = ((config["meme"].get("base_colors") or {}).get("A") or "008000").lstrip("#"),
        c_color    = ((config["meme"].get("base_colors") or {}).get("C") or "0000FF").lstrip("#"),
        g_color    = ((config["meme"].get("base_colors") or {}).get("G") or "FFB300").lstrip("#"),
        t_color    = ((config["meme"].get("base_colors") or {}).get("T") or "FF0000").lstrip("#"),
    threads:
        config["threads"]
    resources:
        mem_mb          = config["resources"]["meme"]["mem_mb"],
        runtime         = config["resources"]["meme"]["runtime"],
    log:
        OUT + "/logs/meme/{sample}.peaks.log"
    shell:
        """
        if [ ! -s {input} ]; then
            touch {output.txt} {output.xml}
        else
            mkdir -p {params.outdir}
            MEME_IN=$(mktemp)
            awk -v MINW={params.minw} '/^>/ {{h=$0; next}} length($0) >= MINW {{print h; print}}' {input} > "$MEME_IN"
            if [ ! -s "$MEME_IN" ]; then
                touch {output.txt} {output.xml}
                rm -f "$MEME_IN"
            else
            if [ -n "{params.has_custom}" ]; then
                printf 'ALPHABET "DNA"\\nA ADENINE {params.a_color} ~ T THYMINE {params.t_color}\\nC CYTOSINE {params.c_color} ~ G GUANINE {params.g_color}\\nEND ALPHABET\\n' \
                    > {params.outdir}/custom_alph.txt
                ALPH="-alph {params.outdir}/custom_alph.txt"
            else
                ALPH="-dna"
            fi
            meme "$MEME_IN" -oc {params.outdir} $ALPH -revcomp \
                -mod {params.mod} -nmotifs {params.nmotifs} \
                -minw {params.minw} -maxw {params.maxw} \
                -markov_order {params.markov_order} \
                -maxsize {params.maxsize} -p {threads} -nostatus {params.extra} \
                2>{log}
            rm -f "$MEME_IN"
            fi
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
        mem_mb          = config["resources"]["meme"]["mem_mb"],
        runtime         = 10,
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


rule narrow_peak_to_fasta_branch:
    input:
        narrowpeak = get_branch_peak_file,
        genome     = config["genome_ref"],
    output:
        OUT + "/fasta/{sample}.{filter_stage}.{peak_type}.fasta",
    params:
        maxpeaks      = config["meme"]["maxpeaks"],
        extend_bp     = lambda wc: config["meme"]["summit_extend"] if wc.peak_type == "summits" else "all",
        fimocoords    = True,
        filter_chroms = config.get("chrom_filter", []),
    resources:
        mem_mb  = config["resources"]["narrow_peak_to_fasta"]["mem_mb"],
        runtime = config["resources"]["narrow_peak_to_fasta"]["runtime"],
    log:
        OUT + "/logs/narrow_peak_to_fasta/{sample}.{filter_stage}.{peak_type}.log"
    script:
        "../scripts/narrow_peak_to_fasta.py"


rule meme_branch:
    input:
        OUT + "/fasta/{sample}.{filter_stage}.{peak_type}.fasta",
    output:
        txt = OUT + "/meme/{sample}/{filter_stage}/{peak_type}/meme.txt",
        xml = OUT + "/meme/{sample}/{filter_stage}/{peak_type}/meme.xml",
    params:
        outdir       = OUT + "/meme/{sample}/{filter_stage}/{peak_type}",
        nmotifs      = config["meme"]["nmotifs"],
        minw         = config["meme"]["minw"],
        maxw         = config["meme"]["maxw"],
        mod          = config["meme"]["mod"],
        markov_order = config["meme"].get("markov_order", 0),
        maxsize      = config["meme"].get("maxsize", 10000000),
        extra        = config["meme"].get("extra", ""),
        has_custom   = "1" if config["meme"].get("base_colors") else "",
        a_color      = ((config["meme"].get("base_colors") or {}).get("A") or "008000").lstrip("#"),
        c_color      = ((config["meme"].get("base_colors") or {}).get("C") or "0000FF").lstrip("#"),
        g_color      = ((config["meme"].get("base_colors") or {}).get("G") or "FFB300").lstrip("#"),
        t_color      = ((config["meme"].get("base_colors") or {}).get("T") or "FF0000").lstrip("#"),
    threads:
        config["threads"]
    resources:
        mem_mb  = config["resources"]["meme"]["mem_mb"],
        runtime = config["resources"]["meme"]["runtime"],
    log:
        OUT + "/logs/meme/{sample}.{filter_stage}.{peak_type}.log"
    shell:
        """
        if [ ! -s {input} ]; then
            touch {output.txt} {output.xml}
        else
            mkdir -p {params.outdir}
            MEME_IN=$(mktemp)
            awk -v MINW={params.minw} '/^>/ {{h=$0; next}} length($0) >= MINW {{print h; print}}' {input} > "$MEME_IN"
            if [ ! -s "$MEME_IN" ]; then
                touch {output.txt} {output.xml}
                rm -f "$MEME_IN"
            else
            if [ -n "{params.has_custom}" ]; then
                printf 'ALPHABET "DNA"\\nA ADENINE {params.a_color} ~ T THYMINE {params.t_color}\\nC CYTOSINE {params.c_color} ~ G GUANINE {params.g_color}\\nEND ALPHABET\\n' \
                    > {params.outdir}/custom_alph.txt
                ALPH="-alph {params.outdir}/custom_alph.txt"
            else
                ALPH="-dna"
            fi
            meme "$MEME_IN" -oc {params.outdir} $ALPH -revcomp \
                -mod {params.mod} -nmotifs {params.nmotifs} \
                -minw {params.minw} -maxw {params.maxw} \
                -markov_order {params.markov_order} \
                -maxsize {params.maxsize} -p {threads} -nostatus {params.extra} \
                2>{log}
            rm -f "$MEME_IN"
            fi
        fi
        """


rule meme_logo_branch:
    input:
        OUT + "/meme/{sample}/{filter_stage}/{peak_type}/meme.txt",
    output:
        logo    = OUT + "/meme/{sample}/{filter_stage}/{peak_type}/logo1.png",
        logo_rc = OUT + "/meme/{sample}/{filter_stage}/{peak_type}/logo_rc1.png",
    params:
        base_colors = config["meme"].get("base_colors") or {},
    resources:
        mem_mb  = config["resources"]["meme"]["mem_mb"],
        runtime = 10,
    log:
        OUT + "/logs/meme/{sample}.{filter_stage}.{peak_type}.logo.log"
    script:
        "../scripts/render_meme_logo.py"


rule fimo_branch:
    input:
        meme_xml = OUT + "/meme/{sample}/{filter_stage}/{peak_type}/meme.xml",
        fasta    = OUT + "/fasta/{sample}.{filter_stage}.{peak_type}.fasta",
    output:
        tsv = OUT + "/fimo/{sample}/{filter_stage}/{peak_type}/fimo.tsv",
    params:
        outdir      = OUT + "/fimo/{sample}/{filter_stage}/{peak_type}",
        extra       = config["fimo"].get("extra", ""),
        fimo_thresh = config["fimo"]["thresh"],
    resources:
        mem_mb  = config["resources"]["fimo"]["mem_mb"],
        runtime = config["resources"]["fimo"]["runtime"],
    log:
        OUT + "/logs/fimo/{sample}.{filter_stage}.{peak_type}.log"
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
