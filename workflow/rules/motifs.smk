rule narrow_peak_to_fasta_summits:
    input:
        narrowpeak = OUT + "/MACS/{sample}_peaks_filt.narrowPeak",
        genome     = config["genome_ref"],
    output:
        OUT + "/fasta/{sample}.summits.fasta",
    params:
        maxpeaks   = config["meme"]["maxpeaks"],
        extend_bp  = config["meme"]["summit_extend"],
        fimocoords = True,
    resources:
        mem_mb          = config["resources"]["narrow_peak_to_fasta"]["mem_mb"],
        runtime         = config["resources"]["narrow_peak_to_fasta"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
    log:
        OUT + "/logs/narrow_peak_to_fasta/{sample}.summits.log"
    script:
        "../scripts/narrow_peak_to_fasta.py"


rule narrow_peak_to_fasta_peaks:
    input:
        narrowpeak = OUT + "/MACS/{sample}_peaks_filt.narrowPeak",
        genome     = config["genome_ref"],
    output:
        OUT + "/fasta/{sample}.peaks.fasta",
    params:
        maxpeaks   = config["meme"]["maxpeaks"],
        extend_bp  = "all",
        fimocoords = True,
    resources:
        mem_mb          = config["resources"]["narrow_peak_to_fasta"]["mem_mb"],
        runtime         = config["resources"]["narrow_peak_to_fasta"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
    log:
        OUT + "/logs/narrow_peak_to_fasta/{sample}.peaks.log"
    script:
        "../scripts/narrow_peak_to_fasta.py"


rule meme_summits:
    input:
        OUT + "/fasta/{sample}.summits.fasta",
    output:
        txt  = OUT + "/meme/{sample}/summits/meme.txt",
        xml  = OUT + "/meme/{sample}/summits/meme.xml",
        logo = OUT + "/meme/{sample}/summits/logo1.png",
    params:
        outdir  = OUT + "/meme/{sample}/summits",
        nmotifs = config["meme"]["nmotifs"],
        minw    = config["meme"]["minw"],
        maxw    = config["meme"]["maxw"],
        mod     = config["meme"]["mod"],
        maxsize = config["meme"].get("maxsize", 10000000),
        extra   = config["meme"].get("extra", ""),
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
            touch {output.txt} {output.xml} {output.logo}
        else
            meme {input} -oc {params.outdir} \
              -dna -revcomp \
              -mod {params.mod} -nmotifs {params.nmotifs} \
              -minw {params.minw} -maxw {params.maxw} \
              -maxsize {params.maxsize} -p {threads} -nostatus {params.extra} 2>{log}
            if [ -f {params.outdir}/logo1.eps ]; then
                gs -dNOPAUSE -dBATCH -sDEVICE=png16m -dEPSCrop -r150 \
                   -sOutputFile={output.logo} {params.outdir}/logo1.eps 2>>{log}
            else
                touch {output.logo}
            fi
        fi
        """


rule meme_peaks:
    input:
        OUT + "/fasta/{sample}.peaks.fasta",
    output:
        txt = OUT + "/meme/{sample}/peaks/meme.txt",
        xml = OUT + "/meme/{sample}/peaks/meme.xml",
    params:
        outdir  = OUT + "/meme/{sample}/peaks",
        nmotifs = config["meme"]["nmotifs"],
        minw    = config["meme"]["minw"],
        maxw    = config["meme"]["maxw"],
        mod     = config["meme"]["mod"],
        maxsize = config["meme"].get("maxsize", 10000000),
        extra   = config["meme"].get("extra", ""),
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
            meme {input} -oc {params.outdir} \
              -dna -revcomp \
              -mod {params.mod} -nmotifs {params.nmotifs} \
              -minw {params.minw} -maxw {params.maxw} \
              -maxsize {params.maxsize} -p {threads} -nostatus {params.extra} 2>{log}
        fi
        """


rule fimo_summits:
    input:
        meme_xml = OUT + "/meme/{sample}/summits/meme.xml",
        fasta    = OUT + "/fasta/{sample}.summits.fasta",
    output:
        tsv = OUT + "/fimo/{sample}/summits/fimo.tsv",
    params:
        outdir = OUT + "/fimo/{sample}/summits",
        extra  = config["fimo"].get("extra", ""),
    resources:
        mem_mb          = config["resources"]["fimo"]["mem_mb"],
        runtime         = config["resources"]["fimo"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
    log:
        OUT + "/logs/fimo/{sample}.summits.log"
    shell:
        """
        if [ ! -s {input.meme_xml} ]; then
            mkdir -p {params.outdir} && touch {output.tsv}
        else
            fimo --parse-genomic-coord --thresh {config[fimo][thresh]} \
              {params.extra} -oc {params.outdir} \
              {input.meme_xml} {input.fasta} 2>{log}
        fi
        """


rule fimo_peaks:
    input:
        meme_xml = OUT + "/meme/{sample}/peaks/meme.xml",
        fasta    = OUT + "/fasta/{sample}.peaks.fasta",
    output:
        tsv = OUT + "/fimo/{sample}/peaks/fimo.tsv",
    params:
        outdir = OUT + "/fimo/{sample}/peaks",
        extra  = config["fimo"].get("extra", ""),
    resources:
        mem_mb          = config["resources"]["fimo"]["mem_mb"],
        runtime         = config["resources"]["fimo"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
    log:
        OUT + "/logs/fimo/{sample}.peaks.log"
    shell:
        """
        if [ ! -s {input.meme_xml} ]; then
            mkdir -p {params.outdir} && touch {output.tsv}
        else
            fimo --parse-genomic-coord --thresh {config[fimo][thresh]} \
              {params.extra} -oc {params.outdir} \
              {input.meme_xml} {input.fasta} 2>{log}
        fi
        """
