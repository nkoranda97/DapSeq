rule narrow_peak_to_fasta_summits:
    input:
        narrowpeak = OUT + "/MACS/{sample}.{bam_type}_peaks_filt.narrowPeak",
        genome     = config["genome_ref"],
    output:
        OUT + "/fasta/{sample}.{bam_type}.summits.fasta",
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
        OUT + "/logs/narrow_peak_to_fasta/{sample}.{bam_type}.summits.log"
    script:
        "../scripts/narrow_peak_to_fasta.py"


rule narrow_peak_to_fasta_peaks:
    input:
        narrowpeak = OUT + "/MACS/{sample}.{bam_type}_peaks_filt.narrowPeak",
        genome     = config["genome_ref"],
    output:
        OUT + "/fasta/{sample}.{bam_type}.peaks.fasta",
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
        OUT + "/logs/narrow_peak_to_fasta/{sample}.{bam_type}.peaks.log"
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
    run:
        import os

        if not (os.path.exists(str(input[0])) and os.path.getsize(str(input[0])) > 0):
            for f in [output.txt, output.xml, output.logo]:
                open(str(f), "w").close()
        else:
            active = {k: v for k, v in params.base_colors.items() if v}
            if active:
                defaults = {"A": "008000", "C": "0000FF", "G": "FFB300", "T": "FF0000"}
                def _hex(k):
                    return active.get(k, defaults[k]).lstrip("#")
                alph_file = os.path.join(params.outdir, "custom_alph.txt")
                os.makedirs(params.outdir, exist_ok=True)
                with open(alph_file, "w") as fh:
                    fh.write('ALPHABET "DNA"\n')
                    fh.write(f'A ADENINE {_hex("A")} ~ T THYMINE {_hex("T")}\n')
                    fh.write(f'C CYTOSINE {_hex("C")} ~ G GUANINE {_hex("G")}\n')
                    fh.write("END ALPHABET\n")
                alph_arg = f"-alph {alph_file}"
            else:
                alph_arg = "-dna"

            shell(
                f"meme {input[0]} -oc {params.outdir} "
                f"{alph_arg} -revcomp "
                f"-mod {params.mod} -nmotifs {params.nmotifs} "
                f"-minw {params.minw} -maxw {params.maxw} "
                f"-maxsize {params.maxsize} -p {threads} -nostatus {params.extra} "
                f"2>{log}"
            )
            if os.path.isfile(os.path.join(params.outdir, "logo1.eps")):
                shell(
                    f"gs -dNOPAUSE -dBATCH -sDEVICE=png16m -r150 "
                    f"-sOutputFile={output.logo} {params.outdir}/logo1.eps "
                    f"2>>{log}"
                )
            else:
                open(str(output.logo), "w").close()


rule meme_peaks:
    input:
        OUT + "/fasta/{sample}.{bam_type}.peaks.fasta",
    output:
        txt = OUT + "/meme/{sample}.{bam_type}/peaks/meme.txt",
        xml = OUT + "/meme/{sample}.{bam_type}/peaks/meme.xml",
    params:
        outdir  = OUT + "/meme/{sample}.{bam_type}/peaks",
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
        OUT + "/logs/meme/{sample}.{bam_type}.peaks.log"
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
        meme_xml = OUT + "/meme/{sample}.{bam_type}/summits/meme.xml",
        fasta    = OUT + "/fasta/{sample}.{bam_type}.summits.fasta",
    output:
        tsv = OUT + "/fimo/{sample}.{bam_type}/summits/fimo.tsv",
    params:
        outdir = OUT + "/fimo/{sample}.{bam_type}/summits",
        extra  = config["fimo"].get("extra", ""),
    resources:
        mem_mb          = config["resources"]["fimo"]["mem_mb"],
        runtime         = config["resources"]["fimo"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
    log:
        OUT + "/logs/fimo/{sample}.{bam_type}.summits.log"
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
        meme_xml = OUT + "/meme/{sample}.{bam_type}/peaks/meme.xml",
        fasta    = OUT + "/fasta/{sample}.{bam_type}.peaks.fasta",
    output:
        tsv = OUT + "/fimo/{sample}.{bam_type}/peaks/fimo.tsv",
    params:
        outdir = OUT + "/fimo/{sample}.{bam_type}/peaks",
        extra  = config["fimo"].get("extra", ""),
    resources:
        mem_mb          = config["resources"]["fimo"]["mem_mb"],
        runtime         = config["resources"]["fimo"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
    log:
        OUT + "/logs/fimo/{sample}.{bam_type}.peaks.log"
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
