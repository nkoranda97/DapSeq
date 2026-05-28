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
    run:
        import os
        import re
        import sys

        sys.path.insert(0, os.path.join(workflow.basedir, "scripts"))
        from logo_utils import _render_logo_with_logomaker

        tf = wildcards.sample.upper()

        # Collect all motif IDs for this TF
        tf_motif_ids = []
        seen = set()
        with open(input.tsv) as fh:
            header = fh.readline().rstrip('\n').split('\t')
            ti = header.index('target')
            ai = header.index('dataset_accession')
            ni = header.index('name')
            for line in fh:
                parts = line.rstrip('\n').split('\t')
                if parts[ti].upper() == tf:
                    mid = f"{parts[ai]}_{parts[ni]}"
                    if mid not in seen:
                        seen.add(mid)
                        tf_motif_ids.append(mid)

        if not tf_motif_ids:
            open(str(output[0]), 'w').close()
        else:
            # Parse catalog for all candidate motifs, pick the one with highest nsites
            tf_id_set    = set(tf_motif_ids)
            motif_ppms   = {}
            motif_nsites = {}
            current      = None
            in_matrix    = False
            with open(input.meme) as fh:
                for line in fh:
                    if line.startswith('MOTIF '):
                        mid       = line.split()[1]
                        current   = mid if mid in tf_id_set else None
                        in_matrix = False
                        if current:
                            motif_ppms[current]   = []
                            motif_nsites[current] = 0
                    elif current is not None:
                        if 'letter-probability matrix' in line:
                            m = re.search(r'nsites= *(\d+)', line)
                            if m:
                                motif_nsites[current] = int(m.group(1))
                            in_matrix = True
                        elif in_matrix:
                            try:
                                vals = [float(x) for x in line.split()]
                                if len(vals) == 4:
                                    motif_ppms[current].append(vals)
                                else:
                                    in_matrix = False
                            except ValueError:
                                in_matrix = False

            motif_id = max(
                (mid for mid in tf_motif_ids if motif_ppms.get(mid)),
                key=lambda mid: motif_nsites.get(mid, 0),
                default=None,
            )

            os.makedirs(os.path.dirname(str(output[0])), exist_ok=True)
            ppm = motif_ppms.get(motif_id) if motif_id else None
            if ppm:
                _render_logo_with_logomaker(ppm, str(output[0]), params.base_colors)
            else:
                open(str(output[0]), 'w').close()


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
            import sys

            sys.path.insert(0, os.path.join(workflow.basedir, "scripts"))
            from logo_utils import _render_logo_with_logomaker, parse_meme_ppm

            ppm = parse_meme_ppm(output.txt) if os.path.getsize(output.txt) > 0 else None
            if ppm:
                _render_logo_with_logomaker(ppm, str(output.logo), params.base_colors)
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
