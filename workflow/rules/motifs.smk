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
        "../scripts/dedup_fasta.py"


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
        "../scripts/tandem_filter.py"


rule meme:
    input:
        OUT + "/fasta/{sample}.fasta.filtered.fasta"
    output:
        txt = OUT + "/meme/{sample}/meme.txt",
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
        "bash {SCRIPTS}/run_meme.sh {input} {params.outdir} {threads} {log} 'filtered FASTA (step 8)' {params.nmotifs} {params.minw} {params.maxw} {params.mod}"


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
        "../scripts/fimo_to_bed.py"


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


rule meme_intersection:
    input:
        OUT + "/fasta/{sample}.motif.fasta"
    output:
        txt = OUT + "/meme/{sample}-intersection/meme.txt",
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
        "bash {SCRIPTS}/run_meme.sh {input} {params.outdir} {threads} {log} 'intersection FASTA (step 12)' {params.nmotifs} {params.minw} {params.maxw} {params.mod}"
