if config.get("gene_annotation"):
    rule homer_annotate:
        input:
            peaks  = OUT + "/MACS/{sample}.{bam_type}_peaks_filt.narrowPeak",
            genome = config["genome_ref"],
            gtf    = config["gene_annotation"],
        output:
            OUT + "/annotations/{sample}.{bam_type}.peak_annotations.txt"
        wildcard_constraints:
            sample   = "|".join(TREATMENT_SAMPLES) if TREATMENT_SAMPLES else "(?!)",
            bam_type = "full|chr",
        params:
            extra = config.get("homer", {}).get("extra", ""),
        resources:
            mem_mb          = config["resources"]["homer_annotate"]["mem_mb"],
            runtime         = config["resources"]["homer_annotate"]["runtime"],
            slurm_partition = config["slurm_partition"],
            slurm_account   = config["slurm_account"],
        log:
            OUT + "/logs/homer_annotate/{sample}.{bam_type}.log"
        shell:
            "annotatePeaks.pl {input.peaks} {input.genome} -gtf {input.gtf} {params.extra} > {output} 2>{log}"
