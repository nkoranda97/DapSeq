if config.get("gene_annotation"):
    rule homer_annotate:
        input:
            peaks  = get_final_filtered_peaks,
            genome = config["genome_ref"],
            gtf    = config["gene_annotation"],
        output:
            OUT + "/annotations/{sample}.peak_annotations.txt"
        wildcard_constraints:
            sample = REPORT_SAMPLE_CONSTRAINT,
        params:
            extra = config.get("homer", {}).get("extra", ""),
        resources:
            mem_mb          = config["resources"]["homer_annotate"]["mem_mb"],
            runtime         = config["resources"]["homer_annotate"]["runtime"],
        log:
            OUT + "/logs/homer_annotate/{sample}.log"
        shell:
            "annotatePeaks.pl {input.peaks} {input.genome} -gtf {input.gtf} {params.extra} > {output} 2>{log}"
