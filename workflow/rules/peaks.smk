rule macs3:
    input:
        sample_bam  = OUT + "/bam/{sample}.bam",
        control_bam = (OUT + f"/bam/{CONTROL}.bam" if CONTROL else []),
    output:
        summits    = OUT + "/MACS/{sample}_summits.bed",
        narrowpeak = OUT + "/MACS/{sample}_peaks.narrowPeak",
    params:
        ctrl     = lambda wc, input: f"-c {input.control_bam}" if CONTROL else "",
        outdir   = OUT + "/MACS",
        keep_dup = config["macs3"].get("keep_dup", 1),
        extra    = config["macs3"].get("extra", ""),
        tmpdir   = OUT + "/temp",
    resources:
        mem_mb          = config["resources"]["macs3"]["mem_mb"],
        runtime         = config["resources"]["macs3"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
    log:
        OUT + "/logs/macs3/{sample}.log"
    shell:
        """
        mkdir -p {params.tmpdir}
        export TMPDIR={params.tmpdir}
        macs3 callpeak \
          -t {input.sample_bam} {params.ctrl} \
          -f {config[macs3][format]} --outdir {params.outdir} \
          -g {config[genome_size]} -n {wildcards.sample} \
          --call-summits --keep-dup {params.keep_dup} {params.extra} --verbose=0 2>{log}
        """


rule filter_peaks:
    input:
        OUT + "/MACS/{sample}_peaks.narrowPeak",
    output:
        filt         = OUT + "/MACS/{sample}_peaks_filt.narrowPeak",
        numpeaks     = OUT + "/MACS/{sample}_numpeaks.txt",
        numpeaks_filt= OUT + "/MACS/{sample}_numpeaks_filt.txt",
    params:
        min_foldch = config["macs3"]["min_foldch"],
    resources:
        mem_mb          = config["resources"]["filter_peaks"]["mem_mb"],
        runtime         = config["resources"]["filter_peaks"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
    log:
        OUT + "/logs/filter_peaks/{sample}.log"
    shell:
        """
        awk -v FCH={params.min_foldch} '$7 >= FCH' {input} > {output.filt} 2>{log}
        wc -l < {input}       > {output.numpeaks}
        wc -l < {output.filt} > {output.numpeaks_filt}
        """
