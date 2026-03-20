rule trimmomatic_se:
    wildcard_constraints:
        sample = "|".join(SE_SAMPLES) if SE_SAMPLES else "(?!)",
    input:
        r1 = lambda wc: config["samples"][wc.sample]["r1"],
    output:
        r1 = OUT + "/trimmed/{sample}.R1.fastq.gz",
    params:
        adapters = config["trimmomatic"].get("adapters") or "",
    resources:
        mem_mb          = config["resources"]["trimmomatic"]["mem_mb"],
        runtime         = config["resources"]["trimmomatic"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
    log:
        OUT + "/logs/trimmomatic/{sample}.log"
    shell:
        """
        ADAPTERS={params.adapters:q}
        if [ -z "$ADAPTERS" ]; then
          ADAPTERS=$(find /opt/conda/envs/dapseq/share -name "TruSeq3-SE.fa" | head -1)
        elif [[ "$ADAPTERS" != /* ]]; then
          ADAPTERS=$(find /opt/conda/envs/dapseq/share -name "$ADAPTERS" | head -1)
        fi
        trimmomatic SE -phred33 \
          {input.r1} {output.r1} \
          ILLUMINACLIP:$ADAPTERS:2:30:10 \
          LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:{config[trimmomatic][minlen]} 2>{log}
        """


rule trimmomatic_pe:
    wildcard_constraints:
        sample = "|".join(PE_SAMPLES) if PE_SAMPLES else "(?!)",
    input:
        r1 = lambda wc: config["samples"][wc.sample]["r1"],
        r2 = lambda wc: config["samples"][wc.sample]["r2"],
    output:
        r1          = OUT + "/trimmed/{sample}.R1.fastq.gz",
        r2          = OUT + "/trimmed/{sample}.R2.fastq.gz",
        r1_unpaired = temp(OUT + "/trimmed/{sample}.R1.unpaired.fastq.gz"),
        r2_unpaired = temp(OUT + "/trimmed/{sample}.R2.unpaired.fastq.gz"),
    params:
        adapters = config["trimmomatic"].get("adapters") or "",
    resources:
        mem_mb          = config["resources"]["trimmomatic"]["mem_mb"],
        runtime         = config["resources"]["trimmomatic"]["runtime"],
        slurm_partition = config["slurm_partition"],
        slurm_account   = config["slurm_account"],
    log:
        OUT + "/logs/trimmomatic/{sample}.log"
    shell:
        """
        ADAPTERS={params.adapters:q}
        if [ -z "$ADAPTERS" ]; then
          ADAPTERS=$(find /opt/conda/envs/dapseq/share -name "TruSeq3-PE-2.fa" | head -1)
        elif [[ "$ADAPTERS" != /* ]]; then
          ADAPTERS=$(find /opt/conda/envs/dapseq/share -name "$ADAPTERS" | head -1)
        fi
        trimmomatic PE -phred33 \
          {input.r1} {input.r2} \
          {output.r1} {output.r1_unpaired} \
          {output.r2} {output.r2_unpaired} \
          ILLUMINACLIP:$ADAPTERS:2:30:10:1:true \
          LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:{config[trimmomatic][minlen]} 2>{log}
        """
